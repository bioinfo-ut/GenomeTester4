/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014 University of Tartu
 *
 * Authors: Maarja Lepamets and Lauris Kaplinski
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>

#include "utils.h"
#include "fasta.h"
#include "wordtable.h"
#include "sequence.h"
#include "queue.h"
#include "common.h"
#include "queue.h"

#include "trie.h"

#define MAX_FILES 200

#define DEFAULT_NUM_THREADS 16
#define DEFAULT_TABLE_SIZE 500000000
#define DEFAULT_MAX_TABLES 32

#define BSIZE 10000000

#define MAX_MERGED_TABLES 256

#define TIME_READ 0
#define TIME_SORT 1
#define TIME_MERGE 2
#define TIME_FF 3

unsigned int wsize = 32;

#define MAX_READ_SIZE 100000000

#define READ_BITS 5
#define READ_BLOCKS (1UL << READ_BITS)
#define BLOCK_SIZE 1000000

typedef struct _ReadData ReadData;
struct _ReadData {
        ReadData *next;
        unsigned long long *words[READ_BLOCKS];
        unsigned int word_counts[READ_BLOCKS];
        unsigned int words_total;
        unsigned int wsize;
};

static ReadData *
read_data_new (unsigned int wsize)
{
        ReadData *rd;
        unsigned int i;
        rd = (ReadData *) malloc (sizeof (ReadData));
        memset (rd, 0, sizeof (ReadData));
        rd->words[0] = (unsigned long long *) malloc (READ_BLOCKS * BLOCK_SIZE * 8);
        for (i = 0; i < READ_BLOCKS; i++) {
                rd->words[i] = rd->words[0] + i * BLOCK_SIZE;
        }
        rd->wsize = wsize;
        return rd;
}

static void
read_data_delete (ReadData *rd) {
        free (rd->words[0]);
        free (rd);
}

typedef struct _GLMQueue GLMQueue;
struct _GLMQueue {
        Queue queue;
        /* Files to process */
        unsigned int nfiles;
        TaskFile *files;
        /* Per prefix tables */
        ReadData *free_reads;
        ReadData *full_reads;
        /* Data */
        unsigned int wordlength;
        unsigned int shift;
        Trie *trie;
        unsigned char locks[256];
};

#define TOKEN_READ_TIME 0
#define TOKEN_READ_COUNT 1
#define TOKEN_SORT_TIME 2
#define TOKEN_INSERT_TIME 3

static int process_word (FastaReader *reader, unsigned long long word, void *data);
static void describe_trie (Trie *trie, TrieRef node, unsigned int level, unsigned long long n[]);
static unsigned long long count_trie (Trie *trie, TrieRef ref);

/* Print usage and help menu */
static void print_help (int exitvalue);

int debug = 0;
int ntables = 0;
const char *outputname = "out";

enum {
        N_NODES,
        N_KMER_NODES,
        N_BRANCH_NODES,
        N_COUNT_NODES,
        N_COUNT_SLOTS,
        N_USED_SLOTS,
        N_WORDS,
        N_UNIQUE_WORDS,
        NUM_DESCRIPTORS
};

static void
print_mem (unsigned long long amount, FILE *ofs)
{
        if (amount < 1000) {
                fprintf (ofs, "%lluB", amount);
        } else if (amount < 1000000) {
                fprintf (ofs, "%.1fKB", amount / 1000.0);
        } else if (amount < 1000000000) {
                fprintf (ofs, "%.1fMB", amount / 1000000.0);
        } else {
                fprintf (ofs, "%.1fGB", amount / 1000000000.0);
        }
}

static void
print_time (const char *str, double time, unsigned long long count)
{
        fprintf (stderr, "%s: %llu words in %.1f seconds (%.0f w/s)\n", str, count, time, count / time);
}

/* Main thread loop */
static void *process (void *arg);

unsigned long long nw = 0;

int 
main (int argc, const char *argv[])
{
	unsigned int argidx, i;
        unsigned int nseqs = 0;
        const char *seqnames[1024];
	char *end;
	GLMQueue glmq;
	unsigned int nthreads = DEFAULT_NUM_THREADS;
	pthread_t threads[256];

	/* default values */
	unsigned int wordlength = 16;
	unsigned int cutoff = 1;

	fprintf (stderr, "Branch size is %u\n", (unsigned int) sizeof (TrieNodeBranch));

	/* parsing commandline arguments */
	for (argidx = 1; argidx < argc; argidx++) {

		if (!strcmp (argv[argidx], "-v") || !strcmp (argv[argidx], "--version")) {
			fprintf (stdout, "glistmaker v%d.%d\n", VERSION_MAJOR, VERSION_MINOR);
			return 0;
		} else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
			print_help (0);
		} else if (!strcmp (argv[argidx], "-o") || !strcmp (argv[argidx], "--outputname")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No output name specified!\n");
				argidx += 1;
				continue;
			}
			outputname = argv[argidx + 1];
			argidx += 1;
		} else if (!strcmp (argv[argidx], "-w") || !strcmp (argv[argidx], "--wordlength")) {
			if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No word-length specified! Using the default value: %d.\n", wordlength);
				argidx += 1;
				continue;
			}
			wordlength = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid word-length: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp (argv[argidx], "-c") || !strcmp (argv[argidx], "--cutoff")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No frequency cut-off specified! Using the default value: %d.\n", cutoff);
				argidx += 1;
				continue;
			}
			cutoff = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid frequency cut-off: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
                } else if (!strcmp (argv[argidx], "-t") || !strcmp (argv[argidx], "--num_threads")) {
                        argidx += 1;
                        if (argidx >= argc) {
                                print_help (1);
                        }
                        nthreads = strtol (argv[argidx], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid number of threads: %s! Must be an integer.\n", argv[argidx]);
				print_help (1);
			}
		} else if (!strcmp (argv[argidx], "-D")) {
			debug += 1;
		} else {
		        if (nseqs >= 1024) {
		                fprintf (stderr, "Maximum number of input sequence files is 1024\n");
		                exit (1);
                        }
                        seqnames[nseqs++] = argv[argidx];
		}
	}

	/* checking parameter values */

	if (!nseqs) {
		fprintf (stderr, "Error: No FastA/FastQ file specified!\n");
		print_help (1);
	}
	if (wordlength < 1 || wordlength > 32) {
		fprintf (stderr, "Error: Invalid word-length: %d!\n", wordlength);
		print_help (1);
	}
	if (cutoff < 1) {
		fprintf (stderr, "Error: Invalid frequency cut-off: %d! Must be positive.\n", cutoff);
		print_help (1);
	}
	if (strlen (outputname) > 200) {
		fprintf (stderr, "Error: Output name exceeds the 200 character limit.");
		return 1;
	}

	wsize = 2 * wordlength;
	
	/* Read files */
	/* Initialize queue */
	memset (&glmq, 0, sizeof (glmq));
	glmq.wordlength = wordlength;
	while ((glmq.shift + 8) < (wordlength * 2)) {
	        glmq.shift += 8;
        }

	if (debug > 0) fprintf (stderr, "Initializing trie...");
	i = 1;
	while (i < nthreads) i = i << 1;
	glmq.trie = trie_new (wsize, 28, i);
	if (debug > 0) fprintf (stderr, " done\n");

	for (i = 0; i < nseqs; i++) {
	        TaskFile *tf = task_file_new (seqnames[i]);
	        tf->next = glmq.files;
	        glmq.files = tf;
	        glmq.nfiles += 1;
        }
        for (i = 0; i < 8; i++) {
                ReadData *rd;
                rd = read_data_new (wsize);
                rd->next = glmq.free_reads;
                glmq.free_reads = rd;
        }
        /* Initialize main mutex and cond */
        queue_init (&glmq.queue);
        /* Lock the mutex */
        queue_lock (&glmq.queue);
        glmq.queue.nthreads = 0;
        for (i = 1; i < nthreads; i++){
                int rc;
                if (debug > 1) fprintf (stderr, "Creating thread %u\n", i);
                rc = pthread_create (&threads[i], NULL, process, &glmq);
                if (rc) {
                        fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
                        exit (-1);
                }
        }
        queue_unlock (&glmq.queue);
        process (&glmq);
        queue_lock (&glmq.queue);
        while (glmq.queue.nthreads > 0) {
                queue_wait (&glmq.queue);
        }
        queue_unlock (&glmq.queue);
        queue_finalize (&glmq.queue);
        if (debug > 0) fprintf (stderr, "Finished counting\n");

#if 0
	for (argidx = firstfasta; argidx <= firstfasta + nfasta - 1; argidx++) {
		FastaReader reader;
		const char *cdata;
	        size_t csize;
#else
        if (1) {
#endif
		unsigned int i, v;
		unsigned long long n[64][NUM_DESCRIPTORS], n_words, n_unique_words, n_branches;

#if 0
		cdata = mmap_by_filename (argv[argidx], &csize);
		if (!cdata) {
			fprintf (stderr, "Error: Cannot read file %s!\n", argv[argidx]);
			return 1;
		}
		fasta_reader_init_from_data (&reader, wordlength, 1, (const unsigned char *) cdata, csize);

		/* reading words from FastA/FastQ */
		v = fasta_reader_read_nwords (&reader, 0xffffffffffffffffULL, NULL, NULL, NULL, NULL, process_word, (void *) &trie);
		if (v) return print_error_message (v);
		fasta_reader_release (&reader);

		/* write the final list into a file */
		if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);

		fprintf (stderr, "File %s\n", argv[argidx]);
		fprintf (stderr, "Kmer matches %llu\n", trie.n_kmer_matches);
		fprintf (stderr, "Kmer mismatches %llu\n", trie.n_kmer_mismatches);
		fprintf (stderr, "Allocated %llu nodes (%.2f MB) (%llu kmers %llu branches)\n", trie.n_nodes_allocated, trie.n_bytes_allocated / 1000000.0, trie.n_kmers_allocated, trie.n_nodes_allocated - trie.n_kmers_allocated);
		fprintf (stderr, "Deallocated %llu nodes (%llu kmers %llu branches)\n", trie.n_nodes_deallocated, trie.n_kmers_deallocated, trie.n_nodes_deallocated - trie.n_kmers_deallocated);
		fprintf (stderr, "Remaining %llu nodes (%llu kmers %llu branches)\n", trie.n_nodes_allocated - trie.n_nodes_deallocated, trie.n_kmers_allocated - trie.n_kmers_deallocated, trie.n_nodes_allocated - trie.n_kmers_allocated - trie.n_nodes_deallocated + trie.n_kmers_deallocated);
#endif
		
		memset (n, 0, 32 * NUM_DESCRIPTORS * 8);
		for (i = 0; i < (1 << glmq.trie->nbits_root); i++) {
		        describe_trie (glmq.trie, glmq.trie->roots[i], 0, n[0]);
		}

		n_words = 0;
		n_unique_words = 0;
		n_branches = 0;
		for (i = 0; i < 64; i++) {
		        if (n[i][N_NODES]) {
		                n_words +=  n[i][N_WORDS];
		                n_unique_words += n[i][N_UNIQUE_WORDS];
		                n_branches += n[i][N_BRANCH_NODES];
		        }
		}
                fprintf (stderr, "Words %llu unique %llu\n", n_words, n_unique_words);
                fprintf (stderr, "Branches %llu\n", n_branches);

                print_time ("Read", glmq.queue.tokens[TOKEN_READ_TIME].dval, glmq.queue.tokens[TOKEN_READ_COUNT].lval);
                print_time ("Sort", glmq.queue.tokens[TOKEN_SORT_TIME].dval, glmq.queue.tokens[TOKEN_READ_COUNT].lval);
                print_time ("Insert", glmq.queue.tokens[TOKEN_INSERT_TIME].dval, glmq.queue.tokens[TOKEN_READ_COUNT].lval);
	}

	/*wordtable_delete (temptable);
	wordtable_delete (table);*/

        pthread_exit (NULL);
}

static void *
process (void *arg)
{
        GLMQueue *glmq;
        Queue *queue;
        int idx = -1;
        unsigned int finished;
        double t_s, t_e;
        Trie *trie;
        unsigned int shift;

        glmq = (GLMQueue *) arg;
        queue = &glmq->queue;
        trie = glmq->trie;
        shift = glmq->shift;

        finished = 0;
        
        /* Do work */
        while (!finished) {
                /* Get exclusive lock on queue */
                queue_lock (queue);
                if (idx < 0) {
                        /* Thread is started, increase counter and get idx */
                        idx = queue->nthreads;
                        queue->nthreads += 1;
                        if (debug > 1) fprintf (stderr, "Thread %d started (total %d)\n", idx, queue->nthreads);
                }
                if (glmq->files && glmq->free_reads) {
                        TaskFile *tf;
                        ReadData *rd;
                        unsigned int i;
                        tf = glmq->files;
                        glmq->files = tf->next;
                        rd = glmq->free_reads;
                        glmq->free_reads = rd->next;

                        queue_unlock (queue);
                        /* Read words from file */
                        if (debug > 0) fprintf (stderr, "Thread %d: reading file %s from %llu\n", idx, tf->filename, tf->reader.cpos);
                        t_s = get_time ();
                        if (task_file_read_nwords (tf, MAX_READ_SIZE, glmq->wordlength, NULL, NULL, NULL, NULL, process_word, rd)) {
                                fprintf (stderr, "Cannot create FastaReader fro %s\n", tf->filename);
                                exit (1);
                        }
                        if (debug > 1) fprintf (stderr, "Thread %d: finished reading %s at %llu\n", idx, tf->filename, tf->reader.cpos);
                        t_e = get_time ();
                        queue_lock (queue);

                        queue->tokens[TOKEN_READ_TIME].dval += (t_e - t_s);
                        for (i = 0; i < READ_BLOCKS; i++) {
                                queue->tokens[TOKEN_READ_COUNT].lval += rd->word_counts[i];
                        }
                        if (tf->reader.in_eof) {
                                glmq->nfiles -= 1;
                                if (debug > 0) fprintf (stderr, "Thread %d: Finished reading %s (%u readers remaining)\n", idx, tf->filename, glmq->nfiles);
                                task_file_delete (tf);
                        } else {
                                tf->next = glmq->files;
                                glmq->files = tf;
                        }
                        /* Create table tasks */
                        rd->next = glmq->full_reads;
                        glmq->full_reads = rd;
                        /* Trailer */
                        queue_broadcast (queue);
                        queue_unlock (queue);
                } else if (glmq->full_reads) {
                        ReadData *rd;
                        unsigned int i, nwords;
                        unsigned long long *words;
                        unsigned int rem;
                        rd = glmq->full_reads;
                        glmq->full_reads = rd->next;
                        /* Search for lock */
                        for (i = 0; i < READ_BLOCKS; i++) {
                                if (!glmq->locks[i] && rd->word_counts[i]) {
                                        glmq->locks[i] = 1;
                                        nwords =  rd->word_counts[i];
                                        words = rd->words[i];
                                        rd->word_counts[i] = 0;
                                        rd->words_total -= nwords;
                                        rem = rd->words_total;
                                        break;
                                }
                        }
                        if (i < READ_BLOCKS) {
                                unsigned int is_last, j;
                                double t_sort, t_insert;
                                /* Found something */
                                is_last = (rd->words_total == 0);

                                if (!is_last) {
                                        /* Add it back (but keep lock) */
                                        if (debug > 1) fprintf (stderr, "Thread %d: Adding block back\n", idx);
                                        rd->next = glmq->full_reads;
                                        glmq->full_reads = rd;
                                }
                                if (debug > 1) fprintf (stderr, "Thread %d: inserting block %u (size %u) to trie\n", idx, i, nwords);

                                queue_unlock (queue);
                                t_s = get_time ();
                                hybridInPlaceRadixSort256 (words, words + nwords, NULL, shift);
                                t_e = get_time ();
                                t_sort = t_e - t_s;
                                t_s = t_e;
                                for (j = 0; j < nwords; j++) {
                                        unsigned long long word = words[j];
                                        trie_add_word_with_allocator (trie, word, 1, idx);
                                }
                                t_e = get_time ();
                                t_insert = t_e - t_s;
                                queue_lock (queue);

                                glmq->locks[i] = 0;
                                queue->tokens[TOKEN_SORT_TIME].dval += t_sort;
                                queue->tokens[TOKEN_INSERT_TIME].dval += t_insert;
                                if (debug > 1) fprintf (stderr, "Thread %d: finished block %u (%u remaining)\n", idx, i, rem);
                                if (is_last) {
                                        /* Insert it into free list if last */
                                        if (debug > 1) fprintf (stderr, "Thread %d: Marking block as free\n", idx);
                                        rd->next = glmq->free_reads;
                                        glmq->free_reads = rd;
                                }
                                /* Trailer */
                                queue_broadcast (queue);
                                queue_unlock (queue);
                        } else {
                                /* Add unused block back */
                                if (debug > 1) fprintf (stderr, "Thread %d: Adding unused block back\n", idx);
                                rd->next = glmq->full_reads;
                                glmq->full_reads = rd;
                                /* All tables are in processing, wait */
                                if (debug > 2) fprintf (stderr, "Thread %d: Waiting\n", idx);
                                queue_wait (queue);
                                if (debug > 2) fprintf (stderr, "Thread %d: Woke up\n", idx);
                                queue_unlock (queue);
                        }
                } else if (glmq->nfiles) {
                        /* All tables are in processing, wait */
                        if (debug > 2) fprintf (stderr, "Thread %d: Waiting\n", idx);
                        queue_wait (queue);
                        if (debug > 2) fprintf (stderr, "Thread %d: Woke up\n", idx);
                        queue_unlock (queue);
                } else {
                        /* Nothing to do */
                        if (debug > 2) fprintf (stderr, "Thread %d: Exiting\n", idx);
                        queue_broadcast (queue);
                        queue_unlock (queue);
                        finished = 1;
                }
        }
        /* Exit if egcc glistmaker2.c wordtable.c wordmap.c fasta.c buffer.c sequence.c queue.c common.c trie.c trie.h utils.c -o glistmaker2 -lm -lpthread -O0verything is done */
        queue_lock (queue);
        glmq->queue.nthreads -= 1;
        if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, glmq->queue.nthreads);
        queue_broadcast (queue);
        queue_unlock (queue);

        return 0;
}

static int 
process_word (FastaReader *reader, unsigned long long word, void *data)
{
        ReadData *rd;
        unsigned long long block;
        rd = (ReadData *) data;
        block = word >> (wsize - READ_BITS);
        rd->words[block][rd->word_counts[block]] = word;
        rd->word_counts[block] += 1;
        rd->words_total += 1;
        return rd->word_counts[block] >= BLOCK_SIZE;
}

static void
describe_trie (Trie *trie, TrieRef ref, unsigned int level, unsigned long long *n)
{
        if (REF_IS_EMPTY(ref)) return;
        n[N_NODES] += 1;
        if (REF_IS_KMER (ref)) {
                if (KMER_GET_COUNT (ref) > 0) {
                        n[N_KMER_NODES] += 1;
                        n[N_WORDS] += KMER_GET_COUNT (ref);
                        n[N_UNIQUE_WORDS] += 1;
                        nw += KMER_GET_COUNT (ref);
                } else {
                        fprintf (stderr, ".");
                }
        } else if (REF_IS_BRANCH (ref)) {
                TrieNodeBranch *node = BRANCH_FROM_REF(trie, ref);
                TrieRef *children = BRANCH_GET_CHILDREN (trie, node);
                unsigned int i;
                n[N_BRANCH_NODES] += 1;
                for (i = 0; i < BRANCH_GET_NUM_CHILDREN (trie, node, level); i++) {
                        describe_trie (trie, children[i], level + 1, n + NUM_DESCRIPTORS);
                }
        }
}

static unsigned long long
count_trie (Trie *trie, TrieRef ref)
{
        unsigned long long n = 0;
        if (REF_IS_EMPTY(ref)) return n;
        if (REF_IS_KMER (ref)) {
                n = KMER_GET_COUNT (ref);
        } else {
                TrieNodeBranch *branch;
                unsigned int i;
                branch = BRANCH_FROM_REF(trie, ref);
                for (i = 0; i < (1 << branch->nbits_children); i++) {
                        n += count_trie (trie, branch->children[i]);
                }
        }
        return n;
}

static void 
print_help (int exitvalue)
{
	fprintf (stderr, "Usage: glistmaker <INPUTFILES> [OPTIONS]\n");
	fprintf (stderr, "Options:\n");
	fprintf (stderr, "    -v, --version           - print version information and exit\n");
	fprintf (stderr, "    -h, --help              - print this usage screen and exit\n");
	fprintf (stderr, "    -w, --wordlength NUMBER - specify index wordsize (1-32) (default 16)\n");
	fprintf (stderr, "    -c, --cutoff NUMBER     - specify frequency cut-off (default 1)\n");
	fprintf (stderr, "    -o, --outputname STRING - specify output name (default \"out\")\n");
	fprintf (stderr, "    --num_threads           - number of threads the program is run on (default MIN(8, num_input_files))\n");
	fprintf (stderr, "    --max_tables            - maximum number of temporary tables (default MAX(num_threads, 2))\n");
	fprintf (stderr, "    --table_size            - maximum size of the temporary table (default 500000000)\n");
	fprintf (stderr, "    -D                      - increase debug level\n");
	exit (exitvalue);
}

