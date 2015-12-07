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
#include <sys/stat.h>
#include <pthread.h>
#include <unistd.h>

#include "utils.h"
#include "fasta.h"
#include "wordtable.h"
#include "sequence.h"
#include "queue.h"
#include "common.h"

#define MAX_FILES 200

#define DEFAULT_NUM_THREADS 8
#define DEFAULT_TABLE_SIZE 500000000
#define DEFAULT_MAX_TABLES 32

#define BSIZE 10000000

#define MAX_MERGED_TABLES 256

#define TIME_READ 0
#define TIME_SORT 1
#define TIME_MERGE 2
#define TIME_FF 3

/* Main thread loop */
static void * process (void *arg);

/* Merge tables directly to disk */
static void merge_write_multi (wordtable **t, unsigned int ntables, const char *filename, unsigned int cutoff);

/* */
int process_word (FastaReader *reader, unsigned long long word, void *data);

/* Print usage and help menu */
void print_help (int exitvalue);

int debug = 0;
int ntables = 0;
const char *outputname = "out";

int 
main (int argc, const char *argv[])
{
	int argidx, firstfasta = -1, nfasta = 1;
	char *end;

	/* default values */
	unsigned int wordlength = 16;
	unsigned int cutoff = 1;
	unsigned int nthreads = 0;
	unsigned long long tablesize = 0;

	/* parsing commandline arguments */
	for (argidx = 1; argidx < argc; argidx++) {

		if (!strcmp (argv[argidx], "-v") || !strcmp (argv[argidx], "--version")) {
			fprintf (stdout, "glistmaker v%d.%d\n", VERSION_MAJOR, VERSION_MINOR);
			return 0;

		} else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
			print_help (0);

		} else if (argidx == 1) {
			/* get the locations of the input files */
			if ((argv[argidx][0] == '-') && (argv[argidx][1] != 0)) {
				fprintf (stderr, "Error: No FastA/FastQ file specified!\n");
				print_help (1);
			}
			firstfasta = argidx;
			while (argv[argidx + 1] && argv[argidx + 1][0] != '-') {
				argidx += 1;
				nfasta += 1;
			}
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
		} else if (!strcmp (argv[argidx], "--num_threads")) {
			if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No num-threads specified! Using the default value: %d.\n", DEFAULT_NUM_THREADS);
				argidx += 1;
				continue;
			}
			nthreads = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid num-threads: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp (argv[argidx], "--max_tables")) {
			if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No max_tables specified! Using the default value: %d.\n", DEFAULT_MAX_TABLES);
				argidx += 1;
				continue;
			}
			ntables = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid max_tables: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp (argv[argidx], "--table_size")) {
			if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
				fprintf (stderr, "Warning: No table-size specified! Using the default value: %llu.\n", (unsigned long long) DEFAULT_TABLE_SIZE);
				argidx += 1;
				continue;
			}
			tablesize = strtoll (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid table-size: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp (argv[argidx], "-D")) {
			debug += 1;
		} else {
			fprintf (stderr, "Error: Unknown argument: %s\n", argv[argidx]);
			print_help (1);
		}
	}

	debug_tables = debug;
	
	if (!nthreads) {
		nthreads = DEFAULT_NUM_THREADS;
		if (nthreads > ((3 * nfasta + 1) >> 1)) nthreads = (3 * nfasta + 1) >> 1;
	}
	if (!tablesize) {
		tablesize = DEFAULT_TABLE_SIZE;
	}
	if (!ntables) {
		ntables = DEFAULT_MAX_TABLES;
		if (ntables > ((3 * nfasta + 1) >> 1) + 1) ntables = ((3 * nfasta + 1) >> 1) + 1;
	}
	if (ntables > MAX_TABLES) ntables = MAX_TABLES;

	/* checking parameter values */
	if (firstfasta == -1) {
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
	if (nthreads < 1) nthreads = 1;
	if (nthreads > 256) nthreads = 256;
	for (argidx = firstfasta; argidx < firstfasta + nfasta; argidx += 1) {
		struct stat s;
		if (stat (argv[argidx], &s)) {
			fprintf (stderr, "Error: Cannot stat %s\n", argv[argidx]);
			exit (1);
		}
	}
	
	if (nthreads > 0) {
		/* CASE: SEVERAL THREADS */
		MakerQueue mq;
	        pthread_t threads[256];
	        int rc;
	        unsigned int t;
	        unsigned int finished = 0;

		maker_queue_setup (&mq);

		for (argidx = firstfasta + nfasta - 1; argidx >= firstfasta; argidx--) {
			maker_queue_add_file (&mq, argv[argidx]);
		}		

	        /* Lock the mutex */
	        pthread_mutex_lock (&mq.queue.mutex);
	        mq.queue.nthreads = 0;
	        mq.wordlen = wordlength;
	        mq.tablesize = tablesize;
	        mq.cutoff = cutoff;

		if (debug) {
			fprintf (stderr, "Num threads is %d\n", nthreads);
			fprintf (stderr, "Num tables is %d\n", ntables);
			fprintf (stderr, "Table size is %lld\n", mq.tablesize);
		}
			  
	        for (t = 1; t < nthreads; t++){
	                if (debug > 1) fprintf (stderr, "Creating thread %u\n", t);
	                rc = pthread_create (&threads[t], NULL, process, &mq);
	                if (rc) {
                                fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
                                exit (-1);
                        }
                }

                pthread_mutex_unlock (&mq.queue.mutex);

                process (&mq);

                while (!finished) {
                        pthread_mutex_lock (&mq.queue.mutex);
                        if (mq.queue.nthreads < 1) finished = 1;
                        pthread_mutex_unlock (&mq.queue.mutex);
                        sleep (1);
                }
                if (mq.nsorted > 0) {
                	/* write the final list into a file */
                	if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);
                	wordtable_write_to_file (mq.sorted[0], outputname, cutoff);
		}

                if (debug) {
                	fprintf (stderr, "Read %.2f\n", mq.queue.tokens[TIME_READ].dval);
                	fprintf (stderr, "Sort %.2f\n", mq.queue.tokens[TIME_SORT].dval);
                	fprintf (stderr, "Collate %.2f\n", mq.queue.tokens[TIME_FF].dval);
                	fprintf (stderr, "Merge %.2f\n", mq.queue.tokens[TIME_MERGE].dval);
                }

                maker_queue_release (&mq);
        } else {
		/* CASE: ONE THREAD */
		const char *ff;
		size_t fsize;
		wordtable *table, *temptable;
		int v;
		
		/* creating initial tables */
		table = wordtable_new (wordlength, 20000);
		temptable = wordtable_new (wordlength, 20000);
		
		for (argidx = firstfasta; argidx <= firstfasta + nfasta - 1; argidx++) {
			FastaReader reader;

			temptable->wordlength = wordlength;

			if (!strcmp (argv[argidx], "-")) {
				/* stdin */
				fasta_reader_init_from_file (&reader, wordlength, 1, stdin);
			} else {
				ff = mmap_by_filename (argv[argidx], &fsize);
				if (!ff) {
					fprintf (stderr, "Error: Cannot read file %s!\n", argv[argidx]);
					return 1;
				}
				fasta_reader_init_from_data (&reader, wordlength, 1, (const unsigned char *) ff, fsize);
			}

			/* reading words from FastA/FastQ */
			v = fasta_reader_read_nwords (&reader, 0xffffffffffffffffULL, NULL, NULL, NULL, NULL, process_word, (void *) temptable);
			if (v) return print_error_message (v);
			fasta_reader_release (&reader);

			/* radix sorting */
			wordtable_sort (temptable, 0);
			v = wordtable_find_frequencies (temptable);
			if (v) return print_error_message (v);


			/* merging two tables */
			if (argidx > firstfasta) {
				v = wordtable_merge (table, temptable);
				if (v) return print_error_message (v);
			} else {
				wordtable *t = table;
				table = temptable;
				temptable = t;
			}

			/* empty the temporary table */
			wordtable_empty (temptable);

		}

		/* write the final list into a file */
		if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);
		wordtable_write_to_file (table, outputname, cutoff);
	}

	/*wordtable_delete (temptable);
	wordtable_delete (table);*/

        pthread_exit (NULL);
}

static void *
process (void *arg)
{
        MakerQueue *mq;
        int idx = -1;
        unsigned int finished;
        double s_t, e_t, d_t;

        mq = (MakerQueue *) arg;

        finished = 0;
        
        /* Do work */
        while (!finished) {
        	unsigned int has_files, has_unsorted, has_unmerged, sorted_tables;
        	
                /* Get exclusive lock on queue */
                pthread_mutex_lock (&mq->queue.mutex);
                if (idx < 0) {
                        /* Thread is started, increase counter and get idx */
                        mq->queue.nthreads += 1;
                        idx = mq->queue.nthreads;
                        if (debug > 1) fprintf (stderr, "Thread %d started (total %d)\n", idx, mq->queue.nthreads);
                }

                if (debug > 1) fprintf (stderr, "Thread %d: FileTasks %u Unsorted %u Sorted %u\n", idx, mq->ntasks[TASK_READ], mq->nunsorted, mq->nsorted);
                
                has_files = mq->files || mq->ntasks[TASK_READ];
                has_unsorted = mq->nunsorted || mq->ntasks[TASK_SORT];
                has_unmerged = mq->ntasks[TASK_MERGE];
                sorted_tables = mq->nsorted + mq->ntasks[TASK_MERGE];

                /* If all files have been read and sorted and there is small enough number of sorted files */
                if (!has_files && !has_unsorted && sorted_tables && (sorted_tables <= MAX_MERGED_TABLES)) {
                	if (!has_unmerged) {
                		/* Merge to disk */
                		wordtable *t[MAX_MERGED_TABLES];
                		unsigned int ntables;
                		char c[1024];
                		ntables = 0;
                		while (mq->nsorted) {
                			t[ntables++] = queue_get_sorted (mq);
				}
				wordtable_build_filename (t[0], c, 1024, outputname);
				mq->ntasks[TASK_MERGE] += 1;
				if (debug) {
                			unsigned int i;
                			fprintf (stderr, "Merging %u tables: %s", ntables, t[0]->id);
                			for (i = 1; i < ntables; i++) {
                				fprintf (stderr, ",%s", t[i]->id);
					}
					fprintf (stderr, " to %s\n", c);
				}
				/* Now we can release mutex */
				pthread_mutex_unlock (&mq->queue.mutex);
				/* merge_write (table, other, c, queue->cutoff); */
				merge_write_multi (t, ntables, c, mq->cutoff);
				pthread_mutex_lock (&mq->queue.mutex);
				mq->ntasks[TASK_MERGE] -= 1;
				pthread_cond_broadcast (&mq->queue.cond);
				pthread_mutex_unlock (&mq->queue.mutex);
				finished = 1;
			} else {
				/* Waiting merging to finish */
                        	if (debug > 1) fprintf (stderr, "Thread %d: Waiting merging to finish\n", idx);
                        	pthread_cond_wait (&mq->queue.cond, &mq->queue.mutex);
                        	pthread_mutex_unlock (&mq->queue.mutex);
			}
			continue;
                }
                if (mq->nsorted > 1) {
                	/* Task 1 - merge sorted tables */
                        wordtable *table, *other;
                        int result;

                        other = queue_get_smallest_sorted (mq);
                        table = queue_get_mostavailable_sorted (mq);
                        if (table->nwordslots < other->nwordslots) {
                        	wordtable *t = table;
                        	table = other;
                        	other = t;
                        }
                        mq->ntasks[TASK_MERGE] += 1;
                        /* Now we can release mutex */
                        pthread_mutex_unlock (&mq->queue.mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Merging tables %s (%llu/%llu) + %s (%llu/%llu) -> %s\n", idx, table->id, table->nwords, table->nwordslots, other->id, other->nwords, other->nwordslots, table->id);
                        s_t = get_time ();
			result = wordtable_merge (table, other);
			e_t = get_time ();
			d_t = e_t - s_t;
                        /* fixme: Error processing */
			if (result) {
			        print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&mq->queue.mutex);
                        /* Add merged table to sorted list */
                        mq->sorted[mq->nsorted++] = table;
                        wordtable_empty (other);
                        other->wordlength = mq->wordlen;
                        mq->available[mq->navailable++] = other;
                        mq->queue.tokens[TIME_MERGE].dval += d_t;
                        /* Release mutex */
                        mq->ntasks[TASK_MERGE] -= 1;
                        pthread_cond_broadcast (&mq->queue.cond);
                        pthread_mutex_unlock (&mq->queue.mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished merging %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (mq->nunsorted > 0) {
                        /* Task 2 - sort table */
                        wordtable *table;
                        int result;
                        
                        table = mq->unsorted[--mq->nunsorted];
                        /* Now we can release mutex */
                        mq->ntasks[TASK_SORT] += 1;
                        pthread_mutex_unlock (&mq->queue.mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Sorting table %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                        s_t = get_time ();
                        wordtable_sort (table, 0);
                        e_t = get_time ();
                        d_t = e_t - s_t;
                        s_t = get_time ();
                        result = wordtable_find_frequencies (table);
                        e_t = get_time ();
                        /* fixme: Error processing */
                        if (result) {
                                print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&mq->queue.mutex);
                        /* Add sorted table to sorted list */
                        mq->sorted[mq->nsorted++] = table;
                        mq->queue.tokens[TIME_SORT].dval += d_t;
                        d_t = e_t - s_t;
                        mq->queue.tokens[TIME_FF].dval += d_t;
                        /* Release mutex */
                        mq->ntasks[TASK_SORT] -= 1;
                        pthread_cond_broadcast (&mq->queue.cond);
                        pthread_mutex_unlock (&mq->queue.mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished sorting %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (mq->files && (mq->ntasks[TASK_READ] < MAX_FILES) && (mq->navailable || (mq->ntablescreated < ntables))) {
                        /* Task 3 - read input file */
                        TaskFile *task;
                        wordtable *table;
                        int result;
                        unsigned long long readsize;

                        task = mq->files;
                        mq->files = task->next;
                        mq->ntasks[TASK_READ] += 1;
                        if (mq->navailable > 0) {
                                /* Has to create new word table */
                                table = queue_get_largest_table (mq);
                        } else {
                                table = wordtable_new (mq->wordlen, 10000000);
                                table->wordlength = mq->wordlen;
                                mq->ntablescreated += 1;
                                if (debug > 0) fprintf (stderr, "Thread %d: Created table %s\n", idx, table->id);
                        }
                        /* Now we can release mutex */
                        pthread_mutex_unlock (&mq->queue.mutex);
                        
                        /* Process file reader task */
                        if (debug > 1) fprintf (stderr, "Thread %d: Processign file %s (%llu) -> %s (%llu)\n", idx, task->filename, (unsigned long long) task->reader.cpos, table->id, table->nwordslots);
                        if (!task->reader.wordlength) {
                        	/* Initialize reader */
                                if (debug > 0) fprintf (stderr, "Thread %d: Creating FastaReader for %s -> %s\n", idx, task->filename, table->id);
                        	if (!strcmp (task->filename, "-")) {
                        		fasta_reader_init_from_file (&task->reader, mq->wordlen, 1, stdin);
                        	} else {
                                	const unsigned char *cdata;
                                	size_t csize;
                                	cdata = (const unsigned char *) mmap_by_filename (task->filename, &csize);
                                	if (cdata) {
                                		scout_mmap (cdata, csize);
                                		task->cdata = cdata;
                                		task->csize = csize;
                                		fasta_reader_init_from_data (&task->reader, mq->wordlen, 1, cdata, csize);
					} else {
                                        	/* Cannot mmap file */
                                        	/* fixme: Error procesing */
					}
                                }
                        }
                        readsize = (mq->tablesize < table->nwordslots) ? table->nwordslots : mq->tablesize;
                        if (debug > 0) fprintf (stderr, "Thread %d: Reading %lld bytes from %s, position %llu/%llu\n", idx, readsize, task->filename, (unsigned long long) task->reader.cpos, (unsigned long long) task->csize);
                        s_t = get_time ();
                        result = fasta_reader_read_nwords (&task->reader, readsize,  NULL, NULL, NULL, NULL, process_word, table);
                        e_t = get_time ();
			d_t = e_t - s_t;
                        if (result) {
                                /* fixme: Error processing */
		                print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&mq->queue.mutex);
                        /* Add generated table to unsorted list */
                        mq->unsorted[mq->nunsorted++] = table;
                        if (task->reader.in_eof) {
                                /* Finished this task */
                                if (debug > 0) fprintf (stderr, "Thread %d: FastaReader for %s finished\n", idx, task->filename);
                                if (task->cdata) {
                                	munmap ((void *) task->cdata, task->csize);
				}
                                free (task);
                        } else {
                                /* Reshedule task */
                                task->next = mq->files;
                                mq->files = task;
                        }
                        mq->ntasks[TASK_READ] -= 1;
                        mq->queue.tokens[TIME_READ].dval += d_t;
                        /* Release mutex */
                        pthread_cond_broadcast (&mq->queue.cond);
                        pthread_mutex_unlock (&mq->queue.mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished reading %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (!has_files && !mq->nunsorted && (mq->nsorted < 2)) {
                        /* Nothing to do */
                        /* Release mutex */
                        pthread_cond_broadcast (&mq->queue.cond);
                        pthread_mutex_unlock (&mq->queue.mutex);
                        finished = 1;
                } else {
                        if (debug > 1) fprintf (stderr, "Thread %d: Waiting\n", idx);
                        /* Release mutex */
                        /* pthread_mutex_unlock (&mq->queue.mutex); */
                        /* fixme: Semaphore */
                        /* sleep (1); */
                        pthread_cond_wait (&mq->queue.cond, &mq->queue.mutex);
                        pthread_mutex_unlock (&mq->queue.mutex);
                }
        }

        /* Exit if everything is done */
        pthread_mutex_lock (&mq->queue.mutex);
        mq->queue.nthreads -= 1;
        if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, mq->queue.nthreads);
        pthread_cond_broadcast (&mq->queue.cond);
        pthread_mutex_unlock (&mq->queue.mutex);

        /* pthread_exit (NULL); */
        return 0;
}

int 
process_word (FastaReader *reader, unsigned long long word, void *data)
{
	wordtable *table = (wordtable *) data;
#if 1
	wordtable_add_word_nofreq (table, word, reader->wordlength);
#endif
	return 0;
}

static void
merge_write_multi (wordtable *t[], unsigned int ntables, const char *filename, unsigned int cutoff)
{
	unsigned long long nwords[MAX_MERGED_TABLES];
	unsigned long long i[MAX_MERGED_TABLES];
	unsigned int nfinished;

	header h;

	FILE *ofs;
	char *b;
	unsigned int bp, j;

	unsigned long long word;
	unsigned int freq;
	double t_s, t_e;

	h.code = glistmaker_code_match;
	h.version_major = VERSION_MAJOR;
	h.version_minor = VERSION_MINOR;
	h.wordlength = t[0]->wordlength;
	h.nwords = 0;
	h.totalfreq = 0;
	h.padding = sizeof (header);

	b = (char *) malloc (BSIZE + 12);
	
	ofs = fopen (filename, "w");

	t_s = get_time ();
	fwrite (&h, sizeof (header), 1, ofs);

	bp = 0;
	nfinished = 0;
	for (j = 0; j < ntables; j++) {
		i[j] = 0;
		nwords[j] = t[j]->nwords;
		if (!nwords[j]) nfinished += 1;
	}
		
	memset (i, 0, sizeof (i));
	
	while (nfinished < ntables) {
		word = 0xffffffffffffffff;
		freq = 0;
		/* Find smalles word and total freq */
		for (j = 0; j < ntables; j++) {
			if (i[j] < nwords[j]) {
				/* This table is not finished */
				if (t[j]->words[i[j]] < word) {
					/* This table has smaller word */
					word = t[j]->words[i[j]];
					freq = t[j]->frequencies[i[j]];
				} else if (t[j]->words[i[j]] == word) {
					/* This table has equal word */
					freq += t[j]->frequencies[i[j]];
				}
				__builtin_prefetch (&t[j]->words[i[j]] + 16);
				__builtin_prefetch (&t[j]->frequencies[i[j]] + 16);
			}
		}
		/* Now we have word and freq */
		if (freq >= cutoff) {
			memcpy (b + bp, &word, 8);
			bp += 8;
			memcpy (b + bp, &freq, 4);
			bp += 4;
			if (bp >= BSIZE) {
				fwrite (b, 1, bp, ofs);
				bp = 0;
			}
			h.nwords += 1;
			h.totalfreq += freq;
		}
		/* Update pointers */
		for (j = 0; j < ntables; j++) {
			if (i[j] < nwords[j]) {
				/* This table is not finished */
				if (t[j]->words[i[j]] == word) {
					i[j] += 1;
					if (i[j] >= nwords[j]) {
						nfinished += 1;
					}
				}
			}
		}
	}
	if (bp) {
		fwrite (b, 1, bp, ofs);
	}
	fseek (ofs, 0, SEEK_SET);
	fwrite (&h, sizeof (header), 1, ofs);
	fclose (ofs);
	t_e = get_time ();
	if (debug > 0) fprintf (stderr, "Writing %d tables with merging %.2f\n", ntables, t_e - t_s);
	
	free (b);
}

void 
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

