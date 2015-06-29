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

/* Main thread loop */
static void * process (void *arg);

/* Merge two tables directly to disk */
static void merge_write (wordtable *table, wordtable *other, const char *filename);

/* */
int process_word (FastaReader *reader, void *data);

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
			if (argv[argidx][0] == '-') {
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
	
	if (nthreads > 1) {
		/* CASE: SEVERAL THREADS */
		Queue queue = { 0 };
	        pthread_t threads[256];
	        int rc;
	        unsigned int t;
	        unsigned int finished = 0;
		
		for (argidx = firstfasta + nfasta - 1; argidx >= firstfasta; argidx--) {
			queue_add_file (&queue, argv[argidx]);
		}		

	        /* Initialize main mutex */
	        pthread_mutex_init (&queue.mutex, NULL);
	        /* Lock the mutex */
	        pthread_mutex_lock (&queue.mutex);
	        queue.nthreads = 0;
	        queue.wordlen = wordlength;
	        queue.tablesize = tablesize;

		if (debug) {
			fprintf (stderr, "Num threads is %d\n", nthreads);
			fprintf (stderr, "Num tables is %d\n", ntables);
			fprintf (stderr, "Table size is %lld\n", queue.tablesize);
		}
			  
	        for (t = 1; t < nthreads; t++){
	                if (debug > 1) fprintf (stderr, "Creating thread %u\n", t);
	                rc = pthread_create (&threads[t], NULL, process, &queue);
	                if (rc) {
                                fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
                                exit (-1);
                        }
                }

                pthread_mutex_unlock (&queue.mutex);

                process (&queue);

                while (!finished) {
                        pthread_mutex_lock (&queue.mutex);
                        if (queue.nthreads < 1) finished = 1;
                        pthread_mutex_unlock (&queue.mutex);
                        sleep (1);
                }
                if (queue.nsorted > 0) {
                	/* write the final list into a file */
                	if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);
                	wordtable_write_to_file (queue.sorted[0], outputname, cutoff);
		}

                pthread_mutex_destroy (&queue.mutex);

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

			ff = mmap_by_filename (argv[argidx], &fsize);
			if (!ff) {
				fprintf (stderr, "Error: Cannot read file %s!\n", argv[argidx]);
				return 1;
			}

			temptable->wordlength = wordlength;

			/* reading words from FastA/FastQ */
			v = fasta_reader_read (ff, fsize, wordlength, (void *) temptable, 1, 0, NULL, NULL, process_word);
			if (v) return print_error_message (v);

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
        Queue *queue;
        int idx = -1;
        unsigned int finished;

        queue = (Queue *) arg;

        finished = 0;
        
        /* Do work */
        while (!finished) {
                /* Get exclusive lock on queue */
                pthread_mutex_lock (&queue->mutex);
                if (idx < 0) {
                        /* Thread is started, increase counter and get idx */
                        queue->nthreads += 1;
                        idx = queue->nthreads;
                        if (debug > 1) fprintf (stderr, "Thread %d started (total %d)\n", idx, queue->nthreads);
                }

                if (debug > 1) fprintf (stderr, "Thread %d: FileTasks %u Unsorted %u Sorted %u\n", idx, queue->nfiletasks, queue->nunsorted, queue->nsorted);

                if (queue->nsorted > 1) {
                	/* Task 1 - merge sorted tables */
                        wordtable *table, *other;
                        int result;

                        /* table = queue->sorted[--queue->nsorted];
                        other = queue->sorted[--queue->nsorted];
                        */
                        other = queue_get_smallest_sorted (queue);
                        table = queue_get_mostavailable_sorted (queue);
                        if (table->nwordslots < other->nwordslots) {
                        	wordtable *t = table;
                        	table = other;
                        	other = t;
                        }
                        /* Now we can release mutex */
                        queue->ntasks += 1;
                        pthread_mutex_unlock (&queue->mutex);
                        if (!queue->files && !queue->nfiletasks && !queue->nunsorted && !queue->nsorted && (queue->ntasks == 1)) {
                        	/* Merge to disk */
                        	char c[1024];
                        	wordtable_build_filename (table, c, 1024, outputname);
                        	if (debug > 0) fprintf (stderr, "Thread %d: Doing final merge %s (%llu/%llu) + %s (%llu/%llu) to %s\n", idx, table->id, table->nwords, table->nwordslots, other->id, other->nwords, other->nwordslots, c);
                        	pthread_mutex_unlock (&queue->mutex);
                        	merge_write (table, other, c);
	                        pthread_mutex_lock (&queue->mutex);
	                        queue->ntasks -= 1;
                        	pthread_mutex_unlock (&queue->mutex);
                        	continue;
                        }
                        if (debug > 0) fprintf (stderr, "Thread %d: Merging tables %s (%llu/%llu) + %s (%llu/%llu) -> %s\n", idx, table->id, table->nwords, table->nwordslots, other->id, other->nwords, other->nwordslots, table->id);
			result = wordtable_merge (table, other);
                        /* fixme: Error processing */
			if (result) {
			        print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&queue->mutex);
                        /* Add merged table to sorted list */
                        queue->sorted[queue->nsorted++] = table;
                        wordtable_empty (other);
                        other->wordlength = queue->wordlen;
                        queue->available[queue->navailable++] = other;
                        /* Release mutex */
                        queue->ntasks -= 1;
                        pthread_mutex_unlock (&queue->mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished merging %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (queue->nunsorted > 0) {
                        /* Task 2 - sort table */
                        wordtable *table;
                        int result;
                        
                        table = queue->unsorted[--queue->nunsorted];
                        /* Now we can release mutex */
                        queue->ntasks += 1;
                        pthread_mutex_unlock (&queue->mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Sorting table %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                        wordtable_sort (table, 0);
                        result = wordtable_find_frequencies (table);
                        /* fixme: Error processing */
                        if (result) {
                                print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&queue->mutex);
                        /* Add sorted table to sorted list */
                        queue->sorted[queue->nsorted++] = table;
                        /* Release mutex */
                        queue->ntasks -= 1;
                        pthread_mutex_unlock (&queue->mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished sorting %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (queue->files && (queue->nfiletasks < MAX_FILES) && (queue->navailable || (queue->ntablescreated < ntables))) {
                        /* Task 3 - read input file */
                        TaskFile *task;
                        wordtable *table;
                        int result;
                        unsigned long long readsize;

                        task = queue->files;
                        queue->files = task->next;
                        queue->nfiletasks += 1;
                        if (queue->navailable > 0) {
                                /* Has to create new word table */
                                table = queue_get_largest_table (queue);
                        } else {
                                table = wordtable_new (queue->wordlen, 10000000);
                                table->wordlength = queue->wordlen;
                                queue->ntablescreated += 1;
                                if (debug > 0) fprintf (stderr, "Thread %d: Created table %s\n", idx, table->id);
                        }
                        /* Now we can release mutex */
                        queue->ntasks += 1;
                        pthread_mutex_unlock (&queue->mutex);
                        
                        /* Process file reader task */
                        if (debug > 1) fprintf (stderr, "Thread %d: Processign file %s (%llu) -> %s (%llu)\n", idx, task->filename, (unsigned long long) task->reader.cpos, table->id, table->nwordslots);
                        if (!task->reader.cdata) {
                                const char *cdata;
                                size_t csize;
                                if (debug > 0) fprintf (stderr, "Thread %d: Creating FastaReader for %s -> %s\n", idx, task->filename, table->id);
                                cdata = mmap_by_filename (task->filename, &csize);
                                if (cdata) {
                                        fasta_reader_init (&task->reader, cdata, csize, queue->wordlen, 1, 0);
                                } else {
                                        /* Cannot mmap file */
                                        /* fixme: Error procesing */
                                }
                        }
                        readsize = (queue->tablesize < table->nwordslots) ? table->nwordslots : queue->tablesize;
                        if (debug > 0) fprintf (stderr, "Thread %d: Reading %lld bytes from %s, position %llu/%llu\n", idx, readsize, task->filename, (unsigned long long) task->reader.cpos, (unsigned long long) task->reader.csize);
                        result = fasta_reader_read_nwords (&task->reader, readsize, table, 0, NULL, NULL, process_word);
                        if (result) {
                                /* fixme: Error processing */
		                print_error_message (result);
                        }
                        /* Lock mutex */
                        pthread_mutex_lock (&queue->mutex);
                        /* Add generated table to unsorted list */
                        queue->unsorted[queue->nunsorted++] = table;
                        if (task->reader.cpos >= task->reader.csize) {
                                /* Finished this task */
                                if (debug > 0) fprintf (stderr, "Thread %d: FastaReader for %s finished\n", idx, task->filename);
                                munmap ((void *) task->reader.cdata, task->reader.csize);
                                free (task);
                        } else {
                                /* Reshedule task */
                                task->next = queue->files;
                                queue->files = task;
                        }
                        queue->nfiletasks -= 1;
                        /* Release mutex */
                        queue->ntasks -= 1;
                        pthread_mutex_unlock (&queue->mutex);
                        if (debug > 0) fprintf (stderr, "Thread %d: Finished reading %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
                } else if (!queue->files && !queue->nfiletasks && !queue->nunsorted && (queue->nsorted < 2)) {
                        /* Nothing to do */
                        /* Release mutex */
                        pthread_mutex_unlock (&queue->mutex);
                        finished = 1;
                } else {
                        if (debug > 1) fprintf (stderr, "Thread %d: Waiting\n", idx);
                        /* Release mutex */
                        pthread_mutex_unlock (&queue->mutex);
                        /* fixme: Semaphore */
                        sleep (1);
                }
        }

        /* Exit if everything is done */
        pthread_mutex_lock (&queue->mutex);
        queue->nthreads -= 1;
        if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, queue->nthreads);
        pthread_mutex_unlock (&queue->mutex);

        /* pthread_exit (NULL); */
        return 0;
}

int 
process_word (FastaReader *reader, void *data)
{
	wordtable *table = (wordtable *) data;
	wordtable_add_word_nofreq (table, reader->currentword, reader->wordlength);
	return 0;
}

static void
merge_write (wordtable *t0, wordtable *t1, const char *filename)
{
	unsigned long long i0, i1;
	header h;
	FILE *ofs;
	char *b;
	unsigned long long word;
	unsigned int freq;
	double start, mid, end;
	h.code = glistmaker_code_match;
	h.version_major = VERSION_MAJOR;
	h.version_minor = VERSION_MINOR;
	h.wordlength = t0->wordlength;
	h.nwords = 0;
	h.totalfreq = 0;
	b = malloc (1024 * 1024);
	ofs = fopen (filename, "w");
	setbuffer (ofs, b, 1024 * 1024);
	start = get_time ();
	fwrite (&h, sizeof (header), 1, ofs);
	i0 = 0;
	i1 = 0;
	while ((i0 < t0->nwords) || (i1 < t1->nwords)) {
		if ((i0 < t0->nwords) && (i1 < t1->nwords) && (t0->words[i0] == t1->words[i1])) {
			word = t0->words[i0];
			freq = t0->frequencies[i0] + t1->frequencies[i1];
			i0 += 1;
			i1 += 1;
		} else if ((i1 >= t1->nwords) || ((i0 < t0->nwords) && (t0->words[i0] < t1->words[i1]))) {
			word = t0->words[i0];
			freq = t0->frequencies[i0];
			i0 += 1;
			if (i0 > t0->nwords) {
				//fprintf (stderr, "i0 %lld i1 %lld word %llu\n", i0, i1, word);
			}
		} else {
			word = t1->words[i1];
			freq = t1->frequencies[i1];
			i1 += 1;
		}
		fwrite (&word, sizeof (word), 1, ofs);
		/* memcpy (&w2, &word, sizeof (word)); */
		fwrite (&freq, sizeof (freq), 1, ofs);
		/* memcpy (&f2, &freq, sizeof (freq)); */
		h.nwords += 1;
		h.totalfreq += freq;
	}
	mid = get_time ();
	fseek (ofs, 0, SEEK_SET);
	fwrite (&h, sizeof (header), 1, ofs);
	fclose (ofs);
	end = get_time ();
	if (debug > 0) fprintf (stderr, "Writing array %.2f, writing header %.2f\n", mid - start, end - mid);
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

