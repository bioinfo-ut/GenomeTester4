#ifndef QUEUE_H_
#define QUEUE_H_

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

#include <pthread.h>

#include "wordtable.h"
#include "fasta.h"

#define MAX_TABLES 256

/*
 * Possible tasks are:
 *   - read segment from FastA file and build table
 *   - sort table
 *   - merge two tables
 */

/* Task for parsing FastA file */
/* These form linked list */

typedef struct _TaskFile TaskFile;
struct _TaskFile {
        TaskFile *next;
        const char *filename;
        FastaReader reader;
};

typedef struct _Queue Queue;
struct _Queue {
        /* Number of threads running */
        int nthreads;
        /* Single mutex for queue management */
        pthread_mutex_t mutex;
        /* Parameters */
        unsigned int wordlen;
        unsigned long long tablesize;
        /* Number of worker tasks */
        unsigned int ntasks;
        /* Input files unread or partially read  */
        unsigned int nfiletasks;
        TaskFile *files;
        /* Total number of tables created */
        unsigned int ntablescreated;
        /* Unsorted tables */
        unsigned int nunsorted;
        wordtable *unsorted[MAX_TABLES];
        /* Sorted tables */
        unsigned int nsorted;
        wordtable *sorted[MAX_TABLES];
        /* Available tables */
        unsigned int navailable;
        wordtable *available[MAX_TABLES];
};

/* Add new file task to queue (not thread-safe) */
void queue_add_file (Queue *queue, const char *filename);
/* Get smallest table */
wordtable * queue_get_smallest_table (Queue *queue);
/* Get largest table */
wordtable * queue_get_largest_table (Queue *queue);

wordtable * queue_get_smallest_sorted (Queue *queue);
wordtable * queue_get_mostavailable_sorted (Queue *queue);


#endif /* SEQUENCE_H_ */
