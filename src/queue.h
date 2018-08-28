#ifndef __GT4_QUEUE_H__
#define __GT4_QUEUE_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-18 University of Tartu
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

typedef struct _GT4Queue GT4Queue;
typedef struct _GT4QueueClass GT4QueueClass;
typedef struct _GT4Task GT4Task;

#define GT4_TYPE_QUEUE gt4_queue_get_type()

#include <pthread.h>

#include <az/class.h>

struct _GT4Task {
  GT4Task *next;
  GT4Queue *queue;
  unsigned int type;
  unsigned int priority;
};

struct _GT4Queue {
  /* Threads */
  unsigned int nthreads_total;
  unsigned int nthreads_running;
  pthread_t *threads;
  /* Single mutex and cond for queue management */
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  GT4Task *tasks;
};

struct _GT4QueueClass {
  AZClass klass;
};

unsigned int gt4_queue_get_type (void);

/* Initialize/destroy mutexes and conds */
/* Thread 0 is main, queue will be created with nthreads - 1 entries */
void gt4_queue_setup (GT4Queue *queue, unsigned int nthreads);

/* Return 0 on success */
unsigned int gt4_queue_create_threads (GT4Queue *queue, void (*process) (GT4Queue *, unsigned int, void *), void *data);
void gt4_queue_lock (GT4Queue *queue);
void gt4_queue_unlock (GT4Queue *queue);
/* Should be called with mutex locked, returns mutex locked */
void gt4_queue_wait (GT4Queue *queue);
void gt4_queue_broadcast (GT4Queue *queue);
/* Put task into queue, ordererd by priority */
void gt4_queue_add_task (GT4Queue *queue, GT4Task *task, unsigned int lock);
void gt4_queue_remove_task (GT4Queue *queue, GT4Task *task, unsigned int lock);

#endif
