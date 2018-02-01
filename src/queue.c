#define __GT4_QUEUE_C__

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

unsigned int queue_debug = 0;

#include <malloc.h>
#include <stdio.h>
#include <string.h>

#include "queue.h"

static void queue_class_init (GT4QueueClass *klass);
static void queue_init (GT4QueueClass *klass, GT4Queue *queue);
static void queue_finalize (GT4QueueClass *klass, GT4Queue *queue);

static unsigned int queue_type = 0;

unsigned int
gt4_queue_get_type (void)
{
  if (!queue_type) {
    az_register_type (&queue_type, AZ_TYPE_STRUCT, (const unsigned char *) "GT4Queue",
      sizeof (GT4QueueClass), sizeof (GT4Queue),
      (void (*) (AZClass *)) queue_class_init,
      (void (*) (AZImplementation *, void *)) queue_init,
      (void (*) (AZImplementation *, void *)) queue_finalize);
  }
  return queue_type;
}

static void
queue_class_init (GT4QueueClass *klass)
{
  klass->klass.flags |= AZ_CLASS_ZERO_MEMORY;
}

static void
queue_init (GT4QueueClass *klass, GT4Queue *queue)
{
  memset (queue, 0, sizeof (GT4Queue));
  pthread_mutex_init (&queue->mutex, NULL);
  pthread_cond_init (&queue->cond, NULL);
}

static void
queue_finalize (GT4QueueClass *klass, GT4Queue *queue)
{
  pthread_cond_destroy (&queue->cond);
  pthread_mutex_destroy (&queue->mutex);
  free (queue->threads);
}

void
gt4_queue_setup (GT4Queue *queue, unsigned int nthreads)
{
  queue->nthreads_total = nthreads;
  queue->nthreads_running = 1;
  queue->threads = (pthread_t *) malloc (nthreads * sizeof (pthread_t));
}

struct _QD {
  GT4Queue *queue;
  void (*process) (GT4Queue *queue, unsigned int idx, void *data);
  void *data;
};

static void *
queue_start_thread (void *data)
{
  struct _QD *qd = (struct _QD *) data;
  unsigned int idx;
  gt4_queue_lock (qd->queue);
  idx = qd->queue->nthreads_running++;
  if (queue_debug) fprintf (stderr, "Thread %d started (total %d)\n", idx, qd->queue->nthreads_running);
  gt4_queue_unlock (qd->queue);
  qd->process (qd->queue, idx, qd->data);
  gt4_queue_lock (qd->queue);
  qd->queue->nthreads_running -= 1;
  if (queue_debug) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, qd->queue->nthreads_running);
  gt4_queue_broadcast (qd->queue);
  gt4_queue_unlock (qd->queue);
  return NULL;
}

unsigned int
gt4_queue_create_threads (GT4Queue *queue, void (*process) (GT4Queue *, unsigned int, void *), void *data)
{
  unsigned i;
  /* fixme: */
  static struct _QD qd;
  qd.queue = queue;
  qd.process = process;
  qd.data = data;
  gt4_queue_lock (queue);
  for (i = 1; i < queue->nthreads_total; i++){
    int rc;
    rc = pthread_create (&queue->threads[i], NULL, queue_start_thread, &qd);
    if (rc) {
      fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
      return (rc);
    }
  }
  gt4_queue_unlock (queue);
  return 0;
}

void
gt4_queue_lock (GT4Queue *queue)
{
  pthread_mutex_lock (&queue->mutex);
}

void
gt4_queue_unlock (GT4Queue *queue)
{
  pthread_mutex_unlock (&queue->mutex);
}

void
gt4_queue_wait (GT4Queue *queue)
{
  pthread_cond_wait (&queue->cond, &queue->mutex);
}

void
gt4_queue_broadcast (GT4Queue *queue)
{
  pthread_cond_broadcast (&queue->cond);
}

void
gt4_queue_add_task (GT4Queue *queue, GT4Task *task, unsigned int lock)
{
  if (lock) gt4_queue_lock (queue);
  if (!queue->tasks || (queue->tasks->priority >= task->priority)) {
    task->next = queue->tasks;
    queue->tasks = task;
  } else {
    GT4Task *prev = queue->tasks;
    while (prev->next && (prev->next->priority < task->priority)) prev = prev->next;
    task->next = prev->next;
    prev->next = task;
  }
  if (lock) gt4_queue_unlock (queue);
}

