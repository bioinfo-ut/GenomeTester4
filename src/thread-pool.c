#define __AOSORA_THREAD_POOL_C__

/*
* Abstraction for worker thread manager
*
* Copyright Lauris Kaplinski 2014
*/

#define debug 0

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <errno.h>

#include "thread-pool.h"

typedef struct _Task Task;

static void *worker_main (void *data);
static Task *task_new (AosoraThreadPool *tpool);
#if 0
static void task_delete (Task *task);
#endif

struct _Task {
	AosoraThreadPool *tpool;
	Task *next;
	void (*pre_run) (void *);
	void (*run) (void *);
	void (*post_run) (void *);
	void (*complete) (void *);
	void *data;
};

struct _AosoraThreadPool {
	unsigned int nthreads;

	/* Main mutex for ThreadPool */
	pthread_cond_t cnd_work;
	pthread_cond_t cnd_finish;
	pthread_mutex_t mtx;
	pthread_t *threads;

	/* Global signal for threads */
	unsigned int req_shutdown;
	unsigned int nworkers;

	/* Task queues */
	Task *submitted;
	Task *scheduled;
	Task *finished;
	Task *completed;
	Task *free;
};

AosoraThreadPool *
aosora_thread_pool_new (unsigned int nthreads)
{
	unsigned int i;
	AosoraThreadPool *tpool = (AosoraThreadPool *) malloc (sizeof (AosoraThreadPool));
	memset (tpool, 0, sizeof (AosoraThreadPool));
	tpool->nthreads = nthreads;
	tpool->threads = (pthread_t *) malloc (nthreads * sizeof (pthread_t));
	pthread_cond_init (&tpool->cnd_work, NULL);
	pthread_cond_init (&tpool->cnd_finish, NULL);
	pthread_mutex_init (&tpool->mtx, NULL);
	for (i = 0; i < nthreads; i++) {
		pthread_create (&tpool->threads[i], NULL, worker_main, tpool);
	}

	return tpool;
}

void
aosora_thread_pool_delete (AosoraThreadPool *tpool)
{
	unsigned int i;
	pthread_mutex_lock (&tpool->mtx);
	tpool->req_shutdown = 1;
	pthread_cond_broadcast (&tpool->cnd_work);
	pthread_mutex_unlock (&tpool->mtx);
	for (i = 0; i < tpool->nthreads; i++) {
		void *res;
		pthread_join (tpool->threads[i], &res);
		if (debug) {
			fprintf (stderr, "aosora_thread_pool_delete: Thread %d exited with result %llu\n", i, (unsigned long long) res);
		}
	}
	free (tpool->threads);
	pthread_cond_destroy (&tpool->cnd_work);
	pthread_cond_destroy (&tpool->cnd_finish);
	pthread_mutex_destroy (&tpool->mtx);
	/* fixme: Clear task lists */
	free (tpool);
}

static void *
worker_main (void *data)
{
	AosoraThreadPool *tpool = (AosoraThreadPool *) data;
	pthread_mutex_lock (&tpool->mtx);
	while (!tpool->req_shutdown) {
		if (tpool->scheduled) {
			Task *task = tpool->scheduled;
			tpool->scheduled = task->next;
			tpool->nworkers += 1;
			pthread_mutex_unlock (&tpool->mtx);
			task->run (task->data);
			pthread_mutex_lock (&tpool->mtx);
			tpool->nworkers -= 1;
			task->next = tpool->finished;
			tpool->finished = task;
			/* Signal that one job is complete */
			pthread_cond_broadcast (&tpool->cnd_finish);
		} else {
			// Yield for 100 ms
			struct timespec t;
			clock_gettime (CLOCK_REALTIME, &t);
			t.tv_nsec += 100000000;
			int result = pthread_cond_timedwait (&tpool->cnd_work, &tpool->mtx, &t);
			if ((debug > 1) && result && (result != ETIMEDOUT)) {
				fprintf (stderr, "worker_main: Error in pthread_cond_timedwait: %d\n", result);
			}
		}
	}
	pthread_mutex_unlock (&tpool->mtx);
	return (void *) 666;
}

void
aosora_thread_pool_submit (AosoraThreadPool *tpool, void (*pre_run) (void *), void (*run) (void *), void (*post_run) (void *), void (*complete) (void *), void *data)
{
	Task *task = NULL;
	pthread_mutex_lock (&tpool->mtx);
	if (tpool->free) {
		task = tpool->free;
		tpool->free = task->next;
	} else {
		task = task_new (tpool);
	}
	task->pre_run = pre_run;
	task->run = run;
	task->post_run = post_run;
	task->complete = complete;
	task->data = data;
	if (tpool->scheduled || tpool->nworkers) {
		task->next = tpool->scheduled;
		tpool->scheduled = task;
	} else {
		task->next = tpool->submitted;
		tpool->submitted = task;
	}
	pthread_mutex_unlock (&tpool->mtx);
}

void
aosora_thread_pool_run (AosoraThreadPool *tpool)
{
	pthread_mutex_lock (&tpool->mtx);
	if (tpool->scheduled) {
		fprintf (stderr, "aosora_thread_pool_run: Scheduled list is not empty\n");
	} else {
		Task *task;
		pthread_mutex_unlock (&tpool->mtx);
		for (task = tpool->submitted; task; task = task->next) {
			if (task->pre_run) task->pre_run (task->data);
		}
		pthread_mutex_lock (&tpool->mtx);
	}
	tpool->scheduled = tpool->submitted;
	tpool->submitted = NULL;
	pthread_cond_broadcast (&tpool->cnd_work);
	while (tpool->scheduled || tpool->nworkers || tpool->finished) {
		while (tpool->finished) {
			Task *task = tpool->finished;
			tpool->finished = task->next;
			if (task->post_run) {
				pthread_mutex_unlock (&tpool->mtx);
				task->post_run (task->data);
				pthread_mutex_lock (&tpool->mtx);
			}
			task->next = tpool->completed;
			tpool->completed = task;
		}
		if (tpool->scheduled || tpool->nworkers) {
			// Yield for 100 ms
			struct timespec t;
			clock_gettime (CLOCK_REALTIME, &t);
			t.tv_nsec += 100000000;
			int result = pthread_cond_timedwait (&tpool->cnd_finish, &tpool->mtx, &t);
			if ((debug > 1) && result && (result != ETIMEDOUT)) {
				fprintf (stderr, "aosora_thread_pool_run: Error in pthread_cond_timedwait: %d\n", result);
			}
		}
	}
	pthread_mutex_unlock (&tpool->mtx);
	/* Now worker threads have completed */
	pthread_mutex_lock (&tpool->mtx);
	while (tpool->completed) {
		Task *task = tpool->completed;
		tpool->completed = task->next;
		if (task->complete) {
			pthread_mutex_unlock (&tpool->mtx);
			task->complete (task->data);
			pthread_mutex_lock (&tpool->mtx);
		}
		task->next = tpool->free;
		tpool->free = task;
	}
	pthread_mutex_unlock (&tpool->mtx);
}

static Task *
task_new (AosoraThreadPool *tpool)
{
	Task *task = (Task *) malloc (sizeof (Task));
	memset (task, 0, sizeof (Task));
	task->tpool = tpool;
	task->next = NULL;
	return task;
}

#if 0
static void
task_delete (Task *task)
{
	free (task);
}
#endif
