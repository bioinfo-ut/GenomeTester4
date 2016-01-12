#ifndef __AOSORA_THREAD_POOL_H__
#define __AOSORA_THREAD_POOL_H__

/*
 * Abstraction for worker thread manager
 *
 * Copyright Lauris Kaplinski 2014
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _AosoraThreadPool AosoraThreadPool;

AosoraThreadPool *aosora_thread_pool_new (unsigned int nthreads);
void aosora_thread_pool_delete (AosoraThreadPool *tpool);

/*
 * The semantics
 *
 * pre_run is called from the main thread before any workers is started (the last opportunity to query the state of other objects)
 * post_run is called from the main thread sometimes after the completion of run (other threads may still run)
 * complete is called from the main thread when all workers threads have finished
 * pre_run, post_run and complete may be NULL
 *
 * It is allowed to submit new workers during pre_run, run and post_run but their pre_run will not be called in that case
 *
 */

void aosora_thread_pool_submit (AosoraThreadPool *tpool, void (*pre_run) (void *), void (*run) (void *), void (*post_run) (void *), void (*complete) (void *), void *data);

/*
 * Run all scheduled threads
 *
 * Main thread waits until all workers have finished
 *
 */
void aosora_thread_pool_run (AosoraThreadPool *tpool);

#ifdef __cplusplus
}
#endif

#endif
