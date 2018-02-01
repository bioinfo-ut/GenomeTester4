#define __GT4_SEQUENCE_FILE_C__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 *
 * Copyright (C) 2014-2017 University of Tartu
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

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "utils.h"

#include "sequence-file.h"

GT4SequenceFile *
gt4_sequence_file_new (const char *path, unsigned int lock)
{
  GT4SequenceFile *seqfile;
  struct stat st;
  if (stat (path, &st) < 0) {
    fprintf (stderr, "Cannot stat file %s\n", path);
    //return NULL;
    st.st_size = 0;
  }
  seqfile = (GT4SequenceFile *) malloc (sizeof (GT4SequenceFile));
  memset (seqfile, 0, sizeof (GT4SequenceFile));
  seqfile->refcount = 1;
  if (path) seqfile->path = strdup (path);
  seqfile->size = st.st_size;
  if (lock) {
    pthread_mutexattr_t attr;
    seqfile->lock = 1;
    pthread_mutexattr_init (&attr);
    pthread_mutexattr_settype (&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init (&seqfile->mutex, &attr);
    pthread_mutexattr_destroy (&attr);
  }
  seqfile->subseq_block_size = 1024;

  return seqfile;
}

void
gt4_sequence_file_ref (GT4SequenceFile *seqf)
{
  gt4_sequence_file_lock (seqf);
  seqf->refcount += 1;
  gt4_sequence_file_unlock (seqf);
}

void
gt4_sequence_file_unref (GT4SequenceFile *seqf)
{
  gt4_sequence_file_lock (seqf);
  if (seqf->refcount <= 1) {
    gt4_sequence_file_unlock (seqf);
    if (seqf->cdata) {
      gt4_munmap (seqf->cdata, seqf->csize);
    }
    if (seqf->subseqs) free (seqf->subseqs);
    if (seqf->lock) pthread_mutex_destroy (&seqf->mutex);
    if (seqf->path) free (seqf->path);
    free (seqf);
  } else {
    seqf->refcount -= 1;
    gt4_sequence_file_unlock (seqf);
  }
}

void
gt4_sequence_file_lock (GT4SequenceFile *seqfile)
{
  if (seqfile->lock) pthread_mutex_lock (&seqfile->mutex);
}

void
gt4_sequence_file_unlock (GT4SequenceFile *seqfile)
{
  if (seqfile->lock) pthread_mutex_unlock (&seqfile->mutex);
}

void
gt4_sequence_file_map_sequence (GT4SequenceFile *seqf)
{
  gt4_sequence_file_lock (seqf);
  seqf->cdata = gt4_mmap (seqf->path, &seqf->csize);
  seqf->pos = 0;
  gt4_sequence_file_unlock (seqf);
}

unsigned int
gt4_sequence_file_add_subsequence (GT4SequenceFile *seqfile, unsigned long long name_pos, unsigned int name_len)
{
  unsigned int index;
  gt4_sequence_file_lock (seqfile);
  if (seqfile->n_subseqs >= seqfile->size_subseqs) {
    seqfile->size_subseqs = seqfile->size_subseqs << 1;
    if (seqfile->size_subseqs < seqfile->subseq_block_size) seqfile->size_subseqs = seqfile->subseq_block_size;
    seqfile->subseqs = (GT4SubSequence *) realloc (seqfile->subseqs, seqfile->size_subseqs * sizeof (GT4SubSequence));
  }
  index = seqfile->n_subseqs++;
  seqfile->subseqs[index].name_pos = name_pos;
  seqfile->subseqs[index].name_len = name_len;
  gt4_sequence_file_unlock (seqfile);
  return index;
}
