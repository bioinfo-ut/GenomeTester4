#define __GT4_SEQUENCE_BLOCK_C__

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

#include <stdio.h>
#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include "sequence-block.h"

extern unsigned int debug;

static void sequence_block_class_init (GT4SequenceBlockClass *klass);
static void sequence_block_init (GT4SequenceBlockClass *klass, GT4SequenceBlock *blk);

/* AZObject implementation */
static void sequence_block_shutdown (AZObject *obj);
/* GT4SequenceSource implementation */
static int sequence_block_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
/* GT4SequenceCollection implementation */
static unsigned int collection_get_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx);
static int collection_add_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned long long name_pos, unsigned int name_len);
static unsigned int collection_set_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx, unsigned int sequence_pos, unsigned int sequence_len);

static unsigned int sequence_block_type = 0;
static GT4SequenceBlockClass *sequence_block_class = NULL;

unsigned int
gt4_sequence_block_get_type (void)
{
  if (!sequence_block_type) {
    sequence_block_class = (GT4SequenceBlockClass *) az_register_type (&sequence_block_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4SequenceBlock",
      sizeof (GT4SequenceBlockClass), sizeof (GT4SequenceBlock),
      (void (*) (AZClass *)) sequence_block_class_init,
      (void (*) (AZImplementation *, void *)) sequence_block_init,
      NULL);
  }
  return sequence_block_type;
}

static void
sequence_block_class_init (GT4SequenceBlockClass *klass)
{
  az_class_set_num_interfaces ((AZClass *) klass, 2);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_SEQUENCE_SOURCE, ARIKKEI_OFFSET(GT4SequenceBlockClass,source_impl), ARIKKEI_OFFSET(GT4SequenceBlock,source_inst));
  az_class_declare_interface ((AZClass *) klass, 1, GT4_TYPE_SEQUENCE_COLLECTION, ARIKKEI_OFFSET(GT4SequenceBlockClass,collection_impl), ARIKKEI_OFFSET(GT4SequenceBlock,collection_inst));
  /* AZObject implementation */
  ((AZObjectClass *) klass)->shutdown = sequence_block_shutdown;
  /* GT4SequenceSource implementation */
  klass->source_impl.read = sequence_block_read;
  /* GT4SequenceCollection implementation */
  klass->collection_impl.get_subsequence = collection_get_subsequence;
  klass->collection_impl.add_subsequence = collection_add_subsequence;
  klass->collection_impl.set_subsequence = collection_set_subsequence;
}

static void
sequence_block_init (GT4SequenceBlockClass *klass, GT4SequenceBlock *blk)
{
  blk->subseq_block_size = 1024;
}

static void
sequence_block_shutdown (AZObject *obj)
{
  GT4SequenceBlock *blk = (GT4SequenceBlock *) obj;
  if (blk->subseqs) {
    free (blk->subseqs);
    blk->subseqs = NULL;
  }
  if (blk->parent) {
    az_object_unref (AZ_OBJECT (blk->parent));
    blk->parent = NULL;
  }
}

static int
sequence_block_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceBlock *blk = GT4_SEQUENCE_BLOCK_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  if (blk->pos >= blk->csize) return 0;
  return blk->cdata[blk->pos++];
}

static unsigned int
collection_get_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx)
{
  return 0;
}

static int
collection_add_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned long long name_pos, unsigned int name_len)
{
  GT4SequenceBlock *blk = GT4_SEQUENCE_BLOCK_FROM_SEQUENCE_COLLECTION_INSTANCE(inst);
  unsigned int index;
  if (blk->n_subseqs >= blk->size_subseqs) {
    blk->size_subseqs = blk->size_subseqs << 1;
    if (blk->size_subseqs < blk->subseq_block_size) blk->size_subseqs = blk->subseq_block_size;
    blk->subseqs = (GT4SubSequence *) realloc (blk->subseqs, blk->size_subseqs * sizeof (GT4SubSequence));
  }
  index = blk->n_subseqs++;
  blk->subseqs[index].name_pos = name_pos;
  blk->subseqs[index].name_len = name_len;
  return index;
}

static unsigned int
collection_set_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx, unsigned int sequence_pos, unsigned int sequence_len)
{
  return 0;
}

GT4SequenceBlock *
gt4_sequence_block_new (const unsigned char *cdata, unsigned long long csize, AZObject *parent)
{
  GT4SequenceBlock *blk;
  arikkei_return_val_if_fail (cdata != NULL, NULL);
  blk = (GT4SequenceBlock *) az_object_new (GT4_TYPE_SEQUENCE_BLOCK);
  blk->parent = parent;
  if (blk->parent) az_object_ref (AZ_OBJECT (blk->parent));
  blk->cdata = cdata;
  blk->csize = csize;
  return blk;
}

unsigned int
gt4_sequence_block_split (GT4SequenceBlock *blk, GT4SequenceBlock *child_blocks[], unsigned int n_child_blocks)
{
  unsigned long long remaining_size, current_pos, current_size;
  unsigned int fastq, child_idx;
  child_idx = 0;
  current_pos = 0;
  remaining_size = blk->csize;
  fastq = (blk->cdata[0] == '@');
  while ((child_idx < n_child_blocks) && (current_pos < blk->csize)) {
    unsigned long long split;
    unsigned int n_splits;
    /* Number of splits in remaining block */
    n_splits = n_child_blocks - child_idx;
    if (n_splits == 1) {
      current_size = remaining_size;
    } else {
      /* Find split position */
      split = remaining_size / n_splits;
      while (split < remaining_size) {
        while (split < remaining_size) {
          if (blk->cdata[current_pos + split] == '\n') break;
          split += 1;
        }
        split += 1;
        if (!fastq) {
          /* Break if FastA name */
          if ((split >= remaining_size) || (blk->cdata[current_pos + split] == '>')) break;
        } else {
          /* If FastQ separator search for @...\n...\n+ pattern */
          if ((split < remaining_size) && (blk->cdata[current_pos + split] == '@')) {
            unsigned int s = split + 1;
            /* Skip name */
            while ((s < remaining_size) && (blk->cdata[current_pos + s] != '\n')) s += 1;
            s += 1;
            /* Skip sequence */
            while ((s < remaining_size) && (blk->cdata[current_pos + s] != '\n')) s += 1;
            s += 1;
            /* If next line starts with '+' we have found it */
            if ((s < remaining_size) && (blk->cdata[current_pos + s] == '+')) break;
          }
        }
      }
      current_size = split;
    }
    /* Add new block current_pos - current_size */
    if (debug) {
      unsigned long long i;
      fprintf (stderr, "Block start:");
      for (i = 0; i < 100; i++) fprintf (stderr, "%c", blk->cdata[current_pos + i]);
      fprintf (stderr, "\n");
    }
    child_blocks[child_idx++] = gt4_sequence_block_new (blk->cdata + current_pos, current_size, &blk->object);
    /* Advance */
    current_pos += current_size;
    remaining_size -= current_size;
  }
  return child_idx;
}
