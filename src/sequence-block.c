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

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include "sequence-block.h"

static void sequence_block_class_init (GT4SequenceBlockClass *klass);

/* AZObject implementation */
static void sequence_block_shutdown (AZObject *obj);
/* GT4SequenceSource implementation */
static int sequence_block_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);

static unsigned int sequence_block_type = 0;
static GT4SequenceBlockClass *sequence_block_class = NULL;

unsigned int
gt4_sequence_block_get_type (void)
{
  if (!sequence_block_type) {
    sequence_block_class = (GT4SequenceBlockClass *) az_register_type (&sequence_block_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4SequenceBlock",
      sizeof (GT4SequenceBlockClass), sizeof (GT4SequenceBlock),
      (void (*) (AZClass *)) sequence_block_class_init,
      NULL, NULL);
  }
  return sequence_block_type;
}

static void
sequence_block_class_init (GT4SequenceBlockClass *klass)
{
  az_class_set_num_interfaces ((AZClass *) klass, 1);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_SEQUENCE_SOURCE, ARIKKEI_OFFSET(GT4SequenceBlockClass,source_implementation), ARIKKEI_OFFSET(GT4SequenceBlock,source_instance));
  /* AZObject implementation */
  ((AZObjectClass *) klass)->shutdown = sequence_block_shutdown;
  /* GT4SequenceSource implementation */
  klass->source_implementation.read = sequence_block_read;
}

static void
sequence_block_shutdown (AZObject *obj)
{
  GT4SequenceBlock *blk = (GT4SequenceBlock *) obj;
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
  unsigned int child_idx;
  child_idx = 0;
  current_pos = 0;
  remaining_size = blk->csize;
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
        if ((split >= remaining_size) || ((blk->cdata[current_pos + split] == '@') || (blk->cdata[current_pos + split] == '>'))) break;
      }
      current_size = split;
    }
    /* Add new block current_pos - current_size */
    child_blocks[child_idx++] = gt4_sequence_block_new (blk->cdata + current_pos, current_size, &blk->object);
    /* Advance */
    current_pos += current_size;
    remaining_size -= current_size;
  }
  return child_idx;
}
