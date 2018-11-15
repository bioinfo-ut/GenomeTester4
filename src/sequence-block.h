#ifndef __GT4_SEQUENCE_BLOCK_H__
#define __GT4_SEQUENCE_BLOCK_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2016 University of Tartu
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

/*
 * A block of FastA/FastQ data in memory
 */

typedef struct _GT4SequenceBlock GT4SequenceBlock;
typedef struct _GT4SequenceBlockClass GT4SequenceBlockClass;

#define GT4_TYPE_SEQUENCE_BLOCK (gt4_sequence_block_get_type ())
#define GT4_SEQUENCE_BLOCK(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_SEQUENCE_BLOCK, GT4SequenceBlock))
#define GT4_IS_SEQUENCE_BLOCK(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_SEQUENCE_BLOCK))

#define GT4_SEQUENCE_BLOCK_FROM_SEQUENCE_SOURCE_INSTANCE(i) (GT4SequenceBlock *) AZ_BASE_ADDRESS(GT4SequenceBlock,source_inst,i)
#define GT4_SEQUENCE_BLOCK_SEQUENCE_SOURCE_IMPLEMENTATION(o) &((GT4SequenceBlockClass *) ((AZObject *) (o))->klass)->source_impl
#define GT4_SEQUENCE_BLOCK_FROM_SEQUENCE_COLLECTION_INSTANCE(i) (GT4SequenceBlock *) AZ_BASE_ADDRESS(GT4SequenceBlock,collection_inst,i)
#define GT4_SEQUENCE_BLOCK_SEQUENCE_COLLECTION_IMPLEMENTATION(o) &((GT4SequenceBlockClass *) ((AZObject *) (o))->klass)->collection_impl

#include <az/object.h>

#include "sequence-collection.h"
#include "sequence-source.h"

struct _GT4SequenceBlock {
  AZObject object;
  /* Reference of parent object to manage chained blocks */
  AZObject *parent;
  const unsigned char *cdata;
  unsigned long long csize;
  unsigned long long pos;
  /* GT4SequenceSource instance */
  GT4SequenceSourceInstance source_inst;
  /* GT4SequenceCollection instance */
  GT4SequenceCollectionInstance collection_inst;
  /* Subsequence data */
  unsigned int n_subseqs;
  GT4SubSequence *subseqs;
  /* Bookkeeping */
  unsigned int size_subseqs;
  unsigned int subseq_block_size;
};

struct _GT4SequenceBlockClass {
  AZObjectClass object_class;
  /* GT4SequenceSource implementation */
  GT4SequenceSourceImplementation source_impl;
  /* GT4SequenceCollection implementation */
  GT4SequenceCollectionImplementation collection_impl;
};

unsigned int gt4_sequence_block_get_type (void);

GT4SequenceBlock *gt4_sequence_block_new (const unsigned char *cdata, unsigned long long csize, AZObject *parent);

/* Split block into array of smaller blocks, return the number of child blocks created */
unsigned int gt4_sequence_block_split (GT4SequenceBlock *blk, GT4SequenceBlock *child_blocks[], unsigned int n_child_blocks);

#endif
