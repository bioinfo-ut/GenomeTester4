#ifndef __GT4_SEQUENCE_COLLECTION_H__
#define __GT4_SEQUENCE_COLLECTION_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2018 University of Tartu
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
 * Interface for managing multi-sequence sources
 *
 */

typedef struct _GT4SequenceCollectionInstance GT4SequenceCollectionInstance;
typedef struct _GT4SequenceCollectionImplementation GT4SequenceCollectionImplementation;
typedef struct _GT4SequenceCollectionClass GT4SequenceCollectionClass;
typedef struct _GT4SubSequence GT4SubSequence;

#define GT4_TYPE_SEQUENCE_COLLECTION (gt4_sequence_collection_get_type ())

#define GT4_SEQUENCE_COLLECTION_INVALID_INTERFACE -1

#include <az/interface.h>

struct _GT4SubSequence {
  /* Start of name from the begining of block */
  unsigned long long name_pos;
  unsigned int name_len;
  /* Start of sequence data relative to name */
  unsigned int seq_pos;
  unsigned int seq_len;
};
            
struct _GT4SequenceCollectionInstance {
  unsigned int writable : 1;
  unsigned int n_subseqs;
  GT4SubSequence subseq;
};
      
struct _GT4SequenceCollectionImplementation {
  AZImplementation implementation;
  unsigned int (*get_subsequence) (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx);
  int (*add_subsequence) (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned long long name_pos, unsigned int name_len);
  unsigned int (*set_subsequence) (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx, unsigned int seq_pos, unsigned int seq_len);
};

struct _GT4SequenceCollectionClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_sequence_collection_get_type (void);

/* Return 0 if error */
unsigned int gt4_sequence_collection_get_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx);
/* Returns new subsequence index, negative if error */
int gt4_sequence_collection_add_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned long long name_pos, unsigned int name_len);
/* Return 0 if error */
unsigned int gt4_sequence_collection_set_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx, unsigned int seq_pos, unsigned int seq_len);

#endif
