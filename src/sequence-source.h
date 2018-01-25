#ifndef __GT4_SEQUENCE_SOURCE_H__
#define __GT4_SEQUENCE_SOURCE_H__

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

typedef struct _GT4SequenceSourceInstance GT4SequenceSourceInstance;
typedef struct _GT4SequenceSourceImplementation GT4SequenceSourceImplementation;
typedef struct _GT4SequenceSourceClass GT4SequenceSourceClass;

#define GT4_TYPE_SEQUENCE_SOURCE (gt4_sequence_source_get_type ())

#define GT4_SEQUENCE_SOURCE_INVALID_INTERFACE -1

#include <az/interface.h>

struct _GT4SequenceSourceInstance {
  unsigned long long position;
  unsigned int open : 1;
  unsigned int eof : 1;
  unsigned int error : 1;
};
      
struct _GT4SequenceSourceImplementation {
  AZImplementation implementation;
  /* Make data available for reading and reset position */
  unsigned int (* open) (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
  /* Read one byte of data */
  int (* read) (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
  /* Release resources if needed */
  unsigned int (* close) (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
};

struct _GT4SequenceSourceClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_sequence_source_get_type (void);

/* Return 1 on succes */
unsigned int gt4_sequence_source_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
/* Return 0 on proper EOF, negative value on error */
int gt4_sequence_source_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
/* Return 1 on succes */
unsigned int gt4_sequence_source_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);

#endif
