#ifndef __GT4_FASTA_H__
#define __GT4_FASTA_H__

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

typedef struct _GT4FastaReader GT4FastaReader;
typedef struct _GT4FastaReaderClass GT4FastaReaderClass;

#define GT4_TYPE_FASTA_READER gt4_fasta_reader_get_type()

#include <stdlib.h>
#include <stdio.h>

#include "sequence-source.h"

#define MAX_NAME_SIZE 1023

#define GT4FR_FASTA 1
#define GT4FR_FASTQ 2
#define FASTA_READER_STATE_NONE 0
#define FASTA_READER_STATE_NAME 1
#define FASTA_READER_STATE_SEQUENCE 2
#define FASTA_READER_STATE_QUALITY 3

struct _GT4FastaReader {
  /* Read settings */
  unsigned int wordlength;
  unsigned long long mask;
  unsigned int canonize;

  /* I/O */
  /* 0 - EOF, negative - error */
  GT4SequenceSourceImplementation *impl;
  GT4SequenceSourceInstance *inst;
  unsigned int in_eof;
  
  /* FastQ or FastA */
  unsigned int type;
  /* Current reading state */
  unsigned int state;
  /* Current character, and word number */
  unsigned long long cpos;
  unsigned long long wpos;
  /* Current nucleotide number relative to the start of subsequence */
  unsigned long long seq_npos;
  /* Current name */
  unsigned long long name_pos;
  unsigned int name_length;
  unsigned int name_idx;
  unsigned char name[MAX_NAME_SIZE + 1];
  /* Words */
  unsigned long long wordfw;
  unsigned long long wordrv;
  unsigned int currentlength;
};

struct _GT4FastaReaderClass {
  AZClass klass;
};

unsigned int gt4_fasta_reader_get_type (void);

/* Set up FastA reader structure */
int fasta_reader_init (GT4FastaReader *reader, unsigned int wordlength, unsigned int canonize, GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
int fasta_reader_release (GT4FastaReader *reader);

/* Read maximum of nwords words from FastA or fastQ file starting from position cpos */
/* If any callback returns not 0, reading stops and result is returned */
int fasta_reader_read_nwords (GT4FastaReader *reader, unsigned long long maxwords,
  /* Called as soon as the full sequence name is known */
  int (*start_sequence) (GT4FastaReader *, void *),
  /* Called when the full sequence has been parsed */
  int (*end_sequence) (GT4FastaReader *, void *),
  int (*read_character) (GT4FastaReader *, unsigned int character, void *),
  int (*read_nucleotide) (GT4FastaReader *, unsigned int nucleotide, void *),
  int (*read_word) (GT4FastaReader *, unsigned long long word, void *),
  void *data);

#endif
