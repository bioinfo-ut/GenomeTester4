#ifndef FASTA_H_
#define FASTA_H_

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
#include <stdio.h>

#include "buffer.h"

#define MAX_NAME_SIZE 1000

#define FASTA_READER_STATE_NONE 0
#define FASTA_READER_STATE_NAME 1
#define FASTA_READER_STATE_SEQUENCE 2
#define FASTA_READER_STATE_QUALITY 3

typedef struct _FastaReader {
	/* Read settings */
	unsigned int wordlength;
	unsigned long long mask;
	unsigned int canonize;

	/* I/O */
	/* 0 - EOF, negative - error */
	int (* read) (void *data);
	void (* free_io_data) (void *data);
	void *read_data;
	unsigned int in_eof;
	
	/* State */
	unsigned int state;
	/* Current character, nucleotide and word number */
	unsigned long long cpos;
	unsigned long long npos;
	unsigned long long wpos;
	/* Name */
	unsigned char name[MAX_NAME_SIZE + 1];
	unsigned int name_length;
	/* Words */
	unsigned long long wordfw;
	unsigned long long wordrv;
	unsigned int currentlength;
} FastaReader;

/* Set up FastA reader structure */
int fasta_reader_init (FastaReader *reader, unsigned int wordlength, unsigned int canonize, int (* read) (void *), void *read_data);
int fasta_reader_release (FastaReader *reader);

int fasta_reader_init_from_data (FastaReader *reader, unsigned int wordlength, unsigned int canonize, const unsigned char *cdata, unsigned long long csize);
int fasta_reader_init_from_file (FastaReader *reader, unsigned int wordlength, unsigned int canonize, FILE *ifs);

/* Read maximum of nwords words from FastA or fastQ file starting from position cpos */
int fasta_reader_read_nwords (FastaReader *reader, unsigned long long maxwords,
	/* Called as soon as the full sequence name is known */
	int (*start_sequence) (FastaReader *, void *),
	/* Called when the full sequence has been parsed */
	int (*end_sequence) (FastaReader *, void *),
	int (*read_character) (FastaReader *, unsigned int character, void *),
	int (*read_nucleotide) (FastaReader *, unsigned int nucleotide, void *),
	int (*read_word) (FastaReader *, unsigned long long word, void *),
	void *data);
                
#endif /* FASTA_H_ */
