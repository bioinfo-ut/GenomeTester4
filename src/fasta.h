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
#include "buffer.h"

#define MAX_NAME_SIZE 1000
#define MAX_WORDS_READ_PER_ROUND 1000

typedef struct _FastaReader {
	const char *filename;
	unsigned int wordlength;
	unsigned long long mask;
	int canonize;
	
	/* Sequence data */
	const char *cdata;
	size_t csize;
	
	/* State */
	size_t cpos;
	unsigned long long currentword;
	
	/* Sequence and word info */
	char seqname[MAX_NAME_SIZE];
	unsigned long long seqstart;
	unsigned long long wordfw;
	unsigned long long wordrv;
	unsigned int currentlength;
	
	/* Buffer for sequence characters */
	SequenceBuffer *sbuffer;
	int noword;
	
	/* Counters */
	unsigned long long ncharinseq;
	unsigned long long nwhitespaceinseq; /* white-space characters */
	unsigned long long nmaskedinseq; /* all non-nucleotide and non-white-space characters */
	unsigned long long nnuclinseq;	/* A, C, G, T characters */
} FastaReader;

/* Set up FastA reader structure */
int fasta_reader_init (FastaReader *reader, const char *cdata, size_t csize, unsigned int wordlength, int canonize, int usebuffer);
/* Allocate new fasta reader and initialize it */
FastaReader *fasta_reader_new (const char *cdata, size_t csize, unsigned int wordlength, int canonize, int usebuffer);

/* Read maximum of nwords words from FastA or fastQ file starting from position cpos */
int fasta_reader_read_nwords (FastaReader *reader, unsigned long long maxwords, void *data, int usebuffer,
	int (*start_sequence) (FastaReader *reader, void *data),
	int (*end_sequence) (FastaReader *reader, void *data),
	int (*process_sequence) (FastaReader *reader, void *data));
                
/* Read words from FastA or FastQ to the end of the file, starting from position cpos */
int fasta_reader_read (const char *fasta, size_t fsize, unsigned int wordlength, void *data, int canonize, int usebuffer,
		int (*start_sequence) (FastaReader *reader, void *data), int (*end_sequence) (FastaReader *reader, void *data),
		int (*process_sequence) (FastaReader *reader, void *data));

#endif /* FASTA_H_ */
