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
#include <string.h>
#include <stdio.h>

#include "buffer.h"
#include "fasta.h"
#include "sequence.h"
#include "utils.h"
#include "common.h"

int
fasta_reader_init (FastaReader *reader, const char *cdata, size_t csize, unsigned int wordlength, int canonize, int usebuffer)
{
	memset (reader, 0, sizeof (FastaReader));
	reader->wordlength = wordlength;
	reader->canonize = canonize;
	reader->cdata = cdata;
	reader->csize = csize;
	reader->mask = create_mask (wordlength);
	if (usebuffer) reader->sbuffer = sequence_buffer_new ();
	return 0;
}

FastaReader *
fasta_reader_new (const char *cdata, size_t csize, unsigned int wordlength, int canonize, int usebuffer)
{
	FastaReader *reader;
	int result;
	reader = (FastaReader *) malloc (sizeof (FastaReader));
	result = fasta_reader_init (reader, cdata, csize, wordlength, canonize, usebuffer);
	if (result) {
		free (reader);
		return NULL;
	}
	return reader;
}


int
fasta_reader_read_nwords (FastaReader *reader, unsigned long long maxwords, void *data, int usebuffer,
	int (*start_sequence) (FastaReader *reader, void *data),
	int (*end_sequence) (FastaReader *reader, void *data),
	int (*process_sequence) (FastaReader *reader, void *data))
{
	int isheader = 0, isquality = 0;
	size_t headerloc = 0;
	unsigned long long nuclval, nquality = 0L;
	unsigned long long nwords = 0;
	
	for (; reader->cpos < reader->csize; reader->cpos++) {
		size_t i = reader->cpos;
		
		/* Return as soon as max words is reached */
		/* FastaReader should retain all relevant states */
		/* We do it here so that cpos would be updated from last word */
		if (nwords >= maxwords) {
			return 0;
		}

		/* iteration over quality scores */
		/* this is necessary because FastQ sequences may span over many lines and/or quality scores may
		 * contain @, > or + characters */
		if (!isheader && isquality && (nquality != reader->nnuclinseq + reader->nmaskedinseq)) {
			if (reader->cdata[i] != 10 || reader->cdata[i] != 13) nquality += 1;
			continue;
		}

		/* The beginning of the header of the next sequence */
		if (reader->cdata[i] == '>' || reader->cdata[i] == '@') {
			if (i != 0 && (end_sequence)) {
				int v = (*end_sequence) (reader, data);
				if (v > 0) return GT_FASTA_END_SEQ_ERROR;
			}
			isheader = 1;
			isquality = 0;
			headerloc = i + 1;
			reader->currentlength = 0;
			reader->wordfw = 0L;
			reader->wordrv = 0L;
			reader->currentword = 0L;
			reader->ncharinseq = 0L;
			reader->nwhitespaceinseq = 0L;
			reader->nmaskedinseq = 0L;
			reader->nnuclinseq = 0L;
			nquality = 0;
		}

		/* Header */
		if (isheader) {
			if (reader->cdata[i] == 10 || reader->cdata[i] == 13) {
				if (!isquality) {
					if (i - headerloc > MAX_NAME_SIZE - 1) {
						/* if header is too long only part of it is kept */
						strncpy (reader->seqname, reader->cdata + headerloc, MAX_NAME_SIZE - 1);
					} else {
						strncpy (reader->seqname, reader->cdata + headerloc, i - headerloc);
					}
					reader->seqname[i - headerloc + 1] = 0; /* end of string character */
					reader->seqstart = i + 1;
					if (start_sequence) {
						int v = (*start_sequence) (reader, data);
						if (v > 0) return GT_FASTA_START_SEQ_ERROR;
					}
				}
				isheader = 0;
			}
		/* Sequence part */
		} else {
			reader->ncharinseq += 1;

			/* Already masked nucleotide or quality */
			if (strchr (alphabet, reader->cdata[i]) == NULL && reader->cdata[i] > ' ') {
				if (reader->cdata[i] == '+') {
					isquality = 1;
					/* there might be another header after + character */
					isheader = 1;
					headerloc = i + 1;
				} 
				if (!isquality && usebuffer) sequence_buffer_add_char (reader->sbuffer, reader->cdata[i]);
				reader->nmaskedinseq += 1;
				reader->wordfw = 0L;
				reader->wordrv = 0L;
				reader->currentword = 0L;
				reader->currentlength = 0;
				continue;
			/* Newline character */	
			} else if (reader->cdata[i] <= ' ') {
				if (usebuffer) sequence_buffer_add_char (reader->sbuffer, reader->cdata[i]);
				reader->nwhitespaceinseq += 1;
				continue;
			}

			/* Next not masked nucleotide */
			if (usebuffer) sequence_buffer_add_char (reader->sbuffer, reader->cdata[i]);
			/* add to the counts */
			reader->nnuclinseq += 1;

			/* add next nucleotide to the word */
			reader->wordfw <<= 2;
			nuclval = get_nucl_value (reader->cdata[i]);
			reader->wordfw |= nuclval;

			if (reader->canonize) {
				reader->wordrv >>= 2;
				reader->wordrv |= ((~nuclval & 3) << ((reader->wordlength - 1) * 2));
			}

			reader->currentlength += 1;
			if (reader->currentlength > reader->wordlength) {
				reader->wordfw &= reader->mask;
				reader->currentlength = reader->wordlength;
			}
			if (reader->currentlength == reader->wordlength) {
				/* Update current word */
				reader->currentword = (!reader->canonize || reader->wordfw < reader->wordrv) ? reader->wordfw : reader->wordrv;
				if (process_sequence) {
					int v = (*process_sequence) (reader, data);
					if (v > 0) return GT_FASTA_PROCESS_SEQ_ERROR;
				}
				nwords += 1;
			}
		}
	}
	if (end_sequence) {
		int v = (*end_sequence) (reader, data);
		if (v > 0) return GT_FASTA_END_SEQ_ERROR;
	}	
	return 0;
}

int 
fasta_reader_read (const char *fasta, size_t fsize, unsigned int wordlength, void *data, int canonize, int usebuffer,
		int (*start_sequence) (FastaReader *reader, void *data), int (*end_sequence) (FastaReader *reader, void *data),
		int (*process_sequence) (FastaReader *reader, void *data))
{ 
	FastaReader *reader = fasta_reader_new (fasta, fsize, wordlength, canonize, usebuffer);

	while (reader->cpos < reader->csize) {
		int result;
		result = fasta_reader_read_nwords (reader, MAX_WORDS_READ_PER_ROUND, data, usebuffer, start_sequence, end_sequence, process_sequence);
		if (result) return result;
	}
	return 0;
}



