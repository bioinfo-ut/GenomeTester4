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
fasta_reader_init (FastaReader *reader, unsigned int wordlength, unsigned int canonize, int (* read) (void *), void *read_data)
{
	memset (reader, 0, sizeof (FastaReader));
	reader->wordlength = wordlength;
	reader->mask = create_mask (wordlength);
	reader->canonize = canonize;
	reader->read = read;
	reader->read_data = read_data;
	return 0;
}

int
fasta_reader_release (FastaReader *reader)
{
	if (reader->free_io_data) {
		reader->free_io_data (reader->read_data);
	}
	return 0;
}

struct BufferData {
	const unsigned char *cdata;
	unsigned long long csize;
	unsigned long long cpos;
};

static int
buffer_read (void *data)
{
	struct BufferData *bdata = (struct BufferData *) data;
	if (bdata->cpos >= bdata->csize) return 0;
	return bdata->cdata[bdata->cpos++];
}

static void
buffer_free (void *data)
{
	free (data);
}

int
fasta_reader_init_from_data (FastaReader *reader, unsigned int wordlength, unsigned int canonize, const unsigned char *cdata, unsigned long long csize)
{
	struct BufferData *bdata = (struct BufferData *) malloc (sizeof (struct BufferData));
	int result;
	bdata->cdata = cdata;
	bdata->csize = csize;
	bdata->cpos = 0;
	result = fasta_reader_init (reader, wordlength, canonize, buffer_read, bdata);
	reader->free_io_data = buffer_free;
	return result;
}

static int
file_read (void *data)
{
	FILE *ifs = (FILE *) data;
	int cval = fgetc (ifs);
	if (cval == EOF) return 0;
	return cval;
}

int
fasta_reader_init_from_file (FastaReader *reader, unsigned int wordlength, unsigned int canonize, FILE *ifs)
{
	int result;
	result = fasta_reader_init (reader, wordlength, canonize, file_read, ifs);
	return result;
}

int
fasta_reader_read_nwords (FastaReader *reader, unsigned long long maxwords,
	int (*start_sequence) (FastaReader *, void *),
	int (*end_sequence) (FastaReader *, void *),
	int (*read_character) (FastaReader *, unsigned int character, void *),
	int (*read_nucleotide) (FastaReader *, unsigned int nucleotide, void *),
	int (*read_word) (FastaReader *, unsigned long long word, void *),
	void *data)
{
	unsigned long long nwords = 0;
	static unsigned int *c2n = NULL;
	if (!c2n) {
		c2n = (unsigned int *) malloc (256 * sizeof (unsigned int));
		memset (c2n, 0xff, 256 * sizeof (unsigned int));
		c2n['A'] = c2n['a'] = 0;
		c2n['C'] = c2n['c'] = 1;
		c2n['G'] = c2n['g'] = 2;
		c2n['T'] = c2n['t'] = c2n['U'] = c2n['u'] = 3;
	}
	
	while (!reader->in_eof && (nwords < maxwords)) {
		int cval = reader->read (reader->read_data);
		if (cval < 0) return cval;
		if (cval == 0) {
			reader->in_eof = 1;
			if (reader->state == FASTA_READER_STATE_SEQUENCE) {
				if (end_sequence) {
					int result = end_sequence (reader, data);
					if (result) return result;
				}
			}
			return 0;
		}
		if (read_character) {
			int result = read_character (reader, cval, data);
			if (result) return result;
		}
		reader->cpos += 1;
		/* Proper character */
		switch (reader->state) {
		case FASTA_READER_STATE_NONE:
			if (cval == '>') {
				reader->state = FASTA_READER_STATE_NAME;
				reader->name_length = 0;
			}
			break;
		case FASTA_READER_STATE_NAME:
			if (cval == '\n') {
				/* Name is complete */
				reader->name[reader->name_length] = 0;
				if (start_sequence) {
					int result = start_sequence (reader, data);
					if (result) return result;
				}
				reader->state = FASTA_READER_STATE_SEQUENCE;
				reader->wordfw = 0;
				reader->wordrv = 0;
				reader->currentlength = 0;
			} else {
				/* Append character to name */
				if (reader->name_length < MAX_NAME_SIZE) {
					reader->name[reader->name_length++] = cval;
				}
			}
			break;
		case FASTA_READER_STATE_SEQUENCE:
			if ((cval == '>') || (cval == '+')) {
				/* End of sequence */
				if (end_sequence) {
					int result = end_sequence (reader, data);
					if (result) return result;
				}
				if (cval == '>') {
					reader->state = FASTA_READER_STATE_NAME;
					reader->name_length = 0;
				} else {
					reader->state = FASTA_READER_STATE_QUALITY;
				}
			} else {
				/* Process */
				unsigned int nuclval = c2n[cval];
				if (nuclval != 0xffffffff) {
					/* Is nucleotide */
					if (read_nucleotide) {
						int result = read_nucleotide (reader, nuclval, data);
						if (result) return result;
					}
					reader->npos += 1;
					reader->wordfw <<= 2;
					reader->wordfw |= nuclval;
					if (reader->canonize) {
						reader->wordrv >>= 2;
						reader->wordrv |= ((unsigned long long) (~nuclval & 3) << ((reader->wordlength - 1) * 2));
					}
					reader->currentlength += 1;
					if (reader->currentlength > reader->wordlength) {
						reader->wordfw &= reader->mask;
						reader->currentlength = reader->wordlength;
					}
					if (reader->currentlength == reader->wordlength) {
						/* Update current word */
						unsigned long long word = (!reader->canonize || reader->wordfw < reader->wordrv) ? reader->wordfw : reader->wordrv;
						if (read_word) {
							int result = read_word (reader, word, data);
							if (result) return result;
						}
						reader->wpos += 1;
						nwords += 1;
					}
				} else if (cval >= ' ') {
					reader->wordfw = 0;
					reader->wordrv = 0;
					reader->currentlength = 0;
				}
			}
			break;
		case FASTA_READER_STATE_QUALITY:
			if (cval == '@') {
				reader->state = FASTA_READER_STATE_NAME;
				reader->name_length = 0;
			}
			break;
		}
	}
	return 0;
}
                                                

#if 0
fasta_reader_read_nwords (FastaReader *reader, unsigned long long maxwords, void *data, int usebuffer,
	int (*start_sequence) (FastaReader *reader, void *data),
	int (*end_sequence) (FastaReader *reader, void *data),
	int (*process_sequence) (FastaReader *reader, void *data))
{
	int isheader = 0, isquality = 0;
	size_t headerloc = 0;
	unsigned long long nuclval, nquality = 0L;
	unsigned long long nwords = 0;

	while (!is_eof) {
		int cval = reader->read (reader->input_data);
		if (val < 0) {
			reader->is_eof = 1;
			val = '\n';
		}
		
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
			if (cval != 10 || cval != 13) nquality += 1;
			reader->cpos += 1;
			continue;
		}

		/* The beginning of the header of the next sequence */
		if (cval == '>' || cval == '@') {
			if (reader->cpos != 0 && end_sequence) {
				int v = (*end_sequence) (reader, data);
				if (v > 0) return GT_FASTA_END_SEQ_ERROR;
			}
			isheader = 1;
			isquality = 0;
			headerloc = reader->cpos + 1;
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
			} else {
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
				reader->cpos += 1;
				continue;
			/* Newline character */	
			} else if (reader->cdata[i] <= ' ') {
				if (usebuffer) sequence_buffer_add_char (reader->sbuffer, reader->cdata[i]);
				reader->nwhitespaceinseq += 1;
				reader->cpos += 1;
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
		reader->cpos += 1;
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
#endif



