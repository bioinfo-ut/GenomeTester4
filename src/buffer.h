#ifndef BUFFER_H_
#define BUFFER_H_

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


#define MAX_SEQ_BUFFER_SIZE 100

typedef struct _SequenceBuffer {
	char buffer[MAX_SEQ_BUFFER_SIZE];
	unsigned int wi;
	unsigned int ri;
	unsigned int nwhitespace;
} SequenceBuffer;

/* create new buffer */
SequenceBuffer *sequence_buffer_new ();

/* a method for finding the lag of the read index with respect to write index */
unsigned int sequence_buffer_get_index_gap (SequenceBuffer *sbuffer);

/* add char to buffer */
void sequence_buffer_add_char (SequenceBuffer *sbuffer, char c);

/* print char from buffer */
char sequence_buffer_get_next_char (SequenceBuffer *sbuffer);

/* empty buffer */
void sequence_buffer_flush_all (SequenceBuffer *sbuffer, FILE *out);
void sequence_buffer_flush_until (SequenceBuffer *sbuffer, FILE *out, unsigned int n);

/* shift read index by 1 */
char sequence_buffer_skip_char (SequenceBuffer *sbuffer);

int sequence_buffer_is_full (SequenceBuffer *sbuffer);


#endif /* BUFFER_H_ */ 
