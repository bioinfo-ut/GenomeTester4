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
#include <assert.h>
#include "buffer.h"

SequenceBuffer *
sequence_buffer_new () 
{
	SequenceBuffer *sbuffer;
	sbuffer = (SequenceBuffer *) malloc (sizeof (SequenceBuffer));
	memset (sbuffer, 0, sizeof (SequenceBuffer));
	return sbuffer;  
}

void
sequence_buffer_add_char (SequenceBuffer *sbuffer, char c) 
{
	if (c <= ' ') sbuffer->nwhitespace += 1;
	sbuffer->buffer[sbuffer->wi] = c;
	sbuffer->wi = sbuffer->wi < MAX_SEQ_BUFFER_SIZE - 1 ? sbuffer->wi + 1 : 0;
	//printf ("lisan tähe: %c, wi = %u\n", c, sbuffer->wi);
	return;
}

char
sequence_buffer_get_next_char (SequenceBuffer *sbuffer)
{
	char c = sbuffer->buffer[sbuffer->ri];
	if (c <= ' ') sbuffer->nwhitespace -= 1;
	sbuffer->ri = sbuffer->ri < MAX_SEQ_BUFFER_SIZE - 1 ? sbuffer->ri + 1 : 0;
	//printf ("võtan tähe: %c, ri = %u\n", c, sbuffer->ri);
	return c;
}

unsigned int
sequence_buffer_get_index_gap (SequenceBuffer *sbuffer)
{
	if (sbuffer->wi >= sbuffer->ri) return sbuffer->wi - sbuffer->ri;
	return MAX_SEQ_BUFFER_SIZE - sbuffer->ri + sbuffer->wi;
}

void 
sequence_buffer_flush_until (SequenceBuffer *sbuffer, FILE *out, unsigned int n)
{
	while (sequence_buffer_get_index_gap (sbuffer) > n) {
		fprintf(out, "%c", sequence_buffer_get_next_char (sbuffer));
	}
	return;	
}

void
sequence_buffer_flush_all (SequenceBuffer *sbuffer, FILE *out) 
{
	sequence_buffer_flush_until (sbuffer, out, 0);
	return;
}

char
sequence_buffer_skip_char (SequenceBuffer *sbuffer)
{
	char c = sbuffer->buffer[sbuffer->ri];
	assert (sequence_buffer_get_index_gap (sbuffer) >= 1);
	if (c <= ' ') sbuffer->nwhitespace -= 1;
	sbuffer->ri += 1;
	return c;
}

int
sequence_buffer_is_full (SequenceBuffer *sbuffer)
{
	if (sequence_buffer_get_index_gap (sbuffer) == MAX_SEQ_BUFFER_SIZE - 1) return 1;
	return 0;
}




