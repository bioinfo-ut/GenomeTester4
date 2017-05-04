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

#define __SEQUENCE_CPP__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sequence.h"

/* all possible nucleotides */
const char *alphabet = "ACGTUacgtu";

/* Complementary table */
static char *ct = NULL;

unsigned long long get_nucl_value (char nucl)
{
	static unsigned long long bit1 = 1 << 2;
	static unsigned long long bit2 = 3 << 1;

	if (nucl & bit1) {
		return ((nucl >> 4) | 2) & 3;
	}
	return (nucl & bit2) >> 1;
}

unsigned long long create_mask (unsigned int wordlength)
{
	unsigned int i;
	unsigned long long mask = 0L;

	for (i = 0; i < 2 * wordlength; i++) {
		mask = (mask << 1) | 1;
	}
	return mask;
}

unsigned long long get_reverse_complement (unsigned long long word, unsigned int wordlength)
{
	unsigned int i;
	unsigned long long mask, v, revcompl = 0L;

	word = ~word;
	mask = 3;
	for (i = 0; i < wordlength; i++) {
		v = word & mask;
		revcompl <<= 2;
		revcompl |= v;
		word >>= 2;
	}
	return revcompl;
}

unsigned long long get_canonical_word (unsigned long long word, unsigned int wordlength)
{
	unsigned long long rev_word;
	rev_word = get_reverse_complement (word, wordlength);
	return word < rev_word ? word : rev_word;
}

char *word_to_string (unsigned long long word, unsigned int wordlength)
{
	char *s = (char *) malloc (wordlength + 1);
	unsigned int i, temp;

	for (i = 0; i < wordlength; i++) {
		temp = word & 3;
		s[wordlength - i - 1] = alphabet[temp];
		word >>= 2;
	}
	s[wordlength] = 0;
	return s;
}

unsigned int
word2string (char *b, unsigned long long word, unsigned int wordlength)
{
	unsigned int i, temp;

	for (i = 0; i < wordlength; i++) {
		temp = word & 3;
		b[wordlength - i - 1] = alphabet[temp];
		word >>= 2;
	}
	b[wordlength] = 0;
	return wordlength;
}

unsigned long long string_to_word (const char *s, unsigned int wordlength)
{
	unsigned int i, l;
	unsigned long long word = 0L;

	l = (wordlength < 32) ? wordlength : 32;
	for (i = 0; i < l; i++) {
		if (strchr (alphabet, s[i]) == NULL) {
			fprintf (stderr, "Invalid character %c in string!\n", s[i]);
		}
		word <<= 2;
		word |= get_nucl_value (s[i]);
	}
	return word;
}

static void
init_ct (void)
{
  unsigned int i;
  ct = (char *) malloc (256);
  for (i = 0; i < 256; i++) ct[i] = 'N';
  ct['a'] = 't';
  ct['c'] = 'g';
  ct['g'] = 'c';
  ct['t'] = 'a';
  ct['u'] = 'a';
  ct['A'] = 'T';
  ct['C'] = 'G';
  ct['G'] = 'C';
  ct['T'] = 'A';
  ct['U'] = 'A';
}

void
gt4_string_revcomp (char *d, const char *s, unsigned int length, unsigned int terminate)
{
  unsigned int i;
  if (!ct) init_ct ();
  for (i = 0; i < length; i++) {
    d[i] = ct[(unsigned char) s[length - 1 - i]];
  }
  if (terminate) d[length] = 0;
}

void
gt4_string_revcomp_inplace (char *s, unsigned int length)
{
  unsigned int i;
  if (!ct) init_ct ();
  for (i = 0; i < (length / 2); i++) {
    unsigned int t = ct[(unsigned char) s[i]];
    s[i] = ct[(unsigned char) s[length - 1 - i]];
    s[length - 1 - i] = t;
  }
  if (length & 1) s[length / 2] = ct[(unsigned char) s[length / 2]];
}

void word_to_bitstring (unsigned long long word)
{
	unsigned long long mask = (unsigned long long) 1 << 63;

	while (mask != 0) {
		if ((mask & word) != 0) putchar('1');
		else putchar('0');
		mask = mask >> 1;
	}
}



