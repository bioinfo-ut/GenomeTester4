#ifndef SEQUENCE_H_
#define SEQUENCE_H_

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

#ifndef __SEQUENCE_CPP__
extern const char *alphabet;
#endif

#include "wordtable.h"

unsigned long long get_nucl_value (char nucl);

unsigned long long create_mask (unsigned int wordlength);

unsigned long long get_reverse_complement (unsigned long long word, unsigned int wordlength);

unsigned long long get_canonical_word (unsigned long long word, unsigned int wordlength);

unsigned long long generate_mismatches (wordtable *mmtable, unsigned long long word, unsigned int wordlength,
		unsigned int givenfreq, unsigned int nmm, unsigned int startsite, int usesmallercomplement, int countonly,
		int equalmmonly);

char *word_to_string (unsigned long long word, unsigned int wordlength);
unsigned int word2string (char *b, unsigned long long word, unsigned int wordlength);

unsigned long long string_to_word (const char *s, unsigned int wordlength);

void word_to_bitstring (unsigned long long word);


#endif /* SEQUENCE_H_ */
