#ifndef COMMON_H_
#define COMMON_H_

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Cpyright (C) 2014 University of Tartu
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

#include <stdio.h>

#define VERSION_MAJOR 4
#define VERSION_MINOR 0

/* errors and warnings */
#define GT_INCOMPATIBLE_WORDLENGTH_ERROR 2
#define GT_INCOMPATIBLE_VERSION_WARNING 3
#define GT_NEWER_VERSION_WARNING 4
#define GT_FASTA_START_SEQ_ERROR 5
#define GT_FASTA_END_SEQ_ERROR 6
#define GT_FASTA_PROCESS_SEQ_ERROR 7
#define GT_OUT_OF_MEMORY_ERROR 8

int print_error_message (int e);



#endif /* COMMON_H_ */
