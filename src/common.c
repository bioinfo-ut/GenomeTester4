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
#include <stdlib.h>
#include "common.h"

int print_error_message (int e)
{
	if (e == 1) {
		return 1;
	} else if (e == GT_INCOMPATIBLE_WORDLENGTH_ERROR) {
		fprintf (stderr, "Error: Incompatible wordlengths!\n");
		return 1;
	} else if (e == GT_OUT_OF_MEMORY_ERROR) {
		fprintf (stderr, "Error: Program out of memory!\n");
		return 1;
	} else if (e == GT_INCOMPATIBLE_VERSION_WARNING) {
		fprintf (stderr, "Warning: Two lists are built under different glistmaker versions!\n");
	} else if (e == GT_NEWER_VERSION_WARNING) {
		fprintf (stderr, "Warning: List is built under newer glistmaker version than your current program version!\n");
	} else if (e == GT_FASTA_START_SEQ_ERROR) {
		fprintf (stderr, "Error: Processing the start of the sequence failed!\n");
		return 1;
	} else if (e == GT_FASTA_END_SEQ_ERROR) {
		fprintf (stderr, "Error: Processing the end of the sequence failed!\n");
		return 1;
	} else if (e == GT_FASTA_PROCESS_SEQ_ERROR) {
		fprintf (stderr, "Error: Processing the sequence failed!\n");
		return 1;
	} else {
		fprintf (stderr, "Error: Program failed due to unknown reasons!\n");
		return 1;
	}
	return 0;
}
 
