#ifndef __GT4_SET_OPERATIONS_H__
#define __GT4_SET_OPERATIONS_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2018 University of Tartu
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

#include "word-map.h"

#define GT4_MAX_SETS 4096

/* Combines N lists into one union, returns 0 on success */
/* Arrays have to be AZObjects implementing GT4SortedWArray interface */

unsigned int gt4_write_union (AZObject *arrays[], unsigned int n_arrays, unsigned int cutoff, int ofile, GT4ListHeader *header);

#endif
