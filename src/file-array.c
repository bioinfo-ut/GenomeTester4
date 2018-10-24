#define __GT4_WORD_ARRAY_SORTED_C__

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

/*
 * An interface for accessing words and counts sequentially
 */

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include "file-array.h"

static void file_array_class_init (GT4FileArrayClass *klass);

unsigned int file_array_type = 0;
GT4FileArrayClass *file_array_class = 0;

unsigned int
gt4_file_array_get_type (void)
{
  if (!file_array_type) {
    file_array_class = (GT4FileArrayClass *) az_register_interface_type (&file_array_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4FileArray",
      sizeof (GT4FileArrayClass), sizeof (GT4FileArrayImplementation), sizeof (GT4FileArrayInstance),
      (void (*) (AZClass *)) file_array_class_init,
      NULL, NULL, NULL);
  }
  return file_array_type;
}

static void
file_array_class_init (GT4FileArrayClass *klass)
{
  klass->interface_class.klass.flags = AZ_CLASS_ZERO_MEMORY;
}

unsigned int
gt4_file_array_get_file (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned int idx)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (idx < inst->num_files, 0);
  return impl->get_file (impl, inst, idx);
}

unsigned int
gt4_file_array_get_sequence (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned long long idx)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (idx < inst->n_sequences, 0);
  return impl->get_sequence (impl, inst, idx);
}

