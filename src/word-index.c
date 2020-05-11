#define __GT4_WORD_INDEX_C__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2020 University of Tartu
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
 * An interface for accessing location data for words
 * Word lookup has to be implemented by separate interface
 */

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include "word-index.h"

static void index_class_init (GT4WordIndexClass *klass);

unsigned int index_type = 0;
GT4WordIndexClass *index_class = 0;

unsigned int
gt4_word_index_get_type (void)
{
  if (!index_type) {
    index_class = (GT4WordIndexClass *) az_register_interface_type (&index_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4WordIndex",
      sizeof (GT4WordIndexClass), sizeof (GT4WordIndexImplementation), sizeof (GT4WordIndexInstance),
      (void (*) (AZClass *)) index_class_init,
      NULL, NULL, NULL);
  }
  return index_type;
}

static void
index_class_init (GT4WordIndexClass *klass)
{
  klass->interface_class.klass.flags = AZ_CLASS_ZERO_MEMORY;
}

unsigned int
gt4_word_index_get_location (GT4WordIndexImplementation *impl, GT4WordIndexInstance *inst, unsigned long long idx)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  return impl->get_location (impl, inst, idx);
}
