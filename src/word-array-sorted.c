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

#include "word-array-sorted.h"

static void sarray_class_init (GT4WordSArrayClass *klass);

unsigned int sarray_type = 0;
GT4WordSArrayClass *sarray_class = 0;

unsigned int
gt4_word_sarray_get_type (void)
{
  if (!sarray_type) {
    sarray_class = (GT4WordSArrayClass *) az_register_interface_type (&sarray_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4WordArraySorted",
      sizeof (GT4WordSArrayClass), sizeof (GT4WordSArrayImplementation), sizeof (GT4WordSArrayInstance),
      (void (*) (AZClass *)) sarray_class_init,
      NULL, NULL, NULL);
  }
  return sarray_type;
}

static void
sarray_class_init (GT4WordSArrayClass *klass)
{
  klass->interface_class.klass.flags = AZ_CLASS_ZERO_MEMORY;
}

unsigned int
gt4_word_sarray_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  inst->idx = 0;
  if (!inst->num_words) return 0;
  return impl->get_first_word (impl, inst);
}

unsigned int
gt4_word_sarray_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  if (inst->idx >= inst->num_words) return 0;
  inst->idx += 1;
  if (inst->idx >= inst->num_words) return 0;
  return impl->get_next_word (impl, inst);
}
