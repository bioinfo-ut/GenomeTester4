#define __GT4_WORD_LIST_SORTED_C__

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

#include "word-list-sorted.h"

static void slist_class_init (GT4WordSListClass *klass);

unsigned int slist_type = 0;
GT4WordSListClass *slist_class = 0;

unsigned int
gt4_word_slist_get_type (void)
{
  if (!slist_type) {
    slist_class = (GT4WordSListClass *) az_register_interface_type (&slist_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4WordArraySorted",
      sizeof (GT4WordSListClass), sizeof (GT4WordSListImplementation), sizeof (GT4WordSListInstance),
      (void (*) (AZClass *)) slist_class_init,
      NULL, NULL, NULL);
  }
  return slist_type;
}

static void
slist_class_init (GT4WordSListClass *klass)
{
  klass->interface_class.klass.flags = AZ_CLASS_ZERO_MEMORY;
}

unsigned int
gt4_word_slist_get_first_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  inst->idx = 0;
  if (!inst->num_words) return 0;
  return impl->get_first_word (impl, inst);
}

unsigned int
gt4_word_slist_get_next_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  if (inst->idx >= inst->num_words) return 0;
  inst->idx += 1;
  if (inst->idx >= inst->num_words) return 0;
  return impl->get_next_word (impl, inst);
}
