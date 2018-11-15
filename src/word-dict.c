#define __GT4_WORD_DICT_C__

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

#include "word-dict.h"

static void dict_class_init (GT4WordDictClass *klass);

unsigned int dict_type = 0;
GT4WordDictClass *dict_class = 0;

unsigned int
gt4_word_dict_get_type (void)
{
  if (!dict_type) {
    dict_class = (GT4WordDictClass *) az_register_interface_type (&dict_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4WordDict",
      sizeof (GT4WordDictClass), sizeof (GT4WordDictImplementation), sizeof (GT4WordDictInstance),
      (void (*) (AZClass *)) dict_class_init,
      NULL, NULL, NULL);
  }
  return dict_type;
}

static void
dict_class_init (GT4WordDictClass *klass)
{
  klass->interface_class.klass.flags = AZ_CLASS_ZERO_MEMORY;
}

unsigned int
gt4_word_dict_lookup (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  return impl->lookup (impl, inst, word);
}

