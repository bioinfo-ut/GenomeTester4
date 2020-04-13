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

#include "sequence.h"
#include "word-table.h"

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
gt4_word_dict_lookup (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word, unsigned int canonize)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  if (canonize) {
    unsigned long long rev = get_reverse_complement (word, inst->word_length);
    if (rev < word) word = rev;
  }
  return impl->lookup (impl, inst, word);
}

unsigned int
gt4_word_dict_lookup_mm (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word, unsigned int n_mm, unsigned int pm_3, unsigned int canonize, int print_all, unsigned int equal_mm_only)
{
  static GT4WordTable mm_table = {0};
  unsigned long long i;
  unsigned int count = 0;

  /* If no mismatches (last resursion) */
  if (!n_mm) {
    return gt4_word_dict_lookup (impl, inst, word, canonize);
  }
  if (!mm_table.n_word_slots) {
    gt4_word_table_setup (&mm_table, inst->word_length, 256, 0);
  }
  gt4_word_table_generate_mismatches (&mm_table, word, NULL, n_mm, pm_3, canonize, 0, equal_mm_only);
  for (i = 0; i < mm_table.n_words; i++) {
    if (gt4_word_dict_lookup (impl, inst, mm_table.words[i], 0)) {
      /* Found it */
      count += inst->value;
      if (print_all) {
        fprintf (stdout, "%s\t%u\n", word_to_string (mm_table.words[i], mm_table.wordlength), inst->value);
      }
    }
  }
  gt4_word_table_clear (&mm_table);
  inst->value = count;
  return count != 0;
}
