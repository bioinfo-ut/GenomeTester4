#define __GT4_SEQUENCE_COLLECTION_C__

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libarikkei/arikkei-utils.h>

#include "sequence-collection.h"

void collection_instance_init (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst);

static unsigned int sequence_collection_type = 0;
static GT4SequenceCollectionClass *sequence_collection_class = NULL;

unsigned int
gt4_sequence_collection_get_type (void)
{
  if (!sequence_collection_type) {
    sequence_collection_class = (GT4SequenceCollectionClass *) az_register_interface_type (&sequence_collection_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4SequenceCollection",
    sizeof (GT4SequenceCollectionClass), sizeof (GT4SequenceCollectionImplementation), sizeof (GT4SequenceCollectionInstance),
    NULL, NULL,
    (void (*) (AZImplementation *, void *)) collection_instance_init,
    NULL);
  }
  return sequence_collection_type;
}

void
collection_instance_init (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst)
{
  inst->writable = 0;
  inst->n_subseqs = 0;
  memset (&inst->subseq, 0, sizeof (GT4SubSequence));
}

unsigned int
gt4_sequence_collection_get_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (idx < inst->n_subseqs, 0);
  return impl->get_subsequence (impl, inst, idx);
}

int
gt4_sequence_collection_add_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned long long name_pos, unsigned int name_len)
{
  arikkei_return_val_if_fail (impl != NULL, GT4_SEQUENCE_COLLECTION_INVALID_INTERFACE);
  arikkei_return_val_if_fail (inst != NULL, GT4_SEQUENCE_COLLECTION_INVALID_INTERFACE);
  arikkei_return_val_if_fail (inst->writable, -1);
  return impl->add_subsequence (impl, inst, name_pos, name_len);
}

unsigned int
gt4_sequence_collection_set_subsequence (GT4SequenceCollectionImplementation *impl, GT4SequenceCollectionInstance *inst, unsigned int idx, unsigned int sequence_pos, unsigned int sequence_len)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (idx < inst->n_subseqs, 0);
  return impl->set_subsequence (impl, inst, idx, sequence_pos, sequence_len);
}

