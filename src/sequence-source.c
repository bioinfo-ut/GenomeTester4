#define __GT4_SEQUENCE_SOURCE_C__

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

#include <libarikkei/arikkei-utils.h>

#include "sequence-source.h"

void source_instance_init (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);

static unsigned int sequence_source_type = 0;
static GT4SequenceSourceClass *sequence_source_class = NULL;

unsigned int
gt4_sequence_source_get_type (void)
{
  if (!sequence_source_type) {
    sequence_source_class = (GT4SequenceSourceClass *) az_register_interface_type (&sequence_source_type, AZ_TYPE_INTERFACE, (const unsigned char *) "GT4SequenceSource",
    sizeof (GT4SequenceSourceClass), sizeof (GT4SequenceSourceImplementation), sizeof (GT4SequenceSourceInstance),
    NULL, NULL,
    (void (*) (AZImplementation *, void *)) source_instance_init,
    NULL);
  }
  return sequence_source_type;
}

void
source_instance_init (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  inst->position = 0;
  inst->eof = 0;
}

unsigned int
gt4_sequence_source_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (!inst->open, 0);
  if (inst->error) return 0;
  if (impl->open && !impl->open (impl, inst)) {
    inst->error = 1;
    return 0;
  }
  inst->open = 1;
  inst->eof = 0;
  inst->position = 0;
  return 1;
}

int
gt4_sequence_source_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  int val;
  arikkei_return_val_if_fail (impl != NULL, GT4_SEQUENCE_SOURCE_INVALID_INTERFACE);
  arikkei_return_val_if_fail (inst != NULL, GT4_SEQUENCE_SOURCE_INVALID_INTERFACE);
  if (inst->error) return -1;
  if (inst->eof) return 0;
  val = impl->read (impl, inst);
  if (val < 0) {
    inst->error = 1;
  } else if (val == 0) {
    inst->eof = 1;
  }
  return val;
}

unsigned int
gt4_sequence_source_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  arikkei_return_val_if_fail (impl != NULL, 0);
  arikkei_return_val_if_fail (inst != NULL, 0);
  arikkei_return_val_if_fail (inst->open, 0);
  if (inst->error) return 0;
  if (impl->close && !impl->close (impl, inst)) {
    inst->error = 1;
    return 0;
  }
  inst->open = 0;
  return 1;
}

