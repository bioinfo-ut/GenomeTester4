#define __AZ_LIST_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <az/interface.h>
#include <az/value.h>

#include <az/list.h>

/* AZInterface implementation */
static void list_implementation_init (AZListImplementation *impl);
/* AZCollection implementation */
static unsigned int list_get_iterator (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
static unsigned int list_iterator_next (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
static unsigned int list_get_element (AZCollectionImplementation *collection_impl, void *collection_inst, const AZValue *iterator, AZValue *value);

static unsigned int list_type = 0;
static AZListClass *list_class;

unsigned int
az_list_get_type (void)
{
	if (!list_type) {
		list_class = (AZListClass *) az_register_interface_type (&list_type, AZ_TYPE_COLLECTION, (const unsigned char *) "AZList",
			sizeof (AZListClass), sizeof (AZListImplementation), 0,
			NULL,
			(void (*) (AZImplementation *)) list_implementation_init,
			NULL, NULL);
		((AZClass *) list_class)->flags |= AZ_CLASS_ZERO_MEMORY;
	}
	return list_type;
}

static void
list_implementation_init (AZListImplementation *impl)
{
	AZCollectionImplementation *collection = (AZCollectionImplementation *) impl;
	collection->get_iterator = list_get_iterator;
	collection->iterator_next = list_iterator_next;
	collection->get_element = list_get_element;
}

static unsigned int
list_get_iterator (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator)
{
	az_value_set_unsigned_int (iterator, AZ_TYPE_UINT32, 0);
	return 1;
}

static unsigned int
list_iterator_next (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator)
{
	if (iterator->uvalue >= az_collection_get_size (collection_impl, collection_inst)) return 0;
	iterator->uvalue += 1;
	return 1;
}

static unsigned int
list_get_element (AZCollectionImplementation *collection_impl, void *collection_inst, const AZValue *iterator, AZValue *value)
{
	return az_list_get_element ((AZListImplementation *) collection_impl, collection_inst, iterator->uvalue, value);
}

unsigned int
az_list_get_element (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value)
{
	if (list_impl->get_element) {
		return list_impl->get_element (list_impl, list_inst, index, value);
	}
	return 0;
}
