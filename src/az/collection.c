#define __AZ_COLLECTION_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <az/collection.h>

static unsigned int collection_type = 0;
static AZCollectionClass *collection_class;

unsigned int
az_collection_get_type (void)
{
	if (!collection_type) {
		collection_class = (AZCollectionClass *) az_register_interface_type (&collection_type, AZ_TYPE_INTERFACE, (const unsigned char *) "AZCollection",
			sizeof (AZCollectionClass), sizeof (AZCollectionImplementation), 0,
			NULL,
			NULL,
			NULL, NULL);
	}
	return collection_type;
}

unsigned int
az_collection_get_element_type (AZCollectionImplementation *impl, void *collection_instance)
{
	if (impl->get_element_type) {
		return impl->get_element_type (impl, collection_instance);
	}
	return AZ_TYPE_NONE;
}

unsigned int
az_collection_get_size (AZCollectionImplementation *impl, void *collection_instance)
{
	if (impl->get_size) {
		return impl->get_size (impl, collection_instance);
	}
	return 0;
}

unsigned int
az_collection_get_iterator (AZCollectionImplementation *impl, void *collection_instance, AZValue *iterator)
{
	if (impl->get_iterator) {
		return impl->get_iterator (impl, collection_instance, iterator);
	}
	return 0;
}

unsigned int
az_collection_iterator_next (AZCollectionImplementation *impl, void *collection_instance, AZValue *iterator)
{
	if (impl->iterator_next) {
		return impl->iterator_next (impl, collection_instance, iterator);
	}
	return 0;
}

unsigned int
az_collection_get_element (AZCollectionImplementation *impl, void *collection_instance, const AZValue *iterator, AZValue *value)
{
	if (impl->get_element) {
		return impl->get_element (impl, collection_instance, iterator, value);
	}
	return 0;
}


