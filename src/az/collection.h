#ifndef __AZ_COLLECTION_H__
#define __AZ_COLLECTION_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_COLLECTION (az_collection_get_type ())

typedef struct _AZCollectionImplementation AZCollectionImplementation;
typedef struct _AZCollectionClass AZCollectionClass;

#include <az/interface.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZCollectionImplementation {
	AZImplementation implementation;
	unsigned int (*get_element_type) (AZCollectionImplementation *collection_impl, void *collection_inst);
	unsigned int (*get_size) (AZCollectionImplementation *collection_impl, void *collection_inst);
	unsigned int (*get_iterator) (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
	unsigned int (*iterator_next) (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
	unsigned int (*get_element) (AZCollectionImplementation *collection_impl, void *collection_inst, const AZValue *iterator, AZValue *value);
};

struct _AZCollectionClass {
	AZInterfaceClass interface_class;
};

unsigned int az_collection_get_type (void);

/* Returns base type that all elements are guaranteed to be assignable to */
unsigned int az_collection_get_element_type (AZCollectionImplementation *collection_impl, void *collection_inst);
unsigned int az_collection_get_size (AZCollectionImplementation *collection_impl, void *collection_inst);
/* Return 0 if unusuccessful */
unsigned int az_collection_get_iterator (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
unsigned int az_collection_iterator_next (AZCollectionImplementation *collection_impl, void *collection_inst, AZValue *iterator);
/* Value will contain the proper type of given element */
unsigned int az_collection_get_element (AZCollectionImplementation *collection_impl, void *collection_inst, const AZValue *iterator, AZValue *value);

#ifdef __cplusplus
};
#endif

#endif

