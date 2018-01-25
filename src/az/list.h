#ifndef __AZ_LIST_H__
#define __AZ_LIST_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_LIST (az_list_get_type ())

typedef struct _AZListImplementation AZListImplementation;
typedef struct _AZListClass AZListClass;

#include <az/collection.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
* An hybrid interface that allows sequential access by index
* The iterator of a list is always unsigned integer
*/

struct _AZListImplementation {
	AZCollectionImplementation collection_implementation;
	unsigned int (*get_element) (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value);
};

struct _AZListClass {
	AZCollectionClass collection_class;
};

unsigned int az_list_get_type (void);

unsigned int az_list_get_element (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value);

#ifdef __cplusplus
};
#endif

#endif

