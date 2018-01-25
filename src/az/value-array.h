#ifndef __AZ_VALUE_ARRAY_H__
#define __AZ_VALUE_ARRAY_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_VALUE_ARRAY az_value_array_get_type ()

typedef struct _AZValueArray AZValueArray;
typedef struct _AZValueArrayClass AZValueArrayClass;

#include <az/list.h>
#include <az/reference.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZValueArray {
	AZReference reference;
	unsigned int length;
	AZValue *values;
};

struct _AZValueArrayClass {
	AZReferenceClass reference_klass;
	AZListImplementation list_implementation;
};

unsigned int az_value_array_get_type (void);

AZValueArray *az_value_array_new (unsigned int length);

void az_value_array_set_element (AZValueArray *varray, unsigned int idx, const AZValue *value);

#ifdef __cplusplus
};
#endif

#endif
