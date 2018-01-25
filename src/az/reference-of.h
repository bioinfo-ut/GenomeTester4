#ifndef __AZ_REFERENCE_OF_H__
#define __AZ_REFERENCE_OF_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2017
*/

typedef struct _AZReferenceOf AZReferenceOf;
typedef struct _AZReferenceOfClass AZReferenceOfClass;

#define AZ_TYPE_REFERENCE_OF(t) az_reference_of_get_type (t)

#define AZ_REFERENCE_OF_INSTANCE(r) ((void *) ((char *)(r) + sizeof (AZReferenceOf)))

#include <az/reference.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZReferenceOf {
	AZReference reference;
};

struct _AZReferenceOfClass {
	AZReferenceClass reference_klass;
	unsigned int instance_type;
};

unsigned int az_reference_of_get_type (unsigned int instance_type);

#ifdef __cplusplus
};
#endif

#endif
