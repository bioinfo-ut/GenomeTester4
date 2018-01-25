#ifndef __AZ_ARRAY_OF_H__
#define __AZ_ARRAY_OF_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

typedef struct _AZArrayOf AZArrayOf;
typedef struct _AZArrayOfClass AZArrayOfClass;

#define AZ_TYPE_ARRAY_OF(t) az_array_of_get_type (t)

#include <az/class.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZArrayOf {
	unsigned int size;
};

struct _AZArrayOfClass {
	AZClass klass;
	unsigned int element_type;
};

unsigned int az_array_of_get_type (unsigned int element_type);

unsigned int az_array_of_get_instance_size (unsigned int element_type, unsigned int size);
void *az_array_of_get_element_instance (AZArrayOf *array_of, unsigned int element_type, unsigned int idx);

void az_array_of_init (AZArrayOf *array_of, unsigned int element_type, unsigned int size);
void az_array_of_finalize (AZArrayOf *array_of, unsigned int element_type);
AZArrayOf *az_array_of_new (unsigned int element_type, unsigned int size);
void az_array_of_free (AZArrayOf *array_of, unsigned int element_type);

#ifdef __cplusplus
};
#endif

#endif
