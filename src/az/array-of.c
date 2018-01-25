#define __AZ_ARRAY_OF_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <az/array-of.h>

static void array_of_class_init (AZArrayOfClass *klass, AZClass *element_class);
static void array_of_instance_init (AZArrayOfClass *klass, AZArrayOf *array_of);
static void array_of_instance_finalize (AZArrayOfClass *klass, AZArrayOf *array_of);

unsigned int
az_array_of_get_type (unsigned int element_type)
{
	static unsigned int num_subtypes = 0;
	static unsigned int *subtypes = NULL;
	if (element_type >= num_subtypes) {
		unsigned int new_size = (element_type + 1 + 255) & 0xffffff00;
		subtypes = realloc (subtypes, num_subtypes * sizeof (unsigned int));
		memset (&subtypes[num_subtypes], 0, (new_size - num_subtypes) * sizeof (unsigned int));
		num_subtypes = new_size;
	}
	if (!subtypes[element_type]) {
		AZClass *element_class = az_type_get_class (element_type);
		unsigned int len = (unsigned int) strlen ((const char *) element_class->name);
		unsigned char *name = malloc (len + 8);
		sprintf ((char *) name, "ArrayOf%s", element_class->name);
		az_register_composite_type (&subtypes[element_type], AZ_TYPE_STRUCT, name,
			sizeof (AZArrayOfClass), sizeof (AZArrayOf),
			(void (*) (AZClass *, void *)) array_of_class_init,
			(void (*) (AZImplementation *, void *)) array_of_instance_init,
			(void (*) (AZImplementation *, void *)) array_of_instance_finalize,
			element_class);
		free (name);
	}
	return subtypes[element_type];
}

static void
array_of_class_init (AZArrayOfClass *klass, AZClass *element_class)
{
	klass->element_type = element_class->implementation.type;
}

static void
array_of_instance_init (AZArrayOfClass *klass, AZArrayOf *array_of)
{
}

static void
array_of_instance_finalize (AZArrayOfClass *klass, AZArrayOf *array_of)
{
	unsigned int i;
	for (i = 0; i < array_of->size; i++) {
		void *instance = az_array_of_get_element_instance (array_of, klass->element_type, i);
		az_instance_finalize (instance, klass->element_type);
	}
}

unsigned int az_array_of_get_instance_size (unsigned int element_type, unsigned int size);
void *az_array_of_get_element_instance (AZArrayOf *array_of, unsigned int element_type, unsigned int idx);

void az_array_of_init (AZArrayOf *array_of, unsigned int element_type, unsigned int size);
void az_array_of_finalize (AZArrayOf *array_of, unsigned int element_type);
AZArrayOf *az_array_of_new (unsigned int element_type, unsigned int size);
void az_array_of_free (AZArrayOf *array_of, unsigned int element_type);

