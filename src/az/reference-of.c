#define __AZ_REFERENCE_OF_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <az/reference-of.h>

static void reference_of_class_init (AZReferenceOfClass *klass, AZClass *instance_class);
static void reference_of_instance_init (AZReferenceOfClass *klass, AZReferenceOf *reference_of);
static void reference_of_instance_finalize (AZReferenceOfClass *klass, AZReferenceOf *reference_of);

unsigned int
az_reference_of_get_type (unsigned int instance_type)
{
	static unsigned int num_subtypes = 0;
	static unsigned int *subtypes = NULL;
	if (instance_type >= num_subtypes) {
		unsigned int new_size = (instance_type + 1 + 255) & 0xffffff00;
		subtypes = realloc (subtypes, num_subtypes * sizeof (unsigned int));
		memset (&subtypes[num_subtypes], 0, (new_size - num_subtypes) * sizeof (unsigned int));
		num_subtypes = new_size;
	}
	if (!subtypes[instance_type]) {
		AZClass *instance_class = az_type_get_class (instance_type);
		unsigned int len = (unsigned int) strlen ((const char *) instance_class->name);
		unsigned char *name = malloc (len + 12);
		sprintf ((char *) name, "ReferenceOf%s", instance_class->name);
		az_register_composite_type (&subtypes[instance_type], AZ_TYPE_REFERENCE, name,
			sizeof (AZReferenceOfClass), sizeof (AZReferenceOf) - 1 + instance_class->instance_size,
			(void (*) (AZClass *, void *)) reference_of_class_init,
			(void (*) (AZImplementation *, void *)) reference_of_instance_init,
			(void (*) (AZImplementation *, void *)) reference_of_instance_finalize,
			instance_class);
		free (name);
	}
	return subtypes[instance_type];
}

static void
reference_of_class_init (AZReferenceOfClass *klass, AZClass *instance_class)
{
	klass->instance_type = instance_class->implementation.type;
}

static void
reference_of_instance_init (AZReferenceOfClass *klass, AZReferenceOf *reference_of)
{
	az_instance_init (AZ_REFERENCE_OF_INSTANCE (reference_of), klass->instance_type);
}

static void
reference_of_instance_finalize (AZReferenceOfClass *klass, AZReferenceOf *reference_of)
{
	az_instance_finalize (AZ_REFERENCE_OF_INSTANCE (reference_of), klass->instance_type);
}

