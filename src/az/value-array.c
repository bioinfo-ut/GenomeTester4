#define __AZ_VALUE_ARRAY_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include <az/class.h>
#include <az/value.h>
#include <az/value-array.h>

static void value_array_class_init (AZValueArrayClass *klass);
static void value_array_finalize (AZValueArrayClass *klass, AZValueArray *varray);

/* AZInstance implementation */
static unsigned int value_array_get_property (AZClass *klass, void *instance, unsigned int idx, AZValue *value);
/* AZCollection implementation */
static unsigned int value_array_get_element_type (AZCollectionImplementation *collection_impl, void *collection_inst);
static unsigned int value_array_get_size (AZCollectionImplementation *collection_impl, void *collection_inst);
/* AZList implementation */
static unsigned int value_array_get_element (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value);

enum {
	PROP_LENGTH,
	NUM_PROPERTIES
};

unsigned int
az_value_array_get_type (void)
{
	static unsigned int type = 0;
	if (!type) {
		AZClass *klass;
		klass = az_register_type (&type, AZ_TYPE_REFERENCE, (const unsigned char *) "AZValueArray",
			sizeof (AZValueArrayClass), sizeof (AZValueArray),
			(void (*) (AZClass *)) value_array_class_init,
			NULL,
			(void (*) (AZImplementation *, void *)) value_array_finalize);
		klass->flags |= AZ_CLASS_ZERO_MEMORY;
	}
	return type;
}

static void
value_array_class_init (AZValueArrayClass *klass)
{
	az_class_set_num_interfaces ((AZClass *) klass, 1);
	az_class_declare_interface ((AZClass *) klass, 0, AZ_TYPE_LIST, ARIKKEI_OFFSET (AZValueArrayClass, list_implementation), 0);
	az_class_set_num_properties ((AZClass *) klass, NUM_PROPERTIES);
	az_class_property_setup ((AZClass *) klass, PROP_LENGTH, (const unsigned char *) "length", AZ_TYPE_UINT32, 0, 1, 0, 1, 0, AZ_TYPE_NONE, NULL);
	((AZClass *) klass)->get_property = value_array_get_property;
	klass->list_implementation.collection_implementation.get_element_type = value_array_get_element_type;
	klass->list_implementation.collection_implementation.get_size = value_array_get_size;
	klass->list_implementation.get_element = value_array_get_element;
}

static void
value_array_finalize (AZValueArrayClass *klass, AZValueArray *varray)
{
	unsigned int i;
	for (i = 0; i < varray->length; i++) {
		az_value_clear (&varray->values[i]);
	}
	if (varray->values) free (varray->values);
}

static unsigned int
value_array_get_property (AZClass *klass, void *inst, unsigned int idx, AZValue *value)
{
	AZValueArray *varray = (AZValueArray *) inst;
	switch (idx) {
	case PROP_LENGTH:
		az_value_set_unsigned_int (value, AZ_TYPE_UINT32, varray->length);
		return 1;
	default:
		break;
	}
	return 0;
}

static unsigned int
value_array_get_element_type (AZCollectionImplementation *collection_impl, void *collection_inst)
{
	return AZ_TYPE_ANY;
}

static unsigned int
value_array_get_size (AZCollectionImplementation *collection_impl, void *collection_inst)
{
	AZValueArray *varray = (AZValueArray *) collection_inst;
	return varray->length;
}

static unsigned int
value_array_get_element (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value)
{
	AZValueArray *varray;
	arikkei_return_val_if_fail (list_impl != NULL, 0);
	arikkei_return_val_if_fail (list_inst != NULL, 0);
	arikkei_return_val_if_fail (value != NULL, 0);
	varray = (AZValueArray *) list_inst;
	arikkei_return_val_if_fail (index < varray->length, 0);
	az_value_copy (value, &varray->values[index]);
	return 1;
}

AZValueArray *
az_value_array_new (unsigned int length)
{
	AZValueArray *varray = (AZValueArray *) az_instance_new (AZ_TYPE_VALUE_ARRAY);
	varray->length = length;
	varray->values = (AZValue *) az_instance_new_array (AZ_TYPE_VALUE, varray->length);
	return varray;
}

void
az_value_array_set_element (AZValueArray *varray, unsigned int idx, const AZValue *value)
{
	arikkei_return_if_fail (varray != NULL);
	arikkei_return_if_fail (idx < varray->length);
	arikkei_return_if_fail (value != NULL);
	az_value_copy (&varray->values[idx], value);
}
