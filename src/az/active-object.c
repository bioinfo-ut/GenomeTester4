#define __AZ_ACTIVE_OBJECT_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <malloc.h>
#include <stdio.h>
#include <string.h>

#include <az/string.h>

#include <az/active-object.h>

static void az_active_object_class_init (AZActiveObjectClass *klass);

/* AZObject implementation */
static void az_active_object_shutdown (AZObject *object);

static AZObjectClass *parent_class;

unsigned int
az_active_object_get_type (void)
{
	static unsigned int type = 0;
	if (!type) {
		az_register_type (&type, AZ_TYPE_OBJECT,
			(const unsigned char *) "AZActiveObject",
			sizeof (AZActiveObjectClass),
			sizeof (AZActiveObject),
			(void (*) (AZClass *)) az_active_object_class_init, NULL, NULL);
	}
	return type;
}

static void
az_active_object_class_init (AZActiveObjectClass *klass)
{
	parent_class = (AZObjectClass *) ((AZClass *) klass)->parent;
	/* AZObject implementation */
	((AZObjectClass *) klass)->shutdown = az_active_object_shutdown;
}

static void
az_active_object_shutdown (AZObject *object)
{
	AZActiveObject *aobj = (AZActiveObject *) object;
	if (aobj->callbacks) {
		unsigned int i;
		for (i = 0; i < aobj->callbacks->length; i++) {
			AZObjectListener *listener;
			listener = aobj->callbacks->listeners + i;
			if (listener->vector->dispose) listener->vector->dispose (aobj, listener->data);
		}
		free (aobj->callbacks);
		aobj->callbacks = NULL;
	}
	if (aobj->attributes) {
		unsigned int i;
		for (i = 0; i < aobj->attributes->length; i++) {
			free (aobj->attributes->attribs[i].key);
			az_value_clear (&aobj->attributes->attribs[i].value);
		}
		free (aobj->attributes);
		aobj->attributes = NULL;
	}
}

static AZObjectAttribute *
az_active_object_get_attribute_slot (AZActiveObject *aobj, const unsigned char *key, unsigned int create)
{
	if (aobj->attributes) {
		unsigned int i;
		for (i = 0; i < aobj->attributes->length; i++) {
			if (!strcmp ((const char *) aobj->attributes->attribs[i].key->str, (const char *) key)) {
				return &aobj->attributes->attribs[i];
			}
		}
	}
	if (!create) return NULL;
	if (!aobj->attributes) {
		aobj->attributes = (AZObjectAttributeArray *) malloc (sizeof (AZObjectAttributeArray) + 3 * sizeof (AZObjectAttribute));
		memset (aobj->attributes, 0, sizeof (AZObjectAttributeArray) + 3 * sizeof (AZObjectAttribute));
		aobj->attributes->size = 4;
		aobj->attributes->length = 0;
	} else if (aobj->attributes->length >= aobj->attributes->size) {
		aobj->attributes->size <<= 1;
		aobj->attributes = (AZObjectAttributeArray *) realloc (aobj->attributes, sizeof (AZObjectAttributeArray) + (aobj->attributes->size - 1) * sizeof (AZObjectAttribute));
		memset (&aobj->attributes->attribs[aobj->attributes->length], 0, (aobj->attributes->size - aobj->attributes->length) * sizeof (AZObjectAttribute));
	}
	aobj->attributes->attribs[aobj->attributes->length].key = az_string_new (key);
	return &aobj->attributes->attribs[aobj->attributes->length++];
}

unsigned int
az_active_object_get_attribute (AZActiveObject *aobj, const unsigned char *key, AZValue *val)
{
	AZObjectAttribute *attr;
	arikkei_return_val_if_fail (AZ_IS_ACTIVE_OBJECT (aobj), 0);
	arikkei_return_val_if_fail (key != NULL, 0);
	arikkei_return_val_if_fail (val != NULL, 0);
	attr = az_active_object_get_attribute_slot (aobj, key, 0);
	if (attr) {
		az_value_copy (val, &attr->value);
		return 1;
	}
	return 0;
}

unsigned int
az_active_object_set_attribute (AZActiveObject *aobj, const unsigned char *key, const AZValue *val)
{
	AZObjectAttribute *attr;
	arikkei_return_val_if_fail (AZ_IS_ACTIVE_OBJECT (aobj), 0);
	arikkei_return_val_if_fail (key != NULL, 0);
	if (!val) return az_active_object_clear_attribute (aobj, key);
	attr = az_active_object_get_attribute_slot (aobj, key, 1);
	az_value_copy (&attr->value, val);
	return 1;
}

unsigned int
az_active_object_clear_attribute (AZActiveObject *aobj, const unsigned char *key)
{
	unsigned int i;
	arikkei_return_val_if_fail (AZ_IS_ACTIVE_OBJECT (aobj), 0);
	arikkei_return_val_if_fail (key != NULL, 0);
	if (!aobj->attributes) return 0;
	for (i = 0; i < aobj->attributes->length; i++) {
		if (!strcmp ((const char *) aobj->attributes->attribs[i].key->str, (const char *) key)) {
			az_string_unref (aobj->attributes->attribs[i].key);
			az_value_clear (&aobj->attributes->attribs[i].value);
			if ((i + 1) < aobj->attributes->length) {
				memcpy (&aobj->attributes->attribs[i], &aobj->attributes->attribs[i + 1], (aobj->attributes->length - (i + 1)) * sizeof (AZObjectAttribute));
			}
			aobj->attributes->length -= 1;
			if (!aobj->attributes->length) {
				free (aobj->attributes);
				aobj->attributes = NULL;
			}
			return 1;
		}
	}
	return 0;
}

void
az_active_object_add_listener (AZActiveObject *aobj, const AZObjectEventVector *vector, unsigned int size, void *data)
{
	AZObjectListener *listener;

	if (!aobj->callbacks) {
		aobj->callbacks = (AZObjectCallbackBlock *) malloc (sizeof (AZObjectCallbackBlock));
		aobj->callbacks->size = 1;
		aobj->callbacks->length = 0;
	}
	if (aobj->callbacks->length >= aobj->callbacks->size) {
		int newsize;
		newsize = aobj->callbacks->size << 1;
		aobj->callbacks = (AZObjectCallbackBlock *) realloc (aobj->callbacks, sizeof (AZObjectCallbackBlock) + (newsize - 1) * sizeof (AZObjectListener));
		aobj->callbacks->size = newsize;
	}
	listener = aobj->callbacks->listeners + aobj->callbacks->length;
	listener->vector = vector;
	listener->size = size;
	listener->data = data;
	aobj->callbacks->length += 1;
}

void
az_active_object_remove_listener_by_data (AZActiveObject *aobj, void *data)
{
	if (aobj->callbacks) {
		unsigned int i;
		for (i = 0; i < aobj->callbacks->length; i++) {
			AZObjectListener *listener;
			listener = aobj->callbacks->listeners + i;
			if (listener->data == data) {
				aobj->callbacks->length -= 1;
				if (aobj->callbacks->length < 1) {
					free (aobj->callbacks);
					aobj->callbacks = NULL;
				} else if (aobj->callbacks->length != i) {
					*listener = aobj->callbacks->listeners[aobj->callbacks->length];
				}
				return;
			}
		}
	}
}

unsigned int
az_active_object_set_attribute_i32 (AZActiveObject *aobj, const unsigned char *key, int value)
{
	unsigned int result;
	AZValue val = { 0 };
	az_value_set_int (&val, AZ_TYPE_INT32, value);
	result = az_active_object_set_attribute (aobj, key, &val);
	az_value_clear (&val);
	return result;
}

unsigned int
az_active_object_set_attribute_object (AZActiveObject *aobj, const unsigned char *key, AZObject *value)
{
	unsigned int result;
	AZValue val = { 0 };
	az_value_set_object (&val, value);
	result = az_active_object_set_attribute (aobj, key, &val);
	az_value_clear (&val);
	return result;
}


