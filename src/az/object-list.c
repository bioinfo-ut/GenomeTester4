#define __AZ_OBJECT_LIST_CPP__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>
#include <string.h>

#include <az/active-object.h>
#include <az/class.h>
#include <az/object-list.h>

static void object_list_class_init (AZObjectListClass *klass);
static void object_list_init (AZObjectListClass *klass, AZObjectList *objl);
static void object_list_finalize (AZObjectListClass *klass, AZObjectList *objl);

/* AZClass implementation */
static unsigned int object_list_get_property (AZClass *klass, void *instance, unsigned int idx, AZValue *value);
/* AZCollection implementation */
static unsigned int object_list_get_element_type (AZCollectionImplementation *impl, void *collection_instance);
static unsigned int object_list_get_size (AZCollectionImplementation *impl, void *collection_instance);
/* AZList implementation */
static unsigned int object_list_get_element (AZListImplementation *impl, void *array_instance, unsigned int index, AZValue *value);

/* Method implementations */
static unsigned int object_list_call_lookUp (AZValue *thisval, AZValue *retval, AZValue *args);

static unsigned int object_list_type = 0;

unsigned int
az_object_list_get_type (void)
{
	if (!object_list_type) {
		az_register_type (&object_list_type, AZ_TYPE_STRUCT, (const unsigned char *) "AZObjectList",
			sizeof (AZObjectListClass), sizeof (AZObjectList),
			(void (*) (AZClass *)) object_list_class_init,
			(void (*) (AZImplementation *, void *)) object_list_init,
			(void (*) (AZImplementation *, void *)) object_list_finalize);
	}
	return object_list_type;
}

static void
object_list_class_init (AZObjectListClass *klass)
{
	((AZClass *) klass)->flags |= AZ_CLASS_ZERO_MEMORY;
	/* Interfaces */
	az_class_set_num_interfaces ((AZClass *) klass, 1);
	az_class_declare_interface ((AZClass *) klass, 0, AZ_TYPE_LIST, AZ_OFFSET (AZObjectListClass, list_implementation), 0);
	/* Array implementation */
	klass->list_implementation.collection_implementation.get_element_type = object_list_get_element_type;
	klass->list_implementation.collection_implementation.get_size = object_list_get_size;
	klass->list_implementation.get_element = object_list_get_element;
}

static void
object_list_init (AZObjectListClass *klass, AZObjectList *objl)
{
	objl->size = 16;
	objl->objects = (AZObject **) malloc (objl->size * sizeof (AZObject *));
}

static void
object_list_finalize (AZObjectListClass *klass, AZObjectList *objl)
{
	if (objl->weak) {
		for (unsigned int i = 0; i < objl->length; i++) {
			az_active_object_remove_listener_by_data (AZ_ACTIVE_OBJECT (objl->objects[i]), objl);
		}
	} else {
		for (unsigned int i = 0; i < objl->length; i++) {
			az_object_unref (objl->objects[i]);
		}
	}
	free (objl->objects);
}

static unsigned int
object_list_get_element_type (AZCollectionImplementation *collection_impl, void *collection_inst)
{
	AZObjectList *objl = (AZObjectList *) collection_inst;
	return objl->type;
}

static unsigned int
object_list_get_size (AZCollectionImplementation *collection_impl, void *collection_inst)
{
	AZObjectList *objl = (AZObjectList *) collection_inst;
	return objl->length;
}

static unsigned int
object_list_get_element (AZListImplementation *list_impl, void *list_inst, unsigned int index, AZValue *value)
{
	AZObjectList *objl = (AZObjectList *) list_inst;
	arikkei_return_val_if_fail (index < objl->length, 0);
	az_value_set_object (value, objl->objects[index]);
	return 1;
}

void
az_object_list_setup (AZObjectList *objl, unsigned int type, unsigned int weak)
{
	arikkei_return_if_fail (objl != NULL);
	arikkei_return_if_fail (az_type_is_a (type, AZ_TYPE_INTERFACE) || (weak && az_type_is_a (type, AZ_TYPE_ACTIVE_OBJECT)) || az_type_is_a (type, AZ_TYPE_OBJECT));
	az_instance_init (objl, AZ_TYPE_OBJECT_LIST);
	objl->type = type;
	objl->interface = az_type_is_a (type, AZ_TYPE_INTERFACE);
	objl->weak = weak;
}

void
az_object_list_release (AZObjectList *objl)
{
	arikkei_return_if_fail (objl != NULL);
	az_instance_finalize (objl, AZ_TYPE_OBJECT_LIST);
}

AZObjectList *
az_object_list_new (unsigned int type, unsigned int weak)
{
	arikkei_return_val_if_fail (az_type_is_a (type, AZ_TYPE_INTERFACE) || (weak && az_type_is_a (type, AZ_TYPE_ACTIVE_OBJECT)) || az_type_is_a (type, AZ_TYPE_OBJECT), NULL);
	AZObjectList *objl = (AZObjectList *) malloc (sizeof (AZObjectList));
	az_object_list_setup (objl, type, weak);
	return objl;
}

void
az_object_list_delete (AZObjectList *objl)
{
	arikkei_return_if_fail (objl != NULL);
	az_object_list_release (objl);
	free (objl);
}

static void object_list_object_dispose (AZActiveObject *object, void *data);
static unsigned int object_list_remove_object_internal (AZObjectList *objl, AZObject *object);

AZObjectEventVector object_list_event_vector = {
	object_list_object_dispose
};

void
az_object_list_append_object (AZObjectList *objl, AZObject *obj)
{
	arikkei_return_if_fail (objl != NULL);
	arikkei_return_if_fail ((objl->weak && AZ_IS_ACTIVE_OBJECT (obj)) || AZ_IS_OBJECT (obj));
	arikkei_return_if_fail ((objl->interface && az_object_implements (obj, objl->type)) || az_object_is_a (obj, objl->type));
	if (objl->length >= objl->size) {
		objl->size = objl->size << 1;
		objl->objects = (AZObject **) realloc (objl->objects, objl->size * sizeof (AZObject *));
	}
	objl->objects[objl->length++] = obj;
	if (objl->weak) {
		az_active_object_add_listener (AZ_ACTIVE_OBJECT (obj), &object_list_event_vector, sizeof (AZObjectEventVector), objl);
	} else {
		az_object_ref (obj);
	}
}

void
az_object_list_remove_object (AZObjectList *objl, AZObject *obj)
{
	arikkei_return_if_fail (objl != NULL);
	arikkei_return_if_fail ((objl->weak && AZ_IS_ACTIVE_OBJECT (obj)) || AZ_IS_OBJECT (obj));
	arikkei_return_if_fail ((objl->interface && az_object_implements (obj, objl->type)) || az_object_is_a (obj, objl->type));
	if (object_list_remove_object_internal (objl, obj)) {
		if (objl->weak) {
			az_active_object_remove_listener_by_data (AZ_ACTIVE_OBJECT (obj), objl);
		} else {
			az_object_unref (obj);
		}
	}
}

void
az_object_list_remove_object_by_index (AZObjectList *objl, unsigned int idx)
{
	AZObject *obj;
	arikkei_return_if_fail (objl != NULL);
	arikkei_return_if_fail (idx < objl->length);
	obj = objl->objects[idx];
	if (object_list_remove_object_internal (objl, obj)) {
		if (objl->weak) {
			az_active_object_remove_listener_by_data (AZ_ACTIVE_OBJECT (obj), objl);
		} else {
			az_object_unref (obj);
		}
	}
}

void
az_object_list_clear (AZObjectList *objl)
{
	unsigned int i;
	arikkei_return_if_fail (objl != NULL);
	for (i = 0; i < objl->length; i++) {
		if (objl->weak) {
			az_active_object_remove_listener_by_data (AZ_ACTIVE_OBJECT (objl->objects[i]), objl);
		} else {
			az_object_unref (objl->objects[i]);
		}
	}
	objl->length = 0;
}

unsigned int
az_object_list_contains (AZObjectList *objl, AZObject *obj)
{
	arikkei_return_val_if_fail (objl != NULL, 0);
	arikkei_return_val_if_fail ((objl->weak && AZ_IS_ACTIVE_OBJECT (obj)) || AZ_IS_OBJECT (obj), 0);
	arikkei_return_val_if_fail ((objl->interface && az_object_implements (obj, objl->type)) || az_object_is_a (obj, objl->type), 0);
	for (unsigned int i = 0; i < objl->length; i++) {
		if (objl->objects[i] == obj) return 1;
	}
	return 0;
}

static void
object_list_object_dispose (AZActiveObject *obj, void *data)
{
	AZObjectList *objl = (AZObjectList *) data;
	object_list_remove_object_internal (objl, (AZObject *) obj);
}

static unsigned int
object_list_remove_object_internal (AZObjectList *objl, AZObject *obj)
{
	for (unsigned int i = 0; i < objl->length; i++) {
		if (objl->objects[i] == obj) {
			if (i < (objl->length - 1)) {
				memcpy (&objl->objects[i], &objl->objects[i + 1], (objl->length - 1 - i) * sizeof (AZObject *));
			}
			objl->length -= 1;
			return 1;
		}
	}
	return 0;
}

