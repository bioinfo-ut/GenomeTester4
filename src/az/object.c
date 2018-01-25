#define __AZ_OBJECT_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdio.h>

#include <az/object.h>

/* AZInstance implementation */
static void object_class_init (AZObjectClass *klass);
static void object_init (AZObjectClass *klass, AZObject *object);
/* AZReference implementation */
void object_dispose (AZReferenceClass *klass, AZReference *ref);

static unsigned int object_type = 0;
static AZObjectClass *object_class;

unsigned int
az_object_get_type (void)
{
	if (!object_type) {
		object_class = (AZObjectClass *) az_register_type (&object_type, AZ_TYPE_REFERENCE, (const unsigned char *) "AZObject",
			sizeof (AZObjectClass), sizeof (AZObject),
			(void (*) (AZClass *)) object_class_init,
			(void (*) (AZImplementation *, void *)) object_init,
			NULL);
	}
	return object_type;
}

static void
object_class_init (AZObjectClass *klass)
{
	((AZClass *) klass)->flags |= AZ_CLASS_ZERO_MEMORY;
	klass->reference_klass.dispose = object_dispose;
}

static void
object_init (AZObjectClass *klass, AZObject *object)
{
	object->flags |= AZ_OBJECT_ALIVE;
}

void
object_dispose (AZReferenceClass *klass, AZReference *ref)
{
	AZObject *obj = AZ_OBJECT (ref);
	if (obj->flags & AZ_OBJECT_ALIVE) {
		if (obj->klass->shutdown) {
			obj->klass->shutdown (obj);
		}
	}
}

AZObject *
az_object_new (unsigned int type)
{
	AZObject *object;
	arikkei_return_val_if_fail (az_type_is_a (type, AZ_TYPE_OBJECT), NULL);
	object = (AZObject *) az_instance_new (type);
	object->klass = (AZObjectClass *) az_type_get_class (type);
	return object;
}

void
az_object_shutdown (AZObject *obj)
{
	arikkei_return_if_fail (obj != NULL);
	arikkei_return_if_fail (AZ_IS_OBJECT (obj));
	arikkei_return_if_fail (obj->flags & AZ_OBJECT_ALIVE);
	if (obj->klass->shutdown) {
		obj->klass->shutdown (obj);
	}
	obj->flags &= ~AZ_OBJECT_ALIVE;
	az_reference_unref (&obj->klass->reference_klass, (AZReference *) obj);
}

void *
az_object_check_instance_cast (void *inst, unsigned int type)
{
	if (inst == NULL) {
		fprintf (stderr, "az_object_check_instance_cast: inst == NULL\n");
		return NULL;
	}
	if (!az_type_is_a (((AZObject *) inst)->klass->reference_klass.klass.implementation.type, type)) {
		AZClass *klass = az_type_get_class (type);
		fprintf (stderr, "az_object_check_instance_cast: %s is not %s\n", ((AZObject *) inst)->klass->reference_klass.klass.name, klass->name);
		return NULL;
	}
	return inst;
}

unsigned int
az_object_check_instance_type (void *inst, unsigned int type)
{
	if (inst == NULL) return 0;
	return az_type_is_a (((AZObject *) inst)->klass->reference_klass.klass.implementation.type, type);
}

unsigned int
az_object_implements (AZObject *obj, unsigned int type)
{
	arikkei_return_val_if_fail (obj != NULL, 0);
	arikkei_return_val_if_fail (AZ_IS_OBJECT (obj), 0);
	return az_type_implements (obj->klass->reference_klass.klass.implementation.type, type);
}

AZImplementation *
az_object_get_interface (AZObject *obj, unsigned int type, void **inst)
{
	arikkei_return_val_if_fail (obj != NULL, NULL);
	arikkei_return_val_if_fail (AZ_IS_OBJECT (obj), NULL);
	return az_get_interface (&obj->klass->reference_klass.klass.implementation, obj, type, inst);
}

#ifdef AZ_HAS_VALUE
void
az_value_set_object (AZValue *val, AZObject *obj)
{
	if (!obj) {
		az_value_set_reference (val, AZ_TYPE_OBJECT, NULL);
	} else {
		az_value_set_reference (val, obj->klass->reference_klass.klass.implementation.type, (AZReference *) obj);
	}
}

void
az_value_transfer_object (AZValue *val, AZObject *obj)
{
	if (!obj) {
		az_value_set_reference (val, AZ_TYPE_OBJECT, NULL);
	} else {
		az_value_transfer_reference (val, obj->klass->reference_klass.klass.implementation.type, (AZReference *) obj);
	}
}

#endif
