#define __AZ_REFERENCE_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include <az/class.h>
#include <az/reference.h>

void
az_reference_unref_last (AZReferenceClass *klass, AZReference *ref)
{
	arikkei_return_if_fail (ref->refcount > 0);
	if (ref->refcount == 1) {
		if (!klass->drop || klass->drop (klass, ref)) {
			if (klass->dispose) klass->dispose (klass, ref);
			az_instance_delete (klass->klass.implementation.type, ref);
		}
	} else {
		ref->refcount -= 1;
	}
}

void
az_reference_dispose (AZReferenceClass *klass, AZReference *ref)
{
	arikkei_return_if_fail (ref->refcount > 0);
	if (klass->dispose) klass->dispose (klass, ref);
	ref->refcount -= 1;
	if (!ref->refcount) {
		az_instance_delete (klass->klass.implementation.type, ref);
	}
}

static AZReferenceClass *reference_class = NULL;

static void
reference_instance_init (AZReferenceClass *klass, void *instance)
{
	((AZReference *) instance)->refcount = 1;
}

void
az_init_reference_class (void)
{
	reference_class = (AZReferenceClass *) az_class_new (AZ_TYPE_REFERENCE, AZ_TYPE_BLOCK, sizeof (AZReferenceClass), sizeof (AZReference), AZ_CLASS_IS_ABSTRACT, "reference");
	az_classes[AZ_TYPE_REFERENCE] = (AZClass *) reference_class;
	az_classes[AZ_TYPE_REFERENCE]->instance_init = (void (*) (AZImplementation *, void *)) reference_instance_init;
}
