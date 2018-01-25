#define __AZ_CLASS_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>
#include <string.h>

#include <libarikkei/arikkei-strlib.h>

#include <az/class.h>
#ifdef AZ_HAS_PROPERTIES
#include <az/field.h>
#include <az/function-object.h>
#endif
#ifdef AZ_HAS_STRING
#include <az/string.h>
#endif

static void az_class_init (AZClass *klass, unsigned int type, unsigned int parent, unsigned int class_size, unsigned int instance_size, unsigned int flags, const char *name);

static unsigned int classes_size = 0;

void
az_classes_init (void)
{
	if (az_classes) return;
	classes_size = AZ_NUM_TYPE_PRIMITIVES + 32;
	az_classes = (AZClass **) malloc (classes_size * sizeof (AZClass *));
}

static AZClass *implementation_class = NULL;
static AZClass *class_class = NULL;

static unsigned int
implementation_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	return arikkei_strncpy (buf, len, (const unsigned char *) "Implementation");
}

void
az_init_implementation_class (void)
{
	implementation_class = az_class_new (AZ_TYPE_IMPLEMENTATION, AZ_TYPE_BLOCK, sizeof (AZClass), sizeof (AZImplementation), AZ_CLASS_IS_ABSTRACT, "imlementation");
	az_classes[AZ_TYPE_IMPLEMENTATION] = implementation_class;
	az_classes[AZ_TYPE_IMPLEMENTATION]->to_string = implementation_to_string;
}

static unsigned int
class_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned int pos;
	pos = arikkei_memcpy_str (buf, len, klass->name);
	pos += arikkei_strncpy (buf + pos, len - pos, (const unsigned char *) " Class");
	return pos;
}

void
az_init_class_class (void)
{
	class_class = az_class_new (AZ_TYPE_CLASS, AZ_TYPE_IMPLEMENTATION, sizeof (AZClass), sizeof (AZClass), AZ_CLASS_IS_ABSTRACT, "class");
	az_classes[AZ_TYPE_CLASS] = class_class;
	az_classes[AZ_TYPE_CLASS]->to_string = class_to_string;
}

static void *
allocate_default (AZClass *klass)
{
	return malloc (klass->instance_size);
}

static void *
allocate_array_default (AZClass *klass, unsigned int nelements)
{
	return malloc (nelements * klass->element_size);
}

static void
free_default (AZClass *klass, void *location)
{
	free (location);
}

static void
free_array_default (AZClass *klass, void *location, unsigned int nelements)
{
	free (location);
}

AZClass *
az_class_new (unsigned int type, unsigned int parent, unsigned int class_size, unsigned int instance_size, unsigned int flags, const char *name)
{
	AZClass *klass = (AZClass *) malloc (class_size);
	az_class_init (klass, type, parent, class_size, instance_size, flags, name);
	return klass;
}

static void
az_class_init (AZClass *klass, unsigned int type, unsigned int parent, unsigned int class_size, unsigned int instance_size, unsigned int flags, const char *name)
{
	memset (klass, 0, class_size);
	if (parent) {
		memcpy (klass, az_classes[parent], az_classes[parent]->class_size);
		/* Overwrite values from supertype */
		klass->flags &= ~AZ_CLASS_IS_ABSTRACT;
		klass->parent = az_classes[parent];
		klass->firstinterface = klass->parent->firstinterface + klass->parent->ninterfaces;
		klass->ninterfaces = 0;
		klass->impl_types = NULL;
		klass->impl_offsets = NULL;
		klass->inst_offsets = NULL;
#ifdef AZ_HAS_PROPERTIES
		klass->firstproperty = klass->parent->firstproperty + klass->parent->nproperties;
		klass->nproperties = 0;
		klass->properties = NULL;
#endif
	}
	klass->implementation.type = type;
	klass->flags |= flags;
	klass->name = (const unsigned char *) name;
	klass->class_size = class_size;
	klass->instance_size = instance_size;
	if (klass->flags | AZ_CLASS_IS_VALUE) {
		klass->value_size = instance_size;
	} else {
		klass->value_size = sizeof (void *);
	}
	klass->element_size = klass->value_size;
	if (!parent) {
		/* Default methods */
		/* Memory management */
		klass->allocate = allocate_default;
		klass->allocate_array = allocate_array_default;
		klass->free = free_default;
		klass->free_array = free_array_default;
	}
}

void
az_register_class (AZClass *klass, unsigned int *type, unsigned int parent, const unsigned char *name, unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *))
{
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_if_fail (klass != NULL);
	arikkei_return_if_fail ((parent == AZ_TYPE_NONE) || (class_size >= az_classes[parent]->class_size));
	arikkei_return_if_fail ((parent == AZ_TYPE_NONE) || (instance_size >= az_classes[parent]->instance_size));
#endif
	if (az_num_classes >= classes_size) {
		classes_size += 32;
		az_classes = (AZClass **) realloc (az_classes, classes_size * sizeof (AZClass *));
	}
	*type = az_num_classes++;

	az_classes[*type] = klass;
	az_class_init (klass, *type, parent, class_size, instance_size, 0, (const char *) name);

	/* Constructors and destructors */
	klass->instance_init = instance_init;
	klass->instance_finalize = instance_finalize;

	if (class_init) class_init (klass);
}

void
az_class_set_num_interfaces (AZClass *klass, unsigned int ninterfaces)
{
	if (klass->parent) klass->firstinterface = klass->parent->firstinterface + klass->parent->ninterfaces;
	klass->ninterfaces = ninterfaces;
	klass->impl_types = (unsigned int *) malloc (ninterfaces * 4);
	klass->impl_offsets = (unsigned int *) malloc (ninterfaces * 4);
	klass->inst_offsets = (unsigned int *) malloc (ninterfaces * 4);
}

void
az_class_declare_interface (AZClass *klass, unsigned int idx, unsigned int type, unsigned int impl_offset, unsigned int inst_offset)
{
	arikkei_return_if_fail (klass != NULL);
	arikkei_return_if_fail (idx < klass->ninterfaces);
	arikkei_return_if_fail (az_type_is_a (type, AZ_TYPE_INTERFACE));
	klass->impl_types[idx] = type;
	klass->impl_offsets[idx] = impl_offset;
	klass->inst_offsets[idx] = inst_offset;
	az_implementation_init ((AZImplementation *) ((char *) klass + impl_offset), type);
}

#ifdef AZ_HAS_PROPERTIES
void
az_class_set_num_properties (AZClass *klass, unsigned int nproperties)
{
	if (klass->parent) klass->firstproperty = klass->parent->firstproperty + klass->parent->nproperties;
	klass->nproperties = nproperties;
	klass->properties = (AZField *) malloc (nproperties * sizeof (AZField));
	memset (klass->properties, 0, nproperties * sizeof (AZField));
}

void
az_class_property_setup (AZClass *klass, unsigned int idx, const unsigned char *key, unsigned int type,
unsigned int is_static, unsigned int can_read, unsigned int can_write, unsigned int is_final, unsigned int is_value,
unsigned int value_type, void *value)
{
	arikkei_return_if_fail (klass != NULL);
	arikkei_return_if_fail (idx < klass->nproperties);
	arikkei_return_if_fail (key != NULL);
	arikkei_return_if_fail (type != AZ_TYPE_NONE);
	arikkei_return_if_fail (!(can_write && is_final));
	arikkei_return_if_fail (!is_value || is_static);
	if (!((value_type == AZ_TYPE_NONE) || (az_type_is_assignable_to (value_type, type)))) {
		return;
	}
	arikkei_return_if_fail ((value_type == AZ_TYPE_NONE) || (az_type_is_assignable_to (value_type, type)));

	az_property_setup (klass->properties + idx, key, type, klass->firstproperty + idx, is_static, can_read, can_write, is_final, is_value, value_type, value);
}

AZField *
az_class_lookup_property (AZClass *klass, const unsigned char *key)
{
	unsigned int i;
	arikkei_return_val_if_fail (klass != NULL, NULL);
	arikkei_return_val_if_fail (key != NULL, NULL);
	for (i = 0; i < klass->nproperties; i++) {
		if (!strcmp ((const char *) key, (const char *) klass->properties[i].key->str)) return &klass->properties[i];
	}
	if (klass->parent) {
		return az_class_lookup_property (klass->parent, key);
	}
	return NULL;
}

void
az_class_method_setup (AZClass *klass, unsigned int idx, const unsigned char *key,
unsigned int rettype, unsigned int nargs, const unsigned int argtypes[],
unsigned int (*call) (AZValue *, AZValue *, AZValue *))
{
	AZFunctionObject *fobj;
	fobj = az_function_object_new (klass->implementation.type, rettype, nargs, argtypes, call);
	/* Property is static although function is not static */
	az_class_property_setup (klass, idx, key, AZ_TYPE_FUNCTION, 1, 1, 0, 1, 1, AZ_TYPE_FUNCTION_OBJECT, fobj);
	az_reference_unref ((AZReferenceClass *) fobj->object.klass, (AZReference *) fobj);
}

void
az_class_static_method_setup (AZClass *klass, unsigned int idx, const unsigned char *key,
unsigned int rettype, unsigned int nargs, const unsigned int argtypes[],
unsigned int (*call) (AZValue *, AZValue *, AZValue *))
{
	AZFunctionObject *fobj;
	fobj = az_function_object_new (AZ_TYPE_NONE, rettype, nargs, argtypes, call);
	az_class_property_setup (klass, idx, key, AZ_TYPE_FUNCTION, 1, 1, 0, 1, 1, AZ_TYPE_FUNCTION_OBJECT, fobj);
	az_reference_unref ((AZReferenceClass *) fobj->object.klass, (AZReference *) fobj);
}

#endif
