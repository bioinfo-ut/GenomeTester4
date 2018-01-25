#define __AZ_TYPES_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <libarikkei/arikkei-utils.h>

#include <az/class.h>
#include <az/interface.h>
#include <az/primitives.h>
#include <az/reference.h>
#ifdef AZ_HAS_STRING
#include <az/string.h>
#endif
#ifdef AZ_HAS_VALUE
#include <az/value.h>
#endif
#ifdef AZ_HAS_PROPERTIES
#include <az/field.h>
#endif

#include <az/types.h>

void
az_types_init (void)
{
	if (az_classes) return;
	az_classes_init ();
	az_init_primitive_classes ();
	az_init_implementation_class ();
	az_init_class_class ();
	az_init_interface_class ();
	az_init_reference_class ();
#ifdef AZ_HAS_STRING
	az_init_string_class ();
#endif
#ifdef AZ_HAS_VALUE
	az_init_value_class ();
#endif
	az_num_classes = AZ_NUM_TYPE_PRIMITIVES;
}

unsigned int
az_type_is_a (unsigned int type, unsigned int test)
{
	AZClass *klass;
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail (type < az_num_classes, 0);
	arikkei_return_val_if_fail (test < az_num_classes, 0);
#endif
	if (type == test) return 1;
	klass = az_classes[type]->parent;
	while (klass) {
		if (klass->implementation.type == test) return 1;
		klass = klass->parent;
	}
	return 0;
}

unsigned int
az_type_implements (unsigned int type, unsigned int test)
{
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail (type < az_num_classes, 0);
	arikkei_return_val_if_fail (test < az_num_classes, 0);
#endif
	return az_get_interface (&az_classes[type]->implementation, NULL, test, NULL) != NULL;
}

unsigned int
az_type_is_assignable_to (unsigned int type, unsigned int test)
{
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail (type < az_num_classes, 0);
	arikkei_return_val_if_fail (test < az_num_classes, 0);
#endif
	if (az_type_is_a (test, AZ_TYPE_INTERFACE)) {
		return az_type_implements (type, test);
	} else {
		return az_type_is_a (type, test);
	}
}

unsigned int
az_type_get_parent_primitive (unsigned int type)
{
	AZClass *klass;
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail (type < az_num_classes, 0);
#endif
	if (type < AZ_NUM_TYPE_PRIMITIVES) return type;
	klass = az_classes[type]->parent;
	while (klass) {
		if (klass->implementation.type < AZ_NUM_TYPE_PRIMITIVES) return klass->implementation.type;
		klass = klass->parent;
	}
	return 0;
}

unsigned int
az_instance_serialize (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	AZClass *klass;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (inst != NULL, 0);
#endif
	klass = az_type_get_class (impl->type);
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (klass != NULL, 0);
#endif
	if (klass->serialize) return klass->serialize (impl, inst, d, dlen, ctx);
	return 0;
}

unsigned int
az_instance_deserialize_value (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	AZClass *klass;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (value != NULL, 0);
#endif
	klass = az_type_get_class (impl->type);
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (klass != NULL, 0);
#endif
	if (klass->deserialize) return klass->deserialize (impl, value, s, slen, ctx);
	return 0;
}

/* fixme: Make signature correct */

unsigned int
az_instance_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	if (klass->to_string) {
		return klass->to_string (klass, instance, buf, len);
	}
	if (len) buf[0] = 0;
	return 0;
}

AZClass *
az_register_type (unsigned int *type, unsigned int parent, const unsigned char *name, unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *))
{
	AZClass *klass;
#ifdef AZ_SAFETY_CHECKS
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail ((parent == AZ_TYPE_NONE) || (class_size >= az_classes[parent]->class_size), NULL);
	arikkei_return_val_if_fail ((parent == AZ_TYPE_NONE) || (instance_size >= az_classes[parent]->instance_size), NULL);
#endif
	klass = (AZClass *) malloc (class_size);
	az_register_class (klass, type, parent, name, class_size, instance_size, class_init, instance_init, instance_finalize);
	return klass;
}

AZClass *az_register_composite_type (unsigned int *type, unsigned int parent, const unsigned char *name, unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *, void *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *),
	void *data)
{
	AZClass *klass = az_register_type (type, parent, name, class_size, instance_size, NULL, instance_init, instance_finalize);
	if (class_init) class_init (klass, data);
	return klass;
}

#ifdef AZ_SAFETY_CHECKS
AZClass *
az_type_get_class (unsigned int type)
{
	if (!az_classes) az_types_init ();
	arikkei_return_val_if_fail (type < az_num_classes, NULL);
	return az_classes[type];
}
#endif

static void
interface_init_recursive (AZClass *klass, AZImplementation *impl, void *inst)
{
	unsigned int i;
	/* Every interface has to be subclass of AZInterface */
	if (klass->parent && (klass->parent->implementation.type >= AZ_TYPE_INTERFACE)) {
		interface_init_recursive (klass->parent, impl, inst);
	}
	/* Interfaces */
	for (i = 0; i < klass->ninterfaces; i++) {
		AZImplementation *sub_impl = (AZImplementation *) ((char *) impl + klass->impl_offsets[i]);
		void *sub_inst = (void *) ((char *) inst + klass->inst_offsets[i]);
		AZClass *sub_class = az_classes[sub_impl->type];
		if (sub_class->flags & AZ_CLASS_ZERO_MEMORY) memset (sub_inst, 0, sub_class->instance_size);
		interface_init_recursive (sub_class, sub_impl, sub_inst);
	}
	/* Instance itself */
	if (klass->instance_init) klass->instance_init (impl, inst);
}

static void
interface_finalize_recursive (AZClass *klass, AZImplementation *impl, void *inst)
{
	unsigned int i;
	if (klass->instance_finalize) klass->instance_finalize (impl, inst);
	for (i = 0; i < klass->ninterfaces; i++) {
		AZImplementation *sub_impl = (AZImplementation *) ((char *) impl + klass->impl_offsets[i]);
		void *sub_inst = (void *) ((char *) inst + klass->inst_offsets[i]);
		AZClass *sub_class = az_classes[sub_impl->type];
		interface_finalize_recursive (sub_class, sub_impl, sub_inst);
	}
	if (klass->parent && (klass->parent->implementation.type >= AZ_TYPE_INTERFACE)) {
		interface_finalize_recursive (klass->parent, impl, inst);
	}
}

static void
instance_init_recursive (AZClass *klass, void *inst)
{
	unsigned int i;
	/* Initialize parent instances */
	if (klass->parent && (klass->parent->implementation.type != AZ_TYPE_NONE)) {
		instance_init_recursive (klass->parent, inst);
	}
	/* Interfaces */
	for (i = 0; i < klass->ninterfaces; i++) {
		AZImplementation *sub_impl = (AZImplementation *) ((char *) klass + klass->impl_offsets[i]);
		void *sub_inst = (void *) ((char *) inst + klass->inst_offsets[i]);
		AZClass *sub_class = az_classes[sub_impl->type];
		if (sub_class->flags & AZ_CLASS_ZERO_MEMORY) memset (sub_inst, 0, sub_class->instance_size);
		interface_init_recursive (sub_class, sub_impl, sub_inst);
	}
	/* Instance itself */
	if (klass->instance_init) klass->instance_init (&klass->implementation, inst);
}

static void
instance_finalize_recursive (AZClass *klass, void *inst)
{
	unsigned int i;
	if (klass->instance_finalize) klass->instance_finalize (&klass->implementation, inst);
	for (i = 0; i < klass->ninterfaces; i++) {
		AZImplementation *sub_impl = (AZImplementation *) ((char *) klass + klass->impl_offsets[i]);
		void *sub_inst = (void *) ((char *) inst + klass->inst_offsets[i]);
		AZClass *sub_class = az_classes[sub_impl->type];
		interface_finalize_recursive (sub_class, sub_impl, sub_inst);
	}
	if (klass->parent && (klass->parent->implementation.type != AZ_TYPE_NONE)) {
		instance_finalize_recursive (klass->parent, inst);
	}
}

static void
implementation_init_recursive (AZInterfaceClass *interface_class, AZImplementation *impl)
{
	unsigned int i;
	AZClass *klass = (AZClass *) interface_class;
	/* Init superimplementations */
	if (klass->parent && (klass->parent->implementation.type >= AZ_TYPE_INTERFACE)) {
		implementation_init_recursive ((AZInterfaceClass *) klass->parent, impl);
	}
	/* Init subimplementations */
	for (i = 0; i < klass->ninterfaces; i++) {
		AZImplementation *sub_impl = (AZImplementation *) ((char *) impl + klass->impl_offsets[i]);
		AZInterfaceClass *sub_class = (AZInterfaceClass *) az_classes[klass->impl_types[i]];
		if (sub_class->klass.flags & AZ_CLASS_ZERO_MEMORY) {
			memset (sub_impl, 0, sub_class->implementation_size);
		}
		sub_impl->type = klass->impl_types[i];
		implementation_init_recursive (sub_class, sub_impl);
	}
	/* Implementation itself */
	if (interface_class->implementation_init) interface_class->implementation_init (impl);
}

void
az_implementation_init (AZImplementation *impl, unsigned int type)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (impl != NULL);
	arikkei_return_if_fail (type != 0);
	arikkei_return_if_fail (type < az_num_classes);
	arikkei_return_if_fail (az_type_is_a (type, AZ_TYPE_INTERFACE));
#endif
	AZInterfaceClass *ifclass = (AZInterfaceClass *) az_type_get_class (type);
	impl->type = type;
	implementation_init_recursive (ifclass, impl);
}

void
az_instance_init (void *inst, unsigned int type)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (inst != NULL);
	arikkei_return_if_fail (type != 0);
	arikkei_return_if_fail (type < az_num_classes);
	arikkei_return_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_classes[type];
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (!(klass->flags & AZ_CLASS_IS_ABSTRACT));
#endif
	if (klass->flags & AZ_CLASS_ZERO_MEMORY) memset (inst, 0, klass->instance_size);
	instance_init_recursive (klass, inst);
}

void
az_instance_finalize (void *inst, unsigned int type)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (inst != NULL);
	arikkei_return_if_fail (type != 0);
	arikkei_return_if_fail (type < az_num_classes);
	arikkei_return_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_classes[type];
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (!(klass->flags & AZ_CLASS_IS_ABSTRACT));
#endif
	instance_finalize_recursive (klass, inst);
}

void
az_interface_init (AZImplementation *impl, void *inst)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (impl != NULL);
	arikkei_return_if_fail (inst != NULL);
	arikkei_return_if_fail (impl->type != 0);
	arikkei_return_if_fail (impl->type < az_num_classes);
	arikkei_return_if_fail (az_type_is_a (impl->type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_classes[impl->type];
	if (klass->flags & AZ_CLASS_ZERO_MEMORY) memset (inst, 0, klass->instance_size);
	interface_init_recursive (klass, impl, inst);
}

void
az_interface_finalize (AZImplementation *impl, void *inst)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (impl != NULL);
	arikkei_return_if_fail (inst != NULL);
	arikkei_return_if_fail (impl->type != 0);
	arikkei_return_if_fail (impl->type < az_num_classes);
	arikkei_return_if_fail (az_type_is_a (impl->type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_classes[impl->type];
	interface_finalize_recursive (klass, impl, inst);
}

void *
az_instance_new (unsigned int type)
{
	void *instance;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (type != 0, NULL);
	arikkei_return_val_if_fail (type < az_num_classes, NULL);
	arikkei_return_val_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE), NULL);
#endif
	AZClass *klass = az_type_get_class (type);
	if (klass->allocate) {
		instance = klass->allocate (klass);
	} else {
		instance = malloc (klass->instance_size);
	}
	if (klass->flags & AZ_CLASS_ZERO_MEMORY) memset (instance, 0, klass->instance_size);
	instance_init_recursive (klass, instance);
	return instance;
}

void *
az_instance_new_array (unsigned int type, unsigned int nelements)
{
	void *elements;
	unsigned int i;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (type != 0, NULL);
	arikkei_return_val_if_fail (type < az_num_classes, NULL);
	arikkei_return_val_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE), NULL);
#endif
	AZClass *klass = az_type_get_class (type);
	if (klass->allocate_array) {
		elements = klass->allocate_array (klass, nelements);
	} else {
		elements = malloc (nelements * klass->element_size);
	}
	if (klass->flags & AZ_CLASS_ZERO_MEMORY) memset (elements, 0, nelements * klass->element_size);
	for (i = 0; i < nelements; i++) {
		void *instance = (char *) elements + i * klass->element_size;
		instance_init_recursive (klass, instance);
	}
	return elements;
}

void
az_instance_delete (unsigned int type, void *instance)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (type != 0);
	arikkei_return_if_fail (type < az_num_classes);
	arikkei_return_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_type_get_class (type);
	instance_finalize_recursive (klass, instance);
	if (klass->free) {
		klass->free (klass, instance);
	} else {
		free (instance);
	}
}

void
az_instance_delete_array (unsigned int type, void *elements, unsigned int nelements)
{
	unsigned int i;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (type != 0);
	arikkei_return_if_fail (type < az_num_classes);
	arikkei_return_if_fail (!az_type_is_a (type, AZ_TYPE_INTERFACE));
#endif
	AZClass *klass = az_type_get_class (type);
	for (i = 0; i < nelements; i++) {
		void *instance = (char *) elements + i * klass->element_size;
		instance_finalize_recursive (klass, instance);
	}
	if (klass->free_array) {
		klass->free_array (klass, elements, nelements);
	} else {
		free (elements);
	}
}

AZImplementation *
az_get_interface (AZImplementation *impl, void *inst, unsigned int type, void **sub_inst)
{
	arikkei_return_val_if_fail (impl != NULL, NULL);
	/* Get our class */
	AZClass *klass = az_classes[impl->type];
	while (klass) {
		unsigned int i;
		/* Iterate over interfaces */
		for (i = 0; i < klass->ninterfaces; i++) {
			AZImplementation *sub_impl = (AZImplementation *) ((char *) impl + klass->impl_offsets[i]);
			/* If this interface is of requested type we are done */
			if (az_type_is_a (sub_impl->type, type)) {
				if (sub_inst) *sub_inst = (char *) inst + klass->inst_offsets[i];
				return sub_impl;
			}
			/* Check sub-interfaces of this interface */
			sub_impl = az_get_interface (sub_impl, (char *) inst + klass->inst_offsets[i], type, sub_inst);
			if (sub_impl) return sub_impl;
		}
		klass = klass->parent;
	}
	if (sub_inst) *sub_inst = NULL;
	return NULL;
}

AZImplementation *
az_get_interface_from_type (unsigned int type, void *inst, unsigned int if_type, void **if_inst)
{
	return az_get_interface (&az_classes[type]->implementation, inst, if_type, if_inst);
}

#ifdef AZ_HAS_PROPERTIES
unsigned int
az_instance_set_property (AZClass *klass, void *instance, const unsigned char *key, const AZValue *val)
{
	AZField *prop;
	arikkei_return_val_if_fail (klass != NULL, 0);
	arikkei_return_val_if_fail (key != NULL, 0);
	arikkei_return_val_if_fail (val != NULL, 0);
	prop = az_class_lookup_property (klass, key);
	if (!prop) return 0;
	return az_instance_set_property_by_id (klass, instance, prop->id, val);
}

unsigned int
az_instance_set_property_by_id (AZClass *klass, void *instance, unsigned int id, const AZValue *val)
{
	arikkei_return_val_if_fail (klass != NULL, 0);
	arikkei_return_val_if_fail (val != NULL, 0);
	if (id >= klass->firstproperty) {
		unsigned int idx = id - klass->firstproperty;
		if (!klass->properties[idx].can_write) return 0;
		if (klass->properties[idx].is_final) return 0;
		if (val && !(az_type_is_a (val->impl->type, klass->properties[idx].type) || az_type_implements (val->impl->type, klass->properties[idx].type))) return 0;
		if (klass->set_property) {
			return klass->set_property (klass, instance, idx, val);
		}
		return 0;
	} else {
		if (klass->parent) {
			return az_instance_set_property_by_id (klass->parent, instance, id, val);
		}
		return 0;
	}
}

unsigned int
az_instance_get_property (AZClass *klass, void *instance, const unsigned char *key, AZValue *val)
{
	AZField *prop;
	arikkei_return_val_if_fail (klass != NULL, 0);
	arikkei_return_val_if_fail (key != NULL, 0);
	arikkei_return_val_if_fail (val != NULL, 0);
	prop = az_class_lookup_property (klass, key);
	if (!prop) return 0;
	return az_instance_get_property_by_id (klass, instance, prop->id, val);
}

unsigned int
az_instance_get_property_by_id (AZClass *klass, void *instance, unsigned int id, AZValue *val)
{
	arikkei_return_val_if_fail (klass != NULL, 0);
	arikkei_return_val_if_fail (val != NULL, 0);
	if (id >= klass->firstproperty) {
		unsigned int idx = id - klass->firstproperty;
		if (!klass->properties[idx].can_read) return 0;
		if (klass->properties[idx].is_value) {
			az_value_copy (val, &klass->properties[idx].value);
			return 1;
		} else if (klass->get_property) {
			return klass->get_property (klass, instance, idx, val);
		}
		return 0;
	} else {
		if (klass->parent) {
			return az_instance_get_property_by_id (klass->parent, instance, id, val);
		}
		return 0;
	}
}

#endif
