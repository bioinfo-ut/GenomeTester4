#ifndef __AZ_CLASS_H__
#define __AZ_CLASS_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <az/types.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZImplementation {
	unsigned int type;
};

/* Class flags */

/* No instancing is allowed (this is not propagated to subclasses) */
#define AZ_CLASS_IS_ABSTRACT 1
/* No subclasses */
#define AZ_CLASS_IS_FINAL 2
/* Instance is value */
#define AZ_CLASS_IS_VALUE 4
/* Subclasses should not remove flag set by parent */
#define AZ_CLASS_ZERO_MEMORY 8

struct _AZClass {
	AZImplementation implementation;

	unsigned int flags;

	AZClass *parent;

	/* All implemented interfaces */
	unsigned int firstinterface;
	unsigned int ninterfaces;
	/* Types are stored in class to automatically initialize sub-implementations of implementation */
	unsigned int *impl_types;
	unsigned int *impl_offsets;
	unsigned int *inst_offsets;

#ifdef AZ_HAS_PROPERTIES
	unsigned int firstproperty;
	unsigned int nproperties;
	AZField *properties;
#endif

	const unsigned char *name;
	/* Size of class structure */
	unsigned int class_size;
	/* Size of instance */
	unsigned int instance_size;
	/* Size of value */
	unsigned int value_size;
	/* Size of values in arrays (rounded to 16 bytes for aligned types) */
	unsigned int element_size;

	/* Memory management */
	void *(*allocate) (AZClass *klass);
	void *(*allocate_array) (AZClass *klass, unsigned int nelements);
	void (*free) (AZClass *klass, void *location);
	void (*free_array) (AZClass *klass, void *location, unsigned int nelements);
	/* Constructors and destructors */
	void (*instance_init) (AZImplementation *impl, void *inst);
	void (*instance_finalize) (AZImplementation *impl, void *inst);
	/* Duplicate creates a copy in uninitialized memory */
	void (*duplicate) (AZClass *klass, void *destination, void *instance);
	/* Assign overwrites existing initialized instance */
	void (*assign) (AZClass *klass, void *destination, void *instance);

	/* Serialization is by instance */
	/* Return number of bytes that should have been written (regardless of dlen) */
	/* It is safe to set d to NULL */
	/* Returns the number of bytes that would have been written if there was enough room in destination */
	unsigned int (*serialize) (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx);
	/* Deserialization is by value, i.e. new instances of reference types should be created */
	/* Returns the number of bytes consumed (0 on error) */
	unsigned int (*deserialize) (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx);

	unsigned int (*to_string) (AZClass *klass, void *instance, unsigned char *buf, unsigned int len);
#ifdef AZ_HAS_PROPERTIES
	unsigned int (*get_property) (AZClass *klass, void *instance, unsigned int idx, AZValue *val);
	unsigned int (*set_property) (AZClass *klass, void *instance, unsigned int idx, const AZValue *val);
#endif
};

/* C array of all defined classes */

#ifndef __AZ_CLASS_C__
extern AZClass **az_classes;
extern unsigned int az_num_classes;
#else
AZClass **az_classes = NULL;
unsigned int az_num_classes = 0;
#endif

/* To be called from class constructors */
void az_class_set_num_interfaces (AZClass *klass, unsigned int ninterfaces);
void az_class_declare_interface (AZClass *klass, unsigned int idx, unsigned int type, unsigned int impl_offset, unsigned int inst_offset);

#ifdef AZ_HAS_PROPERTIES
void az_class_set_num_properties (AZClass *klass, unsigned int nproperties);
void az_class_property_setup (AZClass *klass, unsigned int idx, const unsigned char *key, unsigned int type,
	unsigned int is_static, unsigned int can_read, unsigned int can_write, unsigned int is_final, unsigned int is_value,
	unsigned int value_type, void *value);
AZField *az_class_lookup_property (AZClass *klass, const unsigned char *key);
#endif

#ifdef AZ_HAS_PROPERTIES
void az_class_method_setup (AZClass *klass, unsigned int idx, const unsigned char *key,
	unsigned int rettype, unsigned int nargs, const unsigned int argtypes[],
	unsigned int (*call) (AZValue *, AZValue *, AZValue *));
void az_class_static_method_setup (AZClass *klass, unsigned int idx, const unsigned char *key,
	unsigned int rettype, unsigned int nargs, const unsigned int argtypes[],
	unsigned int (*call) (AZValue *, AZValue *, AZValue *));
#endif

/* Library internals */
void az_classes_init (void);
void az_init_implementation_class (void);
void az_init_class_class (void);
/* Creates an intialized class structure with predefined type (does not register it) */
AZClass *az_class_new (unsigned int type, unsigned int parent, unsigned int class_size, unsigned int instance_size, unsigned int flags, const char *name);
/* Registers and initializes a new class */
void az_register_class (AZClass *klass, unsigned int *type, unsigned int parent, const unsigned char *name, unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *));


#ifdef __cplusplus
};
#endif

#endif
