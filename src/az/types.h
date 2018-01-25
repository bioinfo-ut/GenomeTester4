#ifndef __AZ_TYPES_H__
#define __AZ_TYPES_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

/* Execution context */
typedef struct _AZContext AZContext;

/* Polymorphic parts of type */
typedef struct _AZImplementation AZImplementation;
/* Semantics of a type */
typedef struct _AZClass AZClass;

/* Members and properties */
typedef struct _AZField AZField;

#include <az/base.h>

/* Predeclarations */
typedef struct _AZReference AZReference;
#ifdef AZ_HAS_STRING
typedef struct _AZString AZString;
#endif
#ifdef AZ_HAS_VALUE
typedef struct _AZValue AZValue;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Every value is a collection of bits in memory
 * The layout, meaning and basic handling of these bits is described by class
 * Classes are accessible either by class pointer or u32 type of value
 * Polymorphic parts of a value is described either in class or implementation
 * For standalone types class is the implementation (no implementation declared separately)
 * Implementations of interface types can be nested in other classes and implementations
 * Class is itself an instance of standalone type
 * Both classes and implementation can contain sub-implementations of other types
 */

/* Predefined typecodes */

enum {
	/* Invalid or missing type, typecode 0 */
	AZ_TYPE_NONE,
	/* Universal base class */
	AZ_TYPE_ANY,
#if 0
	/* Number base classes */
	AZ_TYPE_NUMBER,
	AZ_TYPE_INTEGER,
	AZ_TYPE_REAL,
	AZ_TYPE_COMPLEX,
#endif
	/* True primitives */
	AZ_TYPE_BOOLEAN,
	AZ_TYPE_INT8,
	AZ_TYPE_UINT8,
	AZ_TYPE_INT16,
	AZ_TYPE_UINT16,
	AZ_TYPE_INT32,
	AZ_TYPE_UINT32,
	AZ_TYPE_INT64,
	AZ_TYPE_UINT64,
	AZ_TYPE_FLOAT,
	AZ_TYPE_DOUBLE,
	/* Complex numbers */
	AZ_TYPE_COMPLEX_FLOAT,
	AZ_TYPE_COMPLEX_DOUBLE,
	/* Simple pointer */
	AZ_TYPE_POINTER,

	/* Struct is base for all composite value types */
	AZ_TYPE_STRUCT,
	/* Block is the base of all composite reference types */
	AZ_TYPE_BLOCK,
	/* Special types */
	AZ_TYPE_IMPLEMENTATION,
	AZ_TYPE_CLASS,
	AZ_TYPE_INTERFACE,
	/* Predefined composite types */
	/* AZ_TYPE_ARRAY, */
	AZ_TYPE_REFERENCE,
#ifdef AZ_HAS_STRING
	AZ_TYPE_STRING,
#endif
#ifdef AZ_HAS_VALUE
	AZ_TYPE_VALUE,
#endif
	/* Count */
	AZ_NUM_TYPE_PRIMITIVES
};

#define AZ_TYPE_IS_ARITHMETIC(t) (((t) >= AZ_TYPE_INT8) && ((t) <= AZ_TYPE_COMPLEX_DOUBLE))
#define AZ_TYPE_IS_INTEGRAL(t) (((t) >= AZ_TYPE_INT8) && ((t) <= AZ_TYPE_UINT64))
#define AZ_TYPE_IS_SIGNED(t) (((t) == AZ_TYPE_INT8) || ((t) == AZ_TYPE_INT16) || ((t) == AZ_TYPE_INT32) || ((t) == AZ_TYPE_INT64) || ((t) == AZ_TYPE_FLOAT) || ((t) == AZ_TYPE_DOUBLE))
#define AZ_TYPE_IS_UNSIGNED(t) (((t) == AZ_TYPE_UINT8) || ((t) == AZ_TYPE_UINT16) || ((t) == AZ_TYPE_UINT32) || ((t) == AZ_TYPE_UINT64))
#define AZ_TYPE_IS_64(t) (((t) == AZ_TYPE_INT64) || ((t) == AZ_TYPE_UINT64))

/* Constrained type */
typedef struct _AZTypeConstraint AZTypeConstraint;

struct _AZTypeConstraint {
	unsigned int is_a;
	unsigned int implements;
};

/*
 * Initialize type system
 * If AZ_SAFETY_CHECKS is set it is called automatically
 * It is safe to call it more than once
 */

void az_types_init (void);

/*
 * Basic type queries
 */

unsigned int az_type_is_a (unsigned int type, unsigned int test);
unsigned int az_type_implements (unsigned int type, unsigned int test);
unsigned int az_type_is_assignable_to (unsigned int type, unsigned int test);
unsigned int az_type_get_parent_primitive (unsigned int type);

unsigned int az_instance_serialize (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx);
unsigned int az_instance_deserialize_value (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx);

unsigned int az_instance_to_string (AZClass *klass, void *inst, unsigned char *buf, unsigned int len);

/*
 * Frontend to az_register_class
 * Get new typecode, allocate and initialize a class structure
 */

AZClass *az_register_type (unsigned int *type, unsigned int parent, const unsigned char *name,
	unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *));

/*
* Pass arguments to class_init (for composite types)
*/

AZClass *az_register_composite_type (unsigned int *type, unsigned int parent, const unsigned char *name,
	unsigned int class_size, unsigned int instance_size,
	void (*class_init) (AZClass *, void *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *),
	void *data);

#ifdef AZ_SAFETY_CHECKS
AZClass *az_type_get_class (unsigned int type);
#else
#ifndef __AZ_CLASS_C__
extern AZClass **az_classes;
#endif
#define az_type_get_class(t) az_classes[t]
#endif

/* Initialize a new implementation for interface */
void az_implementation_init (AZImplementation *impl, unsigned int type);

void az_instance_init (void *inst, unsigned int type);
void az_instance_finalize (void *inst, unsigned int type);
/* These are needed for unregistered interfaces */
void az_interface_init (AZImplementation *impl, void *inst);
void az_interface_finalize (AZImplementation *impl, void *inst);

void *az_instance_new (unsigned int type);
void *az_instance_new_array (unsigned int type, unsigned int nelements);
void az_instance_delete (unsigned int type, void *inst);
void az_instance_delete_array (unsigned int type, void *elements, unsigned int nelements);

AZImplementation *az_get_interface (AZImplementation *impl, void *inst, unsigned int type, void **sub_inst);
AZImplementation *az_get_interface_from_type (unsigned int type, void *inst, unsigned int if_type, void **if_inst);

#ifdef AZ_HAS_PROPERTIES
/* Properties */
unsigned int az_instance_set_property (AZClass *klass, void *instance, const unsigned char *key, const AZValue *val);
unsigned int az_instance_get_property (AZClass *klass, void *instance, const unsigned char *key, AZValue *val);
unsigned int az_instance_set_property_by_id (AZClass *klass, void *instance, unsigned int id, const AZValue *val);
unsigned int az_instance_get_property_by_id (AZClass *klass, void *instance, unsigned int id, AZValue *val);
#endif

#ifdef __cplusplus
};
#endif

#endif
