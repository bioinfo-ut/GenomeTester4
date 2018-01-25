#ifndef __AZ_ACTIVE_OBJECT_H__
#define __AZ_ACTIVE_OBJECT_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_ACTIVE_OBJECT az_active_object_get_type ()
#define AZ_ACTIVE_OBJECT(o) (AZ_CHECK_INSTANCE_CAST ((o), AZ_TYPE_ACTIVE_OBJECT, AZActiveObject))
#define AZ_IS_ACTIVE_OBJECT(o) (AZ_CHECK_INSTANCE_TYPE ((o), AZ_TYPE_ACTIVE_OBJECT))

typedef struct _AZActiveObject AZActiveObject;
typedef struct _AZActiveObjectClass AZActiveObjectClass;

#include <az/object.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _AZObjectAttribute AZObjectAttribute;
typedef struct _AZObjectAttributeArray AZObjectAttributeArray;
typedef struct _AZObjectListener AZObjectListener;
typedef struct _AZObjectCallbackBlock AZObjectCallbackBlock;
typedef struct _AZObjectEventVector AZObjectEventVector;

struct _AZObjectAttribute {
	AZString *key;
	AZValue value;
};

struct _AZObjectAttributeArray {
	unsigned int size;
	unsigned int length;
	AZObjectAttribute attribs[1];
};

struct _AZObjectEventVector {
	void (*dispose) (AZActiveObject *object, void *data);
};

struct _AZObjectListener {
	const AZObjectEventVector *vector;
	unsigned int size;
	void *data;
};

struct _AZObjectCallbackBlock {
	unsigned int size;
	unsigned int length;
	AZObjectListener listeners[1];
};

struct _AZActiveObject {
	AZObject object;
	AZObjectCallbackBlock *callbacks;
	AZObjectAttributeArray *attributes;
};

struct _AZActiveObjectClass {
	AZObjectClass object_class;
};

unsigned int az_active_object_get_type (void);

unsigned int az_active_object_get_attribute (AZActiveObject *aobj, const unsigned char *key, AZValue *val);
unsigned int az_active_object_set_attribute (AZActiveObject *aobj, const unsigned char *key, const AZValue *val);
unsigned int az_active_object_clear_attribute (AZActiveObject *aobj, const unsigned char *key);

void az_active_object_add_listener (AZActiveObject *aobj, const AZObjectEventVector *vector, unsigned int size, void *data);
void az_active_object_remove_listener_by_data (AZActiveObject *aobj, void *data);

/* Helpers */
unsigned int az_active_object_set_attribute_i32 (AZActiveObject *aobj, const unsigned char *key, int value);
unsigned int az_active_object_set_attribute_object (AZActiveObject *aobj, const unsigned char *key, AZObject *value);

#ifdef __cplusplus
};
#endif

#endif

