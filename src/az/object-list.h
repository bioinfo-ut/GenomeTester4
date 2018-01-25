#ifndef __AZ_OBJECT_LIST_H__
#define __AZ_OBJECT_LIST_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

/*
 * An resizable array of AZObject either being or implementing certain type
 * If weakly referenced, elements have to be active objects
 */

typedef struct _AZObjectList AZObjectList;
typedef struct _AZObjectListClass AZObjectListClass;

#define AZ_TYPE_OBJECT_LIST (az_object_list_get_type ())

#include <az/list.h>
#include <az/object.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZObjectList {
	unsigned int type;
	unsigned int size;
	unsigned int length;
	unsigned int interface : 1;
	unsigned int weak : 1;
	AZObject **objects;
};

struct _AZObjectListClass {
	AZClass klass;
	AZListImplementation list_implementation;
};

unsigned int az_object_list_get_type (void);

void az_object_list_setup (AZObjectList *objl, unsigned int type, unsigned int weak);
void az_object_list_release (AZObjectList *objl);

AZObjectList *az_object_list_new (unsigned int type, unsigned int weak);
void az_object_list_delete (AZObjectList *objl);

void az_object_list_append_object (AZObjectList *objl, AZObject *object);
void az_object_list_remove_object (AZObjectList *objl, AZObject *object);
void az_object_list_remove_object_by_index (AZObjectList *objl, unsigned int idx);
void az_object_list_clear (AZObjectList *objl);

unsigned int az_object_list_contains (AZObjectList *objl, AZObject *object);

#ifdef __cplusplus
};
#endif

#endif
