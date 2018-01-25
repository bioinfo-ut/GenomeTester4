#ifndef __AZ_FUNCTION_OBJECT_H__
#define __AZ_FUNCTION_OBJECT_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_FUNCTION_OBJECT az_function_object_get_type ()
#define AZ_FUNCTION_OBJECT(o) (AZ_CHECK_INSTANCE_CAST ((o), AZ_TYPE_FUNCTION_OBJECT, AZFunctionObject))
#define AZ_IS_FUNCTION_OBJECT(o) (AZ_CHECK_INSTANCE_TYPE ((o), AZ_TYPE_FUNCTION_OBJECT))

typedef struct _AZFunctionObject AZFunctionObject;
typedef struct _AZFunctionObjectClass AZFunctionObjectClass;

#include <az/function.h>
#include <az/object.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZFunctionObject {
	AZObject object;
	AZFunctionInstance function_instance;
	unsigned int (*call) (AZValue *, AZValue *, AZValue *);
};

struct _AZFunctionObjectClass {
	AZObjectClass parent_class;
	AZFunctionImplementation function_implementation;
};

unsigned int az_function_object_get_type (void);

AZFunctionObject *az_function_object_new (unsigned int thistype, unsigned int rettype, unsigned int nargs, const unsigned int argtypes[], unsigned int (*call) (AZValue *, AZValue *, AZValue *));

unsigned int az_function_object_invoke (AZFunctionObject *fobj, AZValue *thisval, AZValue *retval, AZValue *args, unsigned int checktypes);

#ifdef __cplusplus
};
#endif

#endif

