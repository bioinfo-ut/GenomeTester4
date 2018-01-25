#define __AZ_FUNCTION_OBJECT_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>
#include <string.h>

#include <az/class.h>
#include <az/function-object.h>

/* AZInstance implementation */
static void function_object_class_init (AZFunctionObjectClass *klass);
/* AZFunction implementation */
static unsigned int function_object_invoke (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *thisval, AZValue *retval, AZValue *args);

static unsigned int function_object_type = 0;
static AZFunctionObjectClass *function_object_class;
static AZObjectClass *parent_class;

unsigned int
az_function_object_get_type (void)
{
	if (!function_object_type) {
		function_object_class = (AZFunctionObjectClass *) az_register_type (&function_object_type, AZ_TYPE_OBJECT,
			(const unsigned char *) "AZFunctionObject",
			sizeof (AZFunctionObjectClass),
			sizeof (AZFunctionObject),
			(void (*) (AZClass *)) function_object_class_init, NULL, NULL);
	}
	return function_object_type;
}

static void
function_object_class_init (AZFunctionObjectClass *klass)
{
	parent_class = (AZObjectClass *) ((AZClass *) klass)->parent;
	az_class_set_num_interfaces ((AZClass *) klass, 1);
	az_class_declare_interface ((AZClass *) klass, 0, AZ_TYPE_FUNCTION, ARIKKEI_OFFSET (AZFunctionObjectClass, function_implementation), ARIKKEI_OFFSET (AZFunctionObject, function_instance));
	klass->function_implementation.invoke = function_object_invoke;
}

static unsigned int
function_object_invoke (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *thisval, AZValue *retval, AZValue *args)
{
	AZFunctionObject *fobj;
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (az_type_is_a (impl->implementation.type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (inst != NULL, 0);
	fobj = (AZFunctionObject *) ARIKKEI_BASE_ADDRESS (AZFunctionObject, function_instance, inst);
	arikkei_return_val_if_fail (az_object_is_a (fobj, AZ_TYPE_FUNCTION_OBJECT), 0);
	if (fobj->call) {
		return fobj->call (thisval, retval, args);
	}
	return 0;
}

AZFunctionObject *
az_function_object_new (unsigned int thistype, unsigned int rettype, unsigned int nargs, const unsigned int argtypes[], unsigned int (*call) (AZValue *, AZValue *, AZValue *))
{
	AZFunctionObject *fobj;
	fobj = (AZFunctionObject *) az_object_new (AZ_TYPE_FUNCTION_OBJECT);
	fobj->function_instance.this_type = thistype;
	fobj->function_instance.ret_type = rettype;
	if (nargs) {
		fobj->function_instance.n_args = nargs;
		fobj->function_instance.arg_types = (unsigned int *) malloc (nargs * sizeof (unsigned int));
		memcpy (fobj->function_instance.arg_types, argtypes, nargs * sizeof (unsigned int));
	}
	fobj->call = call;
	return fobj;
}

unsigned int
az_function_object_invoke (AZFunctionObject *fobj, AZValue *thisval, AZValue *retval, AZValue *args, unsigned int checktypes)
{
	arikkei_return_val_if_fail (fobj != NULL, 0);
	arikkei_return_val_if_fail (AZ_IS_FUNCTION_OBJECT (fobj), 0);
	return az_function_invoke (&((AZFunctionObjectClass *) fobj->object.klass)->function_implementation, &fobj->function_instance, thisval, retval, args, checktypes);
}

