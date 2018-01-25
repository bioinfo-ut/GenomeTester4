#define __AZ_FUNCTION_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <az/value.h>

#include <az/function.h>

/* AZInstance implementation */
static void function_finalize (AZFunctionImplementation *impl, AZFunctionInstance *inst);

static unsigned int function_type = 0;

unsigned int
az_function_get_type (void)
{
	if (!function_type) {
		az_register_interface_type (&function_type, AZ_TYPE_INTERFACE, (const unsigned char *) "AZFunction",
			sizeof (AZFunctionClass), sizeof (AZFunctionImplementation), sizeof (AZFunctionInstance),
			NULL,
			NULL,
			NULL, (void (*) (AZImplementation *, void *)) function_finalize);
	}
	return function_type;
}

static void
function_finalize (AZFunctionImplementation *impl, AZFunctionInstance *func)
{
	if (func->arg_types) free (func->arg_types);
}

unsigned int
az_function_check_arguments (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *thisval, AZValue *args, unsigned int nargs, unsigned int *canconvert)
{
	unsigned int compatible, i;
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (az_type_is_a (impl->implementation.type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (inst != NULL, 0);
	if (nargs != inst->n_args) {
		if (canconvert) *canconvert = 0;
		return 0;
	}
	if (inst->this_type != AZ_TYPE_NONE) {
		if (!az_type_is_a (thisval->impl->type, inst->this_type)) {
			if (!canconvert) return 0;
			if (!az_value_can_convert (inst->this_type, thisval)) {
				*canconvert = 0;
				return 0;
			}
		}
	}
	compatible = 1;
	for (i = 0; i < inst->n_args; i++) {
		if (!az_type_is_a (args[i].impl->type, inst->arg_types[i])) {
			if (!canconvert) return 0;
			if (!az_value_can_convert (inst->arg_types[i], &args[i])) {
				*canconvert = 0;
				return 0;
			}
			compatible = 0;
		}
	}
	*canconvert = 1;
	return compatible;
}

unsigned int
az_function_convert_arguments (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *dst, AZValue *src)
{
	unsigned int i;
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (az_type_is_a (impl->implementation.type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (inst != NULL, 0);
	for (i = 0; i < inst->n_args; i++) {
		if (!az_value_convert (&dst[i], inst->arg_types[i], &src[i])) return 0;
	}
	return 1;
}

unsigned int
az_function_invoke (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *thisval, AZValue *retval, AZValue *args, unsigned int checktypes)
{
	arikkei_return_val_if_fail (impl != NULL, 0);
	arikkei_return_val_if_fail (az_type_is_a (impl->implementation.type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (inst != NULL, 0);
	if (checktypes) {
		unsigned int i;
		if (inst->this_type != AZ_TYPE_NONE) {
			arikkei_return_val_if_fail (az_type_is_a (thisval->impl->type, inst->this_type), 0);
		}
		arikkei_return_val_if_fail ((inst->ret_type == AZ_TYPE_NONE) || (retval != NULL), 0);
		for (i = 0; i < inst->n_args; i++) {
			arikkei_return_val_if_fail (az_type_is_a (args[i].impl->type, inst->arg_types[i]), 0);
		}
	}
	if (impl->invoke) {
		return impl->invoke (impl, inst, thisval, retval, args);
	}
	return 0;
}

unsigned int
az_function_invoke_direct (AZFunctionImplementation *implementation, AZFunctionInstance *instance, AZValue *thisval, AZValue *retval, ...)
{
	AZValue *vals;
	va_list ap;
	unsigned int result, i;
	arikkei_return_val_if_fail (implementation != NULL, 0);
	arikkei_return_val_if_fail (az_type_is_a (implementation->implementation.type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (instance != NULL, 0);
	if (instance->this_type != AZ_TYPE_NONE) {
		arikkei_return_val_if_fail (az_type_is_a (thisval->impl->type, instance->this_type), 0);
	}
	arikkei_return_val_if_fail ((instance->ret_type == AZ_TYPE_NONE) || (retval != NULL), 0);
	vals = (AZValue *) malloc (instance->n_args * sizeof (AZValue));
	memset (vals, 0, instance->n_args * sizeof (AZValue));
	va_start (ap, instance->n_args);
	for (i = 0; i < instance->n_args; i++) {
		vals[i].impl->type = instance->arg_types[i];
		switch (instance->arg_types[i]) {
		case AZ_TYPE_NONE:
			break;
		case AZ_TYPE_BOOLEAN:
			vals[i].bvalue = va_arg (ap, unsigned int);
			break;
		case AZ_TYPE_INT8:
		case AZ_TYPE_UINT8:
		case AZ_TYPE_INT16:
		case AZ_TYPE_UINT16:
		case AZ_TYPE_INT32:
		case AZ_TYPE_UINT32:
			vals[i].bvalue = va_arg (ap, int);
			break;
		case AZ_TYPE_INT64:
		case AZ_TYPE_UINT64:
			vals[i].lvalue = va_arg (ap, long long);
			break;
		case AZ_TYPE_FLOAT:
			vals[i].fvalue = va_arg (ap, float);
			break;
		case AZ_TYPE_DOUBLE:
			vals[i].dvalue = va_arg (ap, double);
			break;
		case AZ_TYPE_POINTER:
			vals[i].pvalue = va_arg (ap, void *);
			break;
		case AZ_TYPE_STRING:
			az_value_set_string (&vals[i], (AZString *) va_arg (ap, void *));
			break;
		default:
			az_value_set_reference (&vals[i], instance->arg_types[i], (AZReference *) va_arg (ap, void *));
			break;
		}
	}
	va_end (ap);
	result = az_function_invoke (implementation, instance, thisval, retval, vals, 0);
	free (vals);
	return result;
}

unsigned int
az_function_invoke_by_type_instance (unsigned int type, void *instance, AZValue *thisval, AZValue *retval, AZValue *args, unsigned int checktypes)
{
	AZClass *klass;
	AZFunctionImplementation *impl;
	AZFunctionInstance *inst;
	arikkei_return_val_if_fail (az_type_implements (type, AZ_TYPE_FUNCTION), 0);
	arikkei_return_val_if_fail (instance != NULL, 0);
	klass = az_type_get_class (type);
	impl = (AZFunctionImplementation *) az_get_interface (&klass->implementation, instance, AZ_TYPE_FUNCTION, (void **) &inst);
	return az_function_invoke (impl, inst, thisval, retval, args, checktypes);
}
