#ifndef __AZ_FUNCTION_H__
#define __AZ_FUNCTION_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#define AZ_TYPE_FUNCTION (az_function_get_type ())

typedef struct _AZFunctionClass AZFunctionClass;
typedef struct _AZFunctionImplementation AZFunctionImplementation;
typedef struct _AZFunctionInstance AZFunctionInstance;

#include <az/interface.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZFunctionInstance {
	/* This object type */
	unsigned int this_type;
	/* Arguments */
	unsigned int n_args;
	unsigned int *arg_types;
	/* Return type */
	unsigned int ret_type;
};

struct _AZFunctionClass {
	AZInterfaceClass parent_class;
};

struct _AZFunctionImplementation {
	AZImplementation implementation;
		/* All values can be NULL if not present in signature */
	unsigned int (*invoke) (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *this_val, AZValue *ret_val, AZValue *args);
};

unsigned int az_function_get_type (void);

/* Return true if arguments can be submitted as is */
unsigned int az_function_check_arguments (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *this_val, AZValue *args, unsigned int n_args, unsigned int *can_convert);
unsigned int az_function_convert_arguments (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *dst, AZValue *src);

unsigned int az_function_invoke (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *this_val, AZValue *ret_val, AZValue *args, unsigned int check_types);
unsigned int az_function_invoke_direct (AZFunctionImplementation *impl, AZFunctionInstance *inst, AZValue *this_val, AZValue *ret_val, ...);

/* Helper */
unsigned int az_function_invoke_by_type_instance (unsigned int type, void *inst, AZValue *this_val, AZValue *ret_val, AZValue *args, unsigned int check_types);

#ifdef __cplusplus
};
#endif

#endif
