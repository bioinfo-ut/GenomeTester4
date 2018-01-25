#ifndef __AZ_INTERFACE_H__
#define __AZ_INTERFACE_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

typedef struct _AZInterfaceClass AZInterfaceClass;

#include <az/class.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZInterfaceClass {
	AZClass klass;
	/* Size of implementation */
	unsigned int implementation_size;
	/* Constructors and destructors */
	void (*implementation_init) (AZImplementation *impl);
};

/* Register new interface type */
AZInterfaceClass *az_register_interface_type (unsigned int *type, unsigned int parent, const unsigned char *name,
	unsigned int class_size, unsigned int implementation_size, unsigned int instance_size,
	void (*class_init) (AZClass *),
	void (*implementation_init) (AZImplementation *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *));

/* Library internal */
void az_init_interface_class (void);

#ifdef __cplusplus
};
#endif

#endif
