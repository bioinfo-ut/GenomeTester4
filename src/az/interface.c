#define __AZ_INTERFACE_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include <az/class.h>
#include <az/interface.h>

AZInterfaceClass *az_register_interface_type (unsigned int *type, unsigned int parent, const unsigned char *name,
	unsigned int class_size, unsigned int impl_size, unsigned int inst_size,
	void (*class_init) (AZClass *),
	void (*implementation_init) (AZImplementation *),
	void (*instance_init) (AZImplementation *, void *),
	void (*instance_finalize) (AZImplementation *, void *))
{
	AZInterfaceClass *if_klass;
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_val_if_fail (az_type_is_a (parent, AZ_TYPE_INTERFACE), NULL);
#endif
	if_klass = (AZInterfaceClass *) az_register_type (type, parent, name, class_size, inst_size, class_init, instance_init, instance_finalize);
	if_klass->implementation_size = impl_size;
	if_klass->implementation_init = implementation_init;
	return if_klass;
}

static AZInterfaceClass *interface_class = NULL;

void
az_init_interface_class (void)
{
	interface_class = (AZInterfaceClass *) az_class_new (AZ_TYPE_INTERFACE, AZ_TYPE_BLOCK, sizeof (AZInterfaceClass), 0, AZ_CLASS_IS_ABSTRACT, "interface");
	az_classes[AZ_TYPE_INTERFACE] = (AZClass *) interface_class;
	interface_class->implementation_size = sizeof (AZImplementation);
}
