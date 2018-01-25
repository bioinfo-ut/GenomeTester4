#define __AZ_FIELD_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <stdlib.h>

#include <libarikkei/arikkei-utils.h>

#include <az/string.h>
#include <az/field.h>

void
az_property_setup (AZField *prop, const unsigned char *key, unsigned int type, unsigned int id,
unsigned int is_static, unsigned int can_read, unsigned int can_write, unsigned int is_final, unsigned int is_value,
unsigned int value_type, void *value)
{
#ifdef AZ_SAFETY_CHECKS
	arikkei_return_if_fail (prop != NULL);
	arikkei_return_if_fail (key != NULL);
	arikkei_return_if_fail (type != AZ_TYPE_NONE);
	arikkei_return_if_fail (!(can_write && is_final));
	arikkei_return_if_fail (!is_value || is_static);
	arikkei_return_if_fail ((value_type == AZ_TYPE_NONE) || (az_type_is_assignable_to (value_type, type)));
#endif

	prop->key = az_string_new (key);
	prop->type = type;
	prop->id = id;
	prop->is_static = is_static;
	prop->can_read = can_read;
	prop->can_write = can_write;
	prop->is_final = is_final;
	prop->is_value = is_value;
	if (value_type != AZ_TYPE_NONE) {
		az_value_set (&prop->value, value_type, value);
	}
}
