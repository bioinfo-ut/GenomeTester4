#define __AZ_VALUE_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <malloc.h>
#include <string.h>

#include <az/class.h>
#include <az/reference.h>
#include <az/string.h>

#include <az/value.h>

static AZClass *value_class = NULL;

void
az_init_value_class (void)
{
	value_class = az_class_new (AZ_TYPE_VALUE, AZ_TYPE_BLOCK, sizeof (AZClass), sizeof (AZValue), AZ_CLASS_IS_FINAL | AZ_CLASS_ZERO_MEMORY, "value");
	az_classes[AZ_TYPE_VALUE] = value_class;
}

void
az_value_clear (AZValue *value)
{
	if (value->impl && az_type_is_a (value->impl->type, AZ_TYPE_REFERENCE)) {
		if (value->reference) az_reference_unref ((AZReferenceClass *) value->impl, value->reference);
	}
	memset (value, 0, sizeof (AZValue));
}

void
az_value_set_reference (AZValue *value, unsigned int type, AZReference *ref)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[type]->implementation;
	value->reference = ref;
	if (ref) az_reference_ref ((AZReferenceClass *) value->impl, ref);
}

void
az_value_copy_indirect (AZValue *dst, const AZValue *src)
{
	if (dst == src) return;
	if (dst->impl && (dst->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (dst);
	if (src->impl && az_type_is_a (src->impl->type, AZ_TYPE_REFERENCE)) {
		az_value_set_reference (dst, src->impl->type, src->reference);
	} else {
		*dst = *src;
	}
}

#define AZ_VALUE_IS_NULL(v) (((v)->impl->type == AZ_TYPE_STRUCT) && ((v)->pvalue == NULL))

unsigned int
az_value_can_convert (unsigned int to, const AZValue *from)
{
	if (!to || !from->impl || !from->impl->type) return 0;
	if (to == from->impl->type) return 1;
	switch (to) {
	case AZ_TYPE_ANY:
		return 1;
	case AZ_TYPE_BOOLEAN:
		return (from->impl->type >= AZ_TYPE_BOOLEAN);
	case AZ_TYPE_INT8:
	case AZ_TYPE_UINT8:
	case AZ_TYPE_INT16:
	case AZ_TYPE_UINT16:
	case AZ_TYPE_INT32:
	case AZ_TYPE_UINT32:
	case AZ_TYPE_INT64:
	case AZ_TYPE_UINT64:
	case AZ_TYPE_FLOAT:
	case AZ_TYPE_DOUBLE:
		return AZ_TYPE_IS_ARITHMETIC (from->impl->type);
	case AZ_TYPE_POINTER:
		return (from->impl->type >= AZ_TYPE_POINTER);
	default:
		if (AZ_VALUE_IS_NULL (from)) return 1;
		if (az_type_is_a (from->impl->type, to)) return 1;
		break;
	}
	return 0;
}

unsigned int
az_value_convert (AZValue *dst, unsigned int type, const AZValue *from)
{
	unsigned int bvalue;
	int ivalue;
	long long lvalue;
	float fvalue;
	double dvalue;
	if (!type || !from->impl || !from->impl->type) return 0;
	switch (type) {
	case AZ_TYPE_NONE:
		az_value_clear (dst);
		return 1;
	case AZ_TYPE_ANY:
		az_value_copy (dst, from);
		return 1;
	case AZ_TYPE_BOOLEAN:
		if (from->impl->type == AZ_TYPE_BOOLEAN) {
			bvalue = from->bvalue;
		} else if ((from->impl->type >= AZ_TYPE_INT8) && (from->impl->type <= AZ_TYPE_UINT32)) {
			bvalue = from->ivalue != 0;
		} else if ((from->impl->type >= AZ_TYPE_INT64) && (from->impl->type <= AZ_TYPE_UINT64)) {
			bvalue = from->lvalue != 0;
		} else if (from->impl->type == AZ_TYPE_FLOAT) {
			bvalue = from->fvalue != 0;
		} else if (from->impl->type == AZ_TYPE_DOUBLE) {
			bvalue = from->dvalue != 0;
		} else {
			ivalue = from->pvalue != NULL;
			break;
		}
		az_value_clear (dst);
		dst->bvalue = bvalue;
		dst->impl = &az_classes[type]->implementation;
		return 1;
	case AZ_TYPE_INT8:
	case AZ_TYPE_UINT8:
	case AZ_TYPE_INT16:
	case AZ_TYPE_UINT16:
	case AZ_TYPE_INT32:
	case AZ_TYPE_UINT32:
		if ((from->impl->type >= AZ_TYPE_INT8) && (from->impl->type <= AZ_TYPE_UINT32)) {
			ivalue = from->ivalue;
		} else if ((from->impl->type >= AZ_TYPE_INT64) && (from->impl->type <= AZ_TYPE_UINT64)) {
			ivalue = (i32) from->lvalue;
		} else if (from->impl->type == AZ_TYPE_FLOAT) {
			ivalue = (i32) from->fvalue;
		} else if (from->impl->type == AZ_TYPE_DOUBLE) {
			ivalue = (i32) from->dvalue;
		} else {
			ivalue = 0;
			break;
		}
		az_value_clear (dst);
		dst->ivalue = ivalue;
		dst->impl = &az_classes[type]->implementation;
		return 1;
	case AZ_TYPE_INT64:
	case AZ_TYPE_UINT64:
		if ((from->impl->type >= AZ_TYPE_INT8) && (from->impl->type <= AZ_TYPE_UINT32)) {
			lvalue = from->ivalue;
		} else if ((from->impl->type >= AZ_TYPE_INT64) && (from->impl->type <= AZ_TYPE_UINT64)) {
			lvalue = (i64) from->lvalue;
		} else if (from->impl->type == AZ_TYPE_FLOAT) {
			lvalue = (i64) from->fvalue;
		} else if (from->impl->type == AZ_TYPE_DOUBLE) {
			lvalue = (i64) from->dvalue;
		} else {
			lvalue = 0;
			break;
		}
		az_value_clear (dst);
		dst->lvalue = lvalue;
		dst->impl = &az_classes[type]->implementation;
		return 1;
	case AZ_TYPE_FLOAT:
		if ((from->impl->type >= AZ_TYPE_INT8) && (from->impl->type <= AZ_TYPE_UINT32)) {
			fvalue = (f32) from->ivalue;
		} else if ((from->impl->type >= AZ_TYPE_INT64) && (from->impl->type <= AZ_TYPE_UINT64)) {
			fvalue = (f32) from->lvalue;
		} else if (from->impl->type == AZ_TYPE_FLOAT) {
			fvalue = from->fvalue;
		} else if (from->impl->type == AZ_TYPE_DOUBLE) {
			fvalue = (f32) from->dvalue;
		} else {
			fvalue = 0;
			break;
		}
		az_value_clear (dst);
		dst->fvalue = fvalue;
		dst->impl = &az_classes[type]->implementation;
		return 1;
	case AZ_TYPE_DOUBLE:
		if ((from->impl->type >= AZ_TYPE_INT8) && (from->impl->type <= AZ_TYPE_UINT32)) {
			dvalue = (f64) from->ivalue;
		} else if ((from->impl->type >= AZ_TYPE_INT64) && (from->impl->type <= AZ_TYPE_UINT64)) {
			dvalue = (f64) from->lvalue;
		} else if (from->impl->type == AZ_TYPE_FLOAT) {
			dvalue = (f64) from->fvalue;
		} else if (from->impl->type == AZ_TYPE_DOUBLE) {
			dvalue = from->dvalue;
		} else {
			dvalue = 0;
			break;
		}
		az_value_clear (dst);
		dst->dvalue = dvalue;
		dst->impl = &az_classes[type]->implementation;
		return 1;
	case AZ_TYPE_POINTER:
		if (from->impl->type >= AZ_TYPE_POINTER) {
			void *pvalue = from->pvalue;
			az_value_clear (dst);
			dst->pvalue = pvalue;
			dst->impl = &az_classes[type]->implementation;
		} else {
			break;
		}
		return 1;
	default:
		if (AZ_VALUE_IS_NULL (from)) {
			dst->impl = &az_classes[type]->implementation;
			dst->pvalue = NULL;
			return 1;
		} else if (az_type_is_a (from->impl->type, type)) {
			az_value_copy (dst, from);
			return 1;
		}
		break;
	}
	return 0;
}

void
az_value_set (AZValue *dst, unsigned int type, void *val)
{
	int ivalue;
	switch (type) {
	case AZ_TYPE_NONE:
		az_value_clear (dst);
		break;
	case AZ_TYPE_BOOLEAN:
		az_value_set_boolean (dst, ARIKKEI_POINTER_TO_INT (val));
		break;
	case AZ_TYPE_INT8:
		ivalue = ARIKKEI_POINTER_TO_INT (val);
		az_value_clear (dst);
		dst->impl = &az_classes[type]->implementation;
		dst->ivalue = (i8) ivalue;
		break;
	case AZ_TYPE_UINT8:
		ivalue = ARIKKEI_POINTER_TO_INT (val);
		az_value_clear (dst);
		dst->impl = &az_classes[type]->implementation;
		dst->ivalue = (u8) ivalue;
		break;
	case AZ_TYPE_INT16:
		ivalue = ARIKKEI_POINTER_TO_INT (val);
		az_value_clear (dst);
		dst->impl = &az_classes[type]->implementation;
		dst->ivalue = (i16) ivalue;
		break;
	case AZ_TYPE_UINT16:
		ivalue = ARIKKEI_POINTER_TO_INT (val);
		az_value_clear (dst);
		dst->impl = &az_classes[type]->implementation;
		dst->ivalue = (u16) ivalue;
		break;
	case AZ_TYPE_INT32:
		az_value_set_int (dst, AZ_TYPE_INT32, (i32) ARIKKEI_POINTER_TO_INT (val));
		break;
	case AZ_TYPE_UINT32:
		az_value_set_unsigned_int (dst, AZ_TYPE_UINT32, (u32) ARIKKEI_POINTER_TO_INT (val));
		break;
	case AZ_TYPE_INT64:
		/* fixme */
		az_value_set_i64 (dst, (i64) ARIKKEI_POINTER_TO_INT (val));
		break;
	case AZ_TYPE_UINT64:
		az_value_set_u64 (dst, (u64) ARIKKEI_POINTER_TO_INT (val));
		break;
	case AZ_TYPE_FLOAT:
		az_value_set_f32 (dst, *((f32 *) val));
		break;
	case AZ_TYPE_DOUBLE:
		az_value_set_f64 (dst, *((f64 *) val));
		break;
	case AZ_TYPE_POINTER:
	case AZ_TYPE_STRUCT:
	case AZ_TYPE_CLASS:
	case AZ_TYPE_INTERFACE:
		az_value_clear (dst);
		dst->impl = &az_classes[type]->implementation;
		dst->pvalue = val;
		break;
	case AZ_TYPE_STRING:
		az_value_set_string (dst, (AZString *) val);
		break;
	default:
		/* fixme: Make object primitive type? reference? */
		if (az_type_is_a (type, AZ_TYPE_REFERENCE)) {
			az_value_set_reference (dst, type, (AZReference *) val);
		} else if (az_type_is_a (type, AZ_TYPE_STRUCT)) {
			az_value_clear (dst);
			dst->impl = &az_classes[type]->implementation;
			dst->pvalue = val;
			break;
		} else {
			type = az_type_get_parent_primitive (type);
			az_value_set (dst, type, val);
		}
		break;
	}
}

void *
az_value_get_instance (AZValue *value)
{
	if (!value->impl || !value->impl->type) return NULL;
	if (az_type_is_a (value->impl->type, AZ_TYPE_REFERENCE)) {
		return value->reference;
	}
	switch (value->impl->type) {
	case AZ_TYPE_NONE:
	case AZ_TYPE_ANY:
		break;
	case AZ_TYPE_INT8:
	case AZ_TYPE_UINT8:
	case AZ_TYPE_INT16:
	case AZ_TYPE_UINT16:
	case AZ_TYPE_INT32:
	case AZ_TYPE_UINT32:
		return &value->ivalue;
	case AZ_TYPE_INT64:
	case AZ_TYPE_UINT64:
		return &value->lvalue;
	case AZ_TYPE_FLOAT:
		return &value->fvalue;
	case AZ_TYPE_DOUBLE:
		return &value->dvalue;
	default:
		return value->pvalue;
	}
	return NULL;
}

void
az_value_set_from_instance (AZValue *value, unsigned int type, const void *instance)
{
	if (value->impl && value->impl->type >= AZ_TYPE_REFERENCE) az_value_clear (value);
	if (!type) {
		value->impl = NULL;
		return;
	}
	value->impl = &az_classes[type]->implementation;
	switch (type) {
	case AZ_TYPE_NONE:
	case AZ_TYPE_ANY:
		break;
	case AZ_TYPE_INT8:
	case AZ_TYPE_UINT8:
		value->ivalue = *((i8 *) instance);
		break;
	case AZ_TYPE_INT16:
	case AZ_TYPE_UINT16:
		value->ivalue = *((i16 *) instance);
		break;
	case AZ_TYPE_INT32:
	case AZ_TYPE_UINT32:
		value->ivalue = *((i32 *) instance);
		break;
	case AZ_TYPE_INT64:
	case AZ_TYPE_UINT64:
		value->lvalue = *((i64 *) instance);
		break;
	case AZ_TYPE_FLOAT:
		value->fvalue = *((float *) instance);
		break;
	case AZ_TYPE_DOUBLE:
		value->dvalue = *((double *) instance);
		break;
	default:
		if (az_type_is_a (type, AZ_TYPE_REFERENCE)) {
			value->reference = (AZReference *) instance;
			az_reference_ref ((AZReferenceClass *) az_type_get_class (type), value->reference);
			break;
		} else if (az_type_is_a (type, AZ_TYPE_STRING)) {
			value->string = (AZString *) instance;
			az_string_ref (value->string);
			break;
		} else {
			value->pvalue = (void *) instance;
		}
	}
}

