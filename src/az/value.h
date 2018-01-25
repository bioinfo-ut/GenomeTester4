#ifndef __AZ_VALUE_H__
#define __AZ_VALUE_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <az/class.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZValue {
	AZImplementation *impl;
	union {
		unsigned int bvalue;
		i32 ivalue;
		u32 uvalue;
		i64 lvalue;
		f32 fvalue;
		f64 dvalue;
		void *pvalue;
		AZReference *reference;
		AZString *string;
	};
};

void az_value_clear (AZValue *value);

AZ_INLINE void
az_value_set_boolean (AZValue *value, unsigned int val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_BOOLEAN]->implementation;
	value->bvalue = val;
}

AZ_INLINE void
az_value_set_int (AZValue *value, unsigned int type, int val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[type]->implementation;
	value->ivalue = val;
}

AZ_INLINE void
az_value_set_unsigned_int (AZValue *value, unsigned int type, unsigned int val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[type]->implementation;
	value->uvalue = val;
}

AZ_INLINE void
az_value_set_i64 (AZValue *value, i64 val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_INT64]->implementation;
	value->lvalue = val;
}

AZ_INLINE void
az_value_set_u64 (AZValue *value, u64 val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_UINT64]->implementation;
	value->lvalue = (i64) val;
}

AZ_INLINE void
az_value_set_f32 (AZValue *value, f32 val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_FLOAT]->implementation;
	value->fvalue = val;
}

ARIKKEI_INLINE void
az_value_set_f64 (AZValue *value, f64 val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_DOUBLE]->implementation;
	value->dvalue = val;
}

ARIKKEI_INLINE void
az_value_set_pointer (AZValue *value, const void *val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_POINTER]->implementation;
	value->pvalue = (void *) val;
}

ARIKKEI_INLINE void
az_value_set_class (AZValue *value, AZClass *val)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_CLASS]->implementation;
	value->pvalue = val;
}

void az_value_set_reference (AZValue *value, unsigned int type, AZReference *ref);
#define az_value_set_string(v, s) az_value_set_reference (v, AZ_TYPE_STRING, (AZReference *) (s))

AZ_INLINE void
az_value_transfer_reference (AZValue *value, unsigned int type, AZReference *ref)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[type]->implementation;
	value->reference = ref;
}

AZ_INLINE void
az_value_transfer_string (AZValue *value, AZString *str)
{
	if (value->impl && (value->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (value);
	value->impl = &az_classes[AZ_TYPE_STRING]->implementation;
	value->string = str;
}

void az_value_copy_indirect (AZValue *dst, const AZValue *src);

AZ_INLINE void
az_value_copy (AZValue *dst, const AZValue *src)
{
	if (dst == src) return;
	if (dst->impl && (dst->impl->type >= AZ_TYPE_REFERENCE)) az_value_clear (dst);
	if (!src->impl || (src->impl->type < AZ_TYPE_REFERENCE)) {
		*dst = *src;
	} else {
		az_value_copy_indirect (dst, src);
	}
}

unsigned int az_value_can_convert (unsigned int to, const AZValue *from);
unsigned int az_value_convert (AZValue *dst, unsigned int type, const AZValue *from);

void az_value_set (AZValue *dst, unsigned int type, void *val);
void *az_value_get_instance (AZValue *value);
void az_value_set_from_instance (AZValue *value, unsigned int type, const void *instance);

AZ_INLINE unsigned int
az_value_is_a (const AZValue *val, unsigned int type)
{
	if (val->impl) return az_type_is_a (val->impl->type, type);
	return !type;
}

AZ_INLINE unsigned int
az_value_implements (const AZValue *val, unsigned int type)
{
	if (val->impl) return az_type_implements (val->impl->type, type);
	return 0;
}

/* Library internal */
void az_init_value_class (void);

#ifdef __cplusplus
};
#endif


#endif
