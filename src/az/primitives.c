#define __AZ_PRIMITIVES_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <assert.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>

#include <az/class.h>
#include <az/serialization.h>
#include <libarikkei/arikkei-strlib.h>

#include <az/primitives.h>

/* 0 None */

/* 1 Any */

static unsigned int
any_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	if (klass->implementation.type == AZ_TYPE_ANY) {
		/* Pure Any */
		return arikkei_strncpy (buf, len, (const unsigned char *) "Any");
	} else {
		/* Subclass that does not implement to_string */
		char c[32];
		unsigned int pos;
		pos = arikkei_memcpy_str (buf, len, (const unsigned char *) "Instance of ");
		pos += arikkei_memcpy_str (buf + pos, len - pos, klass->name);
		pos += arikkei_memcpy_str (buf + pos, len - pos, (const unsigned char *) " (");
		sprintf (c, "%p", instance);
		pos += arikkei_memcpy_str (buf + pos, len - pos, (const unsigned char *) c);
		pos += arikkei_strncpy (buf + pos, len - pos, (const unsigned char *) ")");
		return pos;
	}
}

/* 2 Boolean */

static unsigned int
serialize_boolean (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	unsigned char v = (*((unsigned int *) inst) != 0);
	return az_serialize_int (d, dlen, &v, 1);
}

static unsigned int
deserialize_boolean (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	if (!slen) return 0;
	*((unsigned int *) value) = *s;
	return 1;
}

static unsigned int
boolean_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	return arikkei_strncpy (buf, len, (*((unsigned int *) instance)) ? (const unsigned char *) "True" : (const unsigned char *) "False");
}

/* 3 Int8 */

static unsigned int
serialize_int (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	AZClass *klass = az_type_get_class (impl->type);
	return az_serialize_int (d, dlen, inst, klass->instance_size);
}

static unsigned int
deserialize_int (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	AZClass *klass = az_type_get_class (impl->type);
	return az_deserialize_int (value, klass->instance_size, s, slen);
}

static unsigned int
copy_int_to_buffer (unsigned char *dst, unsigned int dlen, unsigned long long value, unsigned int sign)
{
	unsigned char c[32];
	unsigned int clen = 0, len = 0;
	c[clen++] = '0' + value % 10;
	value /= 10;
	while (value) {
		c[clen++] = '0' + value % 10;
		value /= 10;
	}
	if (sign) {
		if (dst && (len < dlen)) dst[len++] = '-';
	}
	if (dst) {
		unsigned int i;
		for (i = 0; i < clen; i++) {
			if (len < dlen) dst[len++] = c[clen - 1 - i];
		}
		if (len < dlen) dst[len++] = 0;
	}
	return clen + sign + 1;
}

static unsigned int
int_to_string_any (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned int size = klass->instance_size;
	unsigned int is_signed = ((klass->implementation.type & 1) != 0);
	unsigned long long value = 0;
	unsigned int sign = 0;
	if (size == 1) {
		value = *((unsigned char *) instance);
	} else if (size == 2) {
		value = *((unsigned short *) instance);
	}
	memcpy (&value, instance, size);
	if (is_signed) {
		if (value & (1ULL << ((8 * size) - 1))) {
			sign = 1;
			value = (1ULL << (8 * size)) - value;
		}
	}
	return copy_int_to_buffer (buf, len, value, sign);
}

/* 4 Uint8 */

/* 5 Int16 */

/* 6 Uint16 */

/* 7 Int32 */

/* 8 Uint32 */

/* 9 Int64 */

/* 10 Uint64 */

/* 11 Float */

static unsigned int
float_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned char c[32];
	unsigned int clen = arikkei_dtoa_exp (c, 32, *((float *) instance), 5, 0);
	if (buf) {
		memcpy (buf, c, clen);
		if (clen < len) buf[clen] = 0;
	}
	return clen + 1;
}

/* 12 Double */

static unsigned int
double_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned char c[32];
	unsigned int clen = arikkei_dtoa_exp (c, 32, *((double *) instance), 8, 0);
	if (buf) {
		memcpy (buf, c, clen);
		if (clen < len) buf[clen] = 0;
	}
	return clen + 1;
}

/* 13 Complex float */

static unsigned int
complex_float_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned char c[64];
	unsigned int clen = arikkei_dtoa_exp (c, 32, *((float *) instance), 5, 0);
	c[clen++] = '+';
	clen += arikkei_dtoa_exp (c, 32, *((float *) instance + 1), 5, 0);
	c[clen++] = 'i';
	if (buf) {
		memcpy (buf, c, clen);
		if (clen < len) buf[clen] = 0;
	}
	return clen + 1;
}

static unsigned int
serialize_complex_float (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	if (d && (dlen >= 8)) {
		az_serialize_int (d, dlen, inst, 4);
		az_serialize_int (d + 4, dlen - 4, (float *) inst + 1, 4);
	}
	return 8;
}

static unsigned int
deserialize_complex_float (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	if (slen < 8) return 0;
	az_deserialize_int (value, 4, s, slen);
	az_deserialize_int ((float *) value + 1, 4, s + 4, slen - 4);
	return 8;
}

/* 14 Complex double */

static unsigned int
complex_double_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned char c[64];
	unsigned int clen = arikkei_dtoa_exp (c, 32, *((double *) instance), 8, 0);
	c[clen++] = '+';
	clen += arikkei_dtoa_exp (c, 32, *((double *) instance + 1), 8, 0);
	c[clen++] = 'i';
	if (buf) {
		memcpy (buf, c, clen);
		if (clen < len) buf[clen] = 0;
	}
	return clen + 1;
}

static unsigned int
serialize_complex_double (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	if (d && (dlen >= 16)) {
		az_serialize_int (d, dlen, inst, 8);
		az_serialize_int (d + 8, dlen - 8, (double *) inst + 1, 8);
	}
	return 16;
}

static unsigned int
deserialize_complex_double (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	if (slen < 16) return 0;
	az_deserialize_int (value, 8, s, slen);
	az_deserialize_int ((double *) value + 1, 8, s + 8, slen - 8);
	return 16;
}

/* 15 Pointer */

static unsigned int
pointer_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned int i;
	const char *t;
	char b[32];
	if (instance) {
		static const char *c = "0123456789abcdef";
		unsigned long long v = (unsigned long long) instance;
		unsigned int l = 2 * sizeof (void *);
		for (i = 0; i < l; i++) {
			b[l - 1 - i] = c[v & 0xf];
			v = v >> 4;
		}
		b[i] = 0;
		t = b;
	} else {
		t = "Null";
	}
	return arikkei_strncpy (buf, len, (const unsigned char *) t);
}

struct _PrimitiveDef {
	unsigned int type;
	unsigned int parent;
	unsigned int is_abstract;
	unsigned int is_final;
	unsigned int is_value;
	unsigned int instance_size;
	const char *name;
	unsigned int (*serialize) (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx);
	unsigned int (*deserialize) (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx);
	unsigned int (*to_string) (AZClass *klass, void *instance, unsigned char *buf, unsigned int len);
};

struct _PrimitiveDef defs[] = {
	{ AZ_TYPE_NONE, AZ_TYPE_NONE, 0, 0, 0, 0, NULL, NULL, NULL, NULL },
	{ AZ_TYPE_ANY, AZ_TYPE_NONE, 1, 0, 0, 0, "any", NULL, NULL, any_to_string },
	{ AZ_TYPE_BOOLEAN, AZ_TYPE_ANY, 0, 1, 1, 4, "boolean", serialize_boolean, deserialize_boolean, boolean_to_string },
	{ AZ_TYPE_INT8, AZ_TYPE_ANY, 0, 1, 1, 1, "int8", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_UINT8, AZ_TYPE_ANY, 0, 1, 1, 1, "uint8", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_INT16, AZ_TYPE_ANY, 0, 1, 1, 2, "int16", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_UINT16, AZ_TYPE_ANY, 0, 1, 1, 2, "uint16", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_INT32, AZ_TYPE_ANY, 0, 1, 1, 4, "int32", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_UINT32, AZ_TYPE_ANY, 0, 1, 1, 4, "uint32", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_INT64, AZ_TYPE_ANY, 0, 1, 1, 8, "int64", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_UINT64, AZ_TYPE_ANY, 0, 1, 1, 8, "uint64", serialize_int, deserialize_int, int_to_string_any },
	{ AZ_TYPE_FLOAT, AZ_TYPE_ANY, 0, 1, 1, 4, "float", serialize_int, deserialize_int, float_to_string },
	{ AZ_TYPE_DOUBLE, AZ_TYPE_ANY, 0, 1, 1, 8, "double", serialize_int, deserialize_int, double_to_string },
	{ AZ_TYPE_COMPLEX_FLOAT, AZ_TYPE_ANY, 0, 1, 1, 8, "complex float", serialize_complex_float, deserialize_complex_float, complex_float_to_string },
	{ AZ_TYPE_COMPLEX_DOUBLE, AZ_TYPE_ANY, 0, 1, 1, 16, "complex double", serialize_complex_double, deserialize_complex_double, complex_double_to_string },
	{ AZ_TYPE_POINTER, AZ_TYPE_ANY, 0, 0, 1, 8, "pointer", serialize_int, deserialize_int, pointer_to_string },
	{ AZ_TYPE_STRUCT, AZ_TYPE_ANY, 1, 0, 1, 0, "struct", NULL, NULL, NULL },
	{ AZ_TYPE_BLOCK, AZ_TYPE_ANY, 1, 0, 0, 0, "block", NULL, NULL, NULL }
};

void
az_init_primitive_classes (void)
{
	unsigned int i;
	az_classes[AZ_TYPE_NONE] = NULL;
	for (i = AZ_TYPE_ANY; i <= AZ_TYPE_BLOCK; i++) {
		unsigned int flags = 0;
		assert (defs[i].type == i);
		if (defs[i].is_abstract) flags |= AZ_CLASS_IS_ABSTRACT;
		if (defs[i].is_final) flags |= AZ_CLASS_IS_FINAL;
		if (defs[i].is_value) flags |= AZ_CLASS_IS_VALUE;
		az_classes[i] = az_class_new (defs[i].type, defs[i].parent, sizeof (AZClass), defs[i].instance_size, flags, defs[i].name);
		az_classes[i]->serialize = defs[i].serialize;
		az_classes[i]->deserialize = defs[i].deserialize;
		az_classes[i]->to_string = defs[i].to_string;
	}
}
