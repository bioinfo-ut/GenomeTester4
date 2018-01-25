#define __AZ_STRING_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <az/class.h>
#include <az/serialization.h>
#include <az/string.h>

static unsigned int
string_to_string (AZClass *klass, void *instance, unsigned char *buf, unsigned int len)
{
	unsigned int pos = 0;
	if (instance) {
		AZString *str = (AZString *) instance;
		unsigned int slen = (str->length > len) ? len : str->length;
		memcpy (buf + pos, str->str, slen);
		pos += slen;
	}
	if (pos < len) buf[pos] = 0;
	return pos;
}

static unsigned int
serialize_string (AZImplementation *impl, void *inst, unsigned char *d, unsigned int dlen, AZContext *ctx)
{
	/* fixme: */
	AZString *str = (AZString *) inst;
	if (!str) {
		unsigned int len = 0;
		return az_serialize_int (d, dlen, &len, 4);
	} else {
		if ((4 + str->length) <= dlen) {
			az_serialize_int (d, dlen, &str->length, 4);
			az_serialize_block (d + 4, dlen - 4, str->str, str->length);
		}
		return 4 + str->length;
	}
}

static unsigned int
deserialize_string (AZImplementation *impl, void *value, const unsigned char *s, unsigned int slen, AZContext *ctx)
{
	AZString **str = (AZString **) value;
	unsigned int len;
	if (slen < 4) {
		*str = NULL;
		return 0;
	}
	az_deserialize_int (&len, 4, s, slen);
	if (!len || ((4 + len) > slen)) {
		*str = NULL;
		return 0;
	}
	*str = az_string_new_length (s + 4, len);
	return 4 + len;
}

void
az_init_string_class (void)
{
	string_class = (AZStringClass *) az_class_new (AZ_TYPE_STRING, AZ_TYPE_REFERENCE, sizeof (AZStringClass), 0, AZ_CLASS_IS_FINAL, "string");
	az_classes[AZ_TYPE_STRING] = (AZClass *) string_class;
	az_classes[AZ_TYPE_STRING]->serialize = serialize_string;
	az_classes[AZ_TYPE_STRING]->deserialize = deserialize_string;
	az_classes[AZ_TYPE_STRING]->to_string = string_to_string;
}

AZString *
az_string_new (const unsigned char *str)
{
	if (!str) return NULL;
	return az_string_new_length (str, (unsigned int) strlen ((const char *) str));
}

AZString *
az_string_new_length (const unsigned char *str, unsigned int length)
{
	AZString *astr = (AZString *) malloc (sizeof (AZString) + length);
	az_instance_init (astr, AZ_TYPE_STRING);
	astr->length = length;
	memcpy ((unsigned char *) astr->str, str, length);
	((unsigned char *) astr->str)[length] = 0;
	return astr;
}

AZString *
az_string_concat (AZString *lhs, AZString *rhs)
{
	AZString *astr;
	if (!lhs) return rhs;
	if (!rhs) return lhs;
	astr = (AZString *) malloc (sizeof (AZString) + lhs->length + rhs->length);
	az_instance_init (astr, AZ_TYPE_STRING);
	astr->length = lhs->length + rhs->length;
	if (lhs->length) memcpy ((unsigned char *) astr->str, lhs->str, lhs->length);
	if (rhs->length) memcpy ((unsigned char *) astr->str + lhs->length, rhs->str, rhs->length);
	((unsigned char *) astr->str)[lhs->length + rhs->length] = 0;
	return astr;
}
