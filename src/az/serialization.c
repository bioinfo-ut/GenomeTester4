#define __AZ_SERIALIZATION_C__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2017
*/

#include <stdio.h>
#include <string.h>

#include <az/serialization.h>

unsigned int
az_serialize_block (unsigned char *d, unsigned int dlen, const void *s, unsigned int slen)
{
	if ((slen <= dlen) && d) memcpy (d, s, slen);
	return slen;
}

unsigned int
az_deserialize_block (void *d, unsigned int dlen, const unsigned char *s, unsigned int slen)
{
	if (dlen > slen) return 0;
	memcpy (d, s, dlen);
	return dlen;
}

unsigned int
az_serialize_int (unsigned char *d, unsigned int dlen, const void *inst, unsigned int size)
{
	if ((dlen >= size) && d) {
		unsigned int i;
		for (i = 0; i < size; i++) {
			d[size - 1 - i] = *((unsigned char *) inst + i);
		}
	}
	return size;
}

unsigned int
az_deserialize_int (void *value, unsigned int size, const unsigned char *s, unsigned int slen)
{
	unsigned int i;
	if (slen < size) return 0;
	for (i = 0; i < size; i++) {
		*((unsigned char *) value + i) = s[size - 1 - i];
	}
	return size;
}
