#ifndef __AZ_SERIALIZATION_H__
#define __AZ_SERIALIZATION_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2017
*/

#include <az/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* All serialization methods accept NULL as destination */
/* Returns the total number of bytes that would have been written (regardless of destination length) */
unsigned int az_serialize_block (unsigned char *d, unsigned int dlen, const void *s, unsigned int slen);
/* Return the number of bytes consumed (0 on error) */
unsigned int az_deserialize_block (void *d, unsigned int dlen, const unsigned char *s, unsigned int slen);

/* Size is in bytes */
unsigned int az_serialize_int (unsigned char *d, unsigned int dlen, const void *inst, unsigned int size);
unsigned int az_deserialize_int (void *value, unsigned int size, const unsigned char *s, unsigned int slen);


#ifdef __cplusplus
};
#endif

#endif
