#ifndef __AZ_STRING_H__
#define __AZ_STRING_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

#include <az/reference.h>

typedef struct _AZStringClass AZStringClass;

#ifdef __cplusplus
extern "C" {
#endif

struct _AZString {
	AZReference reference;
	unsigned int length;
	const unsigned char str[1];
};

struct _AZStringClass {
	AZReferenceClass reference_class;
};

#ifndef __AZ_STRING_C__
extern AZStringClass *string_class;
#else
AZStringClass *string_class = NULL;
#endif

AZString *az_string_new (const unsigned char *str);
AZString *az_string_new_length (const unsigned char *str, unsigned int length);

AZ_INLINE void
az_string_ref (AZString *astr)
{
	az_reference_ref (&string_class->reference_class, &astr->reference);
}

AZ_INLINE void
az_string_unref (AZString *astr)
{
	az_reference_unref (&string_class->reference_class, &astr->reference);
}

AZString *az_string_concat (AZString *lhs, AZString *rhs);

/* Library internal */
void az_init_string_class (void);

#ifdef __cplusplus
};
#endif

#endif
