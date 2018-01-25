#ifndef __AZ_REFERENCE_H__
#define __AZ_REFERENCE_H__

/*
* A run-time type library
*
* Copyright (C) Lauris Kaplinski 2016
*/

typedef struct _AZReferenceClass AZReferenceClass;

#include <az/class.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _AZReference {
	unsigned int refcount;
};

struct _AZReferenceClass {
	AZClass klass;
	/* Invoked when the last references to it is about to be dropped */
	/* Returning false means that reference is transferred to another holder */
	unsigned int (*drop) (AZReferenceClass *klass, AZReference *ref);
	/* Dispose should release resources and references hold by this object */
	/* It is called automatically by unref before finalizing object */
	/* Subclasses must ensure that it is safe to call it multiple times */
	void (*dispose) (AZReferenceClass *klass, AZReference *ref);
};

/* For library internal use */
void az_reference_unref_last (AZReferenceClass *klass, AZReference *ref);

AZ_INLINE void
az_reference_ref (AZReferenceClass *klass, AZReference *ref)
{
	ref->refcount += 1;
}

AZ_INLINE void
az_reference_unref (AZReferenceClass *klass, AZReference *ref)
{
	if (ref->refcount > 1) {
		ref->refcount -= 1;
	} else {
		az_reference_unref_last (klass, ref);
	}
}

/* Dispose drops a reference unconditionally (does not allow drop to grab the last one) */
void az_reference_dispose (AZReferenceClass *klass, AZReference *ref);

/* Library internal */
void az_init_reference_class (void);

#ifdef __cplusplus
};
#endif

#endif
