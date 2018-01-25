#ifndef __AZ_BASE_H__
#define __AZ_BASE_H__

/*
 * A run-time type library
 *
 * Copyright (C) Lauris Kaplinski 2016
 */

#include <libarikkei/arikkei-utils.h>

#ifdef __cplusplus
extern "C" {
#endif
	
/*
 * Turn on runtime safety checks in library
 * Well-behaving implementation may disable these for speed
 */

#ifndef AZ_NO_SAFETY_CHECKS
#define AZ_SAFETY_CHECKS 1
#endif

#define _AZ_NO_STRING
#define _AZ_NO_VALUE
#define _AZ_NO_PROPERTIES

#ifndef AZ_NO_STRING
#define AZ_HAS_STRING
#endif

#ifndef AZ_NO_VALUE
#define AZ_HAS_VALUE
#endif

#ifndef AZ_NO_PROPERTIES
#define AZ_HAS_PROPERTIES
#endif

/* Basic primitives */

#ifndef AZ_NO_BASIC_PRIMITIVES
typedef char i8;
typedef unsigned char u8;
typedef short i16;
typedef unsigned short u16;
typedef int i32;
typedef unsigned int u32;
typedef long long i64;
typedef unsigned long long u64;
typedef float f32;
typedef double f64;
#endif

typedef char az_i8;
typedef unsigned char az_u8;
typedef short az_i16;
typedef unsigned short az_u16;
typedef int az_i32;
typedef unsigned int az_u32;
typedef long long az_i64;
typedef unsigned long long az_u64;
typedef float az_f32;
typedef double az_f64;

#define AZ_INLINE ARIKKEI_INLINE

/* Alignment */

#define AZ_ALIGN_16 ARIKKEI_ALIGN_16

/* Pointer arithmetic */
#define AZ_OFFSET ARIKKEI_OFFSET
#define AZ_INT_TO_POINTER ARIKKEI_INT_TO_POINTER
#define AZ_POINTER_TO_INT ARIKKEI_POINTER_TO_INT

#define AZ_BASE_ADDRESS ARIKKEI_BASE_ADDRESS
#define AZ_MEMBER_ADDRESS ARIKKEI_MEMBER_ADDRESS

#ifdef __cplusplus
};
#endif

#endif
