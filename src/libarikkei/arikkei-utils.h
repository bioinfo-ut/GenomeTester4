#ifndef __ARIKKEI_UTILS_H__
#define __ARIKKEI_UTILS_H__

/*
 * Miscellaneous utilities
 *
 * Author:
 *   Lauris Kaplinski <lauris@kaplinski.com>
 *
 * This code is in public domain
 */

/* Basic compatibility */

#ifdef _WIN32
#define ARIKKEI_INLINE __inline
#else
#define ARIKKEI_INLINE static __inline__
#endif

/* Alignment */

#ifdef _WIN32
#define ARIKKEI_ALIGN_16 __declspec(align(16))
#else
#define ARIKKEI_ALIGN_16 __attribute__ ((aligned (16)))
#endif

/* Pointer arithmetic */
#define ARIKKEI_OFFSET(b,m) ((char *) &((b *) 0)->m - (char *) 0)
#define ARIKKEI_INT_TO_POINTER(v) (void *) ((char *) 0 + (v))
#define ARIKKEI_POINTER_TO_INT(p) ((int) ((char *) p - (char *) 0))

#define ARIKKEI_BASE_ADDRESS(klass,member,addr) ((char *) (addr) - ARIKKEI_OFFSET (klass, member))
#define ARIKKEI_MEMBER_ADDRESS(klass,member,addr) ((char *) (addr) + ARIKKEI_OFFSET (klass, member))

#ifdef __cplusplus
extern "C" {
#endif

#define arikkei_return_if_fail(expr) if (!(expr) && arikkei_emit_fail_warning ((const unsigned char *) __FILE__, __LINE__, (const unsigned char *) "?", (const unsigned char *) #expr)) return
#define arikkei_return_val_if_fail(expr,val) if (!(expr) && arikkei_emit_fail_warning ((const unsigned char *) __FILE__, __LINE__, (const unsigned char *) "?", (const unsigned char *) #expr)) return (val)

unsigned int arikkei_emit_fail_warning (const unsigned char *file, unsigned int line, const unsigned char *method, const unsigned char *expr);

#ifdef __cplusplus
};
#endif

#endif

