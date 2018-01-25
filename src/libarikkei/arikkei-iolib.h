#ifndef __ARIKKEI_IOLIB_H__
#define __ARIKKEI_IOLIB_H__

/*
 * Arikkei
 *
 * Basic datatypes and code snippets
 *
 * Author:
 *   Lauris Kaplinski <lauris@kaplinski.com>
 *
 * This code is in public domain
 *
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// All file names are UTF-8

const unsigned char *arikkei_mmap (const unsigned char *filename, size_t *size, const unsigned char *name);
void arikkei_munmap (const unsigned char *buffer, size_t size);

FILE *arikkei_fopen (const unsigned char *filename, const unsigned char *mode);

#ifdef __cplusplus
}
#endif

#endif
