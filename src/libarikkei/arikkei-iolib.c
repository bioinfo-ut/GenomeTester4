#define __ARIKKEI_IOLIB_C__

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

// Disable VS2005 nagging
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#ifndef _CRT_NON_CONFORMING_SWPRINTFS
#define _CRT_NON_CONFORMING_SWPRINTFS 1
#endif

#ifndef _WIN32
#include <unistd.h>
#include <sys/mman.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#include <stdio.h>
#endif

#include "arikkei-strlib.h"
#include "arikkei-dict.h"

#include "arikkei-iolib.h"

#ifdef _DEBUG
//#define PRINT_MAPSIZE 1
#endif

#ifdef PRINT_MAPSIZE
static size_t total = 0;
#endif

#ifdef _WIN32

const unsigned char *
arikkei_mmap (const unsigned char *filename, size_t *size, const unsigned char *mapname)
{
	unsigned short *ucs2filename, *ucs2mapname;
	unsigned char *cdata;
	struct _stat st;
	HANDLE fh, mh;

	if (!filename || !*filename) return NULL;

	ucs2filename = arikkei_utf8_ucs2_strdup (filename);

	ucs2mapname = (mapname && *mapname) ? arikkei_utf8_ucs2_strdup (mapname) : NULL;

	if (_wstat (ucs2filename, &st)) {
		/* No such file */
		/* fprintf (stderr, "arikkei_mmap: File %s not found or not regular file\n", filename); */
		free (ucs2filename);
		if (ucs2mapname) free (ucs2mapname);
		return NULL;
	}

	fh = CreateFile (ucs2filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (fh == INVALID_HANDLE_VALUE) {
		/* Cannot open */
		/* fprintf (stderr, "arikkei_mmap: File %s cannot be opened for reading\n", filename); */
		free (ucs2filename);
		if (ucs2mapname) free (ucs2mapname);
		return NULL;
	}

	mh = CreateFileMapping (fh, NULL, PAGE_READONLY, 0, 0, ucs2mapname);
	if (mh == NULL) {
		/* Mapping failed */
		/* fprintf (stderr, "arikkei_mmap: File %s cannot be mapped as %s\n", filename, mapname); */
		DWORD ecode = GetLastError ();
		fprintf (stderr, "arikkei_mmap: File %s cannot be mapped as %s (Error %d)\n", filename, mapname, ecode);
		CloseHandle (fh);
		free (ucs2filename);
		if (ucs2mapname) free (ucs2mapname);
		return NULL;
	}
	/* Get a pointer to the file-mapped shared memory. */
	cdata = (unsigned char *) MapViewOfFile (mh, FILE_MAP_READ, 0, 0, 0);

#ifdef PRINT_MAPSIZE
	if (!cdata) {
		DWORD ecode = GetLastError ();
		fprintf (stderr, "arikkei_mmap: Error %d\n", ecode);
	} else {
		total += st.st_size;
		fprintf (stderr, "MMap size+: %x\n", (unsigned int) total);
	}
#endif
	CloseHandle (mh);
	CloseHandle (fh);

	free (ucs2filename);
	if (ucs2mapname) free (ucs2mapname);

	*size = st.st_size;

	return cdata;
}

void
arikkei_munmap (const unsigned char *cdata, size_t size)
{
	/* Release data */
	UnmapViewOfFile (cdata);

#ifdef PRINT_MAPSIZE
	total -= size;
	fprintf (stderr, "MMap size-: %x\n", (unsigned int) total);
#endif
}

#else

const unsigned char *
arikkei_mmap (const unsigned char *filename, size_t *size, const unsigned char *name)
{
	unsigned char *cdata;
	struct stat st;
	cdata = NULL;
	if (!stat ((const char *) filename, &st) && S_ISREG (st.st_mode) && (st.st_size > 8)) {
		int fd;
		fd = open ((const char *) filename, O_RDONLY);
		if (fd < 0) return NULL;
		cdata = mmap (NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
		close (fd);
		if ((!cdata) || (cdata == (unsigned char *) -1)) return NULL;
	}

	*size = st.st_size;

	return cdata;
}

void
arikkei_munmap (const unsigned char *cdata, size_t size)
{
	munmap ((void *) cdata, size);

#ifdef DEBUG
	total -= size;
#endif
}

#endif

FILE *
arikkei_fopen (const unsigned char *filename, const unsigned char *mode)
{
#ifdef WIN32
	unsigned short *ucs2filename, *ucs2mode;
	FILE *fs;
	if (!filename || !*filename) return NULL;
	ucs2filename = arikkei_utf8_ucs2_strdup (filename);
	ucs2mode = arikkei_utf8_ucs2_strdup (mode);
	fs = _wfopen (ucs2filename, ucs2mode);
	free (ucs2filename);
	free (ucs2mode);
	return fs;
#else
	return fopen ((const char *) filename, (const char *) mode);
#endif
}

