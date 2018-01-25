#ifndef __ARIKKEI_URL_H__
#define __ARIKKEI_URL_H__

/*
 * URL Handling
 *
 * Copyright Lauris Kaplinski 2010
 */

typedef struct _ArikkeiURL ArikkeiURL;

#ifdef __cplusplus
extern "C" {
#endif

/* Location */

/* [protocol]:[domain][[/]directory/][file][#reference][&arguments] */

struct _ArikkeiURL {
	/* protocol:[domain][[/]directory/][file][#reference][&arguments] */
	unsigned char *address;

	unsigned char *protocol;
	unsigned char *domain;
	unsigned char *directory;
	unsigned char *filename;
	unsigned char *reference;
	unsigned char *arguments;

	/* protocol:domain[[/]directory/] */
	unsigned char *base;
	/* protocol:domain[[/]directory/][file] */
	unsigned char *path;
};

unsigned int arikkei_url_setup (ArikkeiURL *url, const unsigned char *address, const unsigned char *defaultprotocol);
void arikkei_url_release (ArikkeiURL *url);

/* Build canonical file URL - i.e. "file:" prepended and separators set to "/" */
unsigned char *arikkei_build_file_url (const unsigned char *path);
/* Get url relative to parent */
unsigned char *arikkei_build_relative_url (const unsigned char *parent, const unsigned char *path);
/* Get absolute url from parent and relative path */
unsigned char *arikkei_build_absolute_url (const unsigned char *parent, const unsigned char *path);

#ifdef __cplusplus
}
#endif

#endif
