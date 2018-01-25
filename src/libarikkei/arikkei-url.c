#define __ARIKKEI_URL_C__

/*
 * URL Handling
 *
 * Copyright Laurls Kaplinski 2010
 */

#include <string.h>
#include <malloc.h>
#include <ctype.h>

#include "arikkei-url.h"

static unsigned char *
strdup_substr (const unsigned char *str, int start, int end)
{
	unsigned char *d;
	if (!str || (start < 0) || (end < 0) || (start > end)) return NULL;
	d = (unsigned char *) malloc (end - start + 1);
	if (end > start) memcpy (d, str + start, end - start);
	d[end - start] = 0;
	return d;
}

unsigned int
arikkei_url_setup (ArikkeiURL *url, const unsigned char *address, const unsigned char *defaultprotocol)
{
	int next, prot_s, prot_e, dom_s, dom_e, dir_s, dir_e, file_s, file_e, ref_s, ref_e, arg_s;
	int i;

	/* Clear everything */
	memset (url, 0, sizeof (ArikkeiURL));

	if (!address || !*address) return 0;
	if (!defaultprotocol && !strchr ((const char *) address, ':')) return 0;

	next = 0;

	prot_s = -1;
	prot_e = -1;
	/* Protocol - from start to the first ':' anly alphanumeric characters, '-' and '_' are allowed */
	for (i = next; address[i]; i++) {
		if (address[i] == ':') {
			prot_s = next;
			prot_e = i;
			next = i + 1;
			break;
		}
		if (!isalnum (address[i]) && (address[i] != '-') && (address[i] != '_')) break;
	}
#ifdef WIN32
	if ((prot_e - prot_s) == 1) {
		/* One letter protocol is drive letter */
		prot_s = -1;
		prot_e = -1;
	}
#endif
	if (!defaultprotocol && (prot_e < 0)) return 0;

	dom_s = -1;
	dom_e = -1;
	/* Domain - from '//' to next '/' */
	/* Domain - from '(' to ')' */
	if ((address[next] == '/') && (address[next + 1] == '/')) {
		dom_s = next + 2;
		for (dom_e = dom_s; address[dom_e]; dom_e++) if (address[dom_e] == '/') break;
		next = (address[dom_e]) ? dom_e + 1 : dom_e;
	} else if (address[next] == '(') {
		dom_s = next + 1;
		for (dom_e = dom_s; address[dom_e]; dom_e++) if (address[dom_e] == ')') break;
		next = (address[dom_e]) ? dom_e + 1 : dom_e;
	}

	dir_s = -1;
	dir_e = -1;
	file_s = -1;
	file_e = -1;
	/* Directory - from next position to the last '/' before '#' or '&' */
	/* File - from last '/' to '#' or '&' */
	for (i = next; address[i]; i++) {
		if (address[i] == '/') {
			dir_s = next;
			dir_e = i + 1;
			file_s = i + 1;
		}
		if (address[i] == '#') break;
		if (address[i] == '?') break;
	}
	if (i > next) {
		if (file_s < 0) file_s = next;
		file_e = i;
	}
	next = i;

	ref_s = -1;
	ref_e = -1;
	/* Reference - from '#' to '?' */
	if (address[next] == '#') {
		ref_s = next + 1;
		for (ref_e = ref_s; address[ref_e]; ref_e++) if (address[ref_e] == '?') break;
		next = (address[ref_e]) ? ref_e + 1 : ref_e;
	}

	arg_s = -1;
	/* Arguments - from '?' to the end of address */
	if (address[next] == '?') {
		arg_s = next + 1;
	}

	url->protocol = (prot_e > 0) ? strdup_substr (address, prot_s, prot_e) : (unsigned char *) strdup ((const char *) defaultprotocol);
	url->domain = strdup_substr (address, dom_s, dom_e);
	url->directory = strdup_substr (address, dir_s, dir_e);
	url->filename = strdup_substr (address, file_s, file_e);
	url->reference = strdup_substr (address, ref_s, ref_e);
	url->arguments = (arg_s >= 0) ? (unsigned char *) strdup ((const char *) address + arg_s) : NULL;

	if (prot_s >= 0) {
		url->address = (unsigned char *) strdup ((const char *) address);
		url->base = strdup_substr (address, prot_s, dir_e);
		url->path = strdup_substr (address, prot_s, file_e);
	} else {
		size_t plen = strlen ((const char *) defaultprotocol);
		size_t alen = strlen ((const char *) address);
		size_t len = plen + 1 + alen + 1;
		url->address = (unsigned char *) malloc (len);
		memcpy (url->address, defaultprotocol, plen);
		url->address[plen] = ':';
		memcpy (url->address + plen + 1, address, alen);
		url->address[len - 1] = 0;
		if (dir_e >= 0) {
			len = plen + 1 + dir_e + 1;
			url->base = (unsigned char *) malloc (len);
			memcpy (url->base, defaultprotocol, plen);
			url->base[plen] = ':';
			memcpy (url->base + plen + 1, address, dir_e);
			url->base[len - 1] = 0;
		}
		if (file_e >= 0) {
			len = plen + 1 + file_e + 1;
			url->path = (unsigned char *) malloc (len);
			memcpy (url->path, defaultprotocol, plen);
			url->path[plen] = ':';
			memcpy (url->path + plen + 1, address, file_e);
			url->path[len - 1] = 0;
		}
	}

	return 1;
}

void
arikkei_url_release (ArikkeiURL *url)
{
	if (url->address) free (url->address);
	if (url->protocol) free (url->protocol);
	if (url->domain) free (url->domain);
	if (url->directory) free (url->directory);
	if (url->filename) free (url->filename);
	if (url->reference) free (url->reference);
	if (url->arguments) free (url->arguments);
	if (url->base) free (url->base);
	if (url->path) free (url->path);

	/* Clear everything */
	memset (url, 0, sizeof (ArikkeiURL));
}

unsigned char *
arikkei_build_file_url (const unsigned char *path)
{
	char *c, *p;

	if (!path) return NULL;

	c = (char *) malloc (strlen ((const char *) path) + 6);
	strcpy (c, "file:");
	strcpy (c + 5, (const char *) path);
	for (p = c + 5; *p; p++) {
		if (*p == '\\') *p = '/';
	}

	return (unsigned char *) c;
}

static unsigned int
path_is_absolute (const unsigned char *path)
{
	if (!path) return 0;
	if (path[0] == '/') return 1;
#ifdef WIN32
	if (path[0] && (path[1] == ':')) return 1;
#endif
	return 0;
}

static unsigned char *
build_url (const unsigned char *protocol, const unsigned char *domain, const unsigned char *directory, const unsigned char *filename, const unsigned char *reference, const unsigned char *arguments)
{
	size_t plen, mlen, dlen, flen, rlen, alen, p;
	unsigned char *url;
	plen = (protocol) ? strlen ((const char *) protocol) + 1 : 0;
	mlen = (domain) ? strlen ((const char *) domain) + 2: 0;
	dlen = (directory) ? strlen ((const char *) directory) : 0;
	flen = (filename) ? strlen ((const char *) filename) : 0;
	rlen = (reference) ? strlen ((const char *) reference) + 1 : 0;
	alen = (arguments) ? strlen ((const char *) arguments) + 1 : 0;
	url = (unsigned char *) malloc (plen + mlen + dlen + flen + rlen + alen + 1);
	p = 0;
	if (plen) {
		memcpy (url + p, protocol, plen - 1);
		url[plen - 1] = ':';
		p += plen;
	}
	if (mlen) {
		url[p] = '(';
		memcpy (url + p + 1, domain, mlen - 2);
		url[p + mlen - 1] = ')';
		p += mlen;
	}
	if (dlen) {
		memcpy (url + p, directory, dlen);
		p += dlen;
	}
	if (flen) {
		memcpy (url + p, filename, flen);
		p += flen;
	}
	if (rlen) {
		url[p] = '#';
		memcpy (url + p + 1, reference, rlen - 1);
		p += rlen;
	}
	if (alen) {
		url[p] = '?';
		memcpy (url + p + 1, arguments, alen - 1);
		p += alen;
	}
	url[p] = 0;
	return url;
}

unsigned char *
arikkei_build_relative_url (const unsigned char *parent, const unsigned char *path)
{
	ArikkeiURL purl, curl;
	unsigned int p, n;
	unsigned char *url;
	if (!path) return NULL;
	if (!parent) return (unsigned char *) strdup ((const char *) path);
	arikkei_url_setup (&purl, parent, NULL);
	arikkei_url_setup (&curl, path, purl.protocol);
	if (!path_is_absolute (curl.directory) || (purl.protocol && strcmp ((const char *) purl.protocol, (const char *) curl.protocol))) {
		/* Child has relative URL or protocols differ */
		url = (unsigned char *) strdup ((const char *) curl.address);
		arikkei_url_release (&purl);
		arikkei_url_release (&curl);
		return url;
	}
	if (!purl.directory) {
		/* Parent or child does not have directory */
		url = (unsigned char *) strdup ((const char *) curl.address);
		arikkei_url_release (&purl);
		arikkei_url_release (&curl);
		return url;
	}
	p = 0;
	n = 0;
	while (purl.directory[p] && curl.directory[p] && (purl.directory[p] == curl.directory[p])) {
		if (purl.directory[p] == '/') n = p + 1;
		p += 1;
	}
	/* Now add protocol + directory[n] + filename + # + reference + ? + arguments */
	url = build_url (purl.protocol, purl.domain, curl.directory + n, curl.filename, curl.reference, curl.arguments);
	arikkei_url_release (&purl);
	arikkei_url_release (&curl);
	return url;
}

unsigned char *
arikkei_build_absolute_url (const unsigned char *parent, const unsigned char *path)
{
	ArikkeiURL purl, curl;
	size_t plen, clen;
	unsigned char *url, *c;
	if (!path) return NULL;
	if (!parent) return (unsigned char *) strdup ((const char *) path);
	arikkei_url_setup (&purl, parent, NULL);
	arikkei_url_setup (&curl, path, purl.protocol);
	if (path_is_absolute (curl.directory) || (purl.protocol && strcmp ((const char *) purl.protocol, (const char *) curl.protocol))) {
		/* Child has absolute URL or protocols differ */
		url = (unsigned char *) strdup ((const char *) curl.address);
		arikkei_url_release (&purl);
		arikkei_url_release (&curl);
		return url;
	}
	/* Append child directory to parent */
	plen = (purl.directory) ? strlen ((const char *) purl.directory) : 0;
	clen = (curl.directory) ? strlen ((const char *) curl.directory) : 0;
	c = (unsigned char *) malloc (plen + clen + 1);
	if (plen) memcpy (c, purl.directory, plen);
	if (clen) memcpy (c + plen, curl.directory, clen);
	c[plen + clen] = 0;
	/* Now add base + pdirectory+cdirectory + filename + # + reference + ? + arguments */
	url = build_url (purl.protocol, purl.domain, c, curl.filename, curl.reference, curl.arguments);
	free (c);
	arikkei_url_release (&purl);
	arikkei_url_release (&curl);
	return url;
}

