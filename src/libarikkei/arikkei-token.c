#define __ARIKKEI_TOKEN_C__

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

#include <malloc.h>
#include <string.h>

#include "arikkei-token.h"

ArikkeiToken *
arikkei_token_set_from_string (ArikkeiToken *token, const unsigned char *cdata)
{
	token->cdata = cdata;
	token->len = (cdata) ? strlen ((const char *) cdata) : 0;
	return token;
}

ArikkeiToken *
arikkei_token_set_from_data (ArikkeiToken *token, const unsigned char *cdata, size_t start, size_t end)
{
	token->cdata = cdata + start;
	token->len = end - start;
	return token;
}

ArikkeiToken *
arikkei_token_set_from_token (ArikkeiToken *token, const ArikkeiToken *src)
{
	token->cdata = src->cdata;
	token->len = src->len;
	return token;
}

unsigned int
arikkei_token_is_equal (const ArikkeiToken *token, const ArikkeiToken *other)
{
	if (!arikkei_token_is_valid (token)) return 0;
	if (!arikkei_token_is_valid (other)) return 0;
	if (!arikkei_token_is_empty (token)) {
		if (!arikkei_token_is_empty (other)) {
			return ((token->len == other->len) && !strncmp ((const char *) token->cdata, (const char *) other->cdata, token->len));
		} else {
			return 0;
		}
	} else {
		if (!arikkei_token_is_empty (other)) {
			return 0;
		} else {
			return 1;
		}
	}
}

unsigned int
arikkei_token_is_equal_string (const ArikkeiToken *token, const unsigned char *str)
{
	if (!arikkei_token_is_valid (token)) return 0;
	return !arikkei_token_strcmp (token, str);
}

unsigned char *
arikkei_token_strdup (const ArikkeiToken *token)
{
	if (arikkei_token_is_valid (token)) {
		unsigned char *b;
		b = malloc (token->len + 1);
		if (token->len) strncpy ((char *) b, (const char *) token->cdata, token->len);
		b[token->len] = 0;
		return b;
	} else {
		return NULL;
	}
}

size_t
arikkei_token_strcpy (const ArikkeiToken *token, unsigned char *b)
{
	if (arikkei_token_is_valid (token)) {
		if (token->len) strncpy ((char *) b, (const char *) token->cdata, token->len);
		b[token->len] = 0;
		return token->len;
	} else {
		b[0] = 0;
		return 0;
	}
}

size_t
arikkei_token_strncpy (const ArikkeiToken *token, unsigned char *b, size_t size)
{
	if (size < 1) return 0;
	if (arikkei_token_is_valid (token) && (size > 1)) {
		size_t len = token->len;
		if (len > (size - 1)) len = size - 1;
		if (len) strncpy ((char *) b, (const char *) token->cdata, len);
		b[len] = 0;
		return len;
	} else {
		b[0] = 0;
		return 0;
	}
}

int
arikkei_token_strcmp (const ArikkeiToken *token, const unsigned char *str)
{
	if (!arikkei_token_is_valid (token)) return -1;
	if (str) {
		return arikkei_token_strncmp (token, str, strlen ((const char *) str));
	} else {
		return arikkei_token_is_empty (token);
	}
}

int
arikkei_token_strncmp (const ArikkeiToken *token, const unsigned char *str, size_t size)
{
	if (!arikkei_token_is_valid (token)) return -1;
	if (!arikkei_token_is_empty (token)) {
		if (size > 0) {
			size_t len, clen;
			int cval;
			len = token->len;
			clen = (len < size) ? len : size;
			cval = strncmp ((const char *) token->cdata, (const char *) str, clen);
			if (cval) return cval;
			if (len < size) return -1;
			if (len > size) return 1;
			return 0;
		} else {
			return 1;
		}
	} else {
		if (size > 0) {
			return -1;
		} else {
			return 0;
		}
	}
}

ArikkeiToken *
arikkei_token_get_first_line (const ArikkeiToken *token, ArikkeiToken *dst)
{
	return arikkei_token_get_line (token, dst, 0);
}

ArikkeiToken *
arikkei_token_get_line (const ArikkeiToken *token, ArikkeiToken *dst, size_t start)
{
	if (!arikkei_token_is_empty (token)) {
		const unsigned char *p;
		size_t e;
		p = token->cdata;
		e = start;
		while ((e < token->len) && ((p[e] >= 32) || (p[e] == 9))) e += 1;
		arikkei_token_set_from_data (dst, p, start, e);
	} else {
		arikkei_token_set_from_data (dst, token->cdata, 0, 0);
	}
	return dst;
}

ArikkeiToken *
arikkei_token_next_line (const ArikkeiToken *token, ArikkeiToken *dst, const ArikkeiToken *line)
{
	if (!arikkei_token_is_empty (token)) {
		if (arikkei_token_is_valid (line)) {
			const unsigned char *p;
			size_t s;
			p = token->cdata;
			s = (line->cdata + line->len) - token->cdata;
			while ((s < token->len) && (p[s] < 32) && (p[s] != 9)) s += 1;
			return arikkei_token_get_line (token, dst, s);
		} else {
			return arikkei_token_get_first_line (token, dst);
		}
	} else {
		arikkei_token_set_from_data (dst, token->cdata, 0, 0);
	}
	return dst;
}

ArikkeiToken *
arikkei_token_get_token (const ArikkeiToken *token, ArikkeiToken *dst, size_t start, unsigned int space_is_separator)
{
	if (!arikkei_token_is_empty (token)) {
		const unsigned char *p;
		p = token->cdata;
		while ((start < token->len) && (p[start] == 32)) start += 1;
		if (start < token->len) {
			size_t e;
			e = start;
			while ((e < token->len) && ((p[e] > 32) || ((p[e] == 32) && !space_is_separator))) e += 1;
			arikkei_token_set_from_data (dst, token->cdata, start, e);
		} else {
			arikkei_token_set_from_data (dst, token->cdata, start, token->len);
		}
	} else {
		arikkei_token_set_from_data (dst, token->cdata, 0, 0);
	}
	return dst;
}

ArikkeiToken *
arikkei_token_next_token (const ArikkeiToken *token, ArikkeiToken *dst, const ArikkeiToken *prev, unsigned int space_is_separator)
{
	if (!arikkei_token_is_empty (token)) {
		if (arikkei_token_is_valid (prev)) {
			size_t s;
			s = (prev->cdata + prev->len) - token->cdata;
			return arikkei_token_get_token (token, dst, s, space_is_separator);
		} else {
			return arikkei_token_get_token (token, dst, 0, space_is_separator);
		}
	} else {
		arikkei_token_set_from_data (dst, token->cdata, 0, 0);
	}
	return dst;
}

int
arikkei_token_tokenize (ArikkeiToken *token, ArikkeiToken *tokens, int maxtokens, unsigned int space_is_separator, unsigned int multi)
{
	const unsigned char *p;
	size_t s;
	int ntokens;
	if (arikkei_token_is_empty (token)) return 0;
	ntokens = 0;
	p = token->cdata;
	s = 0;
	while ((s < token->len) && (ntokens < maxtokens)) {
		size_t e;
		e = s;
		while ((e < token->len) && ((p[e] > 32) || ((p[e] == 32) && !space_is_separator))) e += 1;
		if (ntokens == (maxtokens - 1)) {
			while ((e < token->len) && ((p[e] >= 32) || (p[e] == 9))) e += 1;
		}
		arikkei_token_set_from_data (tokens + ntokens, token->cdata, s, e);
		s = e + 1;
		if (multi) {
			while ((s < token->len) && ((p[s] < 32) || ((p[s] == 32) && space_is_separator))) s += 1;
		}
		ntokens += 1;
	}
	return ntokens;
}

int
arikkei_token_tokenize_ws (ArikkeiToken *token, ArikkeiToken *tokens, int maxtokens, const unsigned char *ws, unsigned int multi)
{
	size_t len, s;
	int ntokens;
	if (arikkei_token_is_empty (token)) return 0;
	len = strlen ((const char *) ws);
	ntokens = 0;
	s = 0;
	while ((s < token->len) && (ntokens < maxtokens)) {
		size_t e;
		if (ntokens != (maxtokens - 1)) {
			e = s;
			while (e < token->len) {
				size_t i;
				for (i = 0; i < len; i++) {
					if (token->cdata[e] == ws[i]) break;
				}
				if (i < len) break;
				e += 1;
			}
		} else {
			e = token->len;
		}
		arikkei_token_set_from_data (tokens + ntokens, token->cdata, s, e);
		s = e + 1;
		if (multi) {
			while (s < token->len) {
				size_t i;
				for (i = 0; i < len; i++) {
					if (token->cdata[s] == ws[i]) break;
				}
				if (i < len) break;
				s += 1;
			}
		}
		ntokens += 1;
	}
	return ntokens;
}

ArikkeiToken *
arikkei_token_strip_start (ArikkeiToken *token, ArikkeiToken *dst)
{
	const unsigned char *p;
	size_t s;
	p = token->cdata;
	s = 0;
	if (p) {
		while ((s < token->len) && (p[s] <= 32)) s += 1;
	}
	arikkei_token_set_from_data (dst, token->cdata, s, token->len);
	return dst;
}

ArikkeiToken *
arikkei_token_strip_start_ws (ArikkeiToken *token, ArikkeiToken *dst, const unsigned char *ws)
{
	size_t len, s;
	len = strlen ((const char *) ws);
	s = 0;
	if (token->cdata) {
		while (s < token->len) {
			size_t i;
			for (i = 0; i < len; i++) {
				if (token->cdata[s] == ws[i]) break;
			}
			if (i >= len) break;
			s += 1;
		}
	}
	arikkei_token_set_from_data (dst, token->cdata, s, token->len);
	return dst;
}

ArikkeiToken *
arikkei_token_strip_end (ArikkeiToken *token, ArikkeiToken *dst)
{
	const unsigned char *p;
	int e;
	p = token->cdata;
	e = (int) token->len - 1;
	if (p) {
		while ((e >= 0) && (p[e] <= 32)) e -= 1;
	}
	arikkei_token_set_from_data (dst, token->cdata, 0, e + 1);
	return dst;
}

ArikkeiToken *
arikkei_token_strip_end_ws (ArikkeiToken *token, ArikkeiToken *dst, const unsigned char *ws)
{
	size_t len;
	int e;
	len = strlen ((const char *) ws);
	e = (int) token->len - 1;
	if (token->cdata) {
		while (e >= 0) {
			size_t i;
			for (i = 0; i < len; i++) {
				if (token->cdata[e] == ws[i]) break;
			}
			if (i >= len) break;
			e -= 1;
		}
	}
	arikkei_token_set_from_data (dst, token->cdata, 0, e + 1);
	return dst;
}

ArikkeiToken *
arikkei_token_strip (ArikkeiToken *token, ArikkeiToken *dst)
{
	const unsigned char *p;
	int s;
	int e;
	p = token->cdata;
	s = 0;
	e = (int) token->len - 1;
	if (p) {
		while ((s < (int) token->len) && (p[s] <= 32)) s += 1;
		while ((e >= s) && (p[e] <= 32)) e -= 1;
	}
	arikkei_token_set_from_data (dst, token->cdata, s, e + 1);
	return dst;
}

ArikkeiToken *
arikkei_token_strip_ws (ArikkeiToken *token, ArikkeiToken *dst, const unsigned char *ws)
{
	size_t len;
	int s, e;
	len = strlen ((const char *) ws);
	s = 0;
	e = (int) token->len - 1;
	if (token->cdata) {
		while (s < (int) token->len) {
			size_t i;
			for (i = 0; i < len; i++) {
				if (token->cdata[s] == ws[i]) break;
			}
			if (i >= len) break;
			s += 1;
		}
		while (e >= s) {
			size_t i;
			for (i = 0; i < len; i++) {
				if (token->cdata[e] == ws[i]) break;
			}
			if (i >= len) break;
			e -= 1;
		}
	}
	arikkei_token_set_from_data (dst, token->cdata, s, e + 1);
	return dst;
}

unsigned char *
arikkei_token_strconcat (const ArikkeiToken *tokens, int size, const unsigned char *separator)
{
	unsigned char *str, *p;
	size_t slen, len;
	int i;
	slen = strlen ((const char *) separator);
	len = 1;
	for (i = 0; i < size; i++) len += tokens[i].len;
	len += (size - 1) * slen;
	str = malloc (len + 1);
	p = str;
	for (i = 0; i < size; i++) {
		if ((i > 0) && (slen > 0)) {
			strncpy ((char *) p, (const char *) separator, slen);
			p += slen;
		}
		p += arikkei_token_strcpy (tokens + i, p);
	}
	return str;
}







