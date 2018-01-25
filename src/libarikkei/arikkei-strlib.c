#define __ARIKKEI_STRLIB_C__

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

#include <math.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdio.h>

#include "arikkei-strlib.h"

unsigned int
arikkei_memcpy (unsigned char *d, unsigned int d_len, const unsigned char *s, unsigned int s_len)
{
	unsigned int len = s_len;
	if (len > d_len) len = d_len;
	if ((len > 0) && d) memcpy (d, s, len);
	return s_len;
}

unsigned int
arikkei_memcpy_str (unsigned char *d, unsigned int d_len, const unsigned char *s)
{
	unsigned int s_len = (unsigned int) strlen ((const char *) s);
	arikkei_memcpy (d, d_len, s, s_len);
	return s_len;
}

unsigned int
arikkei_strncpy (unsigned char *d, unsigned int d_len, const unsigned char *s)
{
	unsigned int len = arikkei_memcpy_str (d, d_len, s);
	if ((len < d_len) && d) d[len] = 0;
	return len + 1;
}

unsigned int
arikkei_strtod_simple (const unsigned char *str, unsigned int len, double *val)
{
	const unsigned char *p;
	double sign, integra, fract;
	unsigned int valid;
	p = str;
	valid = 0;
	/* NULL string */
	if (!p) return 0;
	/* Skip space */
	while (*p && (*p == ' ')) p += 1;
	/* Empty string */
	if (!*p) return 0;
	/* Sign */
	if (*p == '-') {
		sign = -1.0;
		p += 1;
	} else {
		if (*p == '+') p += 1;
		sign = 1.0;
	}
	/* Value */
	/* Integra */
	integra = 0.0;
	fract = 0.0;
	while (*p && (*p >= '0') && (*p <= '9')) {
		integra = integra * 10.0 + (double) (*p - '0');
		valid = 1;
		p += 1;
	}
	if (*p == '.') {
		double denom;
		/* Fract */
		p += 1;
		denom = 1.0;
		while (*p && (*p >= '0') && (*p <= '9')) {
			fract = fract * 10.0 + (double) (*p - '0');
			denom *= 10.0;
			valid = 1;
			p += 1;
		}
		fract = fract / denom;
	}
	if (!valid) return 0;
	*val = sign * (integra + fract);

	/* assert ((*val > -1e17) && (*val < 1e17)); */

	return (unsigned int) (p - str);
}

unsigned int
arikkei_strtod_exp (const unsigned char *str, unsigned int len, double *val)
{
	const unsigned char *p;
	double rval;
	unsigned int slen;
	p = str;
	slen = arikkei_strtod_simple (p, len, &rval);
	if (!slen) return 0;
	p += slen;
	if ((*p == 'e') || (*p == 'E')) {
		double exval;
		unsigned int exlen;
		exlen = arikkei_strtod_simple (p + 1, (unsigned int) ((str + len) - (p + 1)), &exval);
		if (exlen) {
			p += 1;
			p += exlen;
			rval = rval * pow (10.0, exval);
		}
	}
	*val = rval;

	/* assert ((*val > -1e17) && (*val < 1e17)); */

	return (unsigned int) (p - str);
}

unsigned int
arikkei_itoa (unsigned char *buf, unsigned int len, int val)
{
	char c[32];
	int p, i;
	p = 0;
	if (val < 0) {
		buf[p++] = '-';
		val = -val;
	}
	i = 0;
	do {
		c[32 - (++i)] = '0' + (val % 10);
		val /= 10;
	} while (val > 0);
	memcpy (buf + p, &c[32 - i], i);
	p += i;
	buf[p] = 0;
	return p;
}

unsigned int
arikkei_dtoa_simple (unsigned char *buf, unsigned int len, double val,
		     unsigned int tprec, unsigned int fprec, unsigned int padf)
{
	double dival, fval, epsilon;
	int idigits, ival, ilen, i;
	i = 0;
	if (val < 0.0) {
		buf[i++] = '-';
		val = fabs (val);
	}
	/* Determine number of integral digits */
	if (val >= 1.0) {
		idigits = (int) floor (log10 (val));
	} else {
		idigits = 0;
	}

	/* Determine the actual number of fractional digits */
	if ((tprec - idigits) > fprec) fprec = tprec - idigits;
	/* Find epsilon */
	epsilon = 0.5 * pow (10.0, - (double) fprec);
	/* Round value */
	val += epsilon;
	epsilon *= 2;
	/* Extract integral and fractional parts */
	dival = floor (val);
	ival = (int) dival;
	fval = val - dival;
	/* Write integra */
	ilen = arikkei_itoa (buf + i, len - i, ival);
	i += ilen;
	tprec -= ilen;
	if ((fprec > 0) && (padf || (fval > epsilon))) {
		buf[i++] = '.';
		while ((fprec > 0) && (padf || (fval > epsilon))) {
			fval *= 10.0;
			epsilon *= 10.0;
			dival = floor (fval);
			fval -= dival;
			buf[i++] = '0' + (unsigned char) dival;
			fprec -= 1;
		}

	}
	buf[i] = 0;
	return i;
}

unsigned int
arikkei_dtoa_exp (unsigned char *buf, unsigned int len, double val,
		  unsigned int tprec, unsigned int padf)
{
	if ((val == 0.0) || ((fabs (val) >= 0.01) && (fabs(val) < 10000000))) {
		return arikkei_dtoa_simple (buf, len, val, tprec, 0, padf);
	} else {
		double eval;
		int p;
		eval = floor (log10 (fabs (val)));
		val = val / pow (10.0, eval);
		p = arikkei_dtoa_simple (buf, len, val, tprec, 0, padf);
		buf[p++] = 'e';
		p += arikkei_itoa (buf + p, len - p, (int) eval);
		return p;
	}
}

unsigned int
arikkei_unicode_utf8_bytes (unsigned int uval)
{
	if (uval < 0x80) return 1;
	if (uval < 0x800) return 2;
	if (uval < 0x10000) return 3;
	if (uval < 0x200000) return 4;
	if (uval < 0x4000000) return 5;
	return 6;
}

unsigned int
arikkei_utf8_strlen (const unsigned char *str)
{
	const unsigned char *s;
	unsigned int len;
	s = str;
	len = 0;
	while (*s) {
		if ((*s & 0x80) == 0x0) {
			s += 1;
		} else if ((*s & 0xe0) == 0xc0) {
			s += 2;
		} else if ((*s & 0xf0) == 0xe0) {
			s += 3;
		} else if ((*s & 0xf8) == 0xf0) {
			s += 4;
		} else if ((*s & 0xfc) == 0xf8) {
			s += 5;
		} else if ((*s & 0xfe) == 0xfc) {
			s += 6;
		} else {
			return len;
		}
		len += 1;
	}
	return len;
}

unsigned int
arikkei_utf8_ucs2_strcpy (const unsigned char *s, unsigned short *d)
{
	unsigned int dp;
	dp = 0;
	while (*s) {
		if ((*s & 0x80) == 0x0) {
			d[dp++] = s[0];
			s += 1;
		} else if ((*s & 0xe0) == 0xc0) {
			d[dp++] = ((s[0] & 0x1f) << 6) | (s[1] & 0x3f);
			s += 2;
		} else if ((*s & 0xf0) == 0xe0) {
			d[dp++] = ((s[0] & 0x0f) << 12) | ((s[1] & 0x3f) << 6) | (s[2] & 0x3f);
			s += 3;
		} else if ((*s & 0xf8) == 0xf0) {
			d[dp++] = ((s[0] & 0x07) << 18) | ((s[1] & 0x3f) << 12) |
					  ((s[2] & 0x3f) << 6) | (s[3] & 0x3f);
			s += 4;
		} else if ((*s & 0xfc) == 0xf8) {
			d[dp++] = ((s[0] & 0x03) << 24) | ((s[1] & 0x3f) << 18) |
					  ((s[2] & 0x3f) << 12) | ((s[3] & 0x3f) << 6) | (s[4] & 0x3f);
			s += 5;
		} else if ((*s & 0xfe) == 0xfc) {
			d[dp++] = ((s[0] & 0x01) << 30) | ((s[1] & 0x3f) << 24) | ((s[2] & 0x3f) << 18) |
					  ((s[3] & 0x3f) << 12) | ((s[4] & 0x3f) << 6) | (s[5] & 0x3f);
			s += 6;
		} else {
			break;
		}
	}
	d[dp] = 0;
	return dp;
}

int
arikkei_utf8_get_unicode (const unsigned char **str, unsigned int length)
{
	const unsigned char *s = *str;
	if (((*s & 0x80) == 0x0) && (length >= 1)) {
		*str += 1;
		return s[0];
	} else if (((*s & 0xe0) == 0xc0) && (length >= 2)) {
		*str += 2;
		return ((s[0] & 0x1f) << 6) | (s[1] & 0x3f);
	} else if (((*s & 0xf0) == 0xe0) && (length >= 3)) {
		*str += 3;
		return ((s[0] & 0x0f) << 12) | ((s[1] & 0x3f) << 6) | (s[2] & 0x3f);
	} else if (((*s & 0xf8) == 0xf0) && (length >= 4)) {
		*str += 4;
		return ((s[0] & 0x07) << 18) | ((s[1] & 0x3f) << 12) | ((s[2] & 0x3f) << 6) | (s[3] & 0x3f);
	} else if (((*s & 0xfc) == 0xf8) && (length >= 5)) {
		*str += 5;
		return ((s[0] & 0x03) << 24) | ((s[1] & 0x3f) << 18) | ((s[2] & 0x3f) << 12) | ((s[3] & 0x3f) << 6) | (s[4] & 0x3f);
	} else if (((*s & 0xfe) == 0xfc) && (length >= 6)) {
		*str += 6;
		return ((s[0] & 0x01) << 30) | ((s[1] & 0x3f) << 24) | ((s[2] & 0x3f) << 18) | ((s[3] & 0x3f) << 12) | ((s[4] & 0x3f) << 6) | (s[5] & 0x3f);
	}
	return -1;
}

unsigned short *
arikkei_utf8_ucs2_strdup (const unsigned char *s)
{
	unsigned short *d;
	int len;
	len = arikkei_utf8_strlen (s);
	d = malloc ((len + 1) * sizeof (unsigned short));
	arikkei_utf8_ucs2_strcpy (s, d);
	return d;
}

unsigned int
arikkei_ucs2_strlen (const unsigned short *str)
{
	const unsigned short *s;
	s = str;
	while (*s) s += 1;
	return (unsigned int) (s - str);
}

unsigned int
arikkei_ucs2_utf8_strcpy (const unsigned short *s, unsigned char *d)
{
	unsigned int dp;
	dp = 0;
	while (*s) {
		if (*s < 0x80) {
			d[dp++] = (unsigned char) *s;
		} else if (*s < 0x800) {
			d[dp++] = 0xc0 | (unsigned char) (*s >> 6);
			d[dp++] = 0x80 | (*s & 0x3f);
		} else /* if (*s < 0x10000) */ {
			/* As long as these are 16 bit we do not need more */
			/* fixme: Surrogates (Lauris) */
			d[dp++] = 0xe0 | (*s >> 12);
			d[dp++] = 0x80 | ((*s >> 6) & 0x3f);
			d[dp++] = 0x80 | (*s & 0x3f);
#if 0
		} else if (*s < 0x200000) {
			d[dp++] = 0xf0 | (*s >> 18);
			d[dp++] = 0x80 | ((*s >> 12) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 6) & 0x3f);
			d[dp++] = 0x80 | (*s & 0x3f);
		} else if (*s < 0x4000000) {
			d[dp++] = 0xf8 | (*s >> 24);
			d[dp++] = 0x80 | ((*s >> 18) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 12) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 6) & 0x3f);
			d[dp++] = 0x80 | (*s & 0x3f);
		} else {
			d[dp++] = 0xfc | (*s >> 30);
			d[dp++] = 0x80 | ((*s >> 24) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 18) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 12) & 0x3f);
			d[dp++] = 0x80 | ((*s >> 6) & 0x3f);
			d[dp++] = 0x80 | (*s & 0x3f);
#endif
		}
		s += 1;
	}
	d[dp] = 0;
	return dp;
}

unsigned char *
arikkei_ucs2_utf8_strdup (const unsigned short *s)
{
	unsigned char *d;
	int len;
	len = arikkei_ucs2_utf8_bytelen (s);
	d = malloc ((len + 1) * sizeof (unsigned char));
	arikkei_ucs2_utf8_strcpy (s, d);
	return d;
}

unsigned int
arikkei_ucs2_utf8_bytelen (const unsigned short *str)
{
	const unsigned short *s;
	unsigned int len;
	s = str;
	len = 0;
	while (*s) {
		len += arikkei_unicode_utf8_bytes (*s);
		s += 1;
	}
	return len;
}

unsigned short *
arikkei_ucs2_strdup (const unsigned short *s)
{
	unsigned int len;
	unsigned short *d;
	len = arikkei_ucs2_strlen (s);
	d = (unsigned short *) malloc ((len + 1) * sizeof (unsigned short));
	memcpy (d, s, (len + 1) * sizeof (unsigned short));
	return d;
}

unsigned int
arikkei_ucs2_strncpy (const unsigned short *s, unsigned short *d, unsigned int maxlen)
{
	unsigned int di;
	const unsigned short *p;

	di = 0;
	p = s;

	if (maxlen < 0) {
		while (*p) {
			d[di] = *p;
			++p;
			++di;
		}
		d[di] = 0;
	} else {
		while (*p && di < maxlen) {
			d[di] = *p;
			++p;
			++di;
		}
		if (di < maxlen) d[di] = 0;
	}
	
	return di;
}
