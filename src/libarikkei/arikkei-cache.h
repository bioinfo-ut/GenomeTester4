#ifndef __ARIKKEI_CACHE_H__
#define __ARIKKEI_CACHE_H__

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

#include <libarikkei/arikkei-dict.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _ArikkeiCache ArikkeiCache;
typedef struct _ArikkeiCacheEntry ArikkeiCacheEntry;

struct _ArikkeiCache {
	ArikkeiDict dict;
	unsigned int maxsize;
	unsigned int currentsize;
	unsigned int nentries;
	ArikkeiCacheEntry *entries;
	int first;
	int last;
	int free;
	void *(* key_dup) (const void *key);
	void (* key_free) (void *key);
	void (* object_free) (const void *object);
};

void arikkei_cache_setup_full (ArikkeiCache *cache, unsigned int size,
			      unsigned int (* key_hash) (const void *),
			      unsigned int (* key_equal) (const void *, const void *),
				  void *(* key_dup) (const void *key),
				  void (* key_free) (void *key),
				  void (* object_free) (const void *object));

void arikkei_cache_setup_string (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object));
void arikkei_cache_setup_pointer (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object));
void arikkei_cache_setup_int (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object));
void arikkei_cache_release (ArikkeiCache *cache);

void arikkei_cache_insert (ArikkeiCache *cache, const void *key, void *object, unsigned int size);
void arikkei_cache_remove (ArikkeiCache *cache, const void *key);

const void *arikkei_cache_lookup (ArikkeiCache *cache, const void *key);
const void *arikkei_cache_lookup_notouch (ArikkeiCache *cache, const void *key);

#ifdef __cplusplus
}
#endif

#endif
