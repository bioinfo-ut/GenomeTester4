#define __ARIKKEI_CACHE_C__

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

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef WIN32
/* Disable VS2008 nagging */
#define strdup _strdup
#endif

#include <libarikkei/arikkei-dict.h>

#include "arikkei-cache.h"

struct _ArikkeiCacheEntry {
	int prev;
	int next;
	void *key;
	void *object;
	unsigned int size;
};

static void *
arikkei_string_dup (const void *str)
{
	return strdup ((const char *) str);
}

static void
arikkei_string_free (void *str)
{
	free (str);
}

void arikkei_cache_setup_full (ArikkeiCache *cache, unsigned int size,
			      unsigned int (* key_hash) (const void *),
			      unsigned int (* key_equal) (const void *, const void *),
				  void *(* key_dup) (const void *key),
				  void (* key_free) (void *key),
				  void (* object_free) (const void *object))
{
	arikkei_dict_setup_full (&cache->dict, 313, key_hash, key_equal);
	cache->maxsize = size;
	cache->currentsize = 0;
	cache->nentries = 0;
	cache->entries = NULL;
	cache->first = -1;
	cache->last = -1;
	cache->free = -1;
	cache->key_dup = key_dup;
	cache->key_free = key_free;
	cache->object_free = object_free;
}

void
arikkei_cache_setup_string (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object))
{
	arikkei_cache_setup_full (cache, size, arikkei_string_hash, arikkei_string_equal, arikkei_string_dup, arikkei_string_free, object_free);
}

void
arikkei_cache_setup_pointer (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object))
{
	arikkei_cache_setup_full (cache, size, arikkei_pointer_hash, arikkei_pointer_equal, NULL, NULL, object_free);
}

void
arikkei_cache_setup_int (ArikkeiCache *cache, unsigned int size, void (* object_free) (const void *object))
{
	arikkei_cache_setup_full (cache, size, arikkei_int_hash, arikkei_int_equal, NULL, NULL, object_free);
}

void
arikkei_cache_release (ArikkeiCache *cache)
{
	int pos;
	arikkei_dict_release (&cache->dict);
	for (pos = cache->first; pos >= 0; pos = cache->entries[pos].next) {
		if (cache->key_free) cache->key_free (cache->entries[pos].key);
		if (cache->object_free) cache->object_free (cache->entries[pos].object);
	}
	if (cache->entries) free (cache->entries);
}

static void
remove_entry_from_list (ArikkeiCache *cache, int pos)
{
	if (cache->entries[pos].prev >= 0) {
		cache->entries[cache->entries[pos].prev].next = cache->entries[pos].next;
	} else {
		cache->first = cache->entries[pos].next;
	}
	if (cache->entries[pos].next >= 0) {
		cache->entries[cache->entries[pos].next].prev = cache->entries[pos].prev;
	} else {
		cache->last = cache->entries[pos].prev;
	}
	cache->entries[pos].next = cache->free;
	cache->entries[pos].prev = -1;
	cache->free = pos;
}

static void
arikkei_cache_remove_entry (ArikkeiCache *cache, int pos)
{
	arikkei_dict_remove (&cache->dict, cache->entries[pos].key);
	if (cache->key_free) cache->key_free (cache->entries[pos].key);
	if (cache->object_free) cache->object_free (cache->entries[pos].object);
	cache->currentsize -= cache->entries[pos].size;
	remove_entry_from_list (cache, pos);
}

static void
arikkei_cache_ensure_space (ArikkeiCache *cache, unsigned int requested)
{
	while (requested > (cache->maxsize - cache->currentsize)) {
		assert (cache->last >= 0);
		arikkei_cache_remove_entry (cache, cache->last);
	}
}

void
arikkei_cache_insert (ArikkeiCache *cache, const void *key, void *object, unsigned int size)
{
	int pos;
	if (size > cache->maxsize) {
		if (cache->object_free) cache->object_free (object);
		return;
	}
	pos = (int) ((const char *) arikkei_dict_lookup (&cache->dict, key) - (const char *) 0) - 1;
	if (pos >= 0) {
		/* There is existing entry with the same key */
		if (cache->entries[pos].object == object) {
			/* Set entry as first element */
			remove_entry_from_list (cache, pos);
			cache->free = cache->entries[pos].next;
			cache->entries[pos].next = cache->first;
			cache->first = pos;
			if (cache->entries[pos].next >= 0) cache->entries[cache->entries[pos].next].prev = pos;
			cache->entries[pos].prev = -1;
			if (cache->last < 0) cache->last = pos;
			return;
		} else {
			/* Remove existing object */
			remove_entry_from_list (cache, pos);
			arikkei_dict_remove (&cache->dict, cache->entries[pos].key);
			if (cache->key_free) cache->key_free (cache->entries[pos].key);
			if (cache->object_free) cache->object_free (cache->entries[pos].object);
			cache->currentsize -= cache->entries[pos].size;
		}
	} else {
		if (cache->free < 0) {
			unsigned int i, newnentries;
			newnentries = (cache->nentries) ? cache->nentries << 1 : 32;
			cache->entries = (ArikkeiCacheEntry *) realloc (cache->entries, newnentries * sizeof (ArikkeiCacheEntry));
			for (i = cache->nentries; i < newnentries - 1; i++) cache->entries[i].next = i + 1;
			cache->entries[newnentries - 1].next = -1;
			cache->free = cache->nentries;
			cache->nentries = newnentries;
		}
	}
	/* Ensure we have enough space */
	/* It would be nice to pick new position here */
	arikkei_cache_ensure_space (cache, size);
	/* Set up entry */
	pos = cache->free;
	cache->free = cache->entries[pos].next;
	cache->entries[pos].key = (cache->key_dup) ? cache->key_dup (key) : (void *) key;
	cache->entries[pos].object = object;
	cache->entries[pos].size = size;
	/* Attach entry */
	cache->entries[pos].next = cache->first;
	cache->first = pos;
	if (cache->entries[pos].next >= 0) cache->entries[cache->entries[pos].next].prev = pos;
	cache->entries[pos].prev = -1;
	if (cache->last < 0) cache->last = pos;
	arikkei_dict_insert (&cache->dict, cache->entries[pos].key, (const char *) 0 + pos + 1);
	cache->currentsize += size;
}

void
arikkei_cache_remove (ArikkeiCache *cache, const void *key)
{
	int pos;
	pos = (int) ((const char *) arikkei_dict_lookup (&cache->dict, key) - (const char *) 0) - 1;
	if (pos < 0) return;
	arikkei_cache_remove_entry (cache, pos);
}

const void *
arikkei_cache_lookup (ArikkeiCache *cache, const void *key)
{
	int pos;
	pos = (int) ((const char *) arikkei_dict_lookup (&cache->dict, key) - (const char *) 0) - 1;
	if (pos >= 0) {
		remove_entry_from_list (cache, pos);
		cache->free = cache->entries[pos].next;
		cache->entries[pos].next = cache->first;
		cache->first = pos;
		if (cache->entries[pos].next >= 0) cache->entries[cache->entries[pos].next].prev = pos;
		cache->entries[pos].prev = -1;
		if (cache->last < 0) cache->last = pos;
	}
	return (pos >= 0) ? cache->entries[pos].object : NULL;
}

const void *
arikkei_cache_lookup_notouch (ArikkeiCache *cache, const void *key)
{
	int pos;
	pos = (int) ((const char *) arikkei_dict_lookup (&cache->dict, key) - (const char *) 0) - 1;
	return (pos >= 0) ? cache->entries[pos].object : NULL;
}

