#ifndef UCSF_MSG_PRIVDICT_H
#define UCSF_MSG_PRIVDICT_H

/*
 * Implements a dictionary on the standard C heap (malloc, free); the key
 * and data are assumed to be accessible from the same pointer.
 */

struct PrivdictImpl;
typedef struct PrivdictImpl* Privdict;
typedef const struct PrivdictImpl* PrivdictRO;
struct PrivdictNodeImpl;
typedef struct PrivdictNodeImpl* PrivdictIter;
typedef const struct PrivdictNodeImpl* PrivdictIterRO;

/*
 * Returns -1 if arg1 < arg2, 0 if arg1 == arg2, and 1 if arg1 > arg2.
 */
typedef int (*PrivdictComparer)(const void*, const void*);
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Returns NULL if it was not able to create the dictionary.
 */
Privdict privdict_create(void);

/*
 * Removes all entries in the dictionary and deletes the dictionary
 * itselt.  For each entry, destruct_cb is invoked and passed the
 * pointer to the data/key for the entry.
 */
void privdict_destroy(Privdict dict, void (*destruct_cb)(void*));

/*
 * Removes all entries in the dictionary.  For each entry, destruct_cb
 * is invoked and passed the pointer to the data for the entry.  As a
 * side effect, invalidates all iterators pointing to the dictionary.
 */
void privdict_clear(Privdict dict, void (*destruct_cb)(void*));

/*
 * If key already exists in the dictionary, returns an iterator to it
 * and sets *p_found to 1.  If key does not exist, attempts to create
 * an entry for it; sets *p_found to zero in that case and returns
 * an iterator to the new entry (an invalid iterator if the creation
 * failed).
 */
PrivdictIter privdict_probe(
    Privdict dict, const void* key, PrivdictComparer comparer, int* p_found
);

/*
 * Removes key if exists in the dictionary.  Invalidates any iterator
 * pointing to the entry with that key.  Returns the pointer to
 * the data/key or NULL if key was not found.
 */
void* privdict_remove(
    Privdict dict, const void* key, PrivdictComparer comparer
);

/*
 * Searches the dictionary for given key and returns an iterator to it.
 * If the key could not be found, the iterator returned is invalid.
 */
PrivdictIter privdict_search(
    Privdict dict, const void* key, PrivdictComparer comparer
);
PrivdictIterRO privdict_search_const(
    PrivdictRO dict, const void* key, PrivdictComparer comparer
);

/*
 * Returns an iterator to the entry which has the lowest key value. If the
 * dictionary is empty, the returned iterator is equal to the result of
 * privdict_invalid_iter()
 */
PrivdictIter privdict_min(Privdict dict);
PrivdictIterRO privdict_min_const(PrivdictRO dict);

/*
 * Returns an iterator to the entry which has the highest key value. If the
 * dictionary is empty, the returned iterator is equal to the result of
 * privdict_invalid_iter()
 */
PrivdictIter privdict_max(Privdict dict);
PrivdictIterRO privdict_max_const(PrivdictRO dict);

/*
 * Returns a non-dereferenceable iterator to the dictionary.
 */
PrivdictIter privdict_invalid_iter(Privdict dict);
PrivdictIterRO privdict_invalid_iter_const(PrivdictRO dict);

/*
 * Returns true if the iterator is dereferenceable, false if not.
 */
int privdict_iter_valid(PrivdictIterRO iter);

/*
 * Returns true if the dict is not corrupt and false if not.
 */
int privdict_verify_invariants(PrivdictRO dict, PrivdictComparer comparer);

/*
 * If iter is dereferenceable, returns an iterator to the dictionary
 * entry with the next higher key or an iterator equal to the result of
 * privdict_invalid_iter() if iter is equal to the result of
 * privdict_max().  If the iterator is not dereferenceable,
 * returns an iterator equal to the result of privdict_min();
 */
PrivdictIter privdict_iter_next(PrivdictIter iter);
PrivdictIterRO privdict_iter_next_const(PrivdictIterRO iter);

/*
 * If iter is dereferenceable, returns an iterator to the dictionary
 * entry with the next lower key or an iterator equal to the result of
 * privdict_invalid_iter() if iter is equal to the result of
 * privdict_min().  If the iterator is not dereferenceable,
 * returns an iterator equal to the result of privdict_max();
 */
PrivdictIter privdict_iter_prev(PrivdictIter iter);
PrivdictIterRO privdict_iter_prev_const(PrivdictIterRO iter);

/*
 * Sets the data/key pointer for the dictionary entry pointed to by iter
 * (beware: changing the key can invalidate the ordering of the dictionary).
 * Returns the pointer previously associated with the entry.  The
 * result is undefined if iter is not dereferenceable.
 */
void* privdict_iter_set_data(PrivdictIter iter, void* data);

/*
 * Returns a pointer to the key/data for the entry pointed to by the iterator.
 * The result is undefined if the iterator is not dereferenceable.
 */
void* privdict_iter_data(PrivdictIter iter);
const void* privdict_iter_data_const(PrivdictIterRO iter);

/*
 * Remove the entry pointed to by the given iterator and returns the
 * pointer to the entry's key/data.  The result is undefined if the iterator
 * is not dereferenceable.
 */
void* privdict_iter_remove(PrivdictIter iter);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
