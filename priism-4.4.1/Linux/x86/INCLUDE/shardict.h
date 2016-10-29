#ifndef UCSF_MSG_SHARDICT_H
#define UCSF_MSG_SHARDICT_H

/*
 * Implements a dictionary on the IW shared memory segment; the key
 * and data are assumed to be accessible from the same pointer.
 */

struct ShardictImpl;
typedef struct ShardictImpl* Shardict;
typedef const struct ShardictImpl* ShardictRO;
struct ShardictNodeImpl;
typedef struct ShardictNodeImpl* ShardictIter;
typedef const struct ShardictNodeImpl* ShardictIterRO;

/*
 * Returns -1 if arg1 < arg2, 0 if arg1 == arg2, and 1 if arg1 > arg2.
 */
typedef int (*ShardictComparer)(const void*, const void*);
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Returns NULL if it was not able to create the dictionary.
 */
Shardict shardict_create(void);

/*
 * Removes all entries in the dictionary and deletes the dictionary
 * itselt.  For each entry, destruct_cb is invoked and passed the
 * pointer to the data/key for the entry.
 */
void shardict_destroy(Shardict dict, void (*destruct_cb)(void*));

/*
 * Removes all entries in the dictionary.  For each entry, destruct_cb
 * is invoked and passed the pointer to the data for the entry.  As a
 * side effect, invalidates all iterators pointing to the dictionary.
 */
void shardict_clear(Shardict dict, void (*destruct_cb)(void*));

/*
 * If key already exists in the dictionary, returns an iterator to it
 * and sets *p_found to 1.  If key does not exist, attempts to create
 * an entry for it; sets *p_found to zero in that case and returns
 * an iterator to the new entry (an invalid iterator if the creation
 * failed).
 */
ShardictIter shardict_probe(
    Shardict dict, const void* key, ShardictComparer comparer, int* p_found
);

/*
 * Removes key if exists in the dictionary.  Invalidates any iterator
 * pointing to the entry with that key.  Returns the pointer to
 * the data/key or NULL if key was not found.
 */
void* shardict_remove(
    Shardict dict, const void* key, ShardictComparer comparer
);

/*
 * Searches the dictionary for given key and returns an iterator to it.
 * If the key could not be found, the iterator returned is invalid.
 */
ShardictIter shardict_search(
    Shardict dict, const void* key, ShardictComparer comparer
);
ShardictIterRO shardict_search_const(
    ShardictRO dict, const void* key, ShardictComparer comparer
);

/*
 * Returns an iterator to the entry which has the lowest key value. If the
 * dictionary is empty, the returned iterator is equal to the result of
 * shardict_invalid_iter()
 */
ShardictIter shardict_min(Shardict dict);
ShardictIterRO shardict_min_const(ShardictRO dict);

/*
 * Returns an iterator to the entry which has the highest key value. If the
 * dictionary is empty, the returned iterator is equal to the result of
 * shardict_invalid_iter()
 */
ShardictIter shardict_max(Shardict dict);
ShardictIterRO shardict_max_const(ShardictRO dict);

/*
 * Returns a non-dereferenceable iterator to the dictionary.
 */
ShardictIter shardict_invalid_iter(Shardict dict);
ShardictIterRO shardict_invalid_iter_const(ShardictRO dict);

/*
 * Returns true if the iterator is dereferenceable, false if not.
 */
int shardict_iter_valid(ShardictIterRO iter);

/*
 * Returns true if the dict is not corrupt and false if not.
 */
int shardict_verify_invariants(ShardictRO dict, ShardictComparer comparer);

/*
 * If iter is dereferenceable, returns an iterator to the dictionary
 * entry with the next higher key or an iterator equal to the result of
 * shardict_invalid_iter() if iter is equal to the result of
 * shardict_max().  If the iterator is not dereferenceable,
 * returns an iterator equal to the result of shardict_min();
 */
ShardictIter shardict_iter_next(ShardictIter iter);
ShardictIterRO shardict_iter_next_const(ShardictIterRO iter);

/*
 * If iter is dereferenceable, returns an iterator to the dictionary
 * entry with the next lower key or an iterator equal to the result of
 * shardict_invalid_iter() if iter is equal to the result of
 * shardict_min().  If the iterator is not dereferenceable,
 * returns an iterator equal to the result of shardict_max();
 */
ShardictIter shardict_iter_prev(ShardictIter iter);
ShardictIterRO shardict_iter_prev_const(ShardictIterRO iter);

/*
 * Sets the data/key pointer for the dictionary entry pointed to by iter
 * (beware: changing the key can invalidate the ordering of the dictionary).
 * Returns the pointer previously associated with the entry.  The
 * result is undefined if iter is not dereferenceable.
 */
void* shardict_iter_set_data(ShardictIter iter, void* data);

/*
 * Returns a pointer to the key/data for the entry pointed to by the iterator.
 * The result is undefined if the iterator is not dereferenceable.
 */
void* shardict_iter_data(ShardictIter iter);
const void* shardict_iter_data_const(ShardictIterRO iter);

/*
 * Remove the entry pointed to by the given iterator and returns the
 * pointer to the entry's key/data.  The result is undefined if the iterator
 * is not dereferenceable.
 */
void* shardict_iter_remove(ShardictIter iter);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
