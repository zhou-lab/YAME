#ifndef INDEX_H
#define INDEX_H

#include "khash.h"

KHASH_MAP_INIT_STR(index, int64_t)  
#define index_t khash_t(index)
#define destroyIndex(idx) kh_destroy(index, idx)

typedef struct {
    const char* key;
    int64_t value;
} index_pair_t;

/**
 * Extracts the index filename from the given filename of the form "fname.cg".
 *
 * @param fname_cg The base filename from which to extract the index filename.
 * @return A dynamically allocated string containing the extracted index filename.
 *         The caller is responsible for freeing the memory using `free()`.
 *         Returns NULL if the input filename is NULL or does not match the expected format.
 *
 * Example Usage:
 * ```c
 * char *fname_index = get_fname_index("data.cg");
 * printf("Index filename: %s\n", fname_index);
 * free(fname_index);  // Remember to free the dynamically allocated memory
 * ```
 */
char *get_fname_index(const char *fname_cg);

/*************************************
 ** Loads an index from a file.      **
 *************************************
 * The index file should be in the following format:
 *     sample1\tindex1
 *     sample2\tindex2
 *     ...
 *
 * @param fname_index The filename for the index file.
 * @return A pointer to the loaded index hash table (index_t) if successful, or NULL on failure.
 *         The caller is responsible for freeing the memory.
 */
index_t* loadIndex(char* fname_cg);

/*************************************
 ** Retrieves the value for a given   **
 ** sample name from the index.       **
 *************************************
 *
 * @param index The index hash table (index_t).
 * @param sname The sample name to retrieve the index for.
 * @return The index value for the given sample name if found, or -1 if not found.
 */
int64_t getIndex(index_t* index, const char* sname);

/**
 * Retrieves the key-value pairs from the given index.
 *
 * @param idx The index_t instance from which to retrieve index pairs.
 * @param n number of index pairs successfully read
 * @return A pointer to an array of index_pair_t structures representing the key-value pairs in the given index.
 */
index_pair_t *index_pairs(index_t *idx, int *n);


#endif
