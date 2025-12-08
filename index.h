// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INDEX_H
#define INDEX_H

#include "khash.h"
#include "snames.h"

KHASH_INIT(index, char*,int64_t, 1, kh_str_hash_func, kh_str_hash_equal)
#define index_t khash_t(index)

// use cleanIndex if wishes to clean the key char*s
#define freeIndex(idx) kh_destroy(index, idx)

typedef struct {
    char* key;
    int64_t value;
} index_pair_t;

/**
 * Extracts the index filename from the given filename of the form "fname.cx".
 *
 * @param fname_cx The base filename from which to extract the index filename.
 * @return A dynamically allocated string containing the extracted index filename.
 *         The caller is responsible for freeing the memory using `free()`.
 *         Returns NULL if the input filename is NULL or does not match the expected format.
 *
 * Example Usage:
 * ```c
 * char *fname_index = get_fname_index("data.cx");
 * printf("Index filename: %s\n", fname_index);
 * free(fname_index);  // Remember to free the dynamically allocated memory
 * ```
 */
char *get_fname_index(const char *fname_cx);

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
 *         The caller is responsible for freeing the memory, particularly the char* keys
 * Use cleanIndex instead of freeIndex. The function owns the key memories.
 */
index_t* loadIndex(char* fname_index);

/**
 * @brief This function clears the memory occupied by an index table and its keys.
 *
 * The function iterates over all keys in the provided index table. For each existing key,
 * it frees the memory allocated for the key. After all keys are processed, the function
 * frees the memory allocated for the index table itself.
 *
 * @note This function should only be used when the keys in the index table are dynamically
 * allocated strings. Using this function with index tables that contain statically allocated
 * keys (e.g., string literals) will result in undefined behavior.
 *
 * @param idx A pointer to an index_t object that needs to be cleaned up.
 * @return void
 */
static inline void cleanIndex(index_t *idx) {
  if (idx == NULL) return;
  khiter_t k;
  for (k = kh_begin(idx); k != kh_end(idx); ++k) {
    if (kh_exist(idx, k)) {
      free(kh_key(idx, k));  // frees the string key
    }
  }
  freeIndex(idx);
}

/*************************************
 ** Retrieves the value for a given   **
 ** sample name from the index.       **
 *************************************
 *
 * @param index The index hash table (index_t).
 * @param sname The sample name to retrieve the index for.
 * @return The index value for the given sample name if found, or -1 if not found.
 */
int64_t getIndex(index_t* index, char* sname);

/**
 * Retrieves the key-value pairs from the given index.
 *
 * @param idx The index_t instance from which to retrieve index pairs.
 * @param n number of index pairs successfully read
 * @return A pointer to an array of index_pair_t structures representing the key-value pairs in the given index.
 */
index_pair_t *index_pairs(index_t *idx, int *n);
/* load index pairs, instead of the index */
index_pair_t* load_index_pairs(char *fname_cx, int *n);
void clean_index_pairs(index_pair_t *idx_pairs, int n);

/**
 * Writes an index_t instance to a FILE stream.
 * 
 * This function takes a pointer to a FILE stream and an index_t instance,
 * and writes the data from the index_t instance to the FILE stream. 
 * The data is expected to be in a specific format, matching the structure of the index_t type.
 *
 * @param fp A pointer to the FILE stream to write to.
 * @param idx A pointer to the index_t instance to be written.
 */
void writeIndex(FILE *fp, index_t *idx);

/**
 * Inserts a new element into an index_t structure. 
 * 
 * This function takes a pointer to an index_t structure, a string representing 
 * the sample name (sname), and a 64-bit integer representing the address (addr). 
 * It inserts these values as a new element into the provided index_t structure.
 * If the sample name (sname) already exists in the index, the function will terminate 
 * the program with an error message.
 *
 * @param idx A pointer to the index_t structure to which the new element will be added. 
 * If this pointer is NULL, it will cause undefined behavior.
 * @param sname A string representing the sample name. This value is inserted into the index. 
 * If the sample name already exists in the index, the function will print an error message 
 * to stderr and terminate the program.
 * @param addr A 64-bit integer representing the address. This value is inserted into the index.
 *
 * @return Returns a pointer to the updated index_t structure. 
 * If there is a problem with the insertion (e.g., failure to allocate memory), 
 * the function prints an error message to stdout and terminates the program.
 */
index_t *insert_index(index_t *idx, char *sname, int64_t addr);

snames_t loadSampleNamesFromIndex(char *fname);

#endif
