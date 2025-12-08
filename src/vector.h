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

#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief A dynamic array implementation for storing elements.
 */
typedef struct {
  char** data;    /**< Pointer to the underlying data array */
  size_t size;    /**< Number of elements in the vector */
  size_t capacity;    /**< Capacity of the vector */
} vector_t;

/**
 * @brief Initializes a new vector.
 * @return Pointer to the newly initialized vector.
 */
static inline vector_t* vector_init() {
  vector_t* vector = malloc(sizeof(vector_t));
  if (vector == NULL) {
    return NULL;
  }
  vector->data = NULL;
  vector->size = 0;
  vector->capacity = 0;
  return vector;
}

/**
 * @brief Adds an element to the vector.
 * @param vector Pointer to the vector.
 * @param element The element to be added.
 */
static inline void vector_push(vector_t* vector, const char* element) {
  if (vector->size == vector->capacity) {
    size_t new_capacity = (vector->capacity == 0) ? 1 : vector->capacity * 2;
    char** new_data = realloc(vector->data, new_capacity * sizeof(char*));
    if (new_data == NULL) {
      return;
    }
    vector->data = new_data;
    vector->capacity = new_capacity;
  }
  vector->data[vector->size] = strdup(element);
  vector->size++;
}

/**
 * @brief Retrieves an element from the vector at a specific index.
 * @param vector Pointer to the vector.
 * @param index The index of the element to retrieve.
 * @return The element at the specified index, or NULL if the index is out of range.
 */
static inline char* vector_get(vector_t* vector, size_t index) {
  if (index >= vector->size) {
    return NULL;
  }
  return vector->data[index];
}

/**
 * @brief Frees the memory allocated for the vector.
 * @param vector Pointer to the vector to be destroyed.
 */
static inline void vector_destroy(vector_t* vector) {
  for (size_t i = 0; i < vector->size; i++) {
    free(vector->data[i]);
  }
  free(vector->data);
  free(vector);
}

#endif /* VECTOR_H */
