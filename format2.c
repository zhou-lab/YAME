#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "khash.h"

static uint8_t* compressDataToRLE(uint64_t *data, uint64_t n, uint64_t *rle_n) {
  // Calculate the maximum value
  uint64_t max_value = 0;
  for (uint64_t i = 0; i < n; ++i) {
    uint64_t value = data[i];
    if (value > max_value) {
      max_value = value;
    }}
  
  // Determine the number of bytes needed to encode each value
  int value_bytes;
  if (max_value < (1<<8)) value_bytes = 1;
  else if (max_value < (1<<16)) value_bytes = 2;
  else if (max_value < (1<<24)) value_bytes = 3;
  else value_bytes = 8;

  // Create a buffer for the RLE data
  // value_bytes for the value and 2 bytes for the count
  uint8_t *rle = NULL; *rle_n = 0;

  // Write the number of bytes for each value into the RLE data
  rle = realloc(rle, ((*rle_n)+1));
  rle[(*rle_n)++] = value_bytes;

  // Encode the array into the RLE format
  for (uint64_t i = 0; i < n;) {
    // Get the current value and count
    uint64_t value = data[i];
    uint64_t count = 1;
    while (i + count < n && data[i + count] == value && count < ((1<<16)-1)) {
      ++count;
    }

    // Write the value and count into the RLE data
    rle = realloc(rle, ((*rle_n) + value_bytes + 2));
    memcpy(rle + *rle_n, &value, value_bytes);
    (*rle_n) += value_bytes;
    uint8_t count16[2];
    count16[0] = (count & 0xff);
    count16[1] = (count >> 8) & 0xff;
    memcpy(rle + *rle_n, count16, 2);
    (*rle_n) += 2;

    // Move to the next (value, count) pair
    i += count;
  }
  
  return rle;
}

KHASH_MAP_INIT_STR(str2int, uint64_t) // Initialize a hashmap with keys as strings and values as uint64_t

cdata_t* fmt2_read_raw(char *fname, int verbose) {
  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t *data = calloc(1<<10, sizeof(uint64_t));
  uint64_t data_n = 0, data_m = 1<<10;
  khash_t(str2int) *h = kh_init(str2int); // Initialize the hashmap
  khint_t k;
  uint64_t keys_count = 0;
  char **keys = calloc(1<<10, sizeof(char*));
  uint64_t keys_n = 0, keys_m = 1<<10;
  while (gzFile_read_line(fh, &line) > 0) {
    int ret;
    k = kh_put(str2int, h, strdup(line), &ret);
    if (ret) {  // The key didn't exist before
      kh_val(h, k) = keys_count++;
      if (keys_n+1>keys_m) {
        keys_m <<= 1;
        keys = realloc(keys, keys_m * sizeof(char*));
      }
      keys[keys_n++] = strdup(line);
    }
    if (data_n+1>data_m) {
      data_m <<= 1;
      data = realloc(data, data_m * sizeof(uint64_t));
    }
    data[data_n++] = kh_val(h, k);
  }
  free(line);
  wzclose(fh);

  // Calculate key string size
  uint64_t keys_n_bytes = 0;
  for (uint64_t i = 0; i < keys_n; ++i)
    keys_n_bytes += strlen(keys[i]) + 1;

  // Create a new cdata_t structure
  cdata_t *c = calloc(1, sizeof(cdata_t));
  c->compressed = 0;
  c->fmt = '2';
  c->aux = calloc(1, sizeof(f2_aux_t));
  c->n = data_n;
  c->s = calloc(1, keys_n_bytes + data_n*sizeof(uint64_t) + 1);

  // Write the keys to the data
  f2_aux_t *aux = (f2_aux_t*) c->aux;
  aux->nk = keys_n;
  aux->keys = calloc(keys_n, sizeof(char*));
  uint64_t pos = 0;
  for (uint64_t i = 0; i < keys_n; ++i) {
    size_t len = strlen(keys[i]);
    memcpy(c->s + pos, keys[i], len);
    aux->keys[i] = (char*) (c->s + pos);
    pos += len;
    c->s[pos++] = '\0';  // Separate the keys by null characters
  }

  // Write a separator between the keys and the data
  c->s[pos++] = '\0';

  // Write the data after the keys
  memcpy(c->s + pos, data, data_n*sizeof(uint64_t));

  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, data_n);
    fflush(stderr);
  }

  free(data);

  kh_destroy(str2int, h);  // Destroy the hashmap
  for (uint64_t i = 0; i < keys_n; ++i) {
    free(keys[i]);
  }
  free(keys);

  return c;
}

uint64_t fmt2_get_keys_n(cdata_t *c) {
  uint64_t keys_n = 0;
  for (uint64_t i = 0; ; ++i) {
    if (c->s[i] == '\0') {
      // Count null characters to determine the number of keys
      keys_n++;
    }
    if (c->s[i] == '\0' && c->s[i+1] == '\0') {
      // If we encounter the separator (two null characters in a row), we stop
      break;
    }
  }
  // Subtract 1 to exclude the separator null character
  return keys_n;
}

static uint64_t fmt2_get_keys_nbytes(cdata_t *c) {
  uint64_t i;
  for (i = 0; ; ++i) {
    if (c->s[i] == '\0' && c->s[i+1] == '\0') {
      // If we encounter the separator (two null characters in a row), we stop
      break;
    }
  }
  return i+1;                   /* not counting the 2nd \0 */
}

// c->n is the total nbytes only when compressed
static uint64_t fmt2_get_data_nbytes(cdata_t *c) {
  if (!c->compressed) {
    fprintf(stderr, "[%s:%d] Data is uncompressed.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  uint64_t separator_idx = 0;
  for (uint64_t i = 0; i < c->n; ++i) {
    if (c->s[i] == '\0' && c->s[i+1] == '\0') {
      // If we encounter the separator (two null characters in a row), we store its position
      separator_idx = i;
      break;
    }
  }
  // Subtract the position of the separator and 1 for the value byte count from the total size to get the data size
  return c->n - separator_idx - 1 - 1;
}

// assume c is compressed
static uint8_t fmt2_get_value_byte(cdata_t *c) {
  if (!c->compressed) {
    fprintf(stderr, "[%s:%d] Data is uncompressed.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  uint64_t separator_idx = 0;
  for (uint64_t i = 0; i < c->n; ++i) {
    if (c->s[i] == '\0' && c->s[i+1] == '\0') {
      // If we encounter the separator (two null characters in a row), we store its position
      separator_idx = i;
      break;
    }
  }
  // The value byte count is the byte right after the separator
  return c->s[separator_idx + 2];
}

uint8_t* fmt2_get_data(cdata_t *c) {
  uint64_t separator_idx = 0;
  for (uint64_t i = 0; ; ++i) {
    if (c->s[i] == '\0' && c->s[i+1] == '\0') {
      // If we encounter the separator (two null characters in a row), we store its position
      separator_idx = i;
      break;
    }
  }
  // The data starts right after the separator and the value byte count
  return c->s + separator_idx + 1 + 1;
}

void fmt2_compress(cdata_t *c) {
  uint64_t keys_nb = fmt2_get_keys_nbytes(c);
  uint8_t *data = fmt2_get_data(c);
  uint64_t rle_n;
  uint8_t *rle_data = compressDataToRLE((uint64_t*) data, c->n, &rle_n);
  uint8_t *s_out = calloc(keys_nb + rle_n + 1, sizeof(uint8_t));
  memcpy(s_out, c->s, keys_nb + 1);
  memcpy(s_out + keys_nb + 1, rle_data, rle_n);
  c->s = s_out;
  c->n = keys_nb + rle_n + 1;
  c->compressed = 1;
  c->fmt = '2';
}

void fmt2_decompress(cdata_t *c, cdata_t *inflated) {
  uint64_t keys_nb = fmt2_get_keys_nbytes(c);
  uint8_t *keys = c->s;
  uint8_t *data = fmt2_get_data(c) + 1; // skip value byte
  uint64_t data_nbyte = fmt2_get_data_nbytes(c) - 1;
  inflated->unit = fmt2_get_value_byte(c);

  // Iterate over RLE data to calculate total length of decompressed data
  uint64_t dec_data_n = 0;
  for (uint64_t i = 0; i < data_nbyte; ) {
    i += inflated->unit;
    uint64_t length = ((uint64_t) data[i] | (uint64_t) (data[i+1] << 8));
    i += 2;
    dec_data_n += length;
  }

  // Allocate memory directly to inflated->s
  inflated->n = keys_nb + dec_data_n * inflated->unit + 1;
  inflated->s = malloc(inflated->n);
  if (inflated->s == NULL) {
    fprintf(stderr, "Memory allocation failed. Exiting.\n");
    exit(1);
  }

  // Copy keys to inflated->s
  memcpy(inflated->s, keys, keys_nb);
  inflated->s[keys_nb] = '\0';  // NULL separator

  // Have dec_data point to the data part
  uint64_t sum = 0;
  uint8_t *dec_data = inflated->s + keys_nb + 1;
  uint64_t n = 0;
  for (uint64_t i = 0; i < data_nbyte; ) {
    uint64_t value = 0; uint8_t *d = data+i;
    for (uint64_t j = 0; j < inflated->unit; ++j)
      value |= data[i++] << (8*j);
    uint64_t length = ((uint64_t) data[i] | (uint64_t) (data[i+1] << 8));
    sum += length;
    i += 2;
    for ( ; length--; ++n)
      memcpy(dec_data+n*inflated->unit, d, inflated->unit);
  }

  inflated->compressed = 0;
  inflated->fmt = '2';
  inflated->n = n;
}

void fmt2_set_aux(cdata_t *c) {
  if (c->aux != NULL) {
    fprintf(stderr, "[%s:%d] Aux data exists.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  // Create a keys_t object and allocate memory for s
  f2_aux_t *aux = calloc(1, sizeof(f2_aux_t));
  aux->nk = fmt2_get_keys_n(c);
  aux->keys = (char **)malloc(aux->nk * sizeof(char *));

  char *key_start = (char *)c->s;
  char *key_end;
  uint64_t idx = 0;

  // Iterate through the keys
  uint64_t keys_n_bytes = fmt2_get_keys_nbytes(c);
  while (key_start < (char *)(c->s + keys_n_bytes)) {
    key_end = strchr(key_start, '\0');
    aux->keys[idx++] = key_start;
    key_start = key_end + 1;
  }
  aux->data = fmt2_get_data(c);
  c->aux = aux;
}
