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

/**
 * SHA-256 (FIPS 180-4).
 *
 * Self-contained on purpose: OpenSSL would be a whole dependency for one hash
 * function and CommonCrypto is macOS-only, so a couple of hundred lines is
 * cheaper than either and keeps the bioconda surface at htslib + zlib +
 * libcurl. Ported from kycg's src/digest.c, which in turn shares its structure
 * with sesame-cli's src/sha256.c -- three copies of this existed across the
 * suite before it moved here.
 *
 * SHA-256 only. Every channel YAME fetches from publishes a SHA256SUMS per
 * directory, anchored by a digest compiled into the calling binary, so nothing
 * here ever needs a weaker hash. (kycg carried MD5 for as long as the
 * whole-genome sets came from Zenodo, whose record API publishes nothing
 * stronger; that channel now publishes its own SHA256SUMS and the MD5 half was
 * dropped rather than ported.)
 */

#include "assets.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define YAME_DIGEST_BUFSZ 65536

typedef struct {
  uint32_t h[8];
  uint64_t len;
  uint8_t  buf[64];
  size_t   n;
} sha256_ctx;

static const uint32_t SHA256_K[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

#define ROR32(x,n) (((x) >> (n)) | ((x) << (32 - (n))))

static void sha256_block(sha256_ctx *c, const uint8_t *p) {
  uint32_t w[64], a, b, cc, d, e, f, g, h;
  int i;

  for (i = 0; i < 16; ++i)
    w[i] = ((uint32_t)p[i*4] << 24) | ((uint32_t)p[i*4+1] << 16) |
           ((uint32_t)p[i*4+2] << 8) | (uint32_t)p[i*4+3];
  for (; i < 64; ++i) {
    uint32_t s0 = ROR32(w[i-15],7) ^ ROR32(w[i-15],18) ^ (w[i-15] >> 3);
    uint32_t s1 = ROR32(w[i-2],17) ^ ROR32(w[i-2],19) ^ (w[i-2] >> 10);
    w[i] = w[i-16] + s0 + w[i-7] + s1;
  }

  a=c->h[0]; b=c->h[1]; cc=c->h[2]; d=c->h[3];
  e=c->h[4]; f=c->h[5]; g=c->h[6];  h=c->h[7];

  for (i = 0; i < 64; ++i) {
    uint32_t S1 = ROR32(e,6) ^ ROR32(e,11) ^ ROR32(e,25);
    uint32_t ch = (e & f) ^ (~e & g);
    uint32_t t1 = h + S1 + ch + SHA256_K[i] + w[i];
    uint32_t S0 = ROR32(a,2) ^ ROR32(a,13) ^ ROR32(a,22);
    uint32_t mj = (a & b) ^ (a & cc) ^ (b & cc);
    uint32_t t2 = S0 + mj;
    h=g; g=f; f=e; e=d+t1; d=cc; cc=b; b=a; a=t1+t2;
  }

  c->h[0]+=a; c->h[1]+=b; c->h[2]+=cc; c->h[3]+=d;
  c->h[4]+=e; c->h[5]+=f; c->h[6]+=g;  c->h[7]+=h;
}

static void sha256_init(sha256_ctx *c) {
  c->h[0]=0x6a09e667; c->h[1]=0xbb67ae85; c->h[2]=0x3c6ef372; c->h[3]=0xa54ff53a;
  c->h[4]=0x510e527f; c->h[5]=0x9b05688c; c->h[6]=0x1f83d9ab; c->h[7]=0x5be0cd19;
  c->len = 0; c->n = 0;
}

static void sha256_update(sha256_ctx *c, const void *data, size_t len) {
  const uint8_t *p = data;
  c->len += len;
  while (len) {
    size_t take = 64 - c->n;
    if (take > len) take = len;
    memcpy(c->buf + c->n, p, take);
    c->n += take; p += take; len -= take;
    if (c->n == 64) { sha256_block(c, c->buf); c->n = 0; }
  }
}

static void sha256_final(sha256_ctx *c, uint8_t out[32]) {
  uint64_t bits = c->len * 8;
  int i;

  c->buf[c->n++] = 0x80;
  if (c->n > 56) {
    memset(c->buf + c->n, 0, 64 - c->n);
    sha256_block(c, c->buf);
    c->n = 0;
  }
  memset(c->buf + c->n, 0, 56 - c->n);
  for (i = 0; i < 8; ++i) c->buf[56+i] = (uint8_t)(bits >> (56 - 8*i));
  sha256_block(c, c->buf);

  for (i = 0; i < 8; ++i) {
    out[i*4]   = (uint8_t)(c->h[i] >> 24);
    out[i*4+1] = (uint8_t)(c->h[i] >> 16);
    out[i*4+2] = (uint8_t)(c->h[i] >> 8);
    out[i*4+3] = (uint8_t)(c->h[i]);
  }
}

static void hexify(const uint8_t *raw, size_t n, char *out) {
  static const char HEX[] = "0123456789abcdef";
  size_t i;
  for (i = 0; i < n; ++i) {
    out[i*2]   = HEX[raw[i] >> 4];
    out[i*2+1] = HEX[raw[i] & 0xf];
  }
  out[n*2] = '\0';
}

void yame_assets_sha256_buf(const void *data, size_t len, char out[65]) {
  sha256_ctx c;
  uint8_t raw[32];
  sha256_init(&c);
  sha256_update(&c, data, len);
  sha256_final(&c, raw);
  hexify(raw, 32, out);
}

int yame_assets_sha256_file(const char *path, char out[65]) {
  FILE *fp = fopen(path, "rb");
  if (!fp) return -1;

  sha256_ctx c;
  sha256_init(&c);

  static uint8_t buf[YAME_DIGEST_BUFSZ];
  size_t got;
  while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) sha256_update(&c, buf, got);

  int err = ferror(fp);
  fclose(fp);
  if (err) return -1;

  uint8_t raw[32];
  sha256_final(&c, raw);
  hexify(raw, 32, out);
  return 0;
}

int yame_assets_digest_equal(const char *a, const char *b) {
  if (!a || !b) return 0;
  size_t i;
  for (i = 0; a[i] && b[i]; ++i)
    if (tolower((unsigned char)a[i]) != tolower((unsigned char)b[i])) return 0;
  return a[i] == '\0' && b[i] == '\0';
}
