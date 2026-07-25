#!/bin/bash
set -euo pipefail

# `source: path: ..` copies the working tree as it is on disk, build products
# and all. libyame.a, htslib/libhts.a and the .o they were made from are built
# for the author's machine; linking them into a package for another
# architecture produces a binary that fails at load time rather than at build
# time. Drop them so they rebuild for THIS target.
make clean >/dev/null 2>&1 || true
rm -f libyame.a htslib/libhts.a
find . -name '*.o' -delete 2>/dev/null || true

# conda supplies CC plus the sysroot and $PREFIX include/lib flags in CPPFLAGS
# / CFLAGS / LDFLAGS. YAME's Makefile and the vendored htslib's both assign
# CFLAGS rigidly and ignore CPPFLAGS, so bake conda's compile flags into CC --
# a command-line CC propagates to every sub-make. LDFLAGS is expanded by the
# link rule explicitly; make imports the environment into a variable but an
# explicit recipe never expands it by itself, and conda-forge puts
# -L$PREFIX/lib only in LDFLAGS, so without this a linux-64 build fails to
# find -lz.
#
# CURL is named rather than left to the default `auto`. Detection would find
# the curl-config from the libcurl host dependency on PATH, but naming it
# removes the guesswork in a cross build, where PATH may offer the builder's
# curl-config instead of the target's. Since v1.29 a failed detection stops
# the build rather than dropping fetching silently, so a mistake here is loud
# either way -- this just makes it not happen.
make -j"${CPU_COUNT:-1}" \
     CC="${CC} ${CPPFLAGS:-} ${CFLAGS:-}" \
     LDFLAGS="${LDFLAGS:-}" \
     CURL="${PREFIX}/bin/curl-config"

install -d "${PREFIX}/bin"
install -v -m 0755 yame "${PREFIX}/bin"
