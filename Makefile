# Compiler and flags
CC ?= gcc
CFLAGS = -W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result -O3
CLIB = -lpthread -lz -lm
CF_OPTIMIZE = 1

SRC_DIR = src
INCLUDE = include

OS := $(shell uname)
ifeq ($(OS),  Darwin)
	CFLAGS += -Wno-unused-function
else
	CLIB += -lrt
endif

# Program name
PROG = yame

# Static library for downstream linkers (e.g. MethScope2): all objects except
# the CLI dispatcher (src/main.o), which carries the lone int main().
LIB = libyame.a

# HTSlib
LHTSLIB_DIR = htslib
LHTSLIB = $(LHTSLIB_DIR)/libhts.a

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)

CFLAGS += -I$(SRC_DIR) -I$(LHTSLIB_DIR)

# libcurl, for `yame fetch` and the shared asset store (src/assets.c).
#
# A DEPENDENCY, not a bonus. It used to be detected and silently dropped, so a
# build environment without it produced a yame that looked complete and told
# nobody: `yame fetch` -- and with it the catalogue browser, `summary -b`, and
# resolving -R / -m against the store -- reported itself unavailable at run
# time, on a binary already shipped. A packager cannot notice what the build
# never mentioned, so a missing libcurl now stops the build and says what to
# do about it.
#
# Building without it stays possible, but only on purpose:
#
#   make CURL=0                     no libcurl; fetching is compiled out
#   make CURL=/path/to/curl-config  a specific one, if it is not on PATH
#   make CURL_LIBS="$(curl-config --libs)" CURL_CFLAGS="$(curl-config --cflags)"
#                                   spell it out, for a cross build where the
#                                   host's curl-config cannot be run
#
# libyame.a carries objects only, so a downstream link must still add -lcurl
# itself -- kycg and sesame-cli already do.
CURL ?= auto
ifeq ($(CURL),auto)
  CURL_CFLAGS := $(shell curl-config --cflags 2>/dev/null)
  CURL_LIBS   := $(shell curl-config --libs 2>/dev/null)
else ifeq ($(CURL),0)
  CURL_CFLAGS :=
  CURL_LIBS   :=
else
  CURL_CFLAGS := $(shell $(CURL) --cflags 2>/dev/null)
  CURL_LIBS   := $(shell $(CURL) --libs 2>/dev/null)
endif

ifneq ($(CURL_LIBS),)
  CFLAGS += -DYAME_HAVE_CURL $(CURL_CFLAGS)
  CLIB   += $(CURL_LIBS)
else ifneq ($(CURL),0)
  # Not while merely tidying up: `make clean` must work anywhere.
  ifeq ($(filter clean,$(MAKECMDGOALS)),)
    $(error libcurl not found: no curl-config on PATH. \
      yame fetch needs it. Install libcurl development files (conda: add \
      `libcurl` to the recipe host deps; Debian: libcurl4-openssl-dev; \
      RHEL: libcurl-devel), or point at one with `make CURL=/path/to/curl-config`, \
      or build deliberately without fetching using `make CURL=0`)
  endif
endif

.PHONY: all build debug clean lib

all: build

build: exportcf $(PROG) yame-config

# Build just the static library (used by MethScope2).
lib: exportcf $(LIB) yame-config

# What a downstream tool should compile and link against. Generated so the
# answer reflects this build -- in particular whether libcurl was found and
# which one -- instead of each consumer hardcoding its own guess.
yame-config: yame-config.in Makefile src/yame_version.h
	@sed -e 's|@PREFIX@|$(CURDIR)|g' \
	     -e 's|@VERSION@|$(shell sed -n 's/.*YAME_VERSION "\(.*\)".*/\1/p' src/yame_version.h)|g' \
	     -e 's|@CURL_CFLAGS@|$(CURL_CFLAGS)|g' \
	     -e 's|@CURL_LIBS@|$(CURL_LIBS)|g' \
	     -e 's|@CLIB@|$(CLIB)|g' \
	     yame-config.in > $@
	@chmod +x $@

exportcf:
	$(eval export CF_OPTIMIZE)

debug: CF_OPTIMIZE := 0
debug: CFLAGS += -g
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: build

#####################
##### libraries #####
#####################

$(LHTSLIB) :
	$(MAKE) -C $(LHTSLIB_DIR) libhts.a

###################
### compilation ###
###################

# Compile each .c in src/ into a .o in src/.
#
# -MMD -MP writes a .d beside each .o listing the headers that .c actually
# included; the -include below feeds those back to make. Without it a header
# change rebuilt nothing, so editing yame_version.h and running make produced
# a binary still reporting the previous version, and regenerating registry.h
# produced one still carrying the previous catalogue -- both silent, because
# the .c files were untouched and make had no reason to think otherwise.
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) -MMD -MP $< -o $@

-include $(wildcard $(SRC_DIR)/*.d)

###################
###  linking   ####
###################

$(PROG): $(LHTSLIB) $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LHTSLIB) $(CLIB)

# Archive every object except the CLI entry point (main.o has int main()).
LIBOBJECTS = $(filter-out $(SRC_DIR)/main.o, $(OBJECTS))

$(LIB): $(LHTSLIB) $(LIBOBJECTS)
	ar rcs $@ $(LIBOBJECTS)

###################
###   clean    ####
###################
clean :
	rm -f $(SRC_DIR)/*.o $(PROG) $(LIB)
	$(MAKE) -C $(LHTSLIB_DIR) clean
