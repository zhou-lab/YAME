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

.PHONY: all build debug clean lib

all: build

build: exportcf $(PROG)

# Build just the static library (used by MethScope2).
lib: exportcf $(LIB)

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

# Compile each .c in src/ into a .o in src/
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

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
