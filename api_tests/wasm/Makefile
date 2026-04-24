# WASM variant of the public-API tests.
#
# Compiles the existing api_tests/test_*.c files under emcc and runs them
# under node, linking against the static libskbb_wasm.a produced by
# `make -C src wasm`. This proves the skbb_* C surface behaves identically
# to the native build, within documented tolerances.

.PHONY: all wasm_test wasm_clean

REPO_ROOT  := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
SRC_DIR    := $(REPO_ROOT)/src
WASM_LIB   := $(SRC_DIR)/libskbb_wasm.a
EIGEN_INC  := $(REPO_ROOT)/.wasm-cache/eigen/include

WASM_CC         := emcc
WASM_CFLAGS     := -std=c99 -O2 -Wall \
                   -I$(REPO_ROOT)/api_tests/wasm/include
WASM_LDFLAGS    := -sEXIT_RUNTIME=1 -sALLOW_MEMORY_GROWTH=1 \
                   -sENVIRONMENT=node -sNODERAWFS=0

# Stage the skbb public headers under scikit-bio-binaries/ so the .c
# files' `#include "scikit-bio-binaries/..."` lines resolve under emcc.
include/scikit-bio-binaries/util.h: $(SRC_DIR)/extern/util.h
	@mkdir -p $(@D)
	cp $< $@
include/scikit-bio-binaries/distance.h: $(SRC_DIR)/extern/distance.h
	@mkdir -p $(@D)
	cp $< $@
include/scikit-bio-binaries/ordination.h: $(SRC_DIR)/extern/ordination.h
	@mkdir -p $(@D)
	cp $< $@

STAGED_HEADERS := include/scikit-bio-binaries/util.h \
                  include/scikit-bio-binaries/distance.h \
                  include/scikit-bio-binaries/ordination.h

test_distance_wasm.js: ../test_distance.c $(STAGED_HEADERS) $(WASM_LIB)
	$(WASM_CC) $(WASM_CFLAGS) $< $(WASM_LIB) $(WASM_LDFLAGS) -o $@

test_ordination_wasm.js: ../test_ordination.c $(STAGED_HEADERS) $(WASM_LIB)
	$(WASM_CC) $(WASM_CFLAGS) $< $(WASM_LIB) $(WASM_LDFLAGS) -o $@

wasm_test: test_distance_wasm.js test_ordination_wasm.js
	@echo "--- api/distance ---"
	node test_distance_wasm.js
	@echo "--- api/ordination ---"
	node test_ordination_wasm.js

all: wasm_test

wasm_clean:
	rm -f *.js *.wasm
	rm -rf include
