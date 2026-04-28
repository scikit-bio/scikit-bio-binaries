# WebAssembly build rules for scikit-bio-binaries.
#
# Activated by the top-level `wasm` target. Produces libskbb_wasm.a, a
# single-threaded static archive intended to be linked into downstream
# emscripten projects (e.g. unifrac-binaries/duckdb-miint).
#
# Backend: Eigen header-only (fetched via scripts/fetch_eigen.sh).
# No OpenMP, no pthread, no GPU, no dlopen, no CPU-arch dispatch.
#
# Expected toolchain (activated emsdk on PATH):
#   emcc, em++, emar
#
# Per-translation-unit rules for the CPU sources that are SHARED with
# the native build are defined in src/Makefile via the `skbb_cpu_tu`
# canned recipe (one definition emits both the .o and .wasm.o rules).
# This file owns:
#   - the WASM compiler vars (used by skbb_cpu_tu)
#   - the WASM-only TU (Eigen linalg backend)
#   - the WASM_OBJS list and libskbb_wasm.a archive
#   - the WASM test infrastructure
#
# Invocation: this Makefile fragment is included from src/Makefile and
# must be run with src/ as the working directory. All test source paths
# (e.g. `tests/wasm/test_*.cpp`) are relative to src/, and `-I.` from
# WASM_CXXFLAGS resolves `#include "tests/wasm/..."` against src/.

WASM_REPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../)
WASM_EIGEN_INC := $(WASM_REPO_ROOT)/.wasm-cache/eigen/include

# Compiler-independent target flags. The macros override the native Makefile's
# x86 dispatch (we explicitly don't want x86-v3/v4 object files under WASM
# regardless of what the host machine looks like) and select the Eigen
# backend in linalg_backend_eigen.cpp.
WASM_CXX      := em++
WASM_AR       := emar
WASM_CXXFLAGS := -std=c++17 -O3 -Wall -I. \
                 -I$(WASM_EIGEN_INC) \
                 -DSKBB_WASM=1 \
                 -DSKBB_BLAS_BACKEND_EIGEN=1 \
                 -DNOGPU=1 \
                 -fno-exceptions \
                 -Wno-unknown-pragmas

# Object list for the WASM archive. Stems must match the .wasm.o targets
# emitted by skbb_cpu_tu in src/Makefile, plus the Eigen-backed linalg
# implementation defined below.
WASM_OBJS := \
    util_rand.wasm.o \
    skbb_detect_acc.wasm.o \
    skbb_accapi_cpu.wasm.o \
    dist_permanova.wasm.o \
    permanova_cpu.wasm.o \
    ord_pcoa.wasm.o \
    ord_linalg_backend_eigen.wasm.o \
    skbb_extern_util.wasm.o \
    skbb_extern_distance.wasm.o \
    skbb_extern_ordination.wasm.o

# WASM-only TU: the Eigen-backed linalg backend has no native equivalent
# (native uses linalg_backend_lapacke.cpp), so this rule isn't generated
# by the symmetric skbb_cpu_tu macro.
ord_linalg_backend_eigen.wasm.o: ordination/linalg_backend_eigen.cpp ordination/linalg_backend.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

libskbb_wasm.a: $(WASM_OBJS)
	rm -f $@
	$(WASM_AR) rcs $@ $(WASM_OBJS)

# ---- WASM test binaries ----
# Emscripten links a .wasm + .js pair; node runs the .js, which loads
# the sibling .wasm. NODERAWFS=0 is the default — we don't touch disk.
WASM_TEST_LDFLAGS := -sEXIT_RUNTIME=1 -sALLOW_MEMORY_GROWTH=1 \
                     -sENVIRONMENT=node -sNODERAWFS=0

# The public headers refer to themselves as "scikit-bio-binaries/util.h".
# Provide that prefix via a staged include directory.
wasm_test_include_stage:
	mkdir -p .wasm-test-include/scikit-bio-binaries
	cp extern/util.h .wasm-test-include/scikit-bio-binaries/util.h
	cp extern/distance.h .wasm-test-include/scikit-bio-binaries/distance.h
	cp extern/ordination.h .wasm-test-include/scikit-bio-binaries/ordination.h

test_smoke_wasm.js: libskbb_wasm.a tests/wasm/test_smoke.cpp wasm_test_include_stage
	$(WASM_CXX) $(WASM_CXXFLAGS) -I.wasm-test-include \
	  tests/wasm/test_smoke.cpp libskbb_wasm.a \
	  $(WASM_TEST_LDFLAGS) -o $@

test_permanova_wasm.js: libskbb_wasm.a tests/wasm/test_permanova_wasm.cpp \
                       tests/wasm/permanova_inputs.hpp \
                       tests/wasm/expected/permanova_expected.h \
                       wasm_test_include_stage
	$(WASM_CXX) $(WASM_CXXFLAGS) -I.wasm-test-include \
	  tests/wasm/test_permanova_wasm.cpp libskbb_wasm.a \
	  $(WASM_TEST_LDFLAGS) -o $@

test_center_wasm.js: libskbb_wasm.a tests/wasm/test_center_wasm.cpp \
                    wasm_test_include_stage
	$(WASM_CXX) $(WASM_CXXFLAGS) -I.wasm-test-include \
	  tests/wasm/test_center_wasm.cpp libskbb_wasm.a \
	  $(WASM_TEST_LDFLAGS) -o $@

test_pcoa_wasm.js: libskbb_wasm.a tests/wasm/test_pcoa_wasm.cpp \
                  tests/wasm/pcoa_inputs.hpp \
                  tests/wasm/expected/pcoa_expected.h \
                  wasm_test_include_stage
	$(WASM_CXX) $(WASM_CXXFLAGS) -I.wasm-test-include \
	  tests/wasm/test_pcoa_wasm.cpp libskbb_wasm.a \
	  $(WASM_TEST_LDFLAGS) -o $@

wasm_test: test_smoke_wasm.js test_permanova_wasm.js test_center_wasm.js test_pcoa_wasm.js
	@echo "--- smoke ---"
	node test_smoke_wasm.js
	@echo "--- permanova ---"
	node test_permanova_wasm.js
	@echo "--- center ---"
	node test_center_wasm.js
	@echo "--- pcoa ---"
	node test_pcoa_wasm.js

# Native-build helpers that produce the expected-value headers consumed by
# the WASM tests above. These are not committed to the repo (per reviewer
# guidance); make wasm_test depends on them through the `expected/*.h`
# prerequisites and they are regenerated whenever missing or out of date.
# Both generators require the native build (libskbb_cpu.a) plus BLAS/LAPACKE
# for the PCoA generator.
tests/wasm/expected/permanova_expected.h: tests/wasm/generate_permanova_expected.cpp \
                                         tests/wasm/permanova_inputs.hpp \
                                         libskbb_cpu.a
	$(CXX) $(CXXFLAGS) \
	  tests/wasm/generate_permanova_expected.cpp libskbb_cpu.a \
	  $(LDFLAGS) -o generate_permanova_expected.exe
	mkdir -p tests/wasm/expected
	OMP_NUM_THREADS=1 ./generate_permanova_expected.exe > $@

tests/wasm/expected/pcoa_expected.h: tests/wasm/generate_pcoa_expected.cpp \
                                    tests/wasm/pcoa_inputs.hpp \
                                    libskbb_cpu.a
	$(CXX) $(CXXFLAGS) \
	  tests/wasm/generate_pcoa_expected.cpp libskbb_cpu.a \
	  $(LDFLAGS) $(BLASLIB) -o generate_pcoa_expected.exe
	mkdir -p tests/wasm/expected
	OMP_NUM_THREADS=1 ./generate_pcoa_expected.exe > $@

wasm_clean:
	rm -f libskbb_wasm.a *.wasm.o *_wasm.js *_wasm.wasm
	rm -f generate_permanova_expected.exe generate_pcoa_expected.exe
	rm -rf .wasm-test-include tests/wasm/expected

.PHONY: wasm_clean wasm_test wasm_test_include_stage
