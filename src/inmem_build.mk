# Native in-memory static-archive build for scikit-bio-binaries.
#
# Activated by the top-level `inmem_static` target. Produces
# libskbb_inmem.a — the same in-memory subset shipped as libskbb_wasm.a
# (Eigen backend, no LAPACKE/cblas, no GPU, no CPU-arch dispatch),
# but compiled with the host toolchain so OpenMP is enabled. Intended
# for downstream native projects that want to embed the in-memory
# scikit-bio-binaries surface without dragging in BLAS/LAPACK.
#
# Per-translation-unit rules for sources shared with the native and
# WASM builds are defined in src/Makefile via the `skbb_cpu_tu` canned
# recipe (one definition emits .o, .wasm.o, AND .inmem.o rules). This
# file owns:
#   - the inmem compiler vars (used by skbb_cpu_tu)
#   - the inmem-only TU (Eigen linalg backend, native CXX)
#   - the libskbb_inmem.a archive rule
#   - install / clean targets
#
# Invocation: this Makefile fragment is included from src/Makefile and
# must be run with src/ as the working directory.

# Reuse the Eigen drop fetched by scripts/fetch_eigen.sh for the WASM
# build — the headers don't care which compiler is consuming them.
INMEM_REPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../)
INMEM_EIGEN_INC := $(INMEM_REPO_ROOT)/.wasm-cache/eigen/include

INMEM_CXX      ?= $(CXX)
INMEM_AR       ?= ar
INMEM_MPFLAG   ?= -fopenmp
INMEM_CXXFLAGS := -std=c++17 -O3 -Wall -fPIC -I. \
                  -I$(INMEM_EIGEN_INC) \
                  -DSKBB_BLAS_BACKEND_EIGEN=1 \
                  -DEIGEN_DONT_PARALLELIZE \
                  -DNOGPU=1 \
                  $(INMEM_MPFLAG) \
                  -Wno-unknown-pragmas

# Object list for the inmem archive. Stems must match the .inmem.o
# targets emitted by skbb_cpu_tu in src/Makefile, plus the Eigen-backed
# linalg implementation defined below. Mirrors WASM_OBJS — same TU set.
INMEM_OBJS := \
    util_rand.inmem.o \
    skbb_detect_acc.inmem.o \
    skbb_accapi_cpu.inmem.o \
    dist_permanova.inmem.o \
    permanova_cpu.inmem.o \
    ord_pcoa.inmem.o \
    ord_linalg_backend_eigen.inmem.o \
    skbb_extern_util.inmem.o \
    skbb_extern_distance.inmem.o \
    skbb_extern_ordination.inmem.o

# Inmem-only TU: the Eigen-backed linalg backend has no symmetric-macro
# equivalent because the native build uses linalg_backend_lapacke.cpp
# under the same `ord_linalg_backend.o` stem. Same source as the WASM
# rule, just compiled with $(INMEM_CXX) instead of em++.
ord_linalg_backend_eigen.inmem.o: ordination/linalg_backend_eigen.cpp ordination/linalg_backend.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

libskbb_inmem.a: $(INMEM_OBJS)
	rm -f $@
	$(INMEM_AR) rcs $@ $(INMEM_OBJS)

inmem_static: libskbb_inmem.a

# Install (archive + public headers under a stable prefix layout).
# Guards against the empty-PREFIX footgun: the top-level Makefile falls
# back to CONDA_PREFIX, but if both are unset `mkdir -p /lib` would
# silently target the root filesystem.
install_inmem: libskbb_inmem.a
	@test -n "$(PREFIX)" || { echo "ERROR: PREFIX is unset (and CONDA_PREFIX is unset). Pass PREFIX=/path or activate a conda env."; exit 1; }
	mkdir -p "${PREFIX}/lib" "${PREFIX}/include/scikit-bio-binaries"
	rm -f "${PREFIX}/lib/libskbb_inmem.a"; cp libskbb_inmem.a "${PREFIX}/lib/"
	for f in $(SHBB_EXTERN_HS); do rm -f "${PREFIX}/include/scikit-bio-binaries/$${f}"; cp "extern/$${f}" "${PREFIX}/include/scikit-bio-binaries/"; done

inmem_clean:
	rm -f libskbb_inmem.a *.inmem.o

.PHONY: inmem_static install_inmem inmem_clean
