.PHONY: all api install test_bins test clean clean_install wasm wasm_test wasm_clean

all:
	$(MAKE) api
	$(MAKE) test_bins
	$(MAKE) install

PLATFORM := $(shell uname -s)
ARCH := $(shell uname -m)

ifeq ($(PLATFORM),Darwin)
  # no GPUs on MacOS
  NOGPU := 1
endif

ifeq ($(PREFIX),)
	PREFIX := $(CONDA_PREFIX)
endif

ifndef NOGPU
   # NVIDIA GPUs
   CXXFLAGS += -DSKBB_ENABLE_ACC_NV=1
   SKBB_ENABLE_ACC_NV := 1
   # AMD GPUs
   CXXFLAGS += -DSKBB_ENABLE_ACC_AMD=1
   SKBB_ENABLE_ACC_AMD := 1
   LDFLAGS += -ldl
endif

# Conditional assignment so callers can override via env or command line
# (e.g. BLASLIB=-lopenblas on Ubuntu where libopenblas-dev does not ship
# the libcblas.so / liblapacke.so individual link names). Without `?=`,
# this assignment would shadow the env variable inside the Makefile.
BLASLIB ?= -llapacke -lcblas

ifeq ($(PLATFORM),Darwin)
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -dynamiclib -install_name @rpath/libskbb.so
else
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -shared
endif

ifeq ($(ARCH),x86_64)
  ifeq ($(PLATFORM),Darwin)
   # No new x86 on Mac, so no need to make our life complicated
  else
   # create optimized versions for all the higher x86_64 levels
   CXXFLAGS += -DSKBB_ENABLE_CPU_X86_LEVELS=1
   CXXFLAGS += -DSKBB_ENABLE_CPU_X86V3=1 -DSKBB_ENABLE_CPU_X86V4=1
   SKBB_ENABLE_CPU_X86V3 := 1
   SKBB_ENABLE_CPU_X86V4 := 1
   # cpu_x86_v3 is AVX2
   X86V3FLAGS := -march=x86-64-v3 -mtune=broadwell
   # cpu_x86_v4 is AVX512
   X86V4FLAGS := -march=x86-64-v4 -mtune=znver4
  endif
endif

CXXFLAGS += -O3 -ffast-math -fopenmp -Wall -std=c++17 -pedantic -I. $(OPT) -fPIC

ifneq ($(SKBB_ENABLE_ACC_NV),)
  ifneq ($(NV_CUDA),)
    # CUDA version
    SKBB_ENABLE_ACC_NV_BINS := 1
    NV_CXX := nvcc
    NV_CXXFLAGS := -DSKBB_CUDA
    NV_CXXFLAGS += -O3 --use_fast_math -std=c++17 -I. $(OPT) -Xcompiler -fPIC
    NV_CXXFLAGS += -arch=all-major
  else
    ifneq ($(NV_CXX),)
     # OpenACC version
     SKBB_ENABLE_ACC_NV_BINS := 1
     NV_CXXFLAGS += -acc
     NV_CXXFLAGS += -Ofast -std=c++17 -I. $(OPT) -fPIC
     NV_LDFLAGS += -shared -acc -Bstatic_pgi
     NV_CXXFLAGS += -gpu=ccall
     NV_LDFLAGS += -gpu=ccall
  endif

  endif
endif

ifneq ($(SKBB_ENABLE_ACC_AMD),)
  ifneq ($(AMD_HIP),)
    # HIP version
    SKBB_ENABLE_ACC_AMD_BINS := 1
    AMD_CXX := hipcc
    AMD_CXXFLAGS := -DSKBB_HIP
    AMD_CXXFLAGS += -O3 -ffast-math -std=c++17 -I. $(OPT) -fPIC
    AMD_CXXFLAGS += --offload-arch=gfx1100,gfx1101,gfx1102,gfx1103,gfx1030,gfx1031,gfx90a,gfx942
    AMD_SO_LDDFLAGS += -shared -fgpu-rdc --hip-link
  else
   ifneq ($(AMD_CXX),)
    SKBB_ENABLE_ACC_AMD_BINS := 1
    AMD_CXXFLAGS += -fopenmp -fopenmp-offload-mandatory -DOMPGPU=1
    AMD_CXXFLAGS += -O3 -ffast-math -std=c++17 -I. $(OPT) -fPIC
    AMD_LDFLAGS += -shared -fopenmp
    AMD_CXXFLAGS += --offload-arch=gfx1100,gfx1101,gfx1102,gfx1103,gfx1030,gfx1031,gfx90a,gfx942
    AMD_LDFLAGS += --offload-arch=gfx1100,gfx1101,gfx1102,gfx1103,gfx1030,gfx1031,gfx90a,gfx942
   endif
  endif
endif

# WebAssembly build (emsdk required; see scripts/fetch_eigen.sh).
# Provides WASM_CXX / WASM_CXXFLAGS used by the skbb_cpu_tu macro below,
# plus the libskbb_wasm.a archive rule and WASM-only objects/tests.
include wasm/emscripten_build.mk

wasm: libskbb_wasm.a

clean:
	# all source files are in subdirectories
	rm -f *.cpp *.hpp *.h *.cu
	rm -f *.o *.so *.a *.exe

##
# Shared per-translation-unit recipe.
#
# A handful of CPU-only sources are compiled into BOTH the native libskbb
# (via $(CXX) $(CXXFLAGS)) and the WASM libskbb_wasm (via $(WASM_CXX)
# $(WASM_CXXFLAGS)). The two builds previously had parallel per-file rules
# in src/Makefile and src/wasm/emscripten_build.mk; this macro emits both
# from a single point of definition so adding a new TU only needs one line.
#
# Native-only rules (x86 dispatch, GPU dispatch) and WASM-only rules
# (Eigen linalg backend) remain explicit further down — the macro is for
# the symmetric case only.
#
# Args:
#   $(1) — object stem (e.g. `util_rand` => `util_rand.o`, `util_rand.wasm.o`)
#   $(2) — source path
#   $(3) — header dependencies (space separated)
#   $(4) — extra compile flags applied to BOTH builds (e.g. -DSKBB_ACC_NM=skbb_cpu)

define skbb_cpu_tu
$(1).o: $(2) $(3)
	$$(CXX) $$(CXXFLAGS) $(4) -c $$< -o $$@
$(1).wasm.o: $(2) $(3)
	$$(WASM_CXX) $$(WASM_CXXFLAGS) $(4) -c $$< -o $$@
endef

##
# Utility/helper modules

SKBB_OBJS := util_rand.o
$(eval $(call skbb_cpu_tu,util_rand,util/rand.cpp,util/rand.hpp,))

SKBB_OBJS += skbb_detect_acc.o
$(eval $(call skbb_cpu_tu,skbb_detect_acc,util/skbb_detect_acc.cpp,util/skbb_detect_acc.hpp util/skbb_accapi.hpp,))

##
# skbb_accapi dynamic code wrappers generation
#

#TBD acc-specific binary variants

SKBB_OBJS += skbb_accapi_cpu.o

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
skbb_accapi_cpu.cpp: util/skbb_accapi.hpp
	./tools/generate_skbb_accapi.py cpu direct > $@

$(eval $(call skbb_cpu_tu,skbb_accapi_cpu,skbb_accapi_cpu.cpp,util/skbb_accapi.hpp util/skbb_accapi_impl.hpp,-DSKBB_ACC_NM=skbb_cpu))

#
# True accelerated variants will use a separate compiler
# thus separate generic and acc-cpecific files
#

ifdef SKBB_ENABLE_ACC_NV
SKBB_OBJS += skbb_accapi_acc_nv.o

skbb_accapi_dyn_acc_nv.h: util/skbb_accapi_impl.hpp
	./tools/generate_skbb_accapi.py acc_nv api_h > $@
skbb_accapi_dyn_acc_nv.cpp: util/skbb_accapi_impl.hpp skbb_accapi_dyn_acc_nv.h
	./tools/generate_skbb_accapi.py acc_nv api > $@
skbb_accapi_acc_nv.cpp: util/skbb_accapi_impl.hpp skbb_accapi_dyn_acc_nv.h
	./tools/generate_skbb_accapi.py acc_nv indirect > $@

skbb_accapi_acc_nv.o: skbb_accapi_acc_nv.cpp util/skbb_accapi.hpp skbb_accapi_dyn_acc_nv.h
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

skbb_accapi_dyn_acc_nv.o: skbb_accapi_dyn_acc_nv.cpp skbb_accapi_dyn_acc_nv.h skbb_accapi_dyn_acc_nv.h
	$(NV_CXX) $(NV_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

endif

ifdef SKBB_ENABLE_ACC_AMD
SKBB_OBJS += skbb_accapi_acc_amd.o

skbb_accapi_dyn_acc_amd.h: util/skbb_accapi_impl.hpp
	./tools/generate_skbb_accapi.py acc_amd api_h > $@
skbb_accapi_dyn_acc_amd.cpp: util/skbb_accapi_impl.hpp skbb_accapi_dyn_acc_amd.h
	./tools/generate_skbb_accapi.py acc_amd api > $@
skbb_accapi_acc_amd.cpp: util/skbb_accapi_impl.hpp skbb_accapi_dyn_acc_amd.h
	./tools/generate_skbb_accapi.py acc_amd indirect > $@

skbb_accapi_acc_amd.o: skbb_accapi_acc_amd.cpp util/skbb_accapi.hpp skbb_accapi_dyn_acc_amd.h
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

skbb_accapi_dyn_acc_amd.o: skbb_accapi_dyn_acc_amd.cpp skbb_accapi_dyn_acc_amd.h skbb_accapi_dyn_acc_amd.h
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

endif

##
# Distance methods

SKBB_OBJS += dist_permanova.o
$(eval $(call skbb_cpu_tu,dist_permanova,distance/permanova.cpp,distance/permanova.hpp distance/permanova_dyn.hpp util/skbb_accapi.hpp util/skbb_detect_acc.hpp util/rand.hpp util/portable_shuffle.hpp util/skbb_dgb_info.hpp,))

#
# permanova dynamic code wrappers generation
#
SKBB_OBJS += permanova_cpu.o

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
permanova_cpu.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu direct > $@

$(eval $(call skbb_cpu_tu,permanova_cpu,permanova_cpu.cpp,distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp,-DSKBB_ACC_NM=skbb_cpu))

ifeq ($(SKBB_ENABLE_CPU_X86V3),1)
SKBB_OBJS += permanova_cpu_x86_v3.o

permanova_cpu_x86_v3.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu_x86_v3 direct > $@

permanova_cpu_x86_v3.o: permanova_cpu_x86_v3.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(CXX) $(CXXFLAGS) $(X86V3FLAGS) -DSKBB_ACC_NM=skbb_cpu_x86_v3 -c $< -o $@

endif

ifeq ($(SKBB_ENABLE_CPU_X86V4),1)
SKBB_OBJS += permanova_cpu_x86_v4.o

permanova_cpu_x86_v4.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu_x86_v4 direct > $@

permanova_cpu_x86_v4.o: permanova_cpu_x86_v4.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(CXX) $(CXXFLAGS) $(X86V4FLAGS) -DSKBB_ACC_NM=skbb_cpu_x86_v4 -c $< -o $@

endif

#
# True accelerated variants will use a separate compiler
# thus separate generic and acc-cpecific files
#

ifdef SKBB_ENABLE_ACC_NV
SKBB_OBJS += permanova_acc_nv.o


permanova_dyn_acc_nv.h: distance/permanova_dyn_impl.hpp
	./tools/generate_permanova_dyn.py acc_nv api_h > $@
permanova_acc_nv.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_nv.h
	./tools/generate_permanova_dyn.py acc_nv indirect > $@

permanova_acc_nv.o: permanova_acc_nv.cpp permanova_dyn_acc_nv.h distance/permanova_dyn.hpp util/skbb_dl.cpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

ifneq ($(NV_CUDA),)
# CUDA version
permanova_dyn_acc_nv.cu: distance/permanova_dyn_impl.hpp permanova_dyn_acc_nv.h
	./tools/generate_permanova_dyn.py acc_nv api > $@
permanova_dyn_acc_nv.o: permanova_dyn_acc_nv.cu distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(NV_CXX) $(NV_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

else
# OpenACC version
permanova_dyn_acc_nv.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_nv.h
	./tools/generate_permanova_dyn.py acc_nv api > $@
permanova_dyn_acc_nv.o: permanova_dyn_acc_nv.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(NV_CXX) $(NV_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

endif

endif

##

ifdef SKBB_ENABLE_ACC_AMD
SKBB_OBJS += permanova_acc_amd.o

permanova_dyn_acc_amd.h: distance/permanova_dyn_impl.hpp
	./tools/generate_permanova_dyn.py acc_amd api_h > $@
permanova_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd indirect > $@

permanova_acc_amd.o: permanova_acc_amd.cpp permanova_dyn_acc_amd.h distance/permanova_dyn.hpp util/skbb_dl.cpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

ifneq ($(AMD_HIP),)
# HIP version
permanova_dyn_acc_amd.hip: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd api > $@
permanova_dyn_acc_amd.o: permanova_dyn_acc_amd.hip distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

else
# OpenMP Target version
permanova_dyn_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd api > $@
permanova_dyn_acc_amd.o: permanova_dyn_acc_amd.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

endif

endif




##
# Ordination methods

SKBB_OBJS += ord_pcoa.o
$(eval $(call skbb_cpu_tu,ord_pcoa,ordination/principal_coordinate_analysis.cpp,ordination/principal_coordinate_analysis.hpp ordination/linalg_backend.hpp util/skbb_dgb_info.hpp,))

# Linalg backend split: native uses cblas+LAPACKE, WASM uses Eigen.
# These are different source files, so the symmetric skbb_cpu_tu macro
# doesn't apply — keep two explicit rules. The WASM rule lives in
# wasm/emscripten_build.mk to keep the Eigen include path local to that
# file.

SKBB_OBJS += ord_linalg_backend.o

ord_linalg_backend.o: ordination/linalg_backend_lapacke.cpp ordination/linalg_backend.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

##

SHBB_EXTERN_OBJS := skbb_extern_distance.o skbb_extern_ordination.o skbb_extern_util.o
SHBB_EXTERN_HS := distance.h ordination.h util.h

$(eval $(call skbb_cpu_tu,skbb_extern_util,extern/skbb_util.cpp,extern/util.h util/rand.hpp util/skbb_detect_acc.hpp,))
$(eval $(call skbb_cpu_tu,skbb_extern_distance,extern/skbb_distance.cpp,extern/distance.h distance/permanova.hpp,))
$(eval $(call skbb_cpu_tu,skbb_extern_ordination,extern/skbb_ordination.cpp,extern/ordination.h ordination/principal_coordinate_analysis.hpp,))

##

SKBB_SHLIBS := libskbb.so

libskbb.so: $(SHBB_EXTERN_OBJS) libskbb_cpu.a
	$(CXX) $(SO_LDDFLAGS) -o $@ $(SHBB_EXTERN_OBJS) libskbb_cpu.a $(BLASLIB)

libskbb_cpu.a: $(SKBB_OBJS)
	ar -rc $@ $(SKBB_OBJS)

ifdef SKBB_ENABLE_ACC_NV_BINS

SKBB_SHLIBS += libskbb_acc_nv.so

ifneq ($(NV_CUDA),)
#CUDA version
libskbb_acc_nv.so: skbb_accapi_dyn_acc_nv.o permanova_dyn_acc_nv.o
	$(CXX) $(SO_LDDFLAGS) -o $@ skbb_accapi_dyn_acc_nv.o permanova_dyn_acc_nv.o -lcudart_static
else
#OpenACC version
libskbb_acc_nv.so: skbb_accapi_dyn_acc_nv.o permanova_dyn_acc_nv.o
	$(NV_CXX) $(NV_LDFLAGS) $(SO_LDDFLAGS) -o $@ skbb_accapi_dyn_acc_nv.o permanova_dyn_acc_nv.o
endif

endif

ifdef SKBB_ENABLE_ACC_AMD_BINS

SKBB_SHLIBS += libskbb_acc_amd.so

ifneq ($(AMD_HIP),)
#HIP version
libskbb_acc_amd.so: skbb_accapi_dyn_acc_amd.o permanova_dyn_acc_amd.o
	$(AMD_CXX) $(AMD_SO_LDDFLAGS) -o $@ skbb_accapi_dyn_acc_amd.o permanova_dyn_acc_amd.o
else
#OpenMP Target version
libskbb_acc_amd.so: skbb_accapi_dyn_acc_amd.o permanova_dyn_acc_amd.o
	$(AMD_CXX) $(AMD_LDFLAGS) $(SO_LDDFLAGS) -o $@ skbb_accapi_dyn_acc_amd.o permanova_dyn_acc_amd.o
endif

endif

api: $(SKBB_SHLIBS)

##

install:
	mkdir -p "${PREFIX}/lib"
	mkdir -p "${PREFIX}/include/scikit-bio-binaries"
	for f in $(SKBB_SHLIBS); do rm -f "${PREFIX}/lib/$${f}"; cp "$${f}" "${PREFIX}/lib/"; done
	for f in $(SHBB_EXTERN_HS); do rm -f "${PREFIX}/include/scikit-bio-binaries/$${f}"; cp "extern/$${f}" "${PREFIX}/include/scikit-bio-binaries/"; done

clean_install:
	rm -f "${PREFIX}/lib/libskbb.so"
	for f in $(SKBB_SHLIBS); do rm -f "${PREFIX}/lib/$${f}"; done
	for f in $(SHBB_EXTERN_HS); do rm -f "${PREFIX}/include/scikit-bio-binaries/$${f}"; done

##

test_bins: test_pcoa.exe test_permanova.exe


test: test_bins
	./test_pcoa.exe
	./test_permanova.exe

test_pcoa.exe: tests/test_pcoa.cpp libskbb_cpu.a
	$(CXX) $(CXXFLAGS) $< -o $@ libskbb_cpu.a $(LDFLAGS) $(BLASLIB)

test_permanova.exe: tests/test_permanova.cpp libskbb_cpu.a
	$(CXX) $(CXXFLAGS) $< -o $@ libskbb_cpu.a $(LDFLAGS)

