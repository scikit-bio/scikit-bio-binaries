.PHONY: all api install test_bins test clean clean_install

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

BLASLIB=-llapacke -lcblas

ifeq ($(PLATFORM),Darwin)
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -dynamiclib -install_name @rpath/libskbb.so
else
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -shared
endif

ifeq ($(ARCH),x86_64)
   X86V3FLAGS := -march=x86-64-v3 -mtune=broadwell
   X86V4FLAGS := -march=x86-64-v4 -mtune=znver4
endif

CXXFLAGS += -fopenmp -Wall -std=c++17 -pedantic -I. $(OPT) -fPIC
NV_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC
AMD_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC

clean:
	# all source files are in subdirectories
	rm -f *.cpp *.hpp *.h
	rm -f *.o *.so *.a *.exe

##
# Utility/helper modules

SKBB_OBJS := util_rand.o

util_rand.o: util/rand.cpp util/rand.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

SKBB_OBJS += skbb_detect_acc.o

skbb_detect_acc.o: util/skbb_detect_acc.cpp util/skbb_detect_acc.hpp util/skbb_accapi.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

##
# skbb_accapi dynamic code wrappers generation
#

#TBD acc-specific binary variants

SKBB_OBJS += skbb_accapi_cpu.o

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
skbb_accapi_cpu.cpp: util/skbb_accapi.hpp
	./tools/generate_skbb_accapi.py cpu direct > $@

skbb_accapi_cpu.o: skbb_accapi_cpu.cpp util/skbb_accapi.hpp util/skbb_accapi_impl.hpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_cpu -c $< -o $@

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

skbb_accapi_dyn_acc_nv.o: skbb_accapi_dyn_acc_nv.cpp skbb_accapi_dyn_acc_nv.h skbb_task_impl.hpp
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

skbb_accapi_dyn_acc_amd.o: skbb_accapi_dyn_acc_amd.cpp skbb_accapi_dyn_acc_amd.h skbb_task_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

endif

##
# Distance methods

SKBB_OBJS += dist_permanova.o

dist_permanova.o: distance/permanova.cpp distance/permanova.hpp distance/permanova_dyn.hpp util/skbb_accapi.hpp util/skbb_detect_acc.hpp util/rand.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


#
# permanova dynamic code wrappers generation
#
SKBB_OBJS += permanova_cpu.o

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
permanova_cpu.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu direct > $@

permanova_cpu.o: permanova_cpu.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_cpu -c $< -o $@

ifeq ($(ARCH),x86_64)
# create optimized versions for all the higher x86_64 levels
SKBB_OBJS += permanova_cpu_x86_v3.o permanova_cpu_x86_v4.o

# cpu_x86_v3 is AVX2
permanova_cpu_x86_v3.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu_x86_v3 direct > $@

permanova_cpu_x86_v3.o: permanova_cpu_x86_v3.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(CXX) $(CXXFLAGS) $(X86V3FLAGS) -DSKBB_ACC_NM=skbb_cpu_x86_v3 -c $< -o $@

# cpu_x86_v4 is AVX512
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
permanova_dyn_acc_nv.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_nv.h
	./tools/generate_permanova_dyn.py acc_nv api > $@
permanova_acc_nv.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_nv.h
	./tools/generate_permanova_dyn.py acc_nv indirect > $@

permanova_acc_nv.o: permanova_acc_nv.cpp permanova_dyn_acc_nv.h distance/permanova_dyn.hpp util/skbb_dl.cpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

permanova_dyn_acc_nv.o: permanova_dyn_acc_nv.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(NV_CXX) $(NV_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_nv -c $< -o $@

endif

##

ifdef SKBB_ENABLE_ACC_AMD
SKBB_OBJS += permanova_acc_amd.o


permanova_dyn_acc_amd.h: distance/permanova_dyn_impl.hpp
	./tools/generate_permanova_dyn.py acc_amd api_h > $@
permanova_dyn_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd api > $@
permanova_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd indirect > $@

permanova_acc_amd.o: permanova_acc_amd.cpp permanova_dyn_acc_amd.h distance/permanova_dyn.hpp util/skbb_dl.cpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

permanova_dyn_acc_amd.o: permanova_dyn_acc_amd.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSKBB_ACC_NM=skbb_acc_amd -c $< -o $@

endif




##
# Ordination methods

SKBB_OBJS += ord_pcoa.o

ord_pcoa.o: ordination/principal_coordinate_analysis.cpp ordination/principal_coordinate_analysis.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

##

SHBB_EXTERN_OBJS := skbb_extern_distance.o skbb_extern_ordination.o skbb_extern_util.o
SHBB_EXTERN_HS := distance.h ordination.h util.h

skbb_extern_util.o: extern/skbb_util.cpp extern/util.h util/rand.hpp util/skbb_detect_acc.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

skbb_extern_distance.o: extern/skbb_distance.cpp extern/distance.h distance/permanova.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

skbb_extern_ordination.o: extern/skbb_ordination.cpp extern/ordination.h ordination/principal_coordinate_analysis.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

##

api: libskbb.so

libskbb.so: $(SHBB_EXTERN_OBJS) libskbb_cpu.a
	$(CXX) $(SO_LDDFLAGS) -o $@ $(SHBB_EXTERN_OBJS) libskbb_cpu.a $(BLASLIB)

libskbb_cpu.a: $(SKBB_OBJS)
	ar -rc $@ $(SKBB_OBJS)

##

install:
	mkdir -p "${PREFIX}/lib"
	mkdir -p "${PREFIX}/include/scikit-bio-binaries"
	rm -f "${PREFIX}/lib/libskbb.so"; cp "libskbb.so" "${PREFIX}/lib/"
	for f in $(SHBB_EXTERN_HS); do rm -f "${PREFIX}/include/scikit-bio-binaries/$${f}"; cp "extern/$${f}" "${PREFIX}/include/scikit-bio-binaries/"; done

clean_install:
	rm -f "${PREFIX}/lib/libskbb.so"
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

