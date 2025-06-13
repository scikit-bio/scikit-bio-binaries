.PHONY: all api api_acc install clean

all:
	$(MAKE) api
	$(MAKE) install

PLATFORM := $(shell uname -s)

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
endif

BLASLIB=-llapacke -lcblas

ifeq ($(PLATFORM),Darwin)
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -dynamiclib -install_name @rpath/lib$(SSU).so
        NOGPU := 1
else
        SO_LDDFLAGS = -L$(PREFIX)/lib -fopenmp -shared
        BLASLIB += -lgfortran -lquadmath
        LDFLAGS += -ldl
endif

CXXFLAGS += -fopenmp -Wall -std=c++17 -pedantic -I. $(OPT) -fPIC
NV_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC
AMD_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC

clean:
	# all source files are in subdirectories
	rm -f *.cpp *.hpp *.h
	rm -f *.o *.so *.a *.exe

api: libskbb.so

install:
	mkdir -p "${PREFIX}/lib"
	mkdir -p "${PREFIX}/include/scikit-bio-binaries"
	rm -f "${PREFIX}/lib/libskbb.so"; cp "libskbb.so" "${PREFIX}/lib/"
	rm -f "${PREFIX}/include/scikit-bio-binaries/distance.h"; cp extern/distance.h "${PREFIX}/include/scikit-bio-binaries/"

##
# Utility/helper modules

SKBB_OBJS := util_rand.o

util_rand.o: util/rand.cpp util/rand.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifndef NOGPU

SKBB_OBJS += skbb_detect_acc.o

skbb_detect_acc.o: util/skbb_detect_acc.cpp util/skbb_detect_acc.hpp util/skbb_accapi.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

endif

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
SKBB_OBJS += permanova_dyn_cpu.o

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
permanova_dyn_cpu.cpp: distance/permanova_dyn.hpp
	./tools/generate_permanova_dyn.py cpu direct > $@

permanova_dyn_cpu.o: permanova_dyn_cpu.cpp distance/permanova_dyn.hpp distance/permanova_dyn_impl.hpp
	$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_cpu -c $< -o $@

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

skbb_extern.o: extern/skbb.cpp extern/distance.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

##

libskbb.so: skbb_extern.o $(SKBB_OBJS)
	$(CXX) $(SO_LDDFLAGS) -o $@ skbb_extern.o $(SKBB_OBJS) $(BLASLIB)

libskbb_cpu.a: $(SKBB_OBJS)
	ar -rc $@ $(SKBB_OBJS)

##

test_bins: test_pcoa.exe test_permanova.exe


test: test_bins
	./test_pcoa.exe
	./test_permanova.exe

test_pcoa.exe: tests/test_pcoa.cpp libskbb_cpu.a
	$(CXX) $(CXXFLAGS) $< -o $@ libskbb_cpu.a $(LDFLAGS) $(BLASLIB)

test_permanova.exe: tests/test_permanova.cpp libskbb_cpu.a
	$(CXX) $(CXXFLAGS) $< -o $@ libskbb_cpu.a $(LDFLAGS)

