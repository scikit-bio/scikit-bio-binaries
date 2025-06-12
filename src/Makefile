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
endif

CXXFLAGS += -fopenmp -Wall -std=c++17 -pedantic -I. $(OPT) -fPIC
NV_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC
AMD_CXXFLAGS += -std=c++17 -I. $(OPT) -fPIC

clean:
	# all source files are in subdirectories
	rm -f *.cpp *.hpp *.h
	rm -f *.o *.so

api: libskbb.so

install:
	#TBD

##
# Utility/helper modules

SKBB_OBJS := util_rand.o skbb_detect_acc.o

util_rand.o: util/rand.cpp util/rand.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

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
SKBB_OBJS += permanova_dyn_acc_amd.o


permanova_dyn_acc_amd.h: distance/permanova_dyn_impl.hpp
	./tools/generate_permanova_dyn.py acc_amd api_h > $@
permanova_dyn_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd api > $@
permanova_acc_amd.cpp: distance/permanova_dyn_impl.hpp permanova_dyn_acc_amd.h
	./tools/generate_permanova_dyn.py acc_amd indirect > $@

endif


##


##

skbio_alt_acc_amd.o: skbio_alt_acc_amd.cpp skbio_alt_dyn.hpp skbio_alt_dyn_acc_amd.h
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@
skbio_alt_dyn_acc_amd.o: skbio_alt_dyn_acc_amd.cpp skbio_alt_dyn_acc_amd.h skbio_alt_dyn_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

##
# Ordination methods

SKBB_OBJS += ord_pcoa.o

ord_pcoa.o: ordination/principal_coordinate_analysis.cpp ordination/principal_coordinate_analysis.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

##

#SKBB_OBJS += skbio_alt.o skbb_accapi_cpu.o skbio_alt_dyn_cpu.o

libskbb.so: $(SKBB_OBJS)
	$(CXX) $(SO_LDDFLAGS) -o $@ $(SKBB_OBJS) $(BLASLIB)


