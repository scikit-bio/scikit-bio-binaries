

#
# skbio_alt dynamic code wrappers generation
#

# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
skbb_accapi_cpu.cpp: skbb_accapi.hpp
	./tools/generate_skbb_accapi.py cpu direct > $@
skbio_alt_dyn_cpu.cpp: skbio_alt_dyn.hpp
	./tools/generate_skbio_alt_dyn.py cpu direct > $@

##

skbb_accapi_dyn_acc_nv.h: skbb_accapi_impl.hpp
	./tools/generate_skbb_accapi.py acc_nv api_h > $@
skbb_accapi_dyn_acc_nv.cpp: skbb_accapi_impl.hpp skbb_accapi_dyn_acc_nv.h
	./tools/generate_skbb_accapi.py acc_nv api > $@
skbb_accapi_acc_nv.cpp: skbb_accapi_impl.hpp skbb_accapi_dyn_acc_nv.h
	./tools/generate_skbb_accapi.py acc_nv indirect > $@

skbio_alt_dyn_acc_nv.h: skbio_alt_dyn_impl.hpp
	./tools/generate_skbio_alt_dyn.py acc_nv api_h > $@
skbio_alt_dyn_acc_nv.cpp: skbio_alt_dyn_impl.hpp skbio_alt_dyn_acc_nv.h
	./tools/generate_skbio_alt_dyn.py acc_nv api > $@
skbio_alt_acc_nv.cpp: skbio_alt_dyn_impl.hpp skbio_alt_dyn_acc_nv.h
	./tools/generate_skbio_alt_dyn.py acc_nv indirect > $@

##

skbb_accapi_dyn_acc_amd.h: skbb_accapi_impl.hpp
	./tools/generate_skbb_accapi.py acc_amd api_h > $@
skbb_accapi_dyn_acc_amd.cpp: skbb_accapi_impl.hpp skbb_accapi_dyn_acc_amd.h
	./tools/generate_skbb_accapi.py acc_amd api > $@
skbb_accapi_acc_amd.cpp: skbb_accapi_impl.hpp skbb_accapi_dyn_acc_amd.h
	./tools/generate_skbb_accapi.py acc_amd indirect > $@

skbio_alt_dyn_acc_amd.h: skbio_alt_dyn_impl.hpp
	./tools/generate_skbio_alt_dyn.py acc_amd api_h > $@
skbio_alt_dyn_acc_amd.cpp: skbio_alt_dyn_impl.hpp skbio_alt_dyn_acc_amd.h
	./tools/generate_skbio_alt_dyn.py acc_amd api > $@
skbio_alt_acc_amd.cpp: skbio_alt_dyn_impl.hpp skbio_alt_dyn_acc_amd.h
	./tools/generate_skbio_alt_dyn.py acc_amd indirect > $@

##

skbb_accapi_cpu.o: skbb_accapi_cpu.cpp skbb_accapi.hpp skbb_accapi_impl.hpp
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

skbio_alt_dyn_cpu.o: skbio_alt_dyn_cpu.cpp skbio_alt_dyn.hpp skbio_alt_dyn.hpp
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

##

skbb_accapi_acc_nv.o: skbb_accapi_acc_nv.cpp skbb_accapi.hpp skbb_accapi_dyn_acc_nv.h
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@
skbb_accapi_dyn_acc_nv.o: skbb_accapi_dyn_acc_nv.cpp skbb_accapi_dyn_acc_nv.h skbb_task_impl.hpp
	$(NV_CXX) $(NV_CXXFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@

skbio_alt_acc_nv.o: skbio_alt_acc_nv.cpp skbio_alt_dyn.hpp skbio_alt_dyn_acc_nv.h
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@
skbio_alt_dyn_acc_nv.o: skbio_alt_dyn_acc_nv.cpp skbio_alt_dyn_acc_nv.h skbio_alt_dyn_impl.hpp
	$(NV_CXX) $(NV_CXXFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@

##

skbb_accapi_acc_amd.o: skbb_accapi_acc_amd.cpp skbb_accapi.hpp skbb_accapi_dyn_acc_amd.h
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@
skbb_accapi_dyn_acc_amd.o: skbb_accapi_dyn_acc_amd.cpp skbb_accapi_dyn_acc_amd.h skbb_task_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

skbio_alt_acc_amd.o: skbio_alt_acc_amd.cpp skbio_alt_dyn.hpp skbio_alt_dyn_acc_amd.h
	$(CXX) $(CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@
skbio_alt_dyn_acc_amd.o: skbio_alt_dyn_acc_amd.cpp skbio_alt_dyn_acc_amd.h skbio_alt_dyn_impl.hpp
	$(AMD_CXX) $(AMD_CXXFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

##

