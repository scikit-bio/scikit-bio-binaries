.PHONY: all api install test_bins test clean clean_install wasm wasm_clean

# The native api/install paths and the wasm path both depend on generated
# source files under src/ (permanova_cpu.cpp, skbb_accapi_cpu.cpp, produced
# by python scripts). Each sub-make has its own DAG, so two concurrent
# top-level targets (e.g. `make -j all wasm`) could fire the generators
# in parallel and truncate the same output file. Force top-level targets
# to serialize; sub-makes retain their internal parallelism.
.NOTPARALLEL:

all:
	$(MAKE) api
	$(MAKE) install
	$(MAKE) test_bins

wasm:
	scripts/fetch_eigen.sh
	cd src && $(MAKE) wasm

wasm_test:
	scripts/fetch_eigen.sh
	cd src && $(MAKE) wasm_test

wasm_api_test: wasm
	cd api_tests/wasm && $(MAKE) wasm_test

wasm_clean:
	cd src && $(MAKE) wasm_clean
	cd api_tests/wasm && $(MAKE) wasm_clean

api:
	cd src && $(MAKE) api

install:
	cd src && $(MAKE) install

test_bins:
	cd src && $(MAKE) test_bins
	cd api_tests && $(MAKE) test_bins

test:
	cd src && $(MAKE) test
	cd api_tests && $(MAKE) test

clean:
	cd api_tests && $(MAKE) clean
	cd src && $(MAKE) clean

clean_install:
	cd src && $(MAKE) clean_install

