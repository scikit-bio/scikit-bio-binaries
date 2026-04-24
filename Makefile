.PHONY: all api install test_bins test clean clean_install wasm wasm_clean

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

