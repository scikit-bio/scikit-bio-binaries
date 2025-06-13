.PHONY: all api install test_bins test clean clean_install

all:
	$(MAKE) api
	$(MAKE) install
	$(MAKE) test_bins

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

