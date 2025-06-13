.PHONY: test_binaries test

SKBB_TESTS := test_distance.exe test_ordination.exe

test: $(SKBB_TESTS)
	./test_distance.exe
	./test_ordination.exe

test_binaries: (SKBB_TESTS)

test_distance.exe: test_distance.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

test_ordination.exe: test_ordination.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

clean:
	-rm -f *.o *.exe

