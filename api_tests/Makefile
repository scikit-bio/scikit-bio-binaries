.PHONY: test_bins test clean

SKBB_TESTS := test_distance.exe test_ordination.exe

test: $(SKBB_TESTS)
	./test_distance.exe
	./test_ordination.exe

test_bins: $(SKBB_TESTS)

test_distance.exe: test_distance.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

test_ordination.exe: test_ordination.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

clean:
	-rm -f *.o *.exe

