.PHONY: test_binaries test

test: test_permanova.exe
	./test_permanova.exe

test_binaries: test_permanova.exe

test_permanova.exe: test_permanova.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

clean:
	-rm -f *.o *.exe

