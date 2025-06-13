.PHONY: test_binaries test

test: test_distance.exe
	./test_distance.exe

test_binaries: test_distance.exe

test_distance.exe: test_distance.c
	$(CC) $(CFLAGS) -std=c99 -O0 -g $< $(LDFLAGS) -lskbb -o $@

clean:
	-rm -f *.o *.exe

