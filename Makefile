CC=gcc -m64
CFLAGS=-c -Wall -lm -ldl

all: fqhe

fqhe: main.o hilbert.o lanczos.o lebedev.o comb.o file_io.o angular.o
	$(CC) main.o hilbert.o lanczos.o lebedev.o comb.o file_io.o angular.o -o fqhe -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

hilbert.o: hilbert.c
	$(CC) $(CFLAGS) hilbert.c

lanczos.o: lanczos.c
	$(CC) $(CFLAGS) lanczos.c

lebedev.o: lebedev.c
	$(CC) $(CFLAGS) lebedev.c

comb.o: comb.c
	$(CC) $(CFLAGS) comb.c

file_io.o: file_io.c
	$(CC) $(CFLAGS) file_io.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c


clean:
	rm -rf *.o fqhe
