
CC=g++
CFLAGS=-lm -lz -O3 -Wall

all: mkdir ngsCovar ngsSim ngs2dSFS ngsFST GetSubSfs

mkdir:
	mkdir -p bin

ngsCovar: ngsCovar.cpp
	    $(CC) $(CFLAGS) ngsCovar.cpp -o bin/ngsCovar

ngs2dSFS: ngs2dSFS.cpp
	    $(CC) $(CFLAGS) ngs2dSFS.cpp -o bin/ngs2dSFS

ngsSim: ngsSim.cpp rbeta.cpp
	  $(CC) ngsSim.cpp -o bin/ngsSim $(CFLAGS)

ngsFST: ngsFST.cpp
	  $(CC) ngsFST.cpp -o bin/ngsFST -lgsl -lgslcblas $(CFLAGS)

GetSubSfs: GetSubSfs.cpp
	     $(CC) GetSubSfs.cpp -o bin/GetSubSfs.cpp

clean:
	rm -rf *o bin/ngsFST bin/ngsCovar bin/ngsSim bin/ngs2dSFS

