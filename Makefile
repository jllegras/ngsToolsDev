
CC=g++
CFLAGS=-lm -lz -O3 -Wall -g

all: mkdir ngsCovar ngsSim ngs2dSFS ngsFST GetSubSfs GetSubSim GetSubGeno GetSwitchedGeno GetMergedGeno

mkdir:
	mkdir -p bin

ngsCovar: ngsCovar.cpp
	$(CC) $(CFLAGS) ngsCovar.cpp -o bin/ngsCovar

ngs2dSFS: ngs2dSFS.cpp
	$(CC) $(CFLAGS) ngs2dSFS.cpp -o bin/ngs2dSFS

ngsSim: ngsSim.cpp rbeta.cpp
	$(CC) ngsSim.cpp -o bin/ngsSim $(CFLAGS)

ngsFST: ngsFST.cpp
	$(CC) ngsFST.cpp -o bin/ngsFST $(CFLAGS)

GetSubSfs: GetSubSfs.cpp
	$(CC) GetSubSfs.cpp -o bin/GetSubSfs

GetSubSim: GetSubSim.cpp
	$(CC) GetSubSim.cpp -o bin/GetSubSim

GetSubGeno: GetSubGeno.cpp
	$(CC) GetSubGeno.cpp -o bin/GetSubGeno

GetSwitchedGeno: GetSwitchedGeno.cpp
	$(CC) GetSwitchedGeno.cpp -o bin/GetSwitchedGeno

GetMergedGeno: GetMergedGeno.cpp
	$(CC) GetMergedGeno.cpp -o bin/GetMergedGeno

clean:
	rm -rf *o bin/ngsFST bin/ngsCovar bin/ngsSim bin/ngs2dSFS bin/GetSubSfs bin/GetSubSim bin/GetSubGeno bin/GetSwitchedGeno bin/GetMergedGeno

