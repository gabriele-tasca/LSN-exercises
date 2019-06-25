CC = g++
CFLAGS = -Wall -O3 --std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable

mol : MolDyn_NVE.o random.o
	$(CC) random.o MolDyn_NVE.o -o mol
MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o mol seed.out
