# Makefile
# Author: 	Hannes Herrmann
# File:		Simulation


# Compilers and flags:
CC = g++
CFLAGS = -O3 -std=c++14 -m64 -march=native -mtune=native -pthread
LDFLAGS = 
INCLUDES = brain.o

main: main.cpp $(INCLUDES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


#Compile all Includes

%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $<

#Data generation:

.PHONY: clean

clean:
	-rm main $(INCLUDES) 
