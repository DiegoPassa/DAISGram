# *****************************************************
 
CC = g++
CFLAGS = -Wall -pedantic --std=c++11
OPTORDEBUG = -O1
INCL = ./includes
CI = -I$(INCL)

# *****************************************************

all: main

debug: OPTORDEBUG = -g -O0
debug: main
 
main: main.cpp tensor.o tensor_h.o
	$(CC) $(CFLAGS) $(OPTORDEBUG) $^ -o $@

tensor.o: tensor.cpp
	$(CC) $(CFLAGS) $(OPTORDEBUG) -c $< -o $@

tensor_h.o: tensor.h
	$(CC) $(CFLAGS) $(OPTORDEBUG) -c $< -o $@

clean:
	rm main *.o
