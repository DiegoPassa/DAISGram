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
 
main: test.cpp tensor.o DAISGram.o libbmp.o
	$(CC) $(CFLAGS) $(OPTORDEBUG) -I. $^ -o $@

tensor.o: tensor.cpp tensor.h
	$(CC) -c $(CFLAGS) $(OPTORDEBUG) -I. $< -o $@

DAISGram.o: DAISGram.cpp DAISGram.h
	$(CC) -c $(CFLAGS) $(OPTORDEBUG) -I. $< -o $@

libbmp.o: libbmp.cpp libbmp.h
	$(CC) -c $(CFLAGS) $(OPTORDEBUG) -I. $< -o $@

clean:
	rm main *.o
