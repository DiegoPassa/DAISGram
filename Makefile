MAIN = main
FLAGS = -std=c++11 -lm -O3

all: testbmp main main_tensor

debug: FLAGS=-g
debug: main

tensor.o: tensor.cpp
	g++ tensor.cpp -o tensor.o -c $(FLAGS)

libbmp.o: libbmp.cpp
	g++ libbmp.cpp -o libbmp.o -c $(FLAGS)

DAISGram.o: DAISGram.cpp libbmp.o
	g++ DAISGram.cpp -o DAISGram.o -c $(FLAGS)

main: tensor.o DAISGram.o libbmp.o main.cpp 
	g++ libbmp.o tensor.o DAISGram.o main.cpp -o $(MAIN) $(FLAGS)

testbmp: test_bmplib.cpp libbmp.o
	g++ libbmp.o test_bmplib.cpp -o test_bmplib $(FLAGS)

main_tensor: tensor.o main_tensor.cpp 
	g++ tensor.o main_tensor.cpp -o main_tensor $(FLAGS)

clean:
<<<<<<< HEAD
	rm $(MAIN) *.o *.bmp
=======
	rm $(MAIN) *.o 
>>>>>>> 5994379c2e9a2fef0f3b6192386a473179f024f1
