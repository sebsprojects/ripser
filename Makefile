build: ripser


all: ripser ripser-coeff ripser-debug ripser-bare ripser-1


ripser: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser -O3 -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

ripser-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-debug -g

ripser-bare: ripser-bare.cpp
	c++ -std=c++11 -Wall ripser-bare.cpp -o ripser-bare -g

ripser-1: ripser-1.cpp
	c++ -std=c++11 -Wall ripser-1.cpp -o ripser-1 -g

ripser-t1: ripser-t1.cpp
	c++ -std=c++11 -Wall ripser-t1.cpp -o ripser-t1 -g

clean:
	rm -f ripser ripser-coeff ripser-debug ripser-bare ripser-1 ripser-t1
