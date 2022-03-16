build: ripser


all: ripser ripser-coeff ripser-debug ripser-trc ripser-trc-coeff ripser-trc-debug

ripser: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser -O3 -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

ripser-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-debug -g

ripser-trc: ripser_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_tight_representative_cycles.cpp -o ripser-trc -O3 -D NDEBUG

ripser-trc-coeff: ripser_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_tight_representative_cycles.cpp -o ripser-trc-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

ripser-trc-debug: ripser_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_tight_representative_cycles.cpp -o ripser-trc-debug -g

clean:
	rm -f ripser ripser-coeff ripser-debug ripser-trc ripser-trc-coeff ripser-trc-debug
