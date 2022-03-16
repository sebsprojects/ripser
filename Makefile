all: ripser-rel ripser-rel-coeff ripser-rel-debug ripser-rel-trc ripser-rel-trc-coeff ripser-rel-trc-debug

ripser-rel: ripser_rel.cpp
	c++ -std=c++11 -Wall ripser_rel.cpp -o ripser-rel -O3 -D NDEBUG

ripser-rel-coeff: ripser_rel.cpp
	c++ -std=c++11 -Wall ripser_rel.cpp -o ripser-rel-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

ripser-rel-debug: ripser_rel.cpp
	c++ -std=c++11 -Wall ripser_rel.cpp -o ripser-rel-debug -g

ripser-rel-trc: ripser_rel_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_rel_tight_representative_cycles.cpp -o ripser-rel-trc -O3 -D NDEBUG

ripser-rel-trc-coeff: ripser_rel_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_rel_tight_representative_cycles.cpp -o ripser-rel-trc-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

ripser-rel-trc-debug: ripser_rel_tight_representative_cycles.cpp
	c++ -std=c++11 -Wall ripser_rel_tight_representative_cycles.cpp -o ripser_rel-trc-debug -g

clean:
	rm -f ripser-rel ripser-rel-coeff ripser-rel-debug ripser-rel-trc ripser-rel-trc-coeff ripser-rel-trc-debug
