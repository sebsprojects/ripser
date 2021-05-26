all: ripser_cohom ripser_cohom_clearing ripser_cohom_unionfind ripser_cohom_emergent ripser_cohom_apparent ripser_hom ripser_hom_clearing

# COHOMOLOGY

ripser_cohom: ripser_cohom.cpp
	c++ -std=c++11 -Wall ripser_cohom.cpp -o ripser_cohom -g

ripser_cohom_clearing: ripser_cohom_clearing.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing.cpp -o ripser_cohom_clearing -g

ripser_cohom_unionfind: ripser_cohom_unionfind.cpp
	c++ -std=c++11 -Wall ripser_cohom_unionfind.cpp -o ripser_cohom_unionfind -g

ripser_cohom_emergent: ripser_cohom_emergent.cpp
	c++ -std=c++11 -Wall ripser_cohom_emergent.cpp -o ripser_cohom_emergent -g

ripser_cohom_apparent: ripser_cohom_apparent.cpp
	c++ -std=c++11 -Wall ripser_cohom_apparent.cpp -o ripser_cohom_apparent -g


# HOMOLOGY

ripser_hom: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -g

ripser_hom_clearing: ripser_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripser_hom_clearing -g

ripser_hom_cohom: ripser_hom_cohom.cpp
	c++ -std=c++11 -Wall ripser_hom_cohom.cpp -o ripser_hom_cohom -g

clean:
	rm -f ripser_cohom ripser_cohom_clearing ripser_cohom_unionfind ripser_cohom_emergent ripser_cohom_apprent ripser_hom ripser_hom_clearing ripser_hom_cohom
