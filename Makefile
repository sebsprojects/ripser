all: ripser_cohom ripser_cohom_rel

# COHOMOLOGY

ripser_cohom: ripser_cohom.cpp
	c++ -std=c++11 -Wall ripser_cohom.cpp -o ripser_cohom -g

ripser_cohom_rel: ripser_cohom_rel.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel.cpp -o ripser_cohom_rel -g

clean:
	rm -f ripser_cohom ripser_cohom_rel
