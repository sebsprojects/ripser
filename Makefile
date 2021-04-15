all: ripser_cohom ripser_hom

ripser_cohom: ripser_cohom.cpp
	c++ -std=c++11 -Wall ripser_cohom.cpp -o ripser_cohom -g

ripser_hom: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -g

clean:
	rm -f ripser_cohom ripser_hom
