all: ripser_cohom

ripser_cohom: ripser_cohom.cpp
	c++ -std=c++11 -Wall ripser_cohom.cpp -o ripser_cohom -g

clean:
	rm -f ripser_cohom
