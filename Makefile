ripser_hom: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -g

ripser_hom_clearing: ripser_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripser_hom_clearing -g

# OPTIMIZED VERSIONS

ripser_hom_o3: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -O3

ripser_hom_clearing_o3:
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripser_hom_clearing -O3

clean:
	rm -f \
ripser_hom ripser_hom_o3 \
ripser_hom_clearing ripser_hom_clearing_o3
