ripser_hom: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -g

ripser_hom_optimized: ripser_hom_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_optimized.cpp -o ripser_hom_optimized -g

ripser_hom_clearing: ripser_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripser_hom_clearing -g

ripser_hom_clearing_optimized: ripser_hom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing_optimized.cpp -o ripser_hom_clearing_optimized -g

ripser_cohom_rephom: ripser_cohom_rephom.cpp
	c++ -std=c++11 -Wall ripser_cohom_rephom.cpp -o ripser_cohom_rephom -g

clean:
	rm -f \
ripser_hom \
ripser_hom_optimized \
ripser_hom_clearing \
ripser_hom_clearing_optimized \
ripser_cohom_rephom
