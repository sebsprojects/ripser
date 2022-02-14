ripser_hom: ripser_hom.cpp
	c++ -std=c++11 -Wall ripser_hom.cpp -o ripser_hom -g

ripser_hom_optimized: ripser_hom_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_optimized.cpp -o ripser_hom_optimized -g

ripser_hom_clearing: ripser_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripser_hom_clearing -g

ripser_hom_clearing_optimized: ripser_hom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing_optimized.cpp -o ripser_hom_clearing_optimized -g

ripser_cohom: ripser_cohom.cpp
	c++ -std=c++11 -Wall ripser_cohom.cpp -o ripser_cohom -g

ripser_cohom_clearing: ripser_cohom_clearing.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing.cpp -o ripser_cohom_clearing -g

ripser_cohom_clearing_optimized: ripser_cohom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing_optimized.cpp -o ripser_cohom_clearing_optimized -g

ripser_cohom_rephom: ripser_cohom_rephom.cpp
	c++ -std=c++11 -Wall ripser_cohom_rephom.cpp -o ripser_cohom_rephom -g

ripser_cohom_repinverse: ripser_cohom_repinverse.cpp
	c++ -std=c++11 -Wall ripser_cohom_repinverse.cpp -o ripser_cohom_repinverse -g

ripser_cohom_rel: ripser_cohom_rel.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel.cpp -o ripser_cohom_rel -g

ripser_cohom_relrephom: ripser_cohom_relrephom.cpp
	c++ -std=c++11 -Wall ripser_cohom_relrephom.cpp -o ripser_cohom_relrephom -g

# -O3 versions

ripsero3_cohom_clearing_optimized: ripser_cohom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing_optimized.cpp -o ripsero3_cohom_clearing_optimized -O3

ripsero3_hom_clearing_optimized: ripser_hom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing_optimized.cpp -o ripsero3_hom_clearing_optimized -O3

# Misc versions

ripser_dualbasis: ripser_dualbasis.cpp
	c++ -std=c++11 -Wall ripser_dualbasis.cpp -o ripser_dualbasis -g


clean:
	rm -f \
ripser_hom \
ripser_hom_optimized \
ripser_hom_clearing \
ripser_hom_clearing_optimized \
ripser_cohom \
ripser_cohom_clearing \
ripser_cohom_clearing_optimized \
ripser_cohom_rephom \
ripser_cohom_repinverse \
ripsero3_cohom_clearing_optimized \
ripsero3_hom_clearing_optimized \
ripser_dualbasis
