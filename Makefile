# -g versions

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

ripser_cohom_rel_optimized: ripser_cohom_rel_optimized.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel_optimized.cpp -o ripser_cohom_rel_optimized -g

ripser_cohom_rel_rephom: ripser_cohom_rel_rephom.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel_rephom.cpp -o ripser_cohom_rel_rephom -g


# -O3 versions

ripsero3_hom_clearing: ripser_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing.cpp -o ripsero3_hom_clearing -O3

ripsero3_cohom_clearing: ripser_cohom_clearing.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing.cpp -o ripsero3_cohom_clearing -O3

ripsero3_cohom_clearing_optimized: ripser_cohom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_cohom_clearing_optimized.cpp -o ripsero3_cohom_clearing_optimized -O3

ripsero3_hom_clearing_optimized: ripser_hom_clearing_optimized.cpp
	c++ -std=c++11 -Wall ripser_hom_clearing_optimized.cpp -o ripsero3_hom_clearing_optimized -O3

ripsero3_cohom_rel: ripser_cohom_rel.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel.cpp -o ripsero3_cohom_rel -O3

ripsero3_cohom_rel_optimized: ripser_cohom_rel_optimized.cpp
	c++ -std=c++11 -Wall ripser_cohom_rel_optimized.cpp -o ripsero3_cohom_rel_optimized -O3

ripsero3_cohom_rephom: ripser_cohom_rephom.cpp
	c++ -std=c++11 -Wall ripser_cohom_rephom.cpp -o ripsero3_cohom_rephom -O3

ripsero3_cohom_repinverse: ripser_cohom_repinverse.cpp
	c++ -std=c++11 -Wall ripser_cohom_repinverse.cpp -o ripsero3_cohom_repinverse -O3



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
ripser_cohom_rel \
ripser_cohom_rel_optimized \
ripser_cohom_rel_rephom \
ripsero3_cohom_clearing_optimized \
ripsero3_hom_clearing_optimized \
ripsero3_hom_clearing \
ripsero3_cohom_rel \
ripsero3_cohom_clearing \
ripsero3_cohom_rel_optimized \
ripsero3_cohom_rephom \
ripsero3_cohom_repinverse \
ripser_dualbasis

