all: ripser_hrep_cohom_hom ripser_hrep_cohom_inverse_forwardsubs ripser_hrep_cohom_inverse_implicit

ripser_hrep_hom_ref: ripser_hrep_hom_ref.cpp
	c++ -std=c++11 -Wall ripser_hrep_hom_ref.cpp -o ripser_hrep_hom_ref -g

ripser_hrep_hom: ripser_hrep_hom.cpp
	c++ -std=c++11 -Wall ripser_hrep_hom.cpp -o ripser_hrep_hom -g

ripser_hrep_hom_clearing: ripser_hrep_hom_clearing.cpp
	c++ -std=c++11 -Wall ripser_hrep_hom_clearing.cpp -o ripser_hrep_hom_clearing -g

ripser_hrep_cohom_hom: ripser_hrep_cohom_hom.cpp
	c++ -std=c++11 -Wall ripser_hrep_cohom_hom.cpp -o ripser_hrep_cohom_hom -g

ripser_hrep_cohom_inverse: ripser_hrep_cohom_inverse.cpp
	c++ -std=c++11 -Wall ripser_hrep_cohom_inverse.cpp -o ripser_hrep_cohom_inverse -g

ripser_hrep_cohom_inverse_clearing: ripser_hrep_cohom_inverse_clearing.cpp
	c++ -std=c++11 -Wall ripser_hrep_cohom_inverse_clearing.cpp -o ripser_hrep_cohom_inverse_clearing -g

ripser_hrep_cohom_inverse_forwardsubs: ripser_hrep_cohom_inverse_forwardsubs.cpp
	c++ -std=c++11 -Wall ripser_hrep_cohom_inverse_forwardsubs.cpp -o ripser_hrep_cohom_inverse_forwardsubs -g

ripser_hrep_cohom_inverse_implicit: ripser_hrep_cohom_inverse_implicit.cpp
	c++ -std=c++11 -Wall ripser_hrep_cohom_inverse_implicit.cpp -o ripser_hrep_cohom_inverse_implicit -g

ripser_rel: ripser_rel.cpp
	c++ -std=c++11 -Wall ripser_rel.cpp -o ripser_rel -g

# OPTIMIZED VERSIONS

ripser_hrep_cohomhom_o3:
	c++ -std=c++11 -Wall ripser_hrep_cohom_hom.cpp -o ripser_hrep_cohomhom_o3 -O3

ripser_hrep_cohominv_o3:
	c++ -std=c++11 -Wall ripser_hrep_cohom_inverse_implicit.cpp -o ripser_hrep_cohominv_o3 -O3


ripser_rel_o3:
	c++ -std=c++11 -Wall ripser_rel.cpp -o ripser_rel_o3 -O3

clean:
	rm -f ripser_hrep_hom_ref ripser_hrep_hom ripser_hrep_hom_clearing ripser_hrep_cohom_hom ripser_hrep_cohom_inverse ripser_hrep_cohom_inverse_clearing ripser_hrep_cohom_inverse_forwardsubs ripser_hrep_cohom_inverse_implicit ripser_hrep_cohominv_o3 ripser_hrep_cohomhom_o3 ripser_rel ripser_rel_o3
