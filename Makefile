all: ripser_hrep_hom_ref ripser_hrep_hom ripser_hrep_hom_clearing ripser_hrep_hom_cohom ripser_hrep_cohom_inverse ripser_hrep_cohom_inverse_clearing ripser_hrep_cohom_row

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

clean:
	rm -f ripser_hrep_hom_ref ripser_hrep_hom ripser_hrep_hom_clearing ripser_hrep_cohom_hom ripser_hrep_cohom_inverse ripser_hrep_cohom_inverse_clearing ripser_hrep_cohom_inverse_forwardsubs
