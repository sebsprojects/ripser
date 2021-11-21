#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"


typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            reverse_filtration_order_comp> Column;

index_diameter_t init_coboundary_and_get_pivot(ripser &ripser,
                                               const index_diameter_t simplex,
                                               const index_t dim,
                                               Column& working_coboundary) {
	if(dim == ripser.n - 1) {
		return index_diameter_t(-1, -1);
	}
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		if(get_diameter(cofacet) <= ripser.threshold) {
			working_coboundary.push(cofacet);
		}
	}
	return get_pivot(working_coboundary);
}

void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                const index_t dim) {
	std::vector<index_diameter_t> next_simplices;
	simplex_boundary_enumerator facets(ripser);
	for(index_diameter_t& simplex : simplices) {
		facets.set_simplex(simplex, dim + 1);
		while(facets.has_next()) {
			index_diameter_t facet = facets.next();
			// Threshold check
			if(get_diameter(facet) <= ripser.threshold) {
				//TODO: Uniquness check
				if(std::find(next_simplices.begin(), next_simplices.end(), facet) ==
				   next_simplices.end())
				{
					next_simplices.push_back(facet);
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = simplices.size();
}

void update_hom_class(ripser& ripser,
                      index_t dim,
                      index_diameter_t birth,
                      index_diameter_t death,
                      std::vector<index_diameter_t> R_rep) {
	if(dim == ripser.config.dim_max) {
		if(get_diameter(birth) >
		   std::max(0.0f, get_diameter(death)) * ripser.config.ratio) {
			ripser.add_hom_class(dim + 1, birth, death, R_rep);
		}
	} else {
		std::vector<homology_class>& hc = ripser.hom_classes.at(dim + 1);
		for(index_t i = 0; i < (index_t) hc.size(); i++) {
			homology_class& h = hc.at(i);
			if(get_index(h.birth) == get_index(birth)) {
				// dim(birth) > dim(death), need -inf workaround for death
				if(get_diameter(birth) >
				   std::max(0.0f, get_diameter(death)) * ripser.config.ratio) {
					std::cout << "d=" << dim + 1 << " :: ";
					for(auto s : h.representative) {
						std::cout << get_index(s) << " ";
					}
					std::cout << "vs ";
					for(auto s : R_rep) {
						std::cout << get_index(s) << " ";
					}
					std::cout << std::endl;
					h.death = death;
				} else {
					hc.erase(hc.begin() + i);
				}
				break;
			}
		}
	}
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   const index_t dim) {
	compressed_sparse_matrix V;
	entry_hash_map pivot_column_index;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		V.append_column();
		Column R_j;
		Column V_j;
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
		                                                       sigma_j,
		                                                       dim,
		                                                       R_j);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_coboundary(ripser,
				               V,
				               columns_to_reduce,
				               index_column_to_add,
				               dim,
				               V_j,
				               R_j);
				pivot = get_pivot(R_j);
				ripser.infos.at(dim).addition_count++;
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		// Write V_j to V
		std::vector<index_diameter_t> V_rep;
		V_rep.push_back(sigma_j);
		index_diameter_t e = pop_pivot(V_j);
		while(get_index(e) != -1) {
			V.push_back(e);
			V_rep.push_back(e);
			e = pop_pivot(V_j);
		}
		// Update barcode decomp
		if(get_index(pivot) != -1) {
			std::vector<index_diameter_t> R_rep;
			e = pop_pivot(R_j);
			while(get_index(e) != -1) {
				R_rep.push_back(e);
				e = pop_pivot(R_j);
			}
			update_hom_class(ripser, dim, pivot, sigma_j, R_rep);
		} else {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			ripser.infos.at(dim).class_count++;
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	index_t dim_max = std::min(ripser.config.dim_max, ripser.n - 1);
	assemble_all_simplices(ripser, simplices, dim_max);
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	ripser.infos.at(dim_max).simplex_total_count = simplices.size();
	ripser.infos.at(dim_max).simplex_reduction_count = simplices.size();
	for(index_t dim = dim_max; dim >= 0; dim--) {
		if(dim < ripser.config.dim_max) {
			assemble_columns_to_reduce(ripser, simplices, dim);
		}
		compute_pairs(ripser, simplices, dim);
	}
}


/* **************************************************************************
 * Main
 * *************************************************************************/

int main(int argc, char** argv) {
	ripser_config config;
	if(argc == 2) {
		config = read_config(argv[1]);
	} else {
		std::cerr << "error: missing config path" << std::endl;
		exit(-1);
	}
	ripser ripser(config);
	std::cout << std::endl;
	output_config(ripser, std::cout); std::cout << std::endl;
	//output_simplices(ripser, std::cout, total_filtration_order); std::cout << std::endl;
	compute_barcodes(ripser);
	output_barcode(ripser, std::cout, true); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	//write_standard_output(ripser, false, true);
	exit(0);
}

