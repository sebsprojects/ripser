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
                                std::vector<index_diameter_t>& columns_to_reduce,
                                entry_hash_map& pivot_column_index,
                                const index_t dim) {
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			// Threshold check
			if(get_diameter(cofacet) <= ripser.threshold) {
				// Clearing check
				if(pivot_column_index.find(get_index(cofacet)) ==
				   pivot_column_index.end()) {
					columns_to_reduce.push_back(cofacet);
				} else {
					if(dim <= ripser.config.dim_max) {
						ripser.infos.at(dim).clearing_count++;
					}
				}
				next_simplices.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), reverse_filtration_order);
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = simplices.size();
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const index_t dim) {
	compressed_sparse_matrix V;
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
				ripser.infos.at(dim).addition_count++;
				pivot = get_pivot(R_j);
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		// Write V_j to V
		std::vector<index_diameter_t> V_rep;
		std::vector<index_diameter_t> R_rep;
		V_rep.push_back(sigma_j);
		index_diameter_t e = pop_pivot(V_j);
		while(get_index(e) != -1) {
			V.push_back(e);
			V_rep.push_back(e);
			e = pop_pivot(V_j);
		}
		if(get_index(pivot) != -1) {
			value_t death = std::max(0.0f, get_diameter(sigma_j));
			if(get_diameter(pivot) > death * ripser.config.ratio) {
				e = pop_pivot(R_j);
				while(get_index(e) != -1) {
					R_rep.push_back(e);
					e = pop_pivot(R_j);
				}
				ripser.add_hom_class(dim + 1, pivot, sigma_j, R_rep);
				ripser.infos.at(dim + 1).class_count++;
			}
		} else {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			ripser.infos.at(dim).class_count++;
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
	entry_hash_map pivot_column_index;
	assemble_all_simplices(ripser, simplices, 0);
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	index_t dim_max = std::min(ripser.config.dim_max, ripser.n - 1);
	for(index_t dim = 0; dim <= dim_max; dim++) {
		if(dim == 0) {
			columns_to_reduce = simplices;
		} else {
			columns_to_reduce.clear();
			assemble_columns_to_reduce(ripser,
			                           simplices,
			                           columns_to_reduce,
			                           pivot_column_index,
			                           dim);
		}
		pivot_column_index.clear();
		compute_pairs(ripser, columns_to_reduce, pivot_column_index, dim);
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

