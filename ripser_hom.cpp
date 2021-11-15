#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"


// A working column is represented by this type. The entries are ordered with
// respect to reverse filtration order
typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Column;

// Take a dim-simplex as input and assemble the boundary matrix column
// corresponding to that simplex. Return the pivot element of that column
index_diameter_t init_boundary_and_get_pivot(ripser &ripser,
                                             const index_diameter_t simplex,
                                             const index_t dim,
                                             Column& working_boundary) {
	if(dim == 0) {
		return index_diameter_t(-1, -1);
	}
	simplex_boundary_enumerator facets(ripser);
	// std::vector<index_diameter_t> cofacet_entries;
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		working_boundary.push(facet);
	}
	return get_pivot(working_boundary);
}

// Takes a vector of (dim-1)-simplices as input and sets columns_to_reduce
// to contain all cofacets of those simplices
// ordered column indices (indexed by dim-simplices) of the boundary matrix
// Sets simplices to contain all dim-simplices
void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                std::vector<index_diameter_t>& columns_to_reduce,
                                const index_t dim) {
	columns_to_reduce.clear();
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			next_simplices.push_back(cofacet);
			columns_to_reduce.push_back(cofacet);
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), filtration_order);
	if(dim <= ripser.config.dim_max) {
		ripser.infos.at(dim).simplex_total_count = simplices.size();
		ripser.infos.at(dim).simplex_reduction_count = columns_to_reduce.size();
	}
}

void update_hom_class(ripser& ripser,
                      index_t dim,
                      index_diameter_t birth,
                      index_diameter_t death,
                      std::vector<index_diameter_t> R_rep) {
	for(homology_class& h : ripser.hom_classes) {
		if(h.dim == dim - 1 && get_index(h.birth) == get_index(birth)) {
			std::cout << "d=" << dim - 1 << " :: ";
			for(auto s : h.representative) {
				std::cout << get_index(s) << " ";
			}
			std::cout << "vs ";
			for(auto s : R_rep) {
				std::cout << get_index(s) << " ";
			}
			std::cout << std::endl;
			h.death = death;
			break;
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
		index_diameter_t pivot = init_boundary_and_get_pivot(ripser,
		                                                     sigma_j,
		                                                     dim,
		                                                     R_j);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_boundary(ripser,
				             V,
				             columns_to_reduce,
				             index_column_to_add,
				             dim,
				             V_j,
				             R_j);
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
		e = pop_pivot(R_j);
		while(get_index(e) != -1) {
			R_rep.push_back(e);
			e = pop_pivot(R_j);
		}
		// Update barcode decomp
		if(get_index(pivot) != -1) {
			// Non-essential death index (birth in dimension dim - 1)
			update_hom_class(ripser, dim, pivot, sigma_j, R_rep);
		} else {
			// Zero column
			if(dim <= ripser.config.dim_max) {
				ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
				ripser.infos.at(dim).class_count++;
			}
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
	for(index_t i = ripser.n - 1; i >= 0; i--) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
		columns_to_reduce.push_back(index_diameter_t(i, diam));
	}
	ripser.infos.at(0).simplex_total_count = ripser.n;
	ripser.infos.at(0).simplex_reduction_count = ripser.n;
	for(index_t dim = 0; dim <= ripser.config.dim_max + 1; dim++) {
		if(dim > 0) {
			columns_to_reduce.clear();
			assemble_columns_to_reduce(ripser, simplices, columns_to_reduce, dim);
		}
		compute_pairs(ripser, columns_to_reduce, dim);
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
	output_simplices(ripser, std::cout, total_filtration_order); std::cout << std::endl;
	compute_barcodes(ripser);
	output_barcode(ripser, std::cout); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	write_standard_output(ripser, false, true);
	exit(0);
}

