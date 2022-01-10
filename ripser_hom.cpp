#include <iostream>
#include <algorithm>

#include "ripser_core.hpp"
#include "print_utils.hpp"


typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Column;

index_diameter_t init_boundary_and_get_pivot(ripser &ripser,
                                             const index_diameter_t simplex,
                                             const index_t dim,
                                             Column& working_boundary) {
	if(dim == 0) {
		return index_diameter_t(-1, -1);
	}
	simplex_boundary_enumerator facets(ripser);
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		if(get_diameter(facet) <= ripser.threshold) {
			working_boundary.push(facet);
		}
	}
	return get_pivot(working_boundary);
}

void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                const index_t dim) {
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			// Threshold check
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(simplices.begin(), simplices.end(), filtration_order);
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = simplices.size();
}

void update_hom_class(ripser& ripser,
                      index_t dim,
                      index_diameter_t birth,
                      index_diameter_t death,
                      std::vector<index_diameter_t> R_rep) {
	std::vector<homology_class>& hc = ripser.hom_classes.at(dim - 1);
	for(index_t i = 0; i < (index_t) hc.size(); i++) {
		homology_class& h = hc.at(i);
		if(get_index(h.birth) == get_index(birth)) {
			// dim(death) > dim(birth), need -inf workaround for birth
			if(get_diameter(death) >
			   std::max(0.0f, get_diameter(birth)) * ripser.config.ratio) {
				//std::cout << "replacing: d=" << dim - 1 << " :: ";
				for(auto s : h.representative) {
				//	std::cout << get_index(s) << " ";
				}
				//std::cout << "with ";
				for(auto s : R_rep) {
				//	std::cout << get_index(s) << " ";
				}
				//std::cout << std::endl;
				h.death = death;
			} else {
				hc.erase(hc.begin() + i);
			}
			break;
		}
	}
}

std::ofstream os("./mat-r.txt", std::ofstream::trunc);

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
		//os << dim << "-" << get_index(sigma_j) << " ";
		//std::cout << "d=" << dim << " [ ";
		while(get_index(e) != -1) {
			//std::cout << get_index(e) << " ";
			V.push_back(e);
			V_rep.push_back(e);
			//os << dim << "-" << get_index(e) << " ";
			e = pop_pivot(V_j);
		}
		//os << std::endl;
		//std::cout << get_index(sigma_j) << " ]" << std::endl;
		// Update barcode decomp
		if(get_index(pivot) != -1) {
			std::vector<index_diameter_t> R_rep;
			e = pop_pivot(R_j);
			while(get_index(e) != -1) {
				os << dim - 1 << "-" << get_index(e) << " ";
				R_rep.push_back(e);
				e = pop_pivot(R_j);
			}
			os << std::endl;
			update_hom_class(ripser, dim, pivot, sigma_j, R_rep);
		} else {
			os << std::endl;
			//std::cout << "adding: d=" << dim << " :: ";
			//for(auto s : V_rep) {
			//	std::cout << get_index(s) << " ";
			//}
			//std::cout << std::endl;
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			ripser.infos.at(dim).class_count++;
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	assemble_all_simplices(ripser, simplices, 0);
	std::sort(simplices.begin(), simplices.end(), filtration_order);
	ripser.infos.at(0).simplex_total_count = ripser.n;
	ripser.infos.at(0).simplex_reduction_count = ripser.n;
	index_t dim_max = std::min(ripser.config.dim_max, (int) ripser.n - 1);
	for(index_t dim = 0; dim <= dim_max; dim++) {
		if(dim > 0) {
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
	write_standard_output(ripser, true, true);

	//std::ofstream os(ripser.config.output_path + "/mat.txt", std::ofstream::trunc);
	//output_boundary_matrix(ripser, os);

	exit(0);
}
