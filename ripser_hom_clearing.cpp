#include <iostream>

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
	return get_pivot(ripser, working_boundary);
}

void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                std::vector<index_diameter_t>& columns_to_reduce,
                                entry_hash_map& pivot_column_index,
                                const index_t dim) {
	time_point assemble_start = get_time();
	std::vector<index_diameter_t> next_simplices;
	simplex_boundary_enumerator facets(ripser);
	for(index_diameter_t& simplex : simplices) {
		facets.set_simplex(simplex, dim + 1);
		while(facets.has_next()) {
			index_diameter_t facet = facets.next();
			// Threshold check
			if(get_diameter(facet) <= ripser.threshold) {
				// TODO: This should be done more efficiently in the enumerator
				if(std::find(next_simplices.begin(), next_simplices.end(), facet) ==
				   next_simplices.end())
				{
					if(pivot_column_index.find(get_index(facet)) ==
					   pivot_column_index.end()) {
						columns_to_reduce.push_back(facet);
					} else {
						ripser.infos.at(dim).clearing_count++;
					}
					next_simplices.push_back(facet);
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), filtration_order);
	ripser.infos.at(dim).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = columns_to_reduce.size();
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const index_t dim) {
	time_point reduction_start = get_time();
	compressed_sparse_matrix V;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		ripser.add_reduction_record(dim, j, get_time());
		V.append_column();
		Column R_j;
		Column V_j;
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = init_boundary_and_get_pivot(ripser,
		                                                     sigma_j,
		                                                     dim,
		                                                     R_j);
		// The reduction
		index_t add_count = 0;
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
				add_count++;
				ripser.infos.at(dim).addition_count++;
				pivot = get_pivot(ripser, R_j);
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		// Write V_j to V
		std::vector<index_diameter_t> V_rep;
		std::vector<index_diameter_t> R_rep;
		V_rep.push_back(sigma_j);
		index_diameter_t e = pop_pivot(ripser, V_j);
		while(get_index(e) != -1) {
			V.push_back(e);
			V_rep.push_back(e);
			e = pop_pivot(ripser, V_j);
		}
		if(get_index(pivot) != -1) {
			ripser.get_current_reduction_record().to_zero = false;
			value_t birth = std::max(0.0f, get_diameter(pivot));
			if(get_diameter(sigma_j) > birth * ripser.config.ratio) {
				//e = pop_pivot(R_j);
				//while(get_index(e) != -1) {
				//	R_rep.push_back(e);
				//	e = pop_pivot(R_j);
				//}
				ripser.add_hom_class(dim - 1, pivot, sigma_j, R_rep);
				ripser.infos.at(dim).class_count++;
			} else {
				ripser.infos.at(dim).zero_pers_count++;
			}
		} else if(dim == ripser.n - 1 || dim < ripser.config.dim_max) {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			ripser.infos.at(dim).class_count++;
		}
		ripser.complete_reduction_record(get_time(), add_count, 0, -1);
	}
	ripser.infos.at(dim).reduction_dur = get_duration(reduction_start, get_time());
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
	entry_hash_map pivot_column_index;
	index_t dim_max = std::min(ripser.config.dim_max, (int) ripser.n - 1);
	time_point assemble_start = get_time();
	assemble_all_simplices(ripser, simplices, dim_max);
	std::sort(simplices.begin(), simplices.end(), filtration_order);
	ripser.infos.at(dim_max).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(dim_max).simplex_total_count = simplices.size();
	ripser.infos.at(dim_max).simplex_reduction_count = simplices.size();
	for(index_t dim = dim_max; dim >= 0; dim--) {
		if(dim == dim_max) {
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
	output_barcode(ripser, std::cout, false); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	//write_standard_output(ripser, true, false);
	write_analysis_rr(ripser, "_hom_clearing");
	exit(0);
}

