#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"


typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            reverse_filtration_order_comp> Cohom_Column;

typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Hom_Column;

template <typename Column>
void add_relative_simplex_boundary(ripser &ripser,
                                   const index_diameter_t simplex,
                                   const index_t dim,
                                   Column& working_reduction_column,
                                   Column& working_boundary) {
	simplex_boundary_enumerator facets(ripser);
	working_reduction_column.push(simplex);
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		if(get_diameter(facet) <= ripser.threshold &&
		   !ripser.is_relative_simplex(get_index(facet), dim - 1)) {
			working_boundary.push(facet);
		}
	}
}

template <typename Column>
void add_relative_boundary(ripser& ripser,
                           compressed_sparse_matrix& reduction_matrix,
                           const std::vector<index_diameter_t>& columns_to_reduce,
                           const size_t index_column_to_add,
                           const size_t dim,
                           Column& working_reduction_column,
                           Column& working_boundary) {
	//TODO(seb): Do we need the correct diameter here?
	index_diameter_t column_to_add(columns_to_reduce.at(index_column_to_add));
	// Computation of R_j due to implicit reduction
	add_relative_simplex_boundary(ripser,
	                              column_to_add,
	                              dim,
	                              working_reduction_column,
	                              working_boundary);
	for(index_t i = reduction_matrix.column_start(index_column_to_add);
	    i < reduction_matrix.column_end(index_column_to_add);
	    ++i) {
		index_diameter_t simplex = reduction_matrix.get_entry(i);
		add_relative_simplex_boundary(ripser,
		                              simplex,
		                              dim,
		                              working_reduction_column,
		                              working_boundary);
	}
}

index_diameter_t cohom_init_coboundary_and_get_pivot(ripser &ripser,
                                                     const index_diameter_t simplex,
                                                     const index_t dim,
                                                     Cohom_Column& working_coboundary,
                                                     entry_hash_map& pivot_column_index) {
	if(dim == ripser. n - 1) {
		return index_diameter_t(-1, -1);
	}
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	std::vector<index_diameter_t> working_coboundary_buffer;
	bool check_for_emergent_pair = true;
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		// Threshold check
		if(get_diameter(cofacet) <= ripser.threshold) {
			working_coboundary_buffer.push_back(cofacet);
			// Emergent pair check
			if(check_for_emergent_pair &&
			   (get_diameter(simplex) == get_diameter(cofacet)))
			{
				if((pivot_column_index.find(get_index(cofacet)) ==
				    pivot_column_index.end()) &&
				   (get_index(get_zero_apparent_facet(ripser, cofacet, dim + 1)) == -1)) 
				{
					ripser.infos.at(dim).emergent_count++;
					return cofacet;
				}
				check_for_emergent_pair = false;
			}
		}
	}
	for(index_diameter_t cofacet : working_coboundary_buffer) {
		working_coboundary.push(cofacet);
	}
	return get_pivot(ripser, working_coboundary);
}

index_diameter_t hom_init_boundary_and_get_pivot(ripser &ripser,
                                                 const index_diameter_t simplex,
                                                 const index_t dim,
                                                 Hom_Column& working_boundary) {
	if(dim == 0) {
		return index_diameter_t(-1, -1);
	}
	simplex_boundary_enumerator facets(ripser);
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		// Threshold check
		if(get_diameter(facet) <= ripser.threshold) {
			// Relative check
			if(!ripser.is_relative_simplex(get_index(facet), dim - 1)) {
				working_boundary.push(facet);
			}
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
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			// Threshold check
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
				// Relative check
				if(!ripser.is_relative_simplex(get_index(cofacet), dim)) {
					ripser.infos.at(dim).simplex_total_count++;
					// Clearing check
					if(pivot_column_index.find(get_index(cofacet)) ==
					   pivot_column_index.end()) {
						// Apparent Pair check
						if(!is_in_zero_apparent_pair(ripser, cofacet, dim)) {
							columns_to_reduce.push_back(cofacet);
						} else {
							ripser.infos.at(dim).apparent_count++;
						}
					} else {
						ripser.infos.at(dim).clearing_count++;
					}
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
	ripser.infos.at(dim).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(dim).simplex_reduction_count = columns_to_reduce.size();
}

void compute_cohomology(ripser &ripser,
                        const std::vector<index_diameter_t>& columns_to_reduce,
                        entry_hash_map& pivot_column_index,
                        const index_t dim,
                        std::vector<index_diameter_t>& essential_indices,
                        std::vector<index_diameter_t>& non_essential_indices) {
	time_point reduction_start = get_time();
	compressed_sparse_matrix V;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		ripser.add_reduction_record(dim, j, get_time());
		V.append_column();
		Cohom_Column V_j;
		Cohom_Column R_j;
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = cohom_init_coboundary_and_get_pivot(ripser,
		                                                             sigma_j,
		                                                             dim,
		                                                             R_j,
		                                                             pivot_column_index);
		// The reduction
		index_t add_count;
		index_t app_count;
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
				add_count++;
				ripser.infos.at(dim).addition_count++;
				pivot = get_pivot(ripser, R_j);
			} else {
				index_diameter_t e = get_zero_apparent_facet(ripser, pivot, dim + 1);
				if(get_index(e) != -1) {
					add_simplex_coboundary(ripser,
					                       e,
					                       dim,
					                       V_j,
					                       R_j);
					app_count++;
					pivot = get_pivot(ripser, R_j);
				} else {
					pivot_column_index.insert({get_index(pivot), j});
					break;
				}
			}
		}
		// Write V_j to V
		index_diameter_t e = pop_pivot(ripser, V_j);
		while(get_index(e) != -1) {
			ripser.get_current_reduction_record().to_zero = false;
			V.push_back(e);
			e = pop_pivot(ripser, V_j);
		}
		// Determine Persistence Pair
		// Note: We do need zero-persistence indices as well!
		if(get_index(pivot) != -1) {
			non_essential_indices.push_back(pivot);
		} else {
			essential_indices.push_back(sigma_j);
		}
		ripser.complete_reduction_record(get_time(), add_count, 0, -1);
	}
	ripser.infos.at(dim).reduction_dur = get_duration(reduction_start, get_time());
}

void compute_homology(ripser &ripser,
                      const std::vector<index_diameter_t>& columns_to_reduce,
                      entry_hash_map& pivot_column_index,
                      const index_t dim) {
	//time_point reduction_start = get_time();
	compressed_sparse_matrix V;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		V.append_column();
		Hom_Column R_j;
		Hom_Column V_j;
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = hom_init_boundary_and_get_pivot(ripser,
		                                                         sigma_j,
		                                                         dim,
		                                                         R_j);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_relative_boundary(ripser,
				                      V,
				                      columns_to_reduce,
				                      index_column_to_add,
				                      dim,
				                      V_j,
				                      R_j);
				//ripser.infos.at(dim).addition_count++;
				pivot = get_pivot(ripser, R_j);
			} else {
				index_diameter_t e = get_zero_apparent_cofacet(ripser, pivot, dim - 1);
				if(get_index(e) != -1) {
					add_relative_simplex_boundary(ripser,
					                              e,
					                              dim,
					                              V_j,
					                              R_j);
					pivot = get_pivot(ripser, R_j);
				} else {
					pivot_column_index.insert({get_index(pivot), j});
					break;
				}
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
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			value_t birth = std::max(0.0f, get_diameter(pivot));
			if(get_diameter(sigma_j) > birth * ripser.config.ratio) {
				e = pop_pivot(ripser, R_j);
				while(get_index(e) != -1) {
					R_rep.push_back(e);
					e = pop_pivot(ripser, R_j);
				}
				ripser.add_hom_class(dim - 1, pivot, sigma_j, R_rep);
				//ripser.infos.at(dim).class_count++;
			} else {
				//ripser.infos.at(dim).zero_pers_count++;
			}
		} else {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			//ripser.infos.at(dim).class_count++;
		}
	}
	//ripser.infos.at(dim).reduction_dur = get_duration(reduction_start, get_time());
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
	std::vector<index_diameter_t> essential_indices;
	std::vector<index_diameter_t> non_essential_indices_dim;
	std::vector<index_diameter_t> non_essential_indices_dimm1;
	entry_hash_map pivot_column_index;
	entry_hash_map boundary_pivot_column_index;
	time_point assemble_start = get_time();
	assemble_all_simplices(ripser, simplices, 0, false);
	assemble_all_simplices(ripser, columns_to_reduce, 0, true);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), reverse_filtration_order);
	ripser.infos.at(0).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(0).simplex_total_count = columns_to_reduce.size();
	ripser.infos.at(0).simplex_reduction_count = columns_to_reduce.size();
	index_t dim_max = std::min(ripser.config.dim_max, (int) ripser.n - 1);
	// cohom in dim=0
	for(index_t dim = 0; dim <= dim_max + 1; ++dim) {
		if (dim > 0 && dim <= dim_max) {
			columns_to_reduce.clear();
			assemble_columns_to_reduce(ripser,
			                           simplices,
			                           columns_to_reduce,
			                           pivot_column_index,
			                           dim);
		}
		pivot_column_index.clear();
		non_essential_indices_dimm1.swap(non_essential_indices_dim);
		essential_indices.clear();
		non_essential_indices_dim.clear();
		if(dim <= dim_max) {
			compute_cohomology(ripser,
		                       columns_to_reduce,
		                       pivot_column_index,
		                       dim,
		                       essential_indices,
		                       non_essential_indices_dim);
		}
		columns_to_reduce.clear();
		columns_to_reduce = essential_indices;
		columns_to_reduce.insert(columns_to_reduce.end(),
		                         non_essential_indices_dimm1.begin(),
		                         non_essential_indices_dimm1.end());
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), filtration_order);
		boundary_pivot_column_index.clear();
		time_point rep_start = get_time();
		compute_homology(ripser,
		                 columns_to_reduce,
		                 boundary_pivot_column_index,
		                 dim);
		//ripser.infos.at(dim).representative_dur = get_duration(rep_start, get_time());
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
	write_standard_output(ripser, true, false);
	exit(0);
}
