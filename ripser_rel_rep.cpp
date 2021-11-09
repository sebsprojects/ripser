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
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	bool check_for_emergent_pair = true;
	// TODO: Find a more efficient solution to this annoying problem
	std::vector<index_diameter_t> working_coboundary_buffer;
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		// Threshold check
		if(get_diameter(cofacet) <= ripser.threshold) {
			working_coboundary_buffer.push_back(cofacet);
			// Emergent pair candidate check
			if(check_for_emergent_pair &&
			   (get_diameter(simplex) == get_diameter(cofacet))) {
				// Check if the candidate is viable, if not then we certainly
				// do not have an emergent pair
				// TODO: Why is the check for the apparent facet necessary?
				if((pivot_column_index.find(get_index(cofacet)) ==
				    pivot_column_index.end()) &&
				   (get_index(get_zero_apparent_facet(ripser, cofacet, dim + 1)) == -1)) {
					// working_coboundary is the 0-column
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
	return get_pivot(working_coboundary);
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
		if(get_diameter(facet) <= ripser.threshold &&
		   !ripser.is_relative_simplex(get_index(facet), dim - 1)) {
			working_boundary.push(facet);
		}
	}
	return get_pivot(working_boundary);
}

void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                std::vector<index_diameter_t>& columns_to_reduce,
                                entry_hash_map& pivot_column_index,
                                const index_t dim) {
	info& info = ripser.infos.at(dim);
	time_point assemble_start = get_time();
	columns_to_reduce.clear();
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			// Threshold check
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
				// Clearing check
				if(pivot_column_index.find(get_index(cofacet)) ==
				   pivot_column_index.end()) {
					// Apparent Pair check
					if(!is_in_zero_apparent_pair(ripser, cofacet, dim)) {
						columns_to_reduce.push_back(cofacet);
					} else {
						info.apparent_count++;
					}
				} else {
					info.clearing_count++;
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
	time_point assemble_end = get_time();
	info.assemble_dur = get_duration(assemble_start, assemble_end);
	info.simplex_total_count = simplices.size();
	info.simplex_reduction_count = columns_to_reduce.size();
}

void compute_cohomology(ripser &ripser,
                        const std::vector<index_diameter_t>& columns_to_reduce,
                        entry_hash_map& pivot_column_index,
                        const index_t dim) {
	info& info = ripser.infos.at(dim);
	time_point reduction_start = get_time();
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		ripser.add_reduction_record(dim, j, get_time());
		reduction_matrix.append_column();
		Cohom_Column working_reduction_column; // V_j
		Cohom_Column working_coboundary;       // R_j
		// Assemble
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = cohom_init_coboundary_and_get_pivot(ripser,
		                                                             column_to_reduce,
		                                                             dim,
		                                                             working_coboundary,
		                                                             pivot_column_index);
		// The reduction
		index_t add_count = 0;
		index_t app_count = 0;
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_coboundary(ripser,
				               reduction_matrix,
				               columns_to_reduce,
				               index_column_to_add,
				               dim,
				               working_reduction_column,
				               working_coboundary);
				pivot = get_pivot(working_coboundary);
				add_count++;
			} else {
				index_diameter_t e = get_zero_apparent_facet(ripser, pivot, dim + 1);
				if(get_index(e) != -1) {
					add_simplex_coboundary(ripser,
					                       e,
					                       dim,
					                       working_reduction_column,
					                       working_coboundary);
					pivot = get_pivot(working_coboundary);
					app_count++;
				} else {
					pivot_column_index.insert({get_index(pivot), j});
					break;
				}
			}
		}
		// Write V_j to V
		index_diameter_t e = pop_pivot(working_reduction_column);
		while(get_index(e) != -1) {
			reduction_matrix.push_back(e);
			e = pop_pivot(working_reduction_column);
		}
		//TODO(seb): This may be expensive
		index_t red_count = working_coboundary.size();
		//e = pop_pivot(working_coboundary);
		//while(get_index(e) != -1) {
		//	red_count++;
		//	e = pop_pivot(working_coboundary);
		//}
		ripser.complete_reduction_record(dim, get_time(), add_count, app_count, red_count);
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			// Non-essential pair. Get's ignored for output here
		} else {
			// Essential index (since clearing)
			ripser.add_hom_class(dim, column_to_reduce, index_diameter_t(-1, INF));
		}
	}
	time_point reduction_end = get_time();
	info.reduction_dur = get_duration(reduction_start, reduction_end);
}

// dim corresponds to the dim of simplices in columns_to_reduce
void compute_homology(ripser &ripser,
                      const std::vector<index_diameter_t>& columns_to_reduce,
                      entry_hash_map& pivot_column_index,
                      const index_t dim) {
	info& info = ripser.infos.at(dim - 1);
	time_point rep_start = get_time();
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Hom_Column working_reduction_column; // V_j
		Hom_Column working_boundary;         // R_j
		// Assemble
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = hom_init_boundary_and_get_pivot(ripser,
		                                                         column_to_reduce,
		                                                         dim,
		                                                         working_boundary);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_relative_boundary(ripser,
				                      reduction_matrix,
				                      columns_to_reduce,
				                      index_column_to_add,
				                      dim,
				                      working_reduction_column,
				                      working_boundary);
				pivot = get_pivot(working_boundary);
			} else {
				index_diameter_t e = get_zero_apparent_cofacet(ripser, pivot, dim - 1);
				if(get_index(e) != -1) {
					// TODO: Why is the necessary w.r.t. apparent pairs?
					add_relative_simplex_boundary(ripser,
					                              e,
					                              dim,
					                              working_reduction_column,
					                              working_boundary);
					pivot = get_pivot(working_boundary);
				} else {
					pivot_column_index.insert({get_index(pivot), j});
					break;
				}
			}
		}
		// Write V_j to V
		index_diameter_t e = pop_pivot(working_reduction_column);
		while(get_index(e) != -1) {
			reduction_matrix.push_back(e);
			e = pop_pivot(working_reduction_column);
		}
		if(get_index(pivot) != -1) {
			// Non-essential death index (birth in dimension dim - 1)
			value_t birth = get_diameter(pivot);
			value_t death = get_diameter(column_to_reduce);
			if(death > birth * ripser.config.ratio) {
				std::vector<index_diameter_t> rep;
				e = pop_pivot(working_boundary);
				while(get_index(e) != -1) {
					rep.push_back(e);
					e = pop_pivot(working_boundary);
				}
				ripser.add_hom_class(dim - 1, pivot, column_to_reduce, rep);
			}
		} else {
			// This case will not occur since we only do the reduction
			// for known non-essential pairs
			assert(false);
		}
	}
	time_point rep_end = get_time();
	info.representative_dur = get_duration(rep_start, rep_end);
}

void compute_barcodes(ripser& ripser) {
	// Init 0-simplices
	std::vector<index_diameter_t> simplices;
	for(index_t i = 0; i < ripser.n; i++) {
		if(!ripser.is_relative_vertex(i)) {
			value_t diam = ripser.compute_diameter(i, 0);
			simplices.push_back(index_diameter_t(i, diam));
		}
	}
	std::vector<index_diameter_t> columns_to_reduce(simplices);
	std::vector<index_diameter_t> boundary_columns_to_reduce;
	entry_hash_map pivot_column_index;
	entry_hash_map boundary_pivot_column_index;
	// cohom in dim=0
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	compute_cohomology(ripser, columns_to_reduce, pivot_column_index, 0);
	for(index_t dim = 1; dim <= ripser.config.dim_max + 1; ++dim) {
		// collect columns for homology
		boundary_columns_to_reduce.clear();
		for(auto it = pivot_column_index.begin(); it != pivot_column_index.end(); ++it) {
			auto pair = *it;
			index_t pivot_row_index = pair.first;
			value_t pivot_row_diam = ripser.compute_diameter(pivot_row_index, dim);
			boundary_columns_to_reduce.push_back(std::make_pair(pivot_row_index,
			                                     pivot_row_diam));
		}
		std::sort(boundary_columns_to_reduce.begin(), boundary_columns_to_reduce.end(),
		          filtration_order);
		// compute homology in dim
		boundary_pivot_column_index.clear();
		boundary_pivot_column_index.reserve(boundary_columns_to_reduce.size());
		compute_homology(ripser,
		                 boundary_columns_to_reduce,
		                 boundary_pivot_column_index,
		                 dim);
		// compute cohomology in dim+1
		if(dim <= ripser.config.dim_max) {
			assemble_columns_to_reduce(ripser, simplices, columns_to_reduce,
	                                   pivot_column_index, dim);
			pivot_column_index.clear();
			pivot_column_index.reserve(columns_to_reduce.size());
			compute_cohomology(ripser, columns_to_reduce, pivot_column_index, dim);
		}
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
	print_config(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser);
	print_infos(ripser);
	write_barcode(ripser, true);
	exit(0);
}
