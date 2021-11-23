#include <iostream>
#include <vector>

#include "ripser_core.hpp"
#include "print_utils.hpp"

// A working column is represented by this type. The entries are ordered with
// respect to reverse filtration order
typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            reverse_filtration_order_comp> Column;

// Take a dim-simplex as input and assemble to coboundary matrix column
// corresponding to that simplex. Return the pivot element of that column
index_diameter_t init_coboundary_and_get_pivot(ripser &ripser,
                                               const index_diameter_t simplex,
                                               const index_t dim,
                                               Column& working_coboundary,
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
				// relative check
				if(!ripser.is_relative_simplex(get_index(cofacet), dim)) {
					// Clearing check
					info.simplex_total_count++;
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
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
	// TODO: Is this sort necessary?
	//std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	info.assemble_dur = get_duration(assemble_start, get_time());
	info.simplex_reduction_count = columns_to_reduce.size();
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const index_t dim) {
	info& info = ripser.infos.at(dim);
	time_point reduction_start = get_time();
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		ripser.add_reduction_record(dim, j, get_time());
		reduction_matrix.append_column();
		Column working_reduction_column; // V_j
		Column working_coboundary;       // R_j
		// Assemble
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
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
				index_diameter_t e = get_zero_apparent_facet(ripser,
				                                             pivot,
				                                             dim + 1);
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
		ripser.complete_reduction_record(dim, get_time(), add_count,
		                                 app_count, red_count);
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			// Non-essential index
			value_t birth = get_diameter(column_to_reduce);
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.config.ratio) {
				ripser.infos.at(dim).class_count++;
				ripser.add_hom_class(dim, column_to_reduce, pivot);
			} else {
				ripser.infos.at(dim).zero_pers_count++;
			}
		} else {
			// Essential index (since clearing)
			ripser.infos.at(dim).class_count++;
			ripser.add_hom_class(dim, column_to_reduce, index_diameter_t(-1, INF));
		}
	}
	info.reduction_dur = get_duration(reduction_start, get_time());
}

void compute_dim0_pairs(ripser& ripser,
                        std::vector<index_diameter_t>& edges,
                        std::vector<index_diameter_t>& columns_to_reduce)
{
	time_point unionfind_start = get_time();
	union_find dset(ripser.n);
	edges = get_edges(ripser);
	std::sort(edges.begin(), edges.end(), filtration_order);
	std::vector<index_t> vertices_of_edge(2);
	// Pre-link the whole relative part
	index_t first_rel_vertex = ripser.get_first_relative_vertex();
	if(first_rel_vertex != -1) {
		for(index_t j = 0; j < ripser.n; j++) {
			if(ripser.is_relative_vertex(j) && j != first_rel_vertex) {
				dset.link(first_rel_vertex, j);
			}
		}
	}
	for(index_t j = 0; j < (index_t) edges.size(); j++) {
		index_diameter_t e = edges[j];
		//ripser.add_reduction_record(0, j, get_time());
		ripser.get_simplex_vertices(get_index(e), 1, ripser.n, vertices_of_edge);
		index_t u = dset.find(vertices_of_edge[0]);
		index_t v = dset.find(vertices_of_edge[1]);
		bool is_relative_edge = ripser.is_relative_vertex(vertices_of_edge[0]) &&
		                        ripser.is_relative_vertex(vertices_of_edge[1]);
		if(!is_relative_edge) {
			if(ripser.config.dim_max > 0) {
				ripser.infos.at(1).simplex_total_count++;
			}
			if(u != v) {
				// Zero-persistence check
				if(get_diameter(e) != 0) {
					//TODO: Dummy index
					ripser.add_hom_class(0, index_diameter_t(-1, -INF), e);
					if(ripser.config.dim_max > 0) {
						ripser.infos.at(1).clearing_count++;
					}
				}
				dset.link(u, v);
			//} else if(!is_in_zero_apparent_pair(ripser, e, 1)) {
			//TODO: Why does this cause and error if dim_amax=0
			} else if(ripser.config.dim_max > 0 &&
			          get_index(get_zero_apparent_cofacet(ripser, e, 1)) == -1) {
				columns_to_reduce.push_back(e);
			} else {
				if(ripser.config.dim_max > 0) {
					ripser.infos.at(1).apparent_count++;
				}
			}
		}
		//ripser.complete_reduction_record(0, get_time(), 0, 0, 0);
	}
	// We need the columns_to_reduce in reverse_filtration_order for the
	// reduction algorithm
	if(ripser.config.dim_max > 0) {
		ripser.infos.at(1).simplex_reduction_count = columns_to_reduce.size();
	}
	std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());
	for(index_t i = 0; i < ripser.n; ++i) {
		// Essential index
		if(dset.find(i) == i) {
			// TODO: Dummy index
			ripser.add_hom_class(0, index_diameter_t(-1, -INF), index_diameter_t(-1, INF));
		}
	}
	ripser.infos.at(0).reduction_dur = get_duration(unionfind_start, get_time());
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
	entry_hash_map pivot_column_index;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		if(!ripser.is_relative_vertex(i)) {
			columns_to_reduce.push_back(index_diameter_t(i, diam));
		}
		simplices.push_back(index_diameter_t(i, diam));
	}
	//TODO: Add in reverse order to avoid sort
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), reverse_filtration_order);
	ripser.infos.at(0).simplex_total_count = columns_to_reduce.size();
	ripser.infos.at(0).simplex_reduction_count = columns_to_reduce.size();
	for(index_t dim = 0; dim <= ripser.config.dim_max; ++dim) {
		if(dim == 0 && ripser.config.use_union_find) {
			simplices.clear();
			columns_to_reduce.clear();
			compute_dim0_pairs(ripser, simplices, columns_to_reduce);
			continue;
		}
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser,
		              columns_to_reduce,
		              pivot_column_index,
		              dim);
		if(dim < ripser.config.dim_max) {
			assemble_columns_to_reduce(ripser,
			                           simplices,
			                           columns_to_reduce,
			                           pivot_column_index,
			                           dim + 1);
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
	output_config(ripser, std::cout); std::cout << std::endl;
	//output_simplices(ripser, std::cout, total_reverse_filtration_order);
	std::cout << std::endl;
	compute_barcodes(ripser);
	output_barcode(ripser, std::cout); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	write_standard_output(ripser, false, false, total_filtration_order);
	exit(0);
}
