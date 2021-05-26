#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"

typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            reverse_filtration_order_comp> Cohom_Column;

typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Hom_Column;

index_diameter_t cohom_init_coboundary_and_get_pivot(ripser &ripser,
                                                     const index_diameter_t simplex,
                                                     const index_t dim,
                                                     Cohom_Column& working_coboundary) {
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		// Threshold check
		if(get_diameter(cofacet) <= ripser.threshold) {
			working_coboundary.push(cofacet);
		}
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
		if(get_diameter(facet) <= ripser.threshold) {
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
	columns_to_reduce.clear();
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			index_diameter_t cofacet = cofacets.next();
			// Threshold check
			next_simplices.push_back(cofacet);
			if(get_diameter(cofacet) <= ripser.threshold) {
				// Clearing check
				if(pivot_column_index.find(get_index(cofacet)) ==
				   pivot_column_index.end()) {
					columns_to_reduce.push_back(cofacet);
				} else {
					ripser.barcodes.at(dim).clearing_count++;
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
}

void compute_cohomology(ripser &ripser,
                        const std::vector<index_diameter_t>& columns_to_reduce,
                        entry_hash_map& pivot_column_index,
                        const index_t dim) {
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Cohom_Column working_reduction_column; // V_j
		Cohom_Column working_coboundary;       // R_j
		// Assemble
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = cohom_init_coboundary_and_get_pivot(ripser,
		                                                             column_to_reduce,
		                                                             dim,
		                                                             working_coboundary);
		// The reduction
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
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		// Write V_j to V
		index_diameter_t e = pop_pivot(working_reduction_column);
		while(get_index(e) != -1) {
			reduction_matrix.push_back(e);
			e = pop_pivot(working_reduction_column);
		}
		// Determine Persistence Pair
		value_t birth = get_diameter(column_to_reduce);
		if(get_index(pivot) != -1) {
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.ratio) {
				// Non-essential pair
				//ripser.barcodes.at(dim).add_interval(birth, death);
			}
		} else {
			// Zero column
			//ripser.barcodes.at(dim).add_interval(birth, INF);
		}
	}
}
void compute_homology(ripser &ripser,
                      const std::vector<index_diameter_t>& columns_to_reduce,
                      entry_hash_map& pivot_column_index,
                      const index_t dim) {
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Hom_Column working_reduction_column; // V_j
		Hom_Column working_boundary;         // R_j
		// Assemble
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = hom_init_boundary_and_get_pivot(ripser,
		                                                         column_to_reduce,
		                                                         dim + 1,
		                                                         working_boundary);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_boundary(ripser,
				             reduction_matrix,
				             columns_to_reduce,
				             index_column_to_add,
				             dim + 1,
				             working_reduction_column,
				             working_boundary);
				pivot = get_pivot(working_boundary);
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		print_column(ripser, working_boundary, 1);
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
			if(death > birth * ripser.ratio) {
				// Write the representative cycle to a vector and save it
				std::vector<index_t> rep;
				while(!working_boundary.empty()) {
					rep.push_back(get_index(pop_pivot(working_boundary)));
				}
				ripser.add_hom_class(dim, birth, death, rep);
			}
		} else {
			// Zero column
			value_t birth = get_diameter(column_to_reduce);
			if(dim == 1 && birth == -INF) {
				birth = 0;
			}
			//std::cout << "[" << birth << ", )" << std::endl;
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	// Init 0-simplices
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	index_t dim;
	std::vector<index_diameter_t> columns_to_reduce;
	// dim=0 cohomology (gets discarded)
	dim = 0;
	columns_to_reduce = std::vector<index_diameter_t>(simplices);
	compute_cohomology(ripser, columns_to_reduce, pivot_column_index, dim);
	// TODO: dim=0 homology
	// dim=1 cohomology
	dim = 1;
	pivot_column_index.clear();
	pivot_column_index.reserve(columns_to_reduce.size());
	assemble_columns_to_reduce(ripser, simplices, columns_to_reduce,
	                           pivot_column_index, dim);
	compute_cohomology(ripser, columns_to_reduce, pivot_column_index, dim);
	// get pairs from the cohomology computation
	std::vector<index_diameter_t> boundary_columns_to_reduce;
	for(auto it = pivot_column_index.begin(); it != pivot_column_index.end(); ++it) {
		auto pair = *it;
		index_t pivot_row_index = pair.first;
		value_t pivot_row_diam = ripser.compute_diameter(pivot_row_index, dim + 1);
		boundary_columns_to_reduce.push_back(std::make_pair(pivot_row_index,
		                                                    pivot_row_diam));

	}
	std::sort(boundary_columns_to_reduce.begin(), boundary_columns_to_reduce.end(),
	          filtration_order);
	//for(index_diameter_t i : boundary_columns_to_reduce) {
	//	std::cout << get_index(i) << " ";
	//}
	//std::cout << std::endl;
	// dim=1 homology
	entry_hash_map boundary_pivot_column_index;
	boundary_pivot_column_index.reserve(boundary_columns_to_reduce.size());
	compute_homology(ripser, boundary_columns_to_reduce, boundary_pivot_column_index, dim);
}


/* **************************************************************************
 * Main
 * *************************************************************************/

int main(int argc, char** argv) {
	const char* filename = nullptr;
	if(argc == 2) {
		filename = argv[1];
	} else {
		std::cerr << "error: specify path to lower-distance matrix file as only arg"
				  << std::endl;
		exit(-1);
	}
	// Reading the distance matrix from file
	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "error: couldn't open file " << filename << std::endl;
		exit(-1);
	}
	DistanceMatrix dist = read_lower_distance_matrix(file_stream);
	value_t enclosing_radius = compute_enclosing_radius(dist);
	index_t dim_max = 1;
	float ratio = 1;
	ripser ripser(std::move(dist), dim_max, enclosing_radius, ratio);
	list_all_simplices(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser, false);
	exit(0);
}
