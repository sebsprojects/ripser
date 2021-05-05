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
		// Threshold check
		if(get_diameter(facet) <= ripser.threshold) {
			working_boundary.push(facet);
		}
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
	simplex_boundary_enumerator facets(ripser);
	for(index_diameter_t& simplex : simplices) {
		facets.set_simplex(simplex, dim + 1);
		while(facets.has_next()) {
			index_diameter_t facet = facets.next();
			// Threshold check
			if(get_diameter(facet) <= ripser.threshold) {
				// TODO: This should be done more efficiently in the enumerator
				if(std::find(next_simplices.begin(), next_simplices.end(), facet) == next_simplices.end()) {
					next_simplices.push_back(facet);
					columns_to_reduce.push_back(facet);
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), filtration_order);
}

// TODO: Output for essential pairs is shifted by one dimension
void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   entry_hash_map& previous_pivots,
                   const index_t dim) {
	std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
	compressed_sparse_matrix reduction_matrix; // V
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Column working_reduction_column; // V_j
		Column working_boundary;         // R_j
		// Assemble the column j (corresponding to the simplex column_to_reduce)
		// and get the pivot of that column
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = init_boundary_and_get_pivot(ripser,
		                                                     column_to_reduce,
		                                                     dim,
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
				             dim,
				             working_reduction_column,
				             working_boundary);
				pivot = get_pivot(working_boundary);
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
		if(get_index(pivot) != -1) {
			// Non-essential death index (birth in dimension dim - 1)
			value_t birth = get_diameter(pivot);
			if(dim == 1 && birth == -INF) {
				birth = 0;
			}
			value_t death = get_diameter(column_to_reduce);
			if(death > birth * ripser.ratio) {
				std::cout << " [" << birth << "," << death << ")" << std::endl;
			}
		} else {
			// Zero column
			auto pair = previous_pivots.find(get_index(column_to_reduce));
			if(pair == previous_pivots.end()) {
				value_t birth = get_diameter(column_to_reduce);
				if(dim == 0 && birth == -INF) {
					birth = 0;
				}
				// Since in homology (in contrast to cohomology) we need to go
				// up to dim_max+1 to compute non-essential pairs of dimension
				// dim_max we would get essential pairs in too high dimension
				// Thus this clause
				if(dim <= ripser.dim_max) {
					// Essential birth index (birth in dimension dim)
					std::cout << " [" << birth << ", )" << std::endl;
				}
			} else {
				// Birth index of non-essential pair
			}
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> previous_simplices;
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	entry_hash_map previous_pivots;
	// Vertices
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	// Enumerate all (max_dim+1) simplices by iteratively going up from dim=0
	// and enumerating coboundaries
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_t dim = 0; dim < ripser.dim_max + 1; dim++) {
		previous_simplices.swap(simplices);
		simplices.clear();
		for(index_diameter_t simp : previous_simplices) {
			cofacets.set_simplex(simp, dim);
			while(cofacets.has_next(false)) {
				index_diameter_t cofacet = cofacets.next();
				// Threshold check
				if(get_diameter(cofacet) <= ripser.threshold) {
					simplices.push_back(cofacet);
				}
			}
		}
		std::sort(simplices.begin(), simplices.end(), filtration_order);
	}
	for(index_t dim = ripser.dim_max + 1; dim >= 0; --dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		if(dim == ripser.dim_max + 1) {
			columns_to_reduce = std::vector<index_diameter_t>(simplices);
		} else {
			assemble_columns_to_reduce(ripser, simplices, columns_to_reduce, dim);
		}
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser, columns_to_reduce, pivot_column_index, previous_pivots, dim);
		previous_pivots.swap(pivot_column_index);
	}
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
	exit(0);
}
