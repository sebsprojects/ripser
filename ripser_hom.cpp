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
	simplex_boundary_enumerator facets(ripser);
	// std::vector<index_diameter_t> cofacet_entries;
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		//TODO(seb): Check diam <= threshold
		//cofacet_entries.push_back(cofacet);
	//}
	//for(index_diameter_t cofacet : cofacet_entries) {
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
			//TODO(seb): Check diam <= threshold?
			index_diameter_t cofacet = cofacets.next();
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
				columns_to_reduce.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              filtration_order);
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
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
		                                                     dim + 1,
		                                                     working_boundary);
		print_column(ripser, working_boundary, dim);
		print_simplex(ripser, get_index(pivot), dim);
		value_t birth = get_diameter(column_to_reduce);
		if(dim == 0 && birth == -INF) {
			birth = 0;
		}
		//print_column(ripser, working_coboundary, dim);
		// The reduction
		//print_simplex(ripser, get_index(pivot), dim);
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
		// Write V_j to V
		index_diameter_t e = pop_pivot(working_reduction_column);
		while(get_index(e) != -1) {
			reduction_matrix.push_back(e);
			e = pop_pivot(working_reduction_column);
		}
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.ratio) {
				std::cout << " [" << birth << "," << death << ")" << std::endl;
			}
		} else {
			// Zero column, persistent homology "pair"
			std::cout << " [" << birth << ", )" << std::endl;
		}
		//print_mat(reduction_matrix);
		//print_v(reduction_matrix, columns_to_reduce);
		//std::cout << "-------------------------------------" << std::endl;
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	for(index_t dim = 0; dim <= ripser.dim_max; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		entry_hash_map pivot_column_index;
		assemble_columns_to_reduce(ripser, simplices, columns_to_reduce, dim + 1);
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser, columns_to_reduce, pivot_column_index, dim);
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
