#include <iostream>

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
                                               Column& working_coboundary) {
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

// Takes a vector of (dim-1)-simplices as input and sets columns_to_reduce
// to contain all cofacets of those simplices
// ordered column indices (indexed by dim-simplices) of the coboundary matrix
// Sets simplices to contain all dim-simplices
//
// Clearing: We take the pivots from the reduction in dim-1 and ignore columns
// that appeared as a pivot
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

typedef std::unordered_map<index_t,
                           std::vector<index_diameter_t>> column_hash_map;

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   column_hash_map prev_reduced_cols,
                   column_hash_map reduced_cols,
                   const index_t dim) {
	compressed_sparse_matrix reduction_matrix; // V
	compressed_sparse_matrix v_inv;
	std::vector<std::pair<index_t, index_t>> nonessential_red;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Column working_reduction_column; // V_j
		Column working_coboundary;       // R_j
		// Assemble the column j (corresponding to the simplex column_to_reduce)
		// and get the pivot of that column
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
		                                                       column_to_reduce,
		                                                       dim,
		                                                       working_coboundary);
		// The reduction
		// The loop terminates on either of two conditions
		//   1. R_j got reduced to a zero column with pivot index -1 (the loop condiiton)
		//   2. R_j got fully reduced but is not zero (the break clause)
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
		// Write R_j to reduced_cols
		if(get_index(pivot) != -1) {
			std::vector<index_diameter_t> r_j;
			while(get_index(pivot) != -1) {
				r_j.push_back(pivot);
				pivot = pop_pivot(working_coboundary);
			}
			reduced_cols.insert({get_index(column_to_reduce), r_j});
		}
		// Compute the new inverse of V and store it row-order
		v_inv.append_column();
		for(index_t k = 0; k < j; ++k) {
			// Compute the inner product v_inv(k) * V_j
			index_t sum = 0;
			index_diameter_t row_simplex = columns_to_reduce.at(k);
			for(index_t i = 0; i < v_inv.size(); ++i) {
				if(v_inv.search_column(i, get_index(row_simplex))) {
					index_t col_simplex = get_index(columns_to_reduce.at(i));
					sum += reduction_matrix.search_column(j, col_simplex);
				}
			}
			if(sum % 2 == 1) {
				v_inv.push_back(row_simplex);
			}
		}
		v_inv.push_back(column_to_reduce);
		// Determine Persistence Pair
		value_t birth = get_diameter(column_to_reduce);
		if(get_index(pivot) != -1) {
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.ratio) {
				// Non-essential pair
				ripser.add_hom_class(dim, birth, death, std::vector<index_t>());
				// Store which non-essential hom class corresponds to which column
				nonessential_red.push_back(std::make_pair(ripser.barcodes.at(dim).hom_classes.size() - 1, j));
			}
		} else {
			// Zero column
			// Since we use clearing, zero columns that correspond to killing
			// simplices are filtered out by assemble_columns_to_reduce
			// Thus all zero columns are birth indices of an essential pair
			ripser.add_hom_class(dim, birth, INF, std::vector<index_t>());
		}
	}
	// Assign the rows of v_inv correspoding to non-essential pairs to their
	// respective homology class
	for(index_t i = 0; i < (index_t) nonessential_red.size(); ++i) {
		homology_class& h = ripser.barcodes.at(dim).hom_classes.at(nonessential_red.at(i).first);
		index_t row_index = nonessential_red.at(i).second;
		index_t row_simplex = get_index(columns_to_reduce.at(row_index));
		// Go over all columns in v_inv and check if they contain row_simplex
		// If yes, then add the column simplex to the representative
		for(index_t k = 0; k < v_inv.size(); ++k) {
			if(v_inv.search_column(k, row_simplex)) {
				h.representative.push_back(get_index(columns_to_reduce.at(k)));
			}
		}
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	column_hash_map prev_red_cols;
	column_hash_map red_cols;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	for(index_t dim = 0; dim <= ripser.dim_max; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		if(dim == 0) {
			columns_to_reduce = std::vector<index_diameter_t>(simplices);
		} else {
			// Takes the pivots from the previous iteration
			assemble_columns_to_reduce(ripser, simplices, columns_to_reduce,
			                           pivot_column_index, dim);
		}
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser, columns_to_reduce, pivot_column_index, prev_red_cols, red_cols, dim);
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
	print_barcodes(ripser, true);
	exit(0);
}
