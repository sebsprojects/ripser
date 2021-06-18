#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"

/*
 * Compute the persistence barcode and homology representatives in increasing
 * dimension from [0..dim_threshold]
 * Reduction algorithm for homology. Uses optimizations
 *   dim=0 union-find (replaces the reduction for dim=0)
 *   emergent-pairs (sometimes omittes init_boundary_and_get_pivot)
 *   apparent-pairs (sometimes omittes adding to columns_to_reduce)
 */

typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Column;

// dim is the dimension of simplex
index_diameter_t init_boundary_and_get_pivot(ripser &ripser,
                                             const index_diameter_t simplex,
                                             const index_t dim,
                                             Column& working_boundary,
                                             entry_hash_map& pivot_column_index) {
	if(dim == 0) {
		return index_diameter_t(-1, -1);
	}
	simplex_boundary_enumerator facets(ripser);
	facets.set_simplex(simplex, dim);
	std::vector<index_diameter_t> working_boundary_buffer;
	bool check_for_emergent_pair = true;
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		// Threshold check
		if(get_diameter(facet) <= ripser.threshold) {
			working_boundary_buffer.push_back(facet);
			// Emergent pair: Zero persistence check
			if(check_for_emergent_pair &&
			   get_diameter(simplex) == get_diameter(facet)) {
			    // Emergent pair: Pivot check
				// TODO: Why the apparent cofacet check?
				if((pivot_column_index.find(get_index(facet)) == pivot_column_index.end()) &&
				   (get_index(get_zero_apparent_cofacet(ripser, facet, dim - 1)) == -1)) {
					ripser.barcodes.at(dim - 1).emergent_count++;
					return facet;
				}
				check_for_emergent_pair = false;
			}
		}
	}
	for(index_diameter_t facet : working_boundary_buffer) {
		working_boundary.push(facet);
	}
	return get_pivot(working_boundary);
}

// dim - 1 = dimension of simplices
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
			// Threshold check
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
				// Apparent pair check
				if(!is_in_zero_apparent_pair(ripser, cofacet, dim)) {
					columns_to_reduce.push_back(cofacet);
				} else {
					ripser.barcodes.at(dim - 1).apparent_count++;
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), filtration_order);
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   std::unordered_map<index_t, std::vector<index_t>>& prev_zero_column_index,
                   const index_t dim) {
	compressed_sparse_matrix reduction_matrix; // V
	std::unordered_map<index_t, std::vector<index_t>> zero_column_index;
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
		                                                     working_boundary,
		                                                     pivot_column_index);
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
				index_diameter_t e = get_zero_apparent_cofacet(ripser, pivot, dim - 1);
				if(get_index(e) != -1) {
					// TODO: Why is the necessary w.r.t. apparent pairs?
					add_simplex_boundary(ripser,
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
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			// Non-essential death index (birth in dimension dim - 1)
			value_t birth = get_diameter(pivot);
			value_t death = get_diameter(column_to_reduce);
			if(death > birth * ripser.ratio) {
				std::vector<index_t> rep;
				while(!working_boundary.empty()) {
					rep.push_back(get_index(pop_pivot(working_boundary)));
				}
				ripser.add_hom_class(dim - 1, birth, death, rep);
			}
		} else {
			// Zero column
			std::vector<index_t> rep;
			for(index_t i = reduction_matrix.column_start(j);
			    i < reduction_matrix.column_end(j);
			    ++i)
			{
				rep.push_back(get_index(reduction_matrix.get_entry(i)));
			}
			// Add the diagonal-entry of V that is omitted in reduction_matrix
			rep.push_back(get_index(column_to_reduce));
			zero_column_index.insert({get_index(column_to_reduce), rep});
		}
	}
	// Determine essential pairs by iterating all candidates (zero columns
	// of previous dimension). If such a candidate is not a pivot, we have
	// and essential index
	for(auto it = prev_zero_column_index.begin();
	    it != prev_zero_column_index.end();
	    ++it)
	{
		if(pivot_column_index.find(it->first) == pivot_column_index.end()) {
			// Essential index in dim - 1
			value_t birth = ripser.compute_diameter(it->first, dim - 1);
			ripser.add_hom_class(dim - 1, birth, INF, it->second);
		}
	}
	prev_zero_column_index.swap(zero_column_index);
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	std::unordered_map<index_t, std::vector<index_t>> prev_zero_column_index;
	// List all vertices
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
		std::vector<index_t> rep{ i };
		prev_zero_column_index.insert({i, rep});
	}
	index_t last_dim = std::min(ripser.dim_threshold + 1, ripser.dim_max);
	for(index_t dim = 1; dim <= last_dim; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		assemble_columns_to_reduce(ripser, simplices, columns_to_reduce, dim);
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser,
		              columns_to_reduce,
		              pivot_column_index,
		              prev_zero_column_index,
		              dim);
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
	index_t dim_max = 2;
	index_t dim_threshold = 1;
	float ratio = 1;
	ripser ripser(std::move(dist), dim_max, dim_threshold, enclosing_radius, ratio);
	list_all_simplices(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser, true);
	//simplex_boundary_enumerator b(ripser);
	//b.set_simplex({4, 0.6}, 1);
	//while(b.has_next()) {
	//	std::cout << get_index(b.next()) << " ";
	//}
	//std::cout << std::endl;
	//simplex_coboundary_enumerator cb(ripser);
	//cb.set_simplex({3, -INF}, 0);
	//while(cb.has_next()) {
	//	std::cout << get_index(cb.next()) << " ";
	//}
	//std::cout << std::endl;
	exit(0);
}
