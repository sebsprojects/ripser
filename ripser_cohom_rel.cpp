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
				columns_to_reduce.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const entry_hash_map& previous_pivots,
                   const index_t dim) {
	compressed_sparse_matrix reduction_matrix; // V
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
		// Determine Persistence Pair
		value_t birth = get_diameter(column_to_reduce);
		if(get_index(pivot) != -1) {
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.ratio) {
				// Non-essential pair
				ripser.barcodes.at(dim).add_interval(birth, death);
			}
		} else {
			// Zero column
			auto pair = previous_pivots.find(get_index(column_to_reduce));
			if(pair == previous_pivots.end()) {
				// Essential index!
				ripser.barcodes.at(dim).add_interval(birth, INF);
			} else {
				// Killing index non-essential pair
			}
		}
	}
}

void compute_barcodes(ripser& ripser) {
	ripser.barcodes.clear();
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	entry_hash_map previous_pivots;
	// The relative part of the complex is induced by vertices [0, rel_count)
	// Here we just the complement of that set
	for(index_t i = ripser.rel_count; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	for(index_t dim = 0; dim <= ripser.dim_max; ++dim) {
		ripser.barcodes.push_back(barcode(dim));
		std::vector<index_diameter_t> columns_to_reduce;
		if(dim == 0) {
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
	const char* file_all_points = nullptr;
	const char* file_rel_points = nullptr;
	if(argc == 3) {
		file_all_points = argv[1];
		file_rel_points = argv[2];
	} else {
		std::cerr << "error: specify two point-cloud files for relative cohomology"
				  << std::endl;
		exit(-1);
	}
	// Reading the point cloud from file
	std::ifstream file_stream_all(file_all_points);
	std::ifstream file_stream_rel(file_rel_points);
	if ((file_all_points && file_stream_all.fail()) ||
	    (file_rel_points && file_stream_rel.fail())) {
		std::cerr << "error: couldn't open file(s)" << std::endl;
		exit(-1);
	}
	// Read the point clouds
	std::vector<std::vector<value_t>> all_points = read_point_cloud(file_stream_all);
	std::vector<std::vector<value_t>> rel_points = read_point_cloud(file_stream_rel);
	std::vector<std::vector<value_t>> points = rel_points;
	points.reserve(all_points.size());
	// Push all non-rel points of all_points into the new vector points
	for(size_t i = 0; i < all_points.size(); i++) {
		if(std::find(rel_points.begin(), rel_points.end(), all_points.at(i)) == rel_points.end()) {
			points.push_back(all_points.at(i));
		}
	}
	for(auto v : points) {
		std::cout << v.at(0) << "," << v.at(1) << std::endl;
	}
	DistanceMatrix dist = compressed_lower_distance_matrix(point_cloud_to_distance_vector(points));
	value_t enclosing_radius = compute_enclosing_radius(dist);
	std::cout << "Enclosing radius: " << enclosing_radius << std::endl;
	index_t dim_max = 1;
	float ratio = 1;
	ripser ripser(std::move(dist), rel_points.size(), dim_max, enclosing_radius, ratio);
	list_all_simplices(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser.barcodes);
	exit(0);
}
