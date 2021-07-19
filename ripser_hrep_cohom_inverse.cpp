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
				if(pivot_column_index.find(get_index(cofacet)) ==
				   pivot_column_index.end()) {
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
			next_simplices.push_back(cofacet);
			if(get_diameter(cofacet) <= ripser.threshold) {
				columns_to_reduce.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
	info.assemble_dur = get_duration(assemble_start, get_time());
	info.simplex_total_count = simplices.size();
	info.simplex_reduction_count = columns_to_reduce.size();
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const entry_hash_map& previous_pivots,
                   const index_t dim) {
	info& info = ripser.infos.at(dim);
	compressed_sparse_matrix reduction_matrix; // V
	compressed_sparse_matrix v_inv;
	std::vector<std::pair<index_t, index_t>> nonessential_red;
	for(index_t j = 0; j < (index_t) columns_to_reduce.size(); ++j) { // For j in J
		time_point reduction_start = get_time();
		Column working_reduction_column; // V_j
		Column working_coboundary;       // R_j
		// Assemble the column j (corresponding to the simplex column_to_reduce)
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
		                                                       column_to_reduce,
		                                                       dim,
		                                                       working_coboundary,
		                                                       pivot_column_index);
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
		// Write V_j to V. Does NOT contain duplicate entires anymore
		// (V is not lazy as working_reduction_column is)
		reduction_matrix.append_column();
		index_diameter_t e = pop_pivot(working_reduction_column);
		while(get_index(e) != -1) {
			reduction_matrix.push_back(e);
			e = pop_pivot(working_reduction_column);
		}
		info.reduction_dur += get_duration(reduction_start, get_time());
		time_point rep_start = get_time();
		// Compute the new inverse of V and store it column-order
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
		info.representative_dur += get_duration(rep_start, get_time());
		//print_mat_simplices(ripser, v_inv, dim);
		// Determine Persistence Pair
		value_t birth = get_diameter(column_to_reduce);
		if(get_index(pivot) != -1) {
			value_t death = get_diameter(pivot);
			if(death > birth * ripser.ratio) {
				// Non-essential pair
				ripser.add_hom_class(dim, birth, death, std::vector<index_t>());
				// Store which non-essential hom class corresponds to which column
				nonessential_red.push_back({ripser.barcodes.at(dim).hom_classes.size() - 1, j});
			}
		} else {
			// Zero column
			auto pair = previous_pivots.find(get_index(column_to_reduce));
			if(pair == previous_pivots.end()) {
				// Essential index!
				ripser.add_hom_class(dim, birth, INF, std::vector<index_t>());
			} else {
				// Killing index non-essential pair
			}
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
	//std::cout << std::endl << "V and V^-1: dim=" << dim << std::endl;
	//print_v(reduction_matrix, columns_to_reduce);
	//std::cout << std::endl;
	//print_mat_simplices(ripser, v_inv, dim);
	//print_v(reduction_matrix, columns_to_reduce);
	//std::cout << std::endl;
	//print_v(v_inv, columns_to_reduce);
}

void compute_barcodes(ripser& ripser) {
	ripser.barcodes.clear();
	std::vector<index_diameter_t> simplices;
	entry_hash_map pivot_column_index;
	entry_hash_map previous_pivots;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	index_t last_dim = std::min(ripser.dim_threshold, ripser.dim_max);
	for(index_t dim = 0; dim <= last_dim; ++dim) {
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
    std::string filename = "";
	if(argc == 2) {
		filename = argv[1];
	} else {
		std::cerr << "error: specify path to lower-distance matrix file as only arg"
				  << std::endl;
		exit(-1);
	}
	// Reading the distance matrix from file
	std::ifstream file_stream(filename);
	if (filename != "" && file_stream.fail()) {
		std::cerr << "error: couldn't open file " << filename << std::endl;
		exit(-1);
	}
	DistanceMatrix dist = read_lower_distance_matrix(file_stream);
	value_t enclosing_radius = compute_enclosing_radius(dist);
	index_t dim_max = 1;
	index_t dim_threshold = 1;
	float ratio = 1;
	ripser ripser(std::move(dist), dim_max, dim_threshold, enclosing_radius, ratio);
	//list_all_simplices(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser); std::cout << "\n\n";
	print_infos(ripser);
	std::string fn_pre = "./output/cycle_rep_dim1/";
	std::string fn = filename.substr(filename.find_last_of("/") + 1);
	fn = fn.substr(0, fn.find_last_of("."));
	std::string fn_post = "_cohom_inverse.txt";
	//write_dim1_cycles(ripser, fn_pre + fn + fn_post);
	exit(0);
}
