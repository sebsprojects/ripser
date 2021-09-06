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
	// TODO: Is this sort necessary?
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	info.assemble_dur = get_duration(assemble_start, get_time());
	info.simplex_total_count = simplices.size();
	info.simplex_reduction_count = columns_to_reduce.size();
}

typedef std::unordered_map<index_t,
                           std::vector<index_diameter_t>> column_hash_map;

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   column_hash_map& reduced_cols,
                   compressed_sparse_matrix& reduction_matrix,
                   const index_t dim) {
	info& info = ripser.infos.at(dim);
	time_point reduction_start = get_time();
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
		// Write R_j to reduced_cols
		index_t red_eles = -1;
		if(get_index(pivot) != -1 &&
		   dim < std::min(ripser.dim_threshold, ripser.dim_max)) {
			e = index_diameter_t(pivot);
			std::vector<index_diameter_t> r_j;
			// Additional pop pivot to remove 'diagonal' entry
			e = pop_pivot(working_coboundary);
			while(get_index(e) != -1) {
				r_j.push_back(e);
				e = pop_pivot(working_coboundary);
			}
			reduced_cols.insert({get_index(pivot), r_j});
			red_eles = r_j.size();
		}
		ripser.complete_reduction_record(dim, get_time(), add_count, app_count, red_eles);
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			// Non-essential index
			ripser.add_hom_class(dim, column_to_reduce, pivot);
		} else {
			// Essential index (since clearing)
			ripser.add_hom_class(dim, column_to_reduce, index_diameter_t(-1, INF));
		}
	}
	info.reduction_dur = get_duration(reduction_start, get_time());
}

void compute_reps(ripser& ripser,
                  std::vector<index_diameter_t>& simplices,
                  std::vector<index_diameter_t>& columns_to_reduce,
                  compressed_sparse_matrix& reduction_matrix,
                  column_hash_map& reduced_columns,
                  index_t dim) {
	// Compute a lookup tabel that maps the indices of simplices to either
	//   > the column index in reduction matrix
	//   > -1, if the column was either cleared or apprent (we do not know yet)
	time_point rep_start = get_time();
	index_t ind = 0;
	std::vector<index_t> index_lookup;
	for(index_t i = 0; i < (index_t) simplices.size(); ++i) {
		if(ind < (index_t) columns_to_reduce.size() && simplices.at(i) == columns_to_reduce.at(ind)) {
			index_lookup.push_back(ind);
			ind++;
		} else {
			index_lookup.push_back(-1);
		}
	}
	for(index_t i = 0; i < (index_t) ripser.barcodes.at(dim).hom_classes.size(); ++i) {
		index_diameter_t birth = ripser.barcodes.at(dim).hom_classes.at(i).birth;
		// Get the birth simplex index in simplices
		auto it = std::lower_bound(simplices.begin(), simplices.end(), birth,
		                           reverse_filtration_order);
		index_t column_index = std::distance(simplices.begin(), it);
		std::vector<index_diameter_t>& rep = ripser.barcodes.at(dim).hom_classes.at(i).representative;
		rep.push_back(birth); // j = column_index
		//std::cout << column_index << ": ";
		for(index_t j = column_index + 1; j < (index_t) simplices.size(); j++) {
			index_t parity = 0;
			if(index_lookup.at(j) == -1) {
				auto rcit = reduced_columns.find(get_index(simplices.at(j)));
				if(rcit == reduced_columns.end()) {
					// Apparent, does not contribute
				} else {
					// Clearing
					auto lb = rcit->second.begin();
					for(index_t k = 0; k < (index_t) rep.size(); k++) {
						// TODO: r_j could be empty
						index_t rep_index = rep.size() - k - 1;
						lb = std::lower_bound(lb,
						                      rcit->second.end(),
						                      rep.at(rep_index),
						                      filtration_order);
						//std::cout << "r: " << get_index(rep.at(k)) << ", "
						//            << get_index(simplices.at(j)) << std::endl;
						if(lb != rcit->second.end() && *lb == rep.at(rep_index)) {
							parity += 1;
							std::advance(lb, 1);
						}
						// Linear time find:
						//parity += (std::find(rcit->second.begin(), rcit->second.end(), rep.at(k)) != rcit->second.end());
					}
					if(parity % 2 == 1) {
						rep.push_back(simplices.at(j));
					}
				}
			} else {
				// Reduction Matrix
				index_t reduction_column_index = index_lookup.at(j);
				if(reduction_matrix.column_start(reduction_column_index) ==
				   reduction_matrix.column_end(reduction_column_index)) {
					// Column with only diagonal entry
					continue;
				}
				auto lb = reduction_matrix.entries.begin() +
				          reduction_matrix.column_start(reduction_column_index);
				auto ent_end = reduction_matrix.entries.begin() +
				               reduction_matrix.column_end(reduction_column_index);
				for(index_t k = 0; k < (index_t) rep.size(); k++) {
					index_t rep_index = rep.size() - k - 1;
					lb = std::lower_bound(lb,
					                      ent_end,
					                      rep.at(rep_index),
					                      filtration_order);
					if(lb != ent_end && *lb == rep.at(rep_index)) {
						parity += 1;
						std::advance(lb, 1);
					}
					// Linear time find:
					//parity += reduction_matrix.search_column(reduction_column_index, rep.at(k));
				}
				if(parity % 2 == 1) {
					rep.push_back(simplices.at(j));
				}
			}
		}
	}
	ripser.infos.at(dim).representative_dur += get_duration(rep_start, get_time());
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
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	index_t last_dim = std::min(ripser.dim_threshold, ripser.dim_max);
	for(index_t dim = 0; dim <= last_dim; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		compressed_sparse_matrix reduction_matrix;
		if(dim == 0) {
			columns_to_reduce = std::vector<index_diameter_t>(simplices);
		} else {
			// Takes the pivots from the previous iteration
			assemble_columns_to_reduce(ripser,
			                           simplices,
			                           columns_to_reduce,
			                           pivot_column_index,
			                           dim);
		}
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		red_cols.clear();
		red_cols.reserve(columns_to_reduce.size());
		compute_pairs(ripser,
		              columns_to_reduce,
		              pivot_column_index,
		              red_cols,
		              reduction_matrix,
		              dim);
		compute_reps(ripser,
		             simplices,
		             columns_to_reduce,
		             reduction_matrix,
		             prev_red_cols,
		             dim);
		prev_red_cols.swap(red_cols);
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
	DistanceMatrix dist = read_distance_matrix(file_stream);
	value_t enclosing_radius = compute_enclosing_radius(dist);
	index_t dim_max = 2;
	index_t dim_threshold = 1;
	float ratio = 1;
	ripser ripser(std::move(dist), dim_max, dim_threshold, enclosing_radius, ratio);
	//list_all_simplices(ripser);
	compute_barcodes(ripser);
	print_barcodes(ripser); std::cout << "\n\n";
	print_infos(ripser);
	bool output_cycles = true;
	if(output_cycles) {
		std::string fn_pre = "./output/cycle_rep_dim1/";
		std::string fn = filename.substr(filename.find_last_of("/") + 1);
		fn = fn.substr(0, fn.find_last_of("."));
		std::string fn_post = "_cohom_inv_forward.txt";
		write_dim1_cycles(ripser, fn_pre + fn + fn_post);
	}
	exit(0);
}
