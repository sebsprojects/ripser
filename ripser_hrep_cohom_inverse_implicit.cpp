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

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
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
		ripser.complete_reduction_record(dim, get_time(), add_count, app_count, working_coboundary.size());
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

void add_partial_simplex_coboundary(ripser& ripser,
                                    index_diameter_t simplex,
                                    index_diameter_t min_simplex,
                                    Column& coboundary,
                                    index_t dim)
{
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		if(get_diameter(cofacet) <= ripser.threshold
		   //&& !reverse_filtration_order(cofacet, min_simplex)
		   )
		{
			coboundary.push(cofacet);
		}
	}
}

void compute_reps(ripser& ripser,
                  std::vector<index_diameter_t>& simplices,
                  std::vector<index_diameter_t>& ctrdm1,
                  std::vector<index_diameter_t>& ctrd,
                  compressed_sparse_matrix& Vdm1,
                  compressed_sparse_matrix& Vd,
	              entry_hash_map& pcidm1,
	              index_t dim) {
	time_point rep_start = get_time();
	std::vector<homology_class>& hom_classes = ripser.barcodes.at(dim).hom_classes;
	std::vector<index_t> active_hom_classes;
	index_t hom_class_index = 0;
	index_diameter_t min_simplex;
	index_t ctr_index = 0;
	std::sort(hom_classes.begin(), hom_classes.end(), homology_class_order);
	for(index_t j = 0; j < simplices.size(); j++) {
		if(!active_hom_classes.empty()) {
			Column Vt_row;
			if(ctrd.at(ctr_index) == simplices.at(j)) { // Use V column
				// Add diagonal entry
				Vt_row.push(simplices.at(j));
				// Add remaining entries if they exist
				if(Vd.column_start(ctr_index) != Vd.column_end(ctr_index)) {
					for(index_t l = Vd.column_end(ctr_index) - 1;
						l >= Vd.column_start(ctr_index);
						l--)
					{
						index_diameter_t v_ele = Vd.get_entry(l);
						// TODO: We do not need all elements from Vd
						// TODO: We know that Vd is sorted, maybe push can be faster
						//if(!reverse_filtration_order(v_ele, min_simplex)) {
							Vt_row.push(v_ele);
						//}
					}
				}
			} else {
				// V column missing, check if it was cleared, if not it is
				// an apparent pair index and can be ignored
				auto pair = pcidm1.find(get_index(simplices.at(j)));
				if(pair != pcidm1.end()) {
					index_t vdm1_index = pair->second;
					// Add diagonal entry
					add_partial_simplex_coboundary(ripser, ctrdm1.at(vdm1_index), min_simplex, Vt_row, dim - 1);
					// Add remaining entries
					if(Vdm1.column_start(vdm1_index) == Vdm1.column_end(vdm1_index)) {
						for(index_t l = Vdm1.column_end(vdm1_index) - 1;
							l >= Vdm1.column_start(vdm1_index);
							l--) {
							index_diameter_t v_ele = Vd.get_entry(l);
							add_partial_simplex_coboundary(ripser, v_ele, min_simplex, Vt_row, dim - 1);
						}
					}
				}
			}
			for(index_t active_index : active_hom_classes) {
				// Compute the dot-product (representative * j_row)
				homology_class& hom_class = hom_classes.at(active_index);
				index_t parity = 0;
				for(index_t k = ((index_t) hom_class.representative.size()) - 1; k >= 0; k--) {
					Column Vt_row_copy = Vt_row;
					index_diameter_t rep_ele = hom_class.representative.at(k);
					index_diameter_t v_ele = pop_pivot(Vt_row_copy);
					while(get_index(v_ele) != -1) {
						//if(!reverse_filtration_order(rep_ele, v_ele)) {
							if(rep_ele == v_ele) {
								parity++;
						//	}
						//	break;
						}
						v_ele = pop_pivot(Vt_row_copy);
					}
				}
				if(parity % 2 == 1) {
					hom_class.representative.push_back(simplices.at(j));
				}
			}
		}
		if(ctr_index < ((index_t) ctrd.size()) - 1 && ctrd.at(ctr_index) == simplices.at(j)) {
			ctr_index++;
		}
		if(simplices.at(j) == hom_classes.at(hom_class_index).birth) {
			if(active_hom_classes.empty()) {
				min_simplex = simplices.at(j);
			}
			active_hom_classes.push_back(hom_class_index);
			hom_classes.at(hom_class_index).representative.push_back(simplices.at(j));
			if(hom_class_index < ((index_t) hom_classes.size()) - 1) {
			  hom_class_index++;
			}
		}
	}
	ripser.infos.at(dim).representative_dur += get_duration(rep_start, get_time());
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> ctrdm1;
	entry_hash_map pivot_column_index;
	entry_hash_map pcidm1;
	compressed_sparse_matrix Vdm1;
	compressed_sparse_matrix Vd;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	index_t last_dim = std::min(ripser.dim_threshold, ripser.dim_max);
	for(index_t dim = 0; dim <= last_dim; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
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
		pcidm1 = std::move(pivot_column_index);
		pivot_column_index.clear();
		pivot_column_index.reserve(columns_to_reduce.size());
		compute_pairs(ripser,
		              columns_to_reduce,
		              pivot_column_index,
		              Vd,
		              dim);
		compute_reps(ripser,
		             simplices,
		             ctrdm1,
		             columns_to_reduce,
		             Vdm1,
		             Vd,
		             pcidm1,
		             dim);
		// Move the data structures
		Vdm1.bounds = std::move(Vd.bounds);
		Vdm1.entries = std::move(Vdm1.entries);
		ctrdm1 = std::move(columns_to_reduce);
		// Clear the expired data structures
		Vd.bounds.clear();
		Vd.entries.clear();
		columns_to_reduce.clear();
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
