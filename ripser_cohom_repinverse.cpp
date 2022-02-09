#include <iostream>

#include "ripser_core.hpp"
#include "print_utils.hpp"


typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            reverse_filtration_order_comp> Column;

index_diameter_t init_coboundary_and_get_pivot(ripser &ripser,
                                               const index_diameter_t simplex,
                                               const index_t dim,
                                               Column& working_coboundary,
                                               entry_hash_map& pivot_column_index) {
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	bool check_for_emergent_pair = true;
	std::vector<index_diameter_t> working_coboundary_buffer;
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		// Threshold check
		if(get_diameter(cofacet) <= ripser.threshold) {
			working_coboundary_buffer.push_back(cofacet);
			// Emergent pair candidate check
			if(check_for_emergent_pair &&
			   (get_diameter(simplex) == get_diameter(cofacet)) &&
			   (pivot_column_index.find(get_index(cofacet)) == pivot_column_index.end()) &&
			   (get_index(get_zero_apparent_facet(ripser, cofacet, dim + 1)) == -1))
			{
				ripser.infos.at(dim).emergent_count++;
				return cofacet;
			}
			check_for_emergent_pair = false;
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
                                std::vector<index_diameter_t>& cleared_columns,
                                entry_hash_map& pivot_column_index,
                                const index_t dim) {
	time_point assemble_start = get_time();
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
						ripser.infos.at(dim).apparent_count++;
					}
				} else {
					cleared_columns.push_back(cofacet);
					ripser.infos.at(dim).clearing_count++;
				}
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              reverse_filtration_order);
	std::sort(cleared_columns.begin(), cleared_columns.end(),
	          reverse_filtration_order);
	ripser.infos.at(dim).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = columns_to_reduce.size();
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   compressed_sparse_matrix& V,
                   const index_t dim) {
	time_point reduction_start = get_time();
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		//ripser.add_reduction_record(dim, j, get_time());
		V.append_column();
		Column V_j;
		Column R_j;
		// Assemble
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
		                                                       sigma_j,
		                                                       dim,
		                                                       R_j,
		                                                       pivot_column_index);
		// The reduction
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_coboundary(ripser,
				               V,
				               columns_to_reduce,
				               index_column_to_add,
				               dim,
				               V_j,
				               R_j);
				ripser.infos.at(dim).addition_count++;
				pivot = get_pivot(R_j);
			} else {
				index_diameter_t e = get_zero_apparent_facet(ripser, pivot, dim + 1);
				if(get_index(e) != -1) {
					add_simplex_coboundary(ripser,
					                       e,
					                       dim,
					                       V_j,
					                       R_j);
					pivot = get_pivot(R_j);
				} else {
					pivot_column_index.insert({get_index(pivot), j});
					break;
				}
			}
		}
		// Write V_j to V
		index_diameter_t e = pop_pivot(V_j);
		while(get_index(e) != -1) {
			V.push_back(e);
			e = pop_pivot(V_j);
		}
		//ripser.complete_reduction_record(dim, get_time(), add_count, app_count,
		// R_j.size());
		// Determine Persistence Pair
		if(get_index(pivot) != -1) {
			value_t death = std::max(0.0f, get_diameter(sigma_j));
			value_t birth = get_diameter(pivot);
			if(birth > death * ripser.config.ratio) {
				ripser.add_hom_class(dim, sigma_j, pivot);
			}
		} else {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF));
		}
	}
	ripser.infos.at(dim).reduction_dur = get_duration(reduction_start, get_time());
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
		if(get_diameter(cofacet) <= ripser.threshold &&
		   !reverse_filtration_order(cofacet, min_simplex)
		   )
		{
			coboundary.push(cofacet);
		}
	}
}

void compute_reps(ripser& ripser,
                  std::vector<index_diameter_t>& columns_to_reduce_dimm1,
                  std::vector<index_diameter_t>& columns_to_reduce,
                  std::vector<index_diameter_t>& cleared_columns,
                  compressed_sparse_matrix& V_dimm1,
                  compressed_sparse_matrix& V_dim,
                  entry_hash_map& pivot_column_index_dimm1,
                  index_t dim) {
	time_point rep_start = get_time();
	std::vector<homology_class>& hom_classes = ripser.hom_classes.at(dim);
	std::vector<index_t> active_hom_classes;
	std::vector<index_t> hom_parities(hom_classes.size(), 0);
	std::vector<index_t> hom_offs(hom_classes.size(), 0);
	index_t hom_class_index = 0;
	std::sort(hom_classes.begin(), hom_classes.end(), homology_class_order);
	index_t n = columns_to_reduce.size() + cleared_columns.size() - 1;
	index_t cc_index = 0; // cleared_columns
	index_t ctr_index = 0; // columns_to_reduce
	index_diameter_t current_simplex(-1, -1); // current iteration simplex
	index_diameter_t min_simplex(-1, -1); // minimum considered simplex
	Column Vt_row;
	while(ctr_index < (index_t) columns_to_reduce.size() ||
	      cc_index < (index_t) cleared_columns.size())
	{
		//std::cout << "  " << cc_index + ctr_index << "/" << n << std::flush;
		if((cleared_columns.size() == 0) ||
		   ((ctr_index < (index_t) columns_to_reduce.size()) &&
		    (reverse_filtration_order(columns_to_reduce[ctr_index], cleared_columns[cc_index]))))
		{
			// Case V column: Assemble V
			current_simplex = columns_to_reduce.at(ctr_index);
			if(!active_hom_classes.empty()) {
				Vt_row = Column();
				// Add diagonal entry
				Vt_row.push(current_simplex);
				// Add remaining entries if they exist
				//std::cout << " :: cs=" << get_index(current_simplex) << std::flush;
				if(V_dim.column_start(ctr_index) != V_dim.column_end(ctr_index)) {
					for(index_t l = V_dim.column_end(ctr_index) - 1;
						l >= V_dim.column_start(ctr_index);
						l--)
					{
						index_diameter_t v_ele = V_dim.get_entry(l);
						if(!reverse_filtration_order(v_ele, min_simplex)) {
							Vt_row.push(v_ele);
						}
					}
				}
				//std::cout << " (" << ctr_index << ") :: V=" << Vt_row.size() << std::flush;
			} else {
				//std::cout << " :: no hom class" << std::flush;
			}
			ctr_index++;
		} else {
			// Case V column cleared:
			current_simplex = cleared_columns.at(cc_index);
			if(!active_hom_classes.empty()) {
				Vt_row = Column();
				auto pair = pivot_column_index_dimm1.find(get_index(current_simplex));
				if(pair != pivot_column_index_dimm1.end()) {
					index_t pivot_index = pair->second;
					// Add diagonal entry
					add_partial_simplex_coboundary(ripser,
					                               columns_to_reduce_dimm1.at(pivot_index),
					                               min_simplex,
					                               Vt_row,
					                               dim - 1);
					// Add remaining entries
					//std::cout << " :: cs=" << get_index(current_simplex) << std::flush;
					if(V_dimm1.column_start(pivot_index) != V_dimm1.column_end(pivot_index)) {
						for(index_t l = V_dimm1.column_end(pivot_index) - 1;
							l >= V_dimm1.column_start(pivot_index);
							l--)
						{
							index_diameter_t v_ele = V_dimm1.get_entry(l);
							add_partial_simplex_coboundary(ripser, v_ele, min_simplex, Vt_row, dim - 1);
						}
					}
					//std::cout << " (" << cc_index << ") :: R=" << Vt_row.size() << std::flush;
				} else {
				}
			} else {
				//std::cout << " :: no hom class" << std::flush;
			}
			cc_index++;
		}
		// COMPUTE PRODUCT
		index_t num_pops = 0;
		index_diameter_t v_ele = pop_pivot(Vt_row);
		while(get_index(v_ele) != -1) {
			for(index_t k = 0; k < (index_t) active_hom_classes.size(); k++) {
				std::vector<index_diameter_t>& rep = hom_classes[active_hom_classes[k]].representative;
				// Still elements in rep left to work with?
				index_t rep_index = ((index_t) rep.size()) - 1 - hom_offs[k];
				while(rep_index >= 0) {
					index_diameter_t rep_ele = rep[rep_index];
					// check if !(v_ele < rep_ele) equal to (rep_ele >= v_ele)
					if(!reverse_filtration_order(v_ele, rep_ele)) {
						if(rep_ele == v_ele) {
							hom_parities[k] += 1;
						}
						break;
					} else {
						// v_ele < rep_ele
						//hom_offs[k] += 1;
						rep_index--;
					}
				}
			}
			v_ele = pop_pivot(Vt_row);
			num_pops++;
		}
		//std::cout << ":: np=" << num_pops << std::endl;
		for(index_t k = 0; k < (index_t) active_hom_classes.size(); k++) {
			if(hom_parities.at(k) % 2 == 1) {
				hom_classes[active_hom_classes[k]].representative.push_back(current_simplex);
			}
			hom_parities[k] = 0;
			hom_offs[k] = 0;
		}
		if(current_simplex == hom_classes.at(hom_class_index).birth) {
			if(active_hom_classes.empty()) {
				min_simplex = current_simplex;
			}
			active_hom_classes.push_back(hom_class_index);
			hom_classes.at(hom_class_index).representative.push_back(current_simplex);
			if(hom_class_index < ((index_t) hom_classes.size()) - 1) {
				hom_class_index++;
			}
		}
	}
	ripser.infos.at(dim).representative_dur += get_duration(rep_start, get_time());
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> cleared_columns;
	std::vector<index_diameter_t> columns_to_reduce_dim;
	std::vector<index_diameter_t> columns_to_reduce_dimm1;
	entry_hash_map pivot_column_index_dim;
	entry_hash_map pivot_column_index_dimm1;
	compressed_sparse_matrix V_dimm1;
	compressed_sparse_matrix V_dim;
	time_point assemble_start = get_time();
	assemble_all_simplices(ripser, simplices, 0);
	std::sort(simplices.begin(), simplices.end(), reverse_filtration_order);
	ripser.infos.at(0).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(0).simplex_total_count = simplices.size();
	ripser.infos.at(0).simplex_reduction_count = simplices.size();
	index_t dim_max = std::min(ripser.config.dim_max, (int) ripser.n - 1);
	for(index_t dim = 0; dim <= dim_max; ++dim) {
		if(dim == 0) {
			columns_to_reduce_dim = std::vector<index_diameter_t>(simplices);
		} else {
			columns_to_reduce_dim.clear();
			assemble_columns_to_reduce(ripser,
			                           simplices,
			                           columns_to_reduce_dim,
			                           cleared_columns,
			                           pivot_column_index_dim,
			                           dim);
		}
		pivot_column_index_dimm1 = std::move(pivot_column_index_dim);
		pivot_column_index_dim.clear();
		pivot_column_index_dim.reserve(columns_to_reduce_dim.size());
		compute_pairs(ripser,
		              columns_to_reduce_dim,
		              pivot_column_index_dim,
		              V_dim,
		              dim);
		if(columns_to_reduce_dim.size() == 0) {
			continue;
		}
		compute_reps(ripser,
		             columns_to_reduce_dimm1,
		             columns_to_reduce_dim,
		             cleared_columns,
		             V_dimm1,
		             V_dim,
		             pivot_column_index_dimm1,
		             dim);
		// Move the data structures
		V_dimm1.bounds = std::move(V_dim.bounds);
		V_dimm1.entries = std::move(V_dim.entries);
		columns_to_reduce_dimm1 = std::move(columns_to_reduce_dim);
		// Clear the expired data structures
		V_dim.bounds.clear();
		V_dim.entries.clear();
		columns_to_reduce_dim.clear();
		cleared_columns.clear();
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
	std::cout << std::endl;
	output_config(ripser, std::cout); std::cout << std::endl;
	//output_simplices(ripser, std::cout, total_filtration_order); std::cout << std::endl;
	compute_barcodes(ripser);
	output_barcode(ripser, std::cout, true); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	//write_standard_output(ripser, true,itrue);
	exit(0);
}
