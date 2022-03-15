#include <iostream>
#include <algorithm>

#include "ripser_core.hpp"
#include "print_utils.hpp"


typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            filtration_order_comp> Column;

index_diameter_t init_boundary_and_get_pivot(ripser &ripser,
                                             const index_diameter_t simplex,
                                             const index_t dim,
                                             Column& working_boundary)
{
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
	return get_pivot(ripser, working_boundary);
}

void assemble_columns_to_reduce(ripser &ripser,
                                std::vector<index_diameter_t>& simplices,
                                const index_t dim)
{
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
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(simplices.begin(), simplices.end(), filtration_order);
	ripser.infos.at(dim).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(dim).simplex_total_count = simplices.size();
	ripser.infos.at(dim).simplex_reduction_count = simplices.size();
}

// Used when a death index is encountered, adding it to the corresponding
// homology class
void update_hom_class(ripser& ripser,
                      index_t dim,
                      index_diameter_t birth,
                      index_diameter_t death,
                      std::vector<index_diameter_t> R_rep)
{
	std::vector<homology_class>& hc = ripser.hom_classes.at(dim - 1);
	for(index_t i = 0; i < (index_t) hc.size(); i++) {
		homology_class& h = hc.at(i);
		if(get_index(h.birth.second) == get_index(birth)) {
			// Ratio check
			if(get_diameter(death) >
			   std::max(0.0f, get_diameter(birth)) * ripser.config.ratio) {
				h.death.second = death;
			} else {
				ripser.infos.at(dim - 1).class_count--;
				ripser.infos.at(dim - 1).zero_pers_count++;
				hc.erase(hc.begin() + i);
			}
			break;
		}
	}
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   const index_t dim)
{
	time_point reduction_start = get_time();
	compressed_sparse_matrix V;
	entry_hash_map pivot_column_index;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) {
		ripser.add_reduction_record(dim, j, get_time());
		V.append_column();
		Column R_j;
		Column V_j;
		index_diameter_t sigma_j = columns_to_reduce.at(j);
		index_diameter_t pivot = init_boundary_and_get_pivot(ripser,
		                                                     sigma_j,
		                                                     dim,
		                                                     R_j);
		// The reduction
		index_t add_count = 0;
		while(get_index(pivot) != -1) {
			auto pair = pivot_column_index.find(get_index(pivot));
			if(pair != pivot_column_index.end()) {
				size_t index_column_to_add = pair->second;
				add_boundary(ripser,
				             V,
				             columns_to_reduce,
				             index_column_to_add,
				             dim,
				             V_j,
				             R_j);
				pivot = get_pivot(ripser, R_j);
				add_count++;
				ripser.infos.at(dim).addition_count++;
			} else {
				pivot_column_index.insert({get_index(pivot), j});
				break;
			}
		}
		// Write V_j to V
		std::vector<index_diameter_t> V_rep;
		V_rep.push_back(sigma_j);
		index_diameter_t e = pop_pivot(ripser, V_j);
		while(get_index(e) != -1) {
			V.push_back(e);
			V_rep.push_back(e);
			e = pop_pivot(ripser, V_j);
		}
		// Update barcode decomp
		if(get_index(pivot) != -1) {
			ripser.get_current_reduction_record().to_zero = false;
			std::vector<index_diameter_t> R_rep;
			e = pop_pivot(ripser, R_j);
			while(get_index(e) != -1) {
				R_rep.push_back(e);
				e = pop_pivot(ripser, R_j);
			}
			update_hom_class(ripser, dim, pivot, sigma_j, R_rep);
		} else if(dim == ripser.n - 1 || dim < ripser.config.dim_max) {
			ripser.add_hom_class(dim, sigma_j, index_diameter_t(-1, INF), V_rep);
			ripser.infos.at(dim).class_count++;
		}
		ripser.complete_reduction_record(get_time(), add_count, 0, -1);
	}
	ripser.infos.at(dim).reduction_dur = get_duration(reduction_start, get_time());
}

void compute_barcodes(ripser& ripser)
{
	std::vector<index_diameter_t> simplices;
	time_point assemble_start = get_time();
	assemble_all_simplices(ripser, simplices, 0);
	std::sort(simplices.begin(), simplices.end(), filtration_order);
	ripser.infos.at(0).assemble_dur = get_duration(assemble_start, get_time());
	ripser.infos.at(0).simplex_total_count = ripser.n;
	ripser.infos.at(0).simplex_reduction_count = ripser.n;
	index_t dim_max = std::min(ripser.config.dim_max, (int) ripser.n - 1);
	for(index_t dim = 0; dim <= dim_max; dim++) {
		if(dim > 0) {
			assemble_columns_to_reduce(ripser, simplices, dim);
		}
		compute_pairs(ripser, simplices, dim);
	}
}


/* **************************************************************************
 * Main
 * *************************************************************************/

int main(int argc, char** argv)
{
	ripser_config config;
	if(argc == 2) {
		config = read_config(argv[1]);
	} else {
		std::cerr << "error: missing config path" << std::endl;
		exit(-1);
	}
	ripser ripser(config);
	output_config(ripser, std::cout); std::cout << std::endl;
	compute_barcodes(ripser);
	std::cout << std::endl << std::endl;
	output_barcode(ripser, std::cout, false); std::cout << std::endl;
	output_info(ripser, std::cout); std::cout << std::endl;
	exit(0);
}
