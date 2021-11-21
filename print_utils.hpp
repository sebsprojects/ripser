#ifndef PRINT_UTILS
#define PRINT_UTILS

#include "ripser_core.hpp"


std::string simplex_tos(ripser& ripser,
                        index_diameter_t simplex,
                        index_t dim,
                        bool include_index=true,
                        bool include_simplices=true,
                        bool include_diameter=true,
                        bool include_rel=false)
{
	std::stringstream ss;
	if(include_index) {
		ss << get_index(simplex);
	}
	if(include_simplices) {
		if(include_index) {
		  ss << "-";
		}
		std::vector<index_t> vertices(dim+1, -1);
		ripser.get_simplex_vertices(get_index(simplex), dim, ripser.n, vertices);
		for(auto& i : vertices) {
			ss << i;
			if(i != vertices.back()) {
				ss << "'";
			}
		}
	}
	if(include_diameter) {
		if(include_index || include_simplices) {
			ss << "-";
		}
		ss << std::fixed << std::setprecision(3)
		   << get_diameter(simplex);
	}
	if(include_rel) {
		if(ripser.is_relative_simplex(get_index(simplex), dim)) {
			ss << "-r";
		}
	}
	return ss.str();
}

void output_barcode_simplices(ripser& ripser, std::ostream& os,
                              total_simplex_order ord=total_filtration_order,
                              bool pref=false)
{
	std::vector<std::pair<index_t, index_diameter_t>> all_simplices;
	for(auto& hc : ripser.hom_classes) {
		for(auto& h : hc) {
			all_simplices.push_back(std::make_pair(h.dim, h.birth));
			if(get_index(h.death) != -1) {
				all_simplices.push_back(std::make_pair(h.dim + 1, h.death));
			}
		}
	}
	std::sort(all_simplices.begin(), all_simplices.end(), ord);
	for(auto &s : all_simplices) {
		std::string p = pref ? ("#s" + std::to_string(s.first) + "  ") : "";
		os << p
		   << simplex_tos(ripser, s.second, s.first, true, true, true, true)
		   << std::endl;
	}
}

void output_all_simplices(ripser& ripser, std::ostream& os,
                          total_simplex_order ord=total_filtration_order,
                          bool pref=false)
{
	std::vector<index_diameter_t> simpl_prev;
	std::vector<index_diameter_t> simpl_curr;
	std::vector<std::pair<index_t, index_diameter_t>> all_simplices;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simpl_curr.push_back(index_diameter_t(i, diam));
		all_simplices.push_back(std::make_pair(0, index_diameter_t(i, 0)));
	}
	for(index_t i = 0; i <= ripser.config.dim_max; i++) {
		std::swap(simpl_prev, simpl_curr);
		simpl_curr.clear();
		simplex_coboundary_enumerator cofacets(ripser);
		for(auto s : simpl_prev) {
			cofacets.set_simplex(s, i);
			while(cofacets.has_next(false)) {
			  index_diameter_t cofacet = cofacets.next();
			  simpl_curr.push_back(cofacet);
			  all_simplices.push_back(std::make_pair(i + 1, cofacet));
			}
		}
	}
	std::sort(all_simplices.begin(), all_simplices.end(), ord);
	for(auto &s : all_simplices) {
		std::string p = pref ? ("#s" + std::to_string(s.first) + "  ") : "";
		os << p
		   << simplex_tos(ripser, s.second, s.first, true, true, true, true)
		   << std::endl;
	}
}

void output_all_simplices_by_dim(ripser& ripser, std::ostream& os,
                                 simplex_order ord=filtration_order,
                                 index_t dim=-1, bool pref=false)
{
	std::vector<index_diameter_t> simpl_prev;
	std::vector<index_diameter_t> simpl_curr;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simpl_curr.push_back(index_diameter_t(i, diam));
	}
	std::sort(simpl_curr.begin(), simpl_curr.end(), ord);
	for(index_t i = 0; i <= ripser.config.dim_max; i++) {
		std::string p = pref ? ("#s" + std::to_string(i) + " ") : "";
		if(dim == -1 || i == dim) {
			os << (pref ? "# " : "") << "simplices in dim="
			   << std::to_string(i) << std::endl;
			for(auto& s : simpl_curr) {
				os << p << "  "
				   << simplex_tos(ripser, s, i, true, (i != 0), false)
				   << std::endl;
			}
		}
		if(i == dim) {
			break;
		}
		std::swap(simpl_prev, simpl_curr);
		simpl_curr.clear();
		simplex_coboundary_enumerator cofacets(ripser);
		for(auto s : simpl_prev) {
			cofacets.set_simplex(s, i);
			while(cofacets.has_next(false)) {
			  simpl_curr.push_back(cofacets.next());
			}
		}
	}
}

void output_barcode(ripser& ripser, std::ostream& os, bool with_reps=false, index_t dim=-1, bool pref=false)
{
	for(index_t d = 0; d <= (index_t) ripser.config.dim_max + 1; d++) {
		if(!(dim == -1 || d == dim)) {
			continue;
		}
		std::sort(ripser.hom_classes.at(d).begin(), ripser.hom_classes.at(d).end(),
		          homology_class_print_order);
		std::string p = pref ? ("#b" + std::to_string(d) + " ") : "";
		os << (pref ? "# " : "") << "barcode in dim=" << d << std::endl;
		for(auto& hc : ripser.hom_classes.at(d)) {
			value_t birth = std::max(0.0f, get_diameter(hc.birth));
			value_t death = std::max(0.0f, get_diameter(hc.death));
			os << p << "  [" << birth;
			if(death == INF) {
				os << ", )";
			} else {
				os << "," << death << ")";
			}
			os << " ; [" << simplex_tos(ripser, hc.birth, d, 1, 1, 0, 0);
			if(death == INF) {
				os << ", )";
			} else {
				os << "," << simplex_tos(ripser, hc.death, d + 1, 1, 1, 0, 0)
				   << ")";
			}
			if(with_reps) {
				os << " ; ";
				for(index_diameter_t& simp : hc.representative) {
					os << simplex_tos(ripser, simp, d, false, true, false);
					if(simp != hc.representative.back()) {
						os << " :: ";
					}
				}
			}
			os << std::endl;
		}
		if(d + 1 < (index_t) ripser.hom_classes.size()) {
			os << std::endl;
		}
	}
}

void output_info(ripser& ripser, std::ostream& os, index_t dim=-1, bool pref=false)
{
	for(info& info : ripser.infos) {
		std::string p = pref ? ("#i" + std::to_string(info.dim) + " ") : "";
		if(dim == -1 || info.dim == dim) {
			os << (pref ? "# " : "") << "info in dim=" << info.dim << ":"
			   << std::endl;
			os << p << "  total simplex count:     "
			   << info.simplex_total_count << std::endl;
			os << p << "  reduction simplex count: "
			   << info.simplex_reduction_count << std::endl;
			os << p << "  class count:             "
			   << info.class_count << std::endl;
			os << p << "  clearing count:          "
			   << info.clearing_count << std::endl;
			os << p << "  emergent count:          "
			   << info.emergent_count << std::endl;
			os << p << "  apparent count:          "
			   << info.apparent_count << std::endl;
			os << p << "  addition count:          "
			   << info.addition_count << std::endl;
			os << p << "  assemble duration:       "
			   << info.assemble_dur.count() << "s" << std::endl;
			os << p << "  reduction duration:      "
			   << info.reduction_dur.count() << "s" << std::endl;
			os << p << "  representative duration: "
			   << info.representative_dur.count() << "s" << std::endl;
			if(&info != &ripser.infos.back()) {
				os << std::endl;
			}
		}
	}
}

void output_config(ripser& ripser, std::ostream& os, bool pref=false) {
	std::string p = pref ? "#c " : "";
	ripser_config& config = ripser.config;
	os << (pref ? "# " : "") << "config" << std::endl;
	os << p << "  dataset at file path: " << config.file_path << std::endl;
	os << p << "  output at file path:  " << config.output_path << std::endl;
	os << p << "  dataset input type:   " << config.input_type << std::endl;
	os << p << "  dim max:              " << config.dim_max << std::endl;
	os << p << "  ratio:                " << config.ratio << std::endl;
	os << p << "  threshold:            " << ripser.threshold;
	if(config.use_enclosing_threshold) {
		os << " (enclosing threshold)" << std::endl;
	} else if(config.config_threshold < 0) {
		os << " (maximum distance)" << std::endl;
	} else {
		os << std::endl;
	}
	os << p << "  use union find:       "
	   << config.use_union_find << std::endl;
	os << p << "  relative:             ";
	if(config.relative_subcomplex.size() == 0) {
		os << "(empty)";
	} else {
		auto pair = config.relative_subcomplex.at(0);
		os << pair.first << "-" << pair.second;
		for(index_t i = 1; i < (index_t) config.relative_subcomplex.size(); i++) {
			pair = config.relative_subcomplex.at(i);
			os << "," << pair.first << "-" << pair.second;
		}
	}
	os << std::endl;
}

std::ofstream get_writeout_stream(ripser& ripser) {
	std::string ip = ripser.config.file_path;
	std::string dataset_name = ip.substr(ip.find_last_of("/\\") + 1);
	std::string rel_name = "rel";
	for(index_t i = 0; i < (index_t) ripser.config.relative_subcomplex.size(); i++) {
		rel_name += std::to_string(ripser.config.relative_subcomplex.at(i).first);
		rel_name += "-" + std::to_string(ripser.config.relative_subcomplex.at(i).second);
		if(i < (index_t) ripser.config.relative_subcomplex.size() - 1) {
			rel_name += "x";
		}
	}
	if(ripser.config.relative_subcomplex.empty()) {
		rel_name = "relnone";
	}
	std::string thresh_name = "thresh" + std::to_string(ripser.threshold);
	thresh_name.erase(thresh_name.find_last_not_of('0') + 1, std::string::npos);
	std::replace(thresh_name.begin(), thresh_name.end(), '.', 'd');
	std::string filename = ripser.config.output_path + "/"
		+ dataset_name
		+ "_" + thresh_name + "_"
		+ rel_name
		+ ".txt";
	return std::ofstream(filename, std::ofstream::trunc);
}

void write_standard_output(ripser& ripser,
                           bool with_reps=false,
                           bool with_simplex_list=false,
                           total_simplex_order ord=total_filtration_order) {
	std::ofstream ofs = get_writeout_stream(ripser);
	output_config(ripser, ofs, true); ofs << std::endl;
	if(with_simplex_list) {
		output_barcode_simplices(ripser, ofs, ord, true);
		ofs << std::endl;
	}
	output_info(ripser, ofs, -1, true); ofs << std::endl;
	output_barcode(ripser, ofs, with_reps, -1, true); ofs << std::endl;
}

#endif
