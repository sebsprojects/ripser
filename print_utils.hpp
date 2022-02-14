#ifndef PRINT_UTILS
#define PRINT_UTILS

#include "ripser_core.hpp"


// index-a'b'c-diam-r
std::string simplex_tos(ripser& ripser,
                        index_diameter_t simplex,
                        index_t dim,
                        bool include_dim=false,
                        bool include_index=true,
                        bool include_simplices=true,
                        bool include_diameter=true,
                        bool include_rel=false)
{
	std::stringstream ss;
	if(include_dim) {
		ss << dim;
	}
	if(include_index) {
		if(include_dim) {
			ss << "-";
		}
		ss << get_index(simplex);
	}
	if(include_simplices) {
		if(include_index || include_dim) {
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
		if(include_index || include_simplices || include_dim) {
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

void output_boundary_matrix(ripser& ripser, std::ostream& os)
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
	std::sort(all_simplices.begin(), all_simplices.end(), total_filtration_order);
	simplex_boundary_enumerator facets(ripser);
	for(auto &simp : all_simplices) {
		if(simp.first == 0) {
			os << std::endl;
		} else {
			facets.set_simplex(simp.second, simp.first);
			while(facets.has_next()) {
				index_diameter_t facet = facets.next();
				os << simplex_tos(ripser, facet, simp.first - 1, 1, 1, 0, 0, 0);
				os << " ";
			}
			os << std::endl;
		}
	}
	// Output of row-indices that are non-zero (not ordered correctly)
	/*
	simplex_coboundary_enumerator cofacets(ripser);
	for(auto &simp : all_simplices) {
		cofacets.set_simplex(simp.second, simp.first);
		while(cofacets.has_next()) {
			index_diameter_t cofacet = cofacets.next();
			auto it = std::find(all_simplices.begin(), all_simplices.end(), std::make_pair(simp.first + 1, cofacet));
			if(it != all_simplices.end()) {
				os << std::distance(all_simplices.begin(), it) << " ";
			}
		}
		os << std::endl;
	}
	*/
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
		   << simplex_tos(ripser, s.second, s.first, true, true, true, true, true)
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
		   << simplex_tos(ripser, s.second, s.first, true, true, true, true, true)
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
				   << simplex_tos(ripser, s, i, false, true, (i != 0), false)
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
	for(index_t d = 0; d <= (index_t) ripser.config.dim_max; d++) {
		if(!(dim == -1 || d == dim)) {
			continue;
		}
		// TODO: This sorts by print order instead of homology order
		std::sort(ripser.hom_classes.at(d).begin(), ripser.hom_classes.at(d).end(),
		          homology_class_order);
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
			os << " ; [" << simplex_tos(ripser, hc.birth, d, 0, 1, 1, 0, 0);
			if(death == INF) {
				os << ", )";
			} else {
				// TODO: Dimension needs to be switch to +1 (homology) or -1 (cohomology)
				os << "," << simplex_tos(ripser, hc.death, d + 1, 0, 1, 1, 0, 0)
				   << ")";
			}
			if(with_reps) {
				os << " ; ";
				for(index_diameter_t& simp : hc.representative) {
					os << simplex_tos(ripser, simp, d, false, false, true, false);
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
	duration total_assemble_dur;
	duration total_reduction_dur;
	duration total_representative_dur;
	index_t total_reduction_count = 0;
	index_t total_clearing_count = 0;
	index_t total_emergent_count = 0;
	index_t total_apparent_count = 0;
	index_t total_addition_count = 0;
	for(info& info : ripser.infos) {
		std::string p = pref ? ("#i" + std::to_string(info.dim) + " ") : "";
		if(dim == -1 || info.dim == dim) {
			total_assemble_dur += info.assemble_dur;
			total_reduction_dur += info.reduction_dur;
			total_representative_dur += info.representative_dur;
			total_reduction_count += info.addition_count;
			total_clearing_count += info.clearing_count;
			total_emergent_count += info.emergent_count;
			total_apparent_count += info.apparent_count;
			total_addition_count += info.addition_count;
			duration total_dur = info.assemble_dur + info.reduction_dur + info.representative_dur;
			os << (pref ? "# " : "") << "reduction info in dim=" << info.dim << ":"
			   << std::endl;
			os << p << "  total column count:    "
			   << info.simplex_total_count << std::endl;
			os << p << "  ... cleared:           "
			   << info.clearing_count << std::endl;
			os << p << "  ... emergent:          "
			   << info.emergent_count << std::endl;
			os << p << "  ... apparent:          "
			   << info.apparent_count << std::endl;
			os << p << "  total reduction count: "
			   << info.simplex_reduction_count << std::endl;
			os << p << "  ... non-zero:          "
			   << info.class_count << std::endl;
			os << p << "  ... zero:              "
			   << info.zero_pers_count<< std::endl;
			os << p << "  addition count:        "
			   << info.addition_count << std::endl;
			os << p << "  total duration:        "
			   << total_dur.count() << "s" << std::endl;
			os << p << "  ... assembly:          "
			   << info.assemble_dur.count() << "s" << std::endl;
			os << p << "  ... reduction:         "
			   << info.reduction_dur.count() << "s" << std::endl;
			os << p << "  ... representative:    "
			   << info.representative_dur.count() << "s" << std::endl << std::endl;
		}
	}
	std::string p = (pref ? "# " : "");
	os << p << "total durations:" << std::endl;
	os << p << "  assembly:       " << total_assemble_dur.count() << std::endl;
	os << p << "  reduction:      " << total_reduction_dur.count() << std::endl;
	os << p << "  representative: " << total_representative_dur.count() << std::endl;
	os << p << "total counts:" << std::endl;
	os << p << "  cleared:   " << total_clearing_count << std::endl;
	os << p << "  emergent:  " << total_emergent_count << std::endl;
	os << p << "  apparent:  " << total_apparent_count << std::endl;
	os << p << "  reduction: " << total_reduction_count << std::endl;
	os << p << "  addition:  " << total_addition_count << std::endl;

}

void output_config(ripser& ripser, std::ostream& os, bool pref=false) {
	std::string p = pref ? "#c " : "";
	ripser_config& config = ripser.config;
	os << (pref ? "# " : "") << "config" << std::endl;
	os << p << "  input path:           " << config.input_path << std::endl;
	os << p << "  input type:           " << config.input_type << std::endl;
	os << p << "  output path:          " << config.output_path << std::endl;
	os << p << "  dataset vertex count: " << ripser.n << std::endl;
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
	os << p << "  absolute vertices:    ";
	if(config.absolute_subcomplex.size() == 0) {
		os << "0-" << ripser.n - 1 << " (all)";
	} else {
		auto pair = config.absolute_subcomplex.at(0);
		os << pair.first << "-" << pair.second;
		if(pair.first == 0 && pair.second == ripser.n - 1) {
			os << " (all)";
		} else {
			for(index_t i = 1; i < (index_t) config.absolute_subcomplex.size(); i++) {
				pair = config.absolute_subcomplex.at(i);
				os << "," << pair.first << "-" << pair.second;
			}
		}
	}
	os << std::endl;
	os << p << "  relative vertices:    ";
	if(config.relative_subcomplex.size() == 0) {
		os << "(empty)";
	} else {
		auto pair = config.config_relative_subcomplex.at(0);
		os << pair.first << "-" << pair.second;
		for(index_t i = 1; i < (index_t) config.config_relative_subcomplex.size(); i++) {
			pair = config.config_relative_subcomplex.at(i);
			os << "," << pair.first << "-" << pair.second;
		}
		os << " (is ";
		pair = config.relative_subcomplex.at(0);
		os << pair.first << "-" << pair.second;
		for(index_t i = 1; i < (index_t) config.relative_subcomplex.size(); i++) {
			pair = config.relative_subcomplex.at(i);
			os << "," << pair.first << "-" << pair.second;
		}
		os << ")" << std::endl;
	}
	os << std::endl;
}

std::ofstream get_writeout_stream(ripser& ripser, std::string suffix="") {
	std::string ip = ripser.config.input_path;
	std::string dataset_name = ip.substr(ip.find_last_of("/\\") + 1);
	std::stringstream ss;
	ss << ripser.config.output_path << "/";
	ss << dataset_name.substr(0, dataset_name.find_last_of(".")) << "_";
	ss << "d" << ripser.config.dim_max;
	if(ripser.config.config_threshold < 0) {
		ss << "tM";
	} else {
		std::string ts = std::to_string(ripser.config.config_threshold);
		ts = ts.erase(ts.find_last_not_of('0') + 1, std::string::npos);
		std::replace(ts.begin(), ts.end(), '.', 'f');
		ss << "t" << ts;
	}
	std::string cs = std::to_string(ripser.config.ratio);
	cs = cs.erase(cs.find_last_not_of('0') + 1, std::string::npos);
	std::replace(cs.begin(), cs.end(), '.', 'f');
	ss << "r" << cs;
	ss << "u" << ripser.config.use_enclosing_threshold
	   << ripser.config.use_union_find;
	if(!ripser.config.absolute_subcomplex.empty()) {
		ss << "_abs";
		for(auto& iint : ripser.config.absolute_subcomplex) {
			ss << iint.first << "-" << iint.second;
			if(&iint != &ripser.config.absolute_subcomplex.back()) {
				ss << "x";
			}
		}
	}
	if(!ripser.config.config_relative_subcomplex.empty()) {
		ss << "_rel";
		for(auto& iint : ripser.config.config_relative_subcomplex) {
			ss << iint.first << "-" << iint.second;
			if(&iint != &ripser.config.config_relative_subcomplex.back()) {
				ss << "x";
			}
		}
	}
	ss << suffix << ".txt";
	std::cout << ss.str() << std::endl;
	return std::ofstream(ss.str(), std::ofstream::trunc);
}

void write_standard_output(ripser& ripser,
                           bool with_reps=false,
                           bool with_simplex_list=false,
                           total_simplex_order ord=total_filtration_order) {
	std::ofstream ofs = get_writeout_stream(ripser);
	output_config(ripser, ofs, true); ofs << std::endl;
	if(with_simplex_list) {
		output_all_simplices(ripser, ofs, ord, true);
		ofs << std::endl;
	}
	output_info(ripser, ofs, -1, true); ofs << std::endl;
	output_barcode(ripser, ofs, with_reps, -1, true); ofs << std::endl;
}

void output_reduction_record(std::ostream& os, index_t dim, reduction_record& rr) {
	os << std::setw(2) << dim << " ; "
	   << std::setw(5) << rr.j << " ; "
	   << rr.to_zero << " ; "
	   << std::setw(5) << rr.add_simplex_boundary_count << " ; "
	   << std::setw(7) << rr.pop_count << " ; "
	   << std::setw(7) << rr.push_count << std::endl;
}

void write_analysis_rr(ripser& ripser, std::string suffix="") {
	std::ofstream ofs = get_writeout_stream(ripser, "_rr" + suffix);
	for(info& info : ripser.infos) {
		for(reduction_record& rr : info.red_records) {
			output_reduction_record(ofs, info.dim, rr);
		}
	}
}

#endif
