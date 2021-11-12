#ifndef PRINT_UTILS
#define PRINT_UTILS

#include "ripser_core.hpp"


void print_simplices(ripser& ripser, std::vector<index_diameter_t>& simplices, index_t d) {
  std::vector<index_t> vertices(ripser.n, -1);
  for(auto s : simplices) {
		ripser.get_simplex_vertices(get_index(s), d, ripser.n, vertices);
		std::cout << "    (" << get_index(s) << " :: "<< get_diameter(s) << " :: ";
		for(auto i : vertices) {
			if(i >= 0) {
			  std::cout << i << "'";
			}
		}
		std::cout << ")" << std::endl;
		std::fill(vertices.begin(), vertices.end(), -1);
  }
}

void list_all_simplices(ripser& ripser) {
	std::cout << "info: list of all simplices (id :: diam :: vertices)"
	          << "      in ascending filtration order" << std::endl;
	std::vector<index_diameter_t> simpl_prev;
	std::vector<index_diameter_t> simpl_curr;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simpl_prev.push_back(index_diameter_t(i, diam));
	}
	std::cout << "  dim 0" << std::endl;
	std::sort(simpl_prev.begin(), simpl_prev.end(), filtration_order);
	print_simplices(ripser, simpl_prev, 0);
	
	simplex_coboundary_enumerator e(ripser);
	index_t dim = 1;
	for(auto i : simpl_prev) {
		e.set_simplex(i, dim - 1);
		while(e.has_next(false)) {
		  simpl_curr.push_back(e.next());
		}
	}
	std::cout << "  dim 1:" << std::endl;
	std::sort(simpl_curr.begin(), simpl_curr.end(), filtration_order);
	print_simplices(ripser, simpl_curr, dim);
	
	simpl_prev.swap(simpl_curr);
	simpl_curr.clear();
	dim = 2;
	for(auto i : simpl_prev) {
		e.set_simplex(i, dim - 1);
		while(e.has_next(false)) {
		  simpl_curr.push_back(e.next());
		}
	}
	std::cout << "  dim 2:" << std::endl;
	std::sort(simpl_curr.begin(), simpl_curr.end(), filtration_order);
	print_simplices(ripser, simpl_curr, dim);
}

int sprint_element(char* buf, index_t el, int pad) {
	char format[16]; format[0] = '\0';
	sprintf(format, "%s%i%s ", "%", pad, "i");
	return sprintf(buf, format, el);
}

int sprint_pad(char *buf, int pad) {
	return sprintf(buf, "%*s. ", pad - 1, "");
}

int sprint_simplex(char *buf, ripser& ripser, index_diameter_t simplex, index_t dim) {
	std::vector<index_t> vertices(ripser.n, -1);
	ripser.get_simplex_vertices(get_index(simplex), dim, ripser.n, vertices);
	int offs = 0;
	offs += sprintf(offs + buf, "%li-", get_index(simplex));
	for(auto v : vertices) {
	  if(v == -1) continue;
	  offs += sprintf(offs + buf, "%li'", v);
	}
	offs += sprintf(offs + buf, " %f", get_diameter(simplex));
	return offs;
}

void print_simplex(ripser& ripser, index_diameter_t simplex, index_t dim) {
	char buf[1024]; buf[0] ='\0';
	sprint_simplex(buf, ripser, simplex, dim);
	printf("%s", buf);
}

template <typename Column>
void print_column(ripser& ripser, Column &column, index_t dim) {
	Column col(column); // copy
	char buf[1024]; buf[0] = '\0';
	int offs = 0;
	while(!col.empty()) {
		index_diameter_t simp = col.top();
		col.pop();
		offs += sprint_simplex(buf + offs, ripser, simp, dim);
		offs += sprintf(buf + offs, " ");
	}
	printf("%s\n", buf);
}

void print_v(compressed_sparse_matrix& v,
             std::vector<index_diameter_t> columns_to_reduce) {
	int offs = 0;
	int pad = 3;
	char *buf = (char *) malloc(columns_to_reduce.size() * columns_to_reduce.size() *
	                            2 * pad + columns_to_reduce.size() + 10000);
	buf[0] = '\0';
	for(index_t row = 0; row < (index_t) columns_to_reduce.size(); ++row) {
		index_diameter_t row_ele = columns_to_reduce.at(row);
		for(index_t col = 0; col < (index_t) columns_to_reduce.size(); ++col) {
			if(v.search_column(col, row_ele)) {
				offs += sprint_element(buf + offs, 1, pad);
			} else {
				offs += sprint_pad(buf + offs, pad);
			}
		}
		offs += sprintf(offs + buf, "\n");
	}
	printf("%s", buf);
	free(buf);
}

void print_vrow(compressed_sparse_matrix& v,
                std::vector<index_diameter_t> columns_to_reduce) {
	int offs = 0;
	int pad = 3;
	char *buf = (char *) malloc(columns_to_reduce.size() * columns_to_reduce.size() *
	                            pad + columns_to_reduce.size() + 10000);
	buf[0] = '\0';
	for(index_t row = 0; row < (index_t) columns_to_reduce.size(); ++row) {
		for(index_t col = 0; col < (index_t) columns_to_reduce.size(); ++col) {
			index_diameter_t col_ele = columns_to_reduce.at(col);
			if(v.search_column(row, col_ele)) {
				offs += sprint_element(buf + offs, 1, pad);
			} else {
				offs += sprint_pad(buf + offs, pad);
			}
		}
		offs += sprintf(offs + buf, "\n");
	}
	printf("%s", buf);
	free(buf);
}

void print_mat_simplices(ripser& ripser, compressed_sparse_matrix& v, index_t dim) {
	//std::cout << "Bounds: ";
	//for(auto b : v.bounds) {
	//	std::cout << b << " ";
	//}
	std::cout << std::endl;
	for(index_t i = 0; i < (index_t) v.entries.size(); i++) {
		print_simplex(ripser, v.get_entry(i), dim);
		std::cout << " ";
		if(std::find(v.bounds.begin(), v.bounds.end(), i + 1) != v.bounds.end()) {
			std::cout << std::endl;
		}
	}
}

void print_mat(compressed_sparse_matrix& mat) {
	char buf[1024]; buf[0] = '\0';
	int offs = 0;
	int pad = 3;
	index_t max_row_index = 0;
	for(index_t i = 0; i < mat.size(); ++i) {
		max_row_index = std::max(max_row_index,
		                         mat.column_end(i) - mat.column_start(i));
	}
	std::cout << std::endl;
	for(index_t row = 0; row < max_row_index; ++row) {
		for(index_t col = 0; col < mat.size(); ++col) {
			if(row < mat.column_end(col) - mat.column_start(col)) {
				index_t ele = get_index(mat.get_entry(row, col));
				offs += sprint_element(offs + buf, ele, pad);
			} else {
				offs += sprint_pad(offs + buf, pad);
			}
		}
		offs += sprintf(offs + buf, "\n");
	}
	printf("%s", buf);
}

void output_barcode(ripser& ripser, std::ostream& os, bool with_reps=false, index_t dim=-1, bool pref=false)
{
	std::sort(ripser.hom_classes.begin(), ripser.hom_classes.end(),
	          homology_class_print_order);
	index_t current_dim = -1;
	for(auto& hc : ripser.hom_classes) {
		std::string p = pref ? ("#b" + std::to_string(hc.dim) + " ") : "";
		if(hc.dim > current_dim) {
			if(current_dim > -1) {
				os << std::endl;
			}
			os << (pref ? "# " : "") << "barcode in dim=" << hc.dim
			   << std::endl;
			current_dim = hc.dim;
		}
		value_t birth = std::max(0.0f, get_diameter(hc.birth));
		value_t death = get_diameter(hc.death);
		os << p << "  [" << birth;
		if(death == INF) {
			os << ", )";
		} else {
			os << "," << death << ")";
		}
		os << " - [" << get_index(hc.birth);
		if(death == INF) {
			os << ", )";
		} else {
			os << "," << get_index(hc.death) << ")";
		}
		if(with_reps) {
			os << " - ";
			//for(index_diameter_t simp : hc.representative) {
			//	print_simplex(ripser, simp, hc.dim, false);
			//	os << " ";
			//}
		}
		os << std::endl;
	}
}

void output_info(ripser& ripser, std::ostream& os, index_t dim=-1, bool pref=false)
{
	for(auto info : ripser.infos) {
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
			os << p <<"  assemble duration:       "
			   << info.assemble_dur.count() << "s" << std::endl;
			os << p << "  reduction duration:      "
			   << info.reduction_dur.count() << "s" << std::endl;
			os << p << "  representative duration: "
			   << info.representative_dur.count() << "s" << std::endl;
			if(info.dim < ripser.config.dim_max) {
				os << std::endl;
			}
		}
	}
}

void output_config(ripser& ripser, std::ostream& os, bool pref=false) {
	std::string p = pref ? "#c " : "";
	ripser_config& config = ripser.config;
	os << p << "config" << std::endl;
	os << p << "  dataset at file path: " << config.file_path << std::endl;
	os << p << "  output at file path:  " << config.output_path << std::endl;
	os << p << "  dataset input type:   " << config.input_type << std::endl;
	os << p << "  dim max:              " << config.dim_max << std::endl;
	os << p << "  ratio:                " << config.ratio << std::endl;
	os << p << "  threshold:            " << ripser.threshold;
	if(config.use_enclosing_threshold) {
		os << p << " (enclosing threshold)" << std::endl;
	} else if(config.config_threshold < 0) {
		os << p << " (maximum distance)" << std::endl;
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

void write_output(ripser& ripser, bool with_reps=false) {
	std::string ip = ripser.config.file_path;
	std::string dataset_name = ip.substr(ip.find_last_of("/\\") + 1);
	std::string rel_name = "rel";
	for(index_t i = 0; i < ripser.config.relative_subcomplex.size(); i++) {
		rel_name += std::to_string(ripser.config.relative_subcomplex.at(i).first);
		rel_name += "-" + std::to_string(ripser.config.relative_subcomplex.at(i).second);
		if(i < ripser.config.relative_subcomplex.size() - 1) {
			rel_name += "x";
		}
	}
	if(ripser.config.relative_subcomplex.empty()) {
		rel_name = "relnone";
	}
	std::string thresh_name = "thresh" + std::to_string(ripser.threshold);
	thresh_name.erase(thresh_name.find_last_not_of('0') + 1, std::string::npos);
	std::replace(thresh_name.begin(), thresh_name.end(), '.', 'd');
	std::string rep_name = (with_reps ? "_reps" : "");
	std::string filename = ripser.config.output_path + "/"
		+ dataset_name
		+ rep_name + "_"
		+ thresh_name + "_"
		+ rel_name
		+ ".txt";
	std::ofstream ofs(filename, std::ofstream::trunc);
	output_config(ripser, ofs, true); ofs << std::endl;
	output_info(ripser, ofs, -1, true); ofs << std::endl;
	output_barcode(ripser, ofs, with_reps, -1, true);
}

#endif
