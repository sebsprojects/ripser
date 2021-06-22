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

int sprint_simplex(char *buf, ripser& ripser, index_t simplex, index_t dim) {
	std::vector<index_t> vertices(ripser.n, -1);
	ripser.get_simplex_vertices(simplex, dim, ripser.n, vertices);
	int offs = 0;
	offs += sprintf(offs + buf, "%li-", simplex);
	for(auto v : vertices) {
	  if(v == -1) continue;
	  offs += sprintf(offs + buf, "%li'", v);
	}
	return offs;
}

void print_simplex(ripser& ripser, index_t simplex, index_t dim) {
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
		offs += sprint_simplex(buf + offs, ripser, get_index(simp), dim);
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
		index_t row_ele = get_index(columns_to_reduce.at(row));
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
			index_t col_ele = get_index(columns_to_reduce.at(col));
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
		print_simplex(ripser, get_index(v.get_entry(i)), dim);
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

void print_barcode(ripser& ripser, barcode& barcode) {
	std::cout << "persistence intervals in dim " << barcode.dim << ":" << std::endl;
	std::sort(barcode.hom_classes.begin(),
	          barcode.hom_classes.end(),
	          homology_class_order);
	for(auto hc : barcode.hom_classes) {
		value_t birth = std::max(0.0f, hc.birth);
		value_t death = hc.death;
		std::cout << "[" << birth;
		if(death == INF) {
			std::cout << ", ) :: ";
		} else {
			std::cout << "," << death << ") :: ";
		}
		for(index_t simp : hc.representative) {
			print_simplex(ripser, simp, barcode.dim);
			std::cout << " ";
		}
		std::cout << std::endl;
	}
}

void write_dim1_cycles(ripser& ripser, std::string filename) {
	barcode& barcode = ripser.barcodes.at(1);
	std::sort(barcode.hom_classes.begin(),
	          barcode.hom_classes.end(),
	          homology_class_order);
	std::ofstream ofs(filename, std::ofstream::trunc);
	for(auto hc : barcode.hom_classes) {
		ofs << "# " << std::max(0.0f, hc.birth) << " " << hc.death << std::endl;
		std::vector<index_t> vertices(2, -1);
		for(index_t simplex : hc.representative) {
			ripser.get_simplex_vertices(simplex, 1, ripser.n, vertices);
			ofs << vertices.at(0) << " " << vertices.at(1) << std::endl;
		}
		ofs << std::endl;
	}
	ofs.close();
}

void print_barcodes(ripser& ripser) {
	std::sort(ripser.barcodes.begin(), ripser.barcodes.end(), barcode_order);
	for(auto b : ripser.barcodes) {
		print_barcode(ripser, b);
	}
}

void print_info(ripser& ripser, info& info) {
	std::cout << "info in dim " << info.dim << ":" << std::endl;
	std::cout << "  total simplex count:     "
	          << info.simplex_total_count << std::endl;
	std::cout << "  reduction simplex count: "
	          << info.simplex_reduction_count << std::endl << std::endl;

	std::cout << "  clearing count:          "
	          << info.clearing_count << std::endl;
	std::cout << "  emergent count:          "
	          << info.emergent_count << std::endl;
	std::cout << "  apparent count:          "
	          << info.apparent_count << std::endl;

	std::cout << "  assemble duration:       "
	          << info.assemble_dur.count() << "s" << std::endl;
	std::cout << "  reduction duration:      "
	          << info.reduction_dur.count() << "s" << std::endl;
	std::cout << "  representative duration: "
	          << info.representative_dur.count() << "s" << std::endl;
}

void print_infos(ripser& ripser) {
	for(auto i : ripser.infos) {
		print_info(ripser, i);
	}
}

#endif
