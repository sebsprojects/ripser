#ifndef RIPSER_CORE
#define RIPSER_CORE

#include "ripser_config.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <chrono>
#include <iomanip>


/* **************************************************************************
 * Basic Types
 * *************************************************************************/

typedef int64_t index_t;
typedef float value_t;

static const value_t INF = std::numeric_limits<value_t>::infinity();

void check_overflow(index_t i) {
	if(i < 0) {
		throw std::overflow_error("simplex index " +
		                          std::to_string((uint64_t)i) +
		                          " in filtration is out of bounds");
	}
}

// The basic data-type for simplices, given by the index of the simplex
// w.r.t. binomial numbering system and its (cached) diameter
typedef std::pair<index_t, value_t> index_diameter_t;

index_t get_index(const index_diameter_t i) {
	return i.first;
}

value_t get_diameter(const index_diameter_t i) {
	return i.second;
}


/* **************************************************************************
 * Binomial Numbering System and Ordering
 * *************************************************************************/

typedef bool (simplex_order(const index_diameter_t&, const index_diameter_t&));
typedef bool (total_simplex_order(const std::pair<index_t, index_diameter_t>&,
                                  const std::pair<index_t, index_diameter_t>&));

bool total_filtration_order(const std::pair<index_t, index_diameter_t>& a,
                            const std::pair<index_t, index_diameter_t>& b)
{
	return (get_diameter(a.second) < get_diameter(b.second)) ||
	       ((get_diameter(a.second) == get_diameter(b.second)) &&
	        (a.first < b.first)) ||
	       ((get_diameter(a.second) == get_diameter(b.second)) &&
	        (a.first == b.first) &&
	        (get_index(a.second) > get_index(b.second)));
}

bool total_reverse_filtration_order(const std::pair<index_t, index_diameter_t>& a,
                                    const std::pair<index_t, index_diameter_t>& b)
{
	return total_filtration_order(b, a);
}


// a < b iff diam(a) < diam(b), if equal index(a) > index(b)
bool filtration_order(const index_diameter_t& a,
                      const index_diameter_t& b) {
	return (get_diameter(a) < get_diameter(b)) ||
	        ((get_diameter(a) == get_diameter(b)) &&
	         (get_index(a) > get_index(b)));
}

struct filtration_order_comp {
  bool operator()(const index_diameter_t& a, const index_diameter_t& b) {
	return filtration_order(a, b);
  }
};

// TODO(seb): Should be equivalent to filtration_order(b, a)

// a < b iff diam(a) > diam(b), if equal index(a) < index(b)
bool reverse_filtration_order(const index_diameter_t& a,
                              const index_diameter_t& b) {
	return (get_diameter(a) > get_diameter(b)) ||
	       ((get_diameter(a) == get_diameter(b)) &&
	        (get_index(a) < get_index(b)));
}

struct reverse_filtration_order_comp {
  bool operator()(const index_diameter_t& a, const index_diameter_t& b) {
	return reverse_filtration_order(a, b);
  }
};

// Binomial lookup table
struct binomial_coeff_table {
	std::vector<std::vector<index_t>> B;

	binomial_coeff_table(index_t n, index_t k) : B(n + 1) {
		for (index_t i = 0; i <= n; ++i) {
			B[i].resize(k + 1, 0);
			B[i][0] = 1;
			for (index_t j = 1; j < std::min(i, k + 1); ++j)
				B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			if (i <= k) B[i][i] = 1;
			check_overflow(B[i][std::min(i >> 1, k)]);
		}
	}

	index_t operator()(index_t n, index_t k) const {
		// TODO: Unsafe casts?
		assert((size_t) n < B.size() && (size_t) k < B[n].size() && n >= k - 1);
		return B[n][k];
	}
};


/* **************************************************************************
 * Storage and Matrices
 * *************************************************************************/

// Type to store the input data as a distance matrix
struct compressed_lower_distance_matrix {

	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	compressed_lower_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)),
			rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	compressed_lower_distance_matrix(const compressed_lower_distance_matrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();
		for (size_t i = 1; i < size(); ++i)
			for (size_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const {
		return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
	}

	void set_entry(const index_t i, const index_t j, const value_t val) {
		if(i < j) {
			rows[j][i] = val;
		}
		if(i > j) {
			rows[i][j] = val;
		}
	}

	size_t size() const {
	  return rows.size();
	}
	
	void init_rows() {
		value_t* pointer = &distances[0];
		for (size_t i = 1; i < size(); ++i) {
			rows[i] = pointer;
			pointer += i;
		}
	}
};
typedef struct compressed_lower_distance_matrix DistanceMatrix;

struct compressed_sparse_matrix {

	std::vector<index_t> bounds;
	std::vector<index_diameter_t> entries;

	compressed_sparse_matrix() {}

	index_t column_start(const index_t index) const {
		return index == 0 ? 0 : bounds.at(index - 1);
	}

	index_t column_end(const index_t index) const {
		return bounds.at(index);
	}

	index_diameter_t get_entry(const index_t index) {
		return entries.at(index);
	}

	index_diameter_t get_entry(const index_t row_index, const index_t col_index) {
		return entries.at(column_start(col_index) + row_index);
	}

	bool search_column(const index_t col_index, const index_diameter_t el) {
		// Linear time find:
		for(index_t i = column_start(col_index); i < column_end(col_index); ++i) {
			if(get_entry(i) == el) {
				return true;
			}
		}
		return false;
	}

	index_t size() const {
		return (index_t) bounds.size();
	}

	index_t get_total_size() const {
		return (index_t) bounds.at(bounds.size() - 1);
	}

	void append_column() {
		bounds.push_back(entries.size());
	}

	void push_back(const index_diameter_t idx) {
		assert(size() >= 0);
		entries.push_back(idx);
		++bounds.back();
	}
};


/* **************************************************************************
 * Matrix Read-In
 * *************************************************************************/

bool is_absolute_index(ripser_config& config, index_t idx) {
	if(config.absolute_subcomplex.empty()) {
		return true;
	}
	for(index_t i = 0; i < (index_t) config.absolute_subcomplex.size(); i++) {
		auto pair = config.absolute_subcomplex.at(i);
		if(idx >= pair.first && idx <= pair.second) {
			return true;
		}
	}
	return false;
}

compressed_lower_distance_matrix read_distance_matrix(ripser_config& config,
                                                      std::istream& input_stream)
{
	std::vector<value_t> distances;
	std::string line;
	value_t value;
	int offs = config.input_type == "lower_distance_matrix" ? 0 : 0;
	for(int i = 0; std::getline(input_stream, line); ++i) {
		if(!is_absolute_index(config, i)) {
			continue;
		}
		std::istringstream s(line);
		for (int j = 0; j < i + offs; ++j) {
			s >> value;
			if(!is_absolute_index(config, j)) {
				continue;
			}
			distances.push_back(value);
			s.ignore();
		}
	}
	return compressed_lower_distance_matrix(std::move(distances));
}

std::vector<std::vector<value_t>> read_point_cloud(ripser_config& config,
                                                   std::istream& input_stream)
{
	std::vector<std::vector<value_t>> points;
	std::string line;
	value_t value;
	for(int i = 0; std::getline(input_stream, line); ++i) {
		if(!is_absolute_index(config, i)) {
			continue;
		}
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if(!point.empty()) {
			points.push_back(point);
		}
		assert(point.size() == points.front().size());
	}
	return points;
}

struct diff_squared {
	value_t operator()(value_t u, value_t v) {
		return (u - v) * (u - v);
	}
};

std::vector<value_t> point_cloud_to_distance_vector(std::vector<std::vector<value_t>>& points)
{
	std::vector<value_t> distances;
	for(size_t i = 0; i < points.size(); i++) {
		for(size_t j = 0; j < i; j++) {
			value_t d = std::sqrt(std::inner_product(points.at(i).begin(),
			                                         points.at(i).end(),
			                                         points.at(j).begin(),
			                                         0.0,
			                                         std::plus<value_t>(),
			                                         diff_squared()));
			distances.push_back(d);
		}
	}
	return distances;
}

compressed_lower_distance_matrix read_input(ripser_config config) {
	std::ifstream file_stream(config.input_path);
	if(!file_stream) {
		std::cerr << "error: couldn't open matrix file: " << config.input_path
		          << std::endl;
		exit(-1);
	}
	if(config.input_type == "lower_distance_matrix" ||
	   config.input_type == "full_distance_matrix") {
		return read_distance_matrix(config, file_stream);
	} else if(config.input_type == "point_cloud") {
		std::vector<std::vector<value_t>> points =
			read_point_cloud(config, file_stream);
		std::vector<value_t> distances =
			point_cloud_to_distance_vector(points);
		return compressed_lower_distance_matrix(std::move(distances));
	} else {
		std::cerr << "error: invalid matrix type: " << config.input_type
		          << std::endl;
		exit(-1);
	}
}

value_t compute_enclosing_radius(const DistanceMatrix& dist) {
	value_t enclosing_radius = INF;
	for(size_t i = 0; i < dist.size(); ++i) {
		value_t r_i = -INF;
		for (size_t j = 0; j < dist.size(); ++j) {
			r_i = std::max(r_i, dist(i, j));
		}
		enclosing_radius = std::min(enclosing_radius, r_i);
	}
	return enclosing_radius;
}

value_t compute_max_distance(const DistanceMatrix& dist) {
	value_t max = -INF;
	for(auto d : dist.distances) {
		max = std::max(max, d);
	}
	return max;
}


/* **************************************************************************
 * Types for Reduction Algorithm
 * *************************************************************************/


struct index_hash {
	size_t operator()(const index_t& e) const {
		return std::hash<index_t>()(e);
	}
};

struct equal_index {
	bool operator()(const index_t& e, const index_t& f) const {
		return e == f;
	}
};

// hash map with key = pivot index and value = column index
// the column index is relative to the current dimension since we only have
// a slice of V w.r.t to that dimension
typedef std::unordered_map<index_t,
                           size_t,
                           index_hash,
                           equal_index> entry_hash_map;


/* **************************************************************************
 * Ripser
 * *************************************************************************/

struct homology_class {
	index_t dim;
	index_diameter_t birth;
	index_diameter_t death;
	std::vector<index_diameter_t> representative;

	homology_class(index_t _dim, index_diameter_t _birth, index_diameter_t _death, std::vector<index_diameter_t> rep)
		: dim(_dim),
		  birth(_birth),
		  death(_death),
		  representative(rep)
	{ }
};

bool homology_class_print_order(homology_class& a, homology_class& b) {
	return (a.dim < b.dim) ||
	       ((a.dim == b.dim) &&
	        (get_diameter(a.birth) < get_diameter(b.birth))) ||
	       ((a.dim == b.dim) &&
	        (get_diameter(a.birth) == get_diameter(b.birth)) &&
	        (get_diameter(a.death) < get_diameter(b.death)));
}

bool homology_class_order(homology_class& a, homology_class& b) {
	return filtration_order(a.birth, b.birth);
}

bool reverse_homology_class_order(homology_class& a, homology_class& b) {
	return reverse_filtration_order(a.birth, b.birth);
}

typedef std::chrono::steady_clock::time_point time_point;
typedef std::chrono::duration<double> duration;

struct reduction_record {
	index_t j;
	time_point start;
	time_point end;
	uint64_t addition_count;
	uint64_t addition_apparent_count;
	uint64_t coboundary_element_count;

	bool to_zero;

	uint64_t add_simplex_boundary_count;
	uint64_t pop_count;
	uint64_t push_count;
	uint64_t app_facet_count;
	uint64_t app_cofacet_count;
	uint64_t max_mem;

	reduction_record()
		: j(-1), start(), end(), addition_count(0), addition_apparent_count(0),
		  coboundary_element_count(0), to_zero(true), add_simplex_boundary_count(0),
		  pop_count(0), push_count(0), app_facet_count(0), app_cofacet_count(0), max_mem(0)
  { }

	reduction_record(index_t _j, time_point _start)
		: j(_j), start(_start), end(), addition_count(0), addition_apparent_count(0),
		  coboundary_element_count(0), to_zero(true), add_simplex_boundary_count(0),
		  pop_count(0), push_count(0), app_facet_count(0), app_cofacet_count(0), max_mem(0)
	{
	}
};

struct info {
	index_t dim;
	uint64_t clearing_count;
	uint64_t emergent_count;
	uint64_t apparent_count;
	uint64_t simplex_total_count;
	uint64_t simplex_reduction_count;
	uint64_t class_count;
	uint64_t zero_pers_count;
	uint64_t addition_count;
	
	duration assemble_dur;
	duration reduction_dur;
	duration representative_dur;

	std::vector<reduction_record> red_records;

	info(index_t _dim)
		: dim(_dim), clearing_count(0), emergent_count(0), apparent_count(0),
		  simplex_total_count(0), simplex_reduction_count(0), class_count(0),
		  zero_pers_count(0), addition_count(0),
		  assemble_dur(), reduction_dur(), representative_dur()
	{ }
};

time_point get_time() {
	return std::chrono::steady_clock::now();
}

duration get_duration(time_point start, time_point end) {
	return end - start;
}

struct ripser {

	const DistanceMatrix dist;
	const index_t n;
	const binomial_coeff_table binomial_coeff;

	mutable ripser_config config;
	mutable value_t threshold;
	mutable std::vector<std::vector<homology_class>> hom_classes;
	mutable std::vector<info> infos;
	mutable index_t current_reduction_record_dim;
	mutable reduction_record catch_all_record;

	ripser(ripser_config _config)
		: dist(read_input(_config)),
		  n(dist.size()),
		  binomial_coeff(n, _config.dim_max + 2),
		  config(_config),
		  threshold(-1),
		  hom_classes(),
		  infos(),
		  current_reduction_record_dim(-1),
		  catch_all_record()
	{
		if(config.use_enclosing_threshold) {
			threshold = compute_enclosing_radius(dist);
		}
		if(config.config_threshold < 0.0) {
			threshold = compute_max_distance(dist);
		} else {
			threshold = config.config_threshold;
		}
		for(index_t i = 0; i <= config.dim_max + 1; ++i) {
			infos.push_back(info(i));
			hom_classes.push_back(std::vector<homology_class>());
		}
	}
	
	// For a k-simplex idx, return the vertex (of idx) of maximal index
	// n is an upper bound for the search
	index_t get_max_vertex(const index_t idx, const index_t k, const index_t n) const
	{
		index_t top = n;
		index_t bot = k - 1;
		// Find max{ i | binomial_coeff(i,k) <= idx } via binary search
		if(binomial_coeff(top, k) > idx) {
			index_t count = top - bot;
			while(count > 0) {
				index_t half = count >> 1; // divide by 2
				index_t mid = top - half;
				if(binomial_coeff(mid, k) > idx) { // lower half
					top = mid - 1;
					count -= half + 1;
				} else { // upper half
					count = half;
				}
			}
		}
		return top;
	}

	// Write all vertices of the dim-simplex idx into vertices in ascending index order
	// n is an upper bound on indices
	void get_simplex_vertices(index_t idx, const index_t dim, index_t n,
	                          std::vector<index_t>& vertices) const {
		--n;
		for(index_t k = dim + 1; k > 0; --k) {
			n = get_max_vertex(idx, k, n);
			vertices.at(k - 1) = n;
			idx -= binomial_coeff(n, k);
		}
	}

	index_t get_first_relative_vertex() const {
		for(index_t i = 0; i < (index_t) config.relative_subcomplex.size(); i++) {
			return config.relative_subcomplex.at(i).first;
		}
		return -1;
	}

	// TODO: Might be problematic if filtration order instead of
	// reverse_filt_order since indices in relative_subcomplex no
	// longer correspond to vertex indices
	bool is_relative_vertex(index_t idx) const {
		for(index_t i = 0; i < (index_t) config.relative_subcomplex.size(); i++) {
			auto pair = config.relative_subcomplex.at(i);
			if(idx >= pair.first && idx <= pair.second) {
				return true;
			}
		}
		return false;
	}

	bool is_relative_simplex(index_t idx, const index_t dim) const {
		index_t n = this->n - 1;
		for(index_t k = dim + 1; k > 0; --k) {
			n = get_max_vertex(idx, k, n);
			if(!is_relative_vertex(n)) {
				return false;
			}
			idx -= binomial_coeff(n, k);
		}
		return true;
	}

	// Compute the diameter for a k-simplex idx
	// TODO(seb): dim vs k issue
	value_t compute_diameter(const index_t idx, const index_t dim) const {
		value_t diam = -INF;
		std::vector<index_t> vertices(dim + 1, -1); // TODO(seb): Default value?
		get_simplex_vertices(idx, dim, dist.size(), vertices);
		for (index_t i = 0; i <= dim; ++i)
			for (index_t j = 0; j < i; ++j) {
				diam = std::max(diam, dist(vertices[i], vertices[j]));
			}
		return diam;
	}

	void add_hom_class(index_t dim, index_diameter_t birth, index_diameter_t death, std::vector<index_diameter_t> rep = std::vector<index_diameter_t>()) {
		hom_classes.at(dim).push_back(homology_class(dim, birth, death, rep));
	}

	reduction_record& get_current_reduction_record() {
		// Return dummy
		if(current_reduction_record_dim == -1) {
			return catch_all_record;
		}
		return infos.at(current_reduction_record_dim).red_records.back();
	}

	void add_reduction_record(index_t dim, index_t j, time_point start) {
		infos.at(dim).red_records.push_back(reduction_record(j, start));
		current_reduction_record_dim = dim;
	}

	void complete_reduction_record(time_point end,
	                               index_t add_count,
	                               index_t add_app_count,
	                               index_t coboundary_count) {
		//TODO(seb): hacky hacky hacky output handling
		index_t dim = current_reduction_record_dim;
		static index_t thresh_counter = -1;
		static index_t curr_dim = -1;
		static index_t thresh_lim = 50;
		static index_t thresh_linelim = 50;
		reduction_record& rec = infos.at(dim).red_records.back();
		rec.end = end;
		rec.addition_count = add_count;
		rec.addition_apparent_count = add_app_count;
		rec.coboundary_element_count = coboundary_count;
		if(config.print_progress) {
			if(dim != curr_dim) {
				std::cout << std::endl;
				std::cout << "reduction progress in dim=" << dim << std::endl;
				curr_dim = dim;
				thresh_counter = -1;
			}
			if(add_count + add_app_count < thresh_lim) {
				if(thresh_counter == 0) {
				  std::cout << std::endl << "  ";
				}
				if(thresh_counter == -1) {
					thresh_counter++;
					std::cout << "  ";
				}
				std::cout << "." << std::flush;
				thresh_counter = (thresh_counter + 1) % (thresh_linelim - 2);
			} else {
				index_t ms = (index_t)
					(get_duration(rec.start, end).count() * 1000.0);
				if(thresh_counter > 0) {
					std:: cout << std::endl;
				}
				// was previously in add_red_record
				std::cout << "  "
						  << (rec.j+1) << "/"
						  << infos.at(dim).simplex_reduction_count;
				//std::cout << " :: tim=(" << std::setw(6) << ms << "ms)"
				//		  << " :: add=(" << std::setw(5) << add_count << "/"
				//		  << std::setw(5) << add_app_count << ")"
				//		  << " :: cob=" << std::setw(4) << coboundary_count
				//		  << std::endl;
				std::cout << (rec.to_zero ? " z=YES" : " z=NO ")
				          << " :: add=" << std::setw(3) << rec.add_simplex_boundary_count
				          << " :: push=" << std::setw(3) << rec.push_count
						  << " :: pop=" << std::setw(3) << rec.pop_count << std::endl;
				thresh_counter = -1;
			}
		}
	}
};


/* **************************************************************************
 * Boundary and Coboundary Enumerators
 * *************************************************************************/

class simplex_boundary_enumerator {

private:
	const ripser& parent;

	// Simplex to enumerate off given as ((index, diam), dim)
	mutable index_diameter_t simplex;
	mutable index_t dim;

	// Ongoing enumeration
	mutable index_t idx_below;
	mutable index_t idx_above;
	mutable index_t j;
	mutable index_t k;

public:
	simplex_boundary_enumerator(const ripser& _parent) : parent(_parent)
	{ }

	// _simplex should be a _dim-simplex
	void set_simplex(const index_diameter_t _simplex, const index_t _dim) {
		idx_below = get_index(_simplex);
		idx_above = 0;
		j = parent.n - 1;
		k = _dim;
		simplex = _simplex;
		dim = _dim;
	}

	bool has_next(bool all_facets = true) {
		return (k >= 0);
	}

	index_diameter_t next() {
		j = parent.get_max_vertex(idx_below, k + 1, j);
		//std::cout << "j=" << j << " idx_above=" << idx_above
		//          << " binom=" << parent.binomial_coeff(j, k + 1)
		//          << " idx_below=" << idx_below << std::endl;
		index_t face_index = idx_above - parent.binomial_coeff(j, k + 1) + idx_below;
		value_t face_diameter = parent.compute_diameter(face_index, dim - 1);
		// Update values for next call
		idx_below -= parent.binomial_coeff(j, k + 1);
		idx_above += parent.binomial_coeff(j, k);
		--k;
		return index_diameter_t(face_index, face_diameter);
	}
};

class simplex_coboundary_enumerator {

private:
	const ripser& parent;

	// Simplex to enumerate off given as ((index, diam), dim)
	mutable index_diameter_t simplex;
	mutable index_t dim;
	std::vector<index_t> vertices;

	// Ongoing enumeration
	mutable index_t idx_below;
	mutable index_t idx_above;
	mutable index_t j;
	mutable index_t k;

public:
	simplex_coboundary_enumerator(const ripser& _parent) : parent(_parent)
	{ }

	// _simplex should be a _dim-simplex
	void set_simplex(const index_diameter_t _simplex, const index_t _dim) {
		idx_below = get_index(_simplex);
		idx_above = 0;
		j = parent.n - 1;
		k = _dim + 1;
		simplex = _simplex;
		vertices.resize(_dim + 1);
		parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices);
	}

	bool has_next(bool all_cofacets = true) {
		return (j >= k && (all_cofacets || parent.binomial_coeff(j, k) > idx_below));
	}

	index_diameter_t next() {
		while (parent.binomial_coeff(j, k) <= idx_below) {
			idx_below -= parent.binomial_coeff(j, k);
			idx_above += parent.binomial_coeff(j, k + 1);
			--j;
			--k;
			assert(k != -1);
		}
		//TODO(seb): Do we really need the diameter of simplex in this enumerator?
		value_t cofacet_diameter = get_diameter(simplex);
		for (index_t i : vertices) {
		  cofacet_diameter = std::max(cofacet_diameter, parent.dist(j, i));
		}
		index_t cofacet_index = idx_above + parent.binomial_coeff(j--, k + 1) + idx_below;
		return index_diameter_t(cofacet_index, cofacet_diameter);
	}
};

void assemble_all_simplices(ripser& ripser,
                            std::vector<index_diameter_t>& simplices,
                            index_t dim,
                            bool respect_relative=true)
{
	simplices.clear();
	if(dim < ripser.n) {
		std::vector<index_diameter_t> next_simplices;
		for(index_t i = 0; i < ripser.n; i++) {
			if(respect_relative && ripser.is_relative_vertex(i)) {
				continue;
			}
			value_t diam = ripser.compute_diameter(i, 0);
			simplices.push_back(index_diameter_t(i, diam));
		}
		for(index_t i = 1; i <= dim; i++) {
			simplex_coboundary_enumerator cofacets(ripser);
			for(index_diameter_t& simplex : simplices) {
				cofacets.set_simplex(simplex, i - 1);
				while(cofacets.has_next(false)) {
					index_diameter_t cofacet = cofacets.next();
					if(get_diameter(cofacet) <= ripser.threshold) {
						next_simplices.push_back(cofacet);
					}
				}
			}
			simplices.swap(next_simplices);
			next_simplices.clear();
		}
	}
}


/* **************************************************************************
 * Matrix Operations
 * *************************************************************************/

// Take a column (formal sum of simplices) and return the pivot element,
// i.e. the largest simplex that does not cancel out
// If and only if the column is zero, return (-1, -1)
template <typename Column>
index_diameter_t pop_pivot(ripser& ripser, Column& column) {
	index_diameter_t pivot(-1, -1);
	while(!column.empty()) {
		pivot = column.top();
		ripser.get_current_reduction_record().pop_count++;
		column.pop();
		if(column.empty()) {
			return pivot;
		}
		if(get_index(column.top()) != get_index(pivot)) { // different index on top
			return pivot;
		} else { // same index on top
			// The same row index appears twice, they sum up to 0 (in F2), continue
			ripser.get_current_reduction_record().pop_count++;
			column.pop();
		}
	}
	return index_diameter_t(-1, -1);
}

// Note: May reduce to size of column by popping 'canceling' pivots but only
// replacing one
template <typename Column>
index_diameter_t get_pivot(ripser& ripser, Column& column) {
	index_diameter_t pivot = pop_pivot(ripser, column);
	if(get_index(pivot) != -1) {
		ripser.get_current_reduction_record().push_count++;
		column.push(pivot); // push back the popped pivot
	}
	return pivot;
}

// Add the coboundary column of simplex to the coboundary column working_coboundary
// Add simplex to the working_reduction_column (matrix V)
template <typename Column>
void add_simplex_coboundary(ripser& ripser,
                            const index_diameter_t simplex,
                            const index_t dim,
                            Column& working_reduction_column,
                            Column& working_coboundary) {
	simplex_coboundary_enumerator cofacets(ripser);
	working_reduction_column.push(simplex);
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		// Threshold check
		if(get_diameter(cofacet) <= ripser.threshold) {
			ripser.get_current_reduction_record().push_count++;
			working_coboundary.push(cofacet);
		}
	}
}

template <typename Column>
void add_coboundary(ripser& ripser,
                    compressed_sparse_matrix& reduction_matrix,
                    const std::vector<index_diameter_t>& columns_to_reduce,
                    const size_t index_column_to_add,
                    const size_t dim,
                    Column& working_reduction_column,
                    Column& working_coboundary) {
	//TODO(seb): Do we need the correct diameter here?
	index_diameter_t column_to_add(columns_to_reduce.at(index_column_to_add));
	ripser.get_current_reduction_record().add_simplex_boundary_count++;
	add_simplex_coboundary(ripser,
	                       column_to_add,
	                       dim,
	                       working_reduction_column,
	                       working_coboundary);
	// Computation of R_j due to implicit reduction
	for(index_t i = reduction_matrix.column_start(index_column_to_add);
	    i < reduction_matrix.column_end(index_column_to_add);
	    ++i) {
		index_diameter_t simplex = reduction_matrix.get_entry(i);
		ripser.get_current_reduction_record().add_simplex_boundary_count++;
		add_simplex_coboundary(ripser,
		                       simplex,
		                       dim,
		                       working_reduction_column,
		                       working_coboundary);
	}
}

template <typename Column>
void add_simplex_boundary(ripser &ripser,
                          const index_diameter_t simplex,
                          const index_t dim,
                          Column& working_reduction_column,
                          Column& working_boundary) {
	simplex_boundary_enumerator facets(ripser);
	working_reduction_column.push(simplex);
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		index_diameter_t facet = facets.next();
		// Threshold check
		if(get_diameter(facet) <= ripser.threshold) {
			ripser.get_current_reduction_record().push_count++;
			working_boundary.push(facet);
		}
	}
}

template <typename Column>
void add_boundary(ripser& ripser,
                  compressed_sparse_matrix& reduction_matrix,
                  const std::vector<index_diameter_t>& columns_to_reduce,
                  const size_t index_column_to_add,
                  const size_t dim,
                  Column& working_reduction_column,
                  Column& working_boundary) {
	//TODO(seb): Do we need the correct diameter here?
	index_diameter_t column_to_add(columns_to_reduce.at(index_column_to_add));
	// Computation of R_j due to implicit reduction
	ripser.get_current_reduction_record().add_simplex_boundary_count++;
	add_simplex_boundary(ripser,
	                     column_to_add,
	                     dim,
	                     working_reduction_column,
	                     working_boundary);
	for(index_t i = reduction_matrix.column_start(index_column_to_add);
	    i < reduction_matrix.column_end(index_column_to_add);
	    ++i)
	{
		index_diameter_t simplex = reduction_matrix.get_entry(i);
		ripser.get_current_reduction_record().add_simplex_boundary_count++;
		add_simplex_boundary(ripser,
		                     simplex,
		                     dim,
		                     working_reduction_column,
		                     working_boundary);
	}
}


/* **************************************************************************
 * Union Find
 * *************************************************************************/

struct union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;

	union_find(const index_t n) : parent(n), rank(n, 0) {
		for (index_t i = 0; i < n; ++i) {
			parent[i] = i;
		}
	}

	index_t find(index_t x) {
		index_t y = x, z;
		while((z = parent[y]) != y) {
		  y = z;
		}
		while((z = parent[x]) != y) {
			parent[x] = y;
			x = z;
		}
		return z;
	}

	void link(index_t x, index_t y) {
		if((x = find(x)) == (y = find(y))) {
		  return;
		}
		if(rank[x] > rank[y]) {
			parent[y] = x;
		} else {
			parent[x] = y;
			if (rank[x] == rank[y]) {
			  ++rank[y];
			}
		}
	}
};

std::vector<index_diameter_t> get_edges(ripser& ripser) {
	std::vector<index_diameter_t> edges;
	std::vector<index_t> vertices(2);
	for(index_t index = ripser.binomial_coeff(ripser.n, 2); index-- > 0;) {
		ripser.get_simplex_vertices(index, 1, ripser.dist.size(), vertices);
		value_t length = ripser.dist(vertices[0], vertices[1]);
		// Threshold check
		if(length <= ripser.threshold) {
			edges.push_back(std::make_pair(index, length));
		}
	}
	return edges;
}


/* **************************************************************************
 * Apparent Pairs
 * *************************************************************************/

index_diameter_t get_zero_pivot_facet(ripser& ripser, index_diameter_t simplex, index_t dim) {
	simplex_boundary_enumerator facets(ripser);
	facets.set_simplex(simplex, dim);
	while(facets.has_next()) {
		ripser.get_current_reduction_record().app_facet_count++;
		index_diameter_t facet = facets.next();
		if(get_diameter(facet) == get_diameter(simplex)) {
			// Relative check
			if(!ripser.is_relative_simplex(get_index(facet), dim - 1)) {
				return facet;
			}
		}
	}
	return index_diameter_t(-1, 0);
}

index_diameter_t get_zero_pivot_cofacet(ripser& ripser, index_diameter_t simplex, index_t dim) {
	simplex_coboundary_enumerator cofacets(ripser);
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		ripser.get_current_reduction_record().app_cofacet_count++;
		index_diameter_t cofacet = cofacets.next();
		if(get_diameter(cofacet) == get_diameter(simplex)) {
			return cofacet;
		}
	}
	return index_diameter_t(-1, 0);
}

index_diameter_t get_zero_apparent_facet(ripser& ripser, index_diameter_t simplex, index_t dim) {
	index_diameter_t facet = get_zero_pivot_facet(ripser, simplex, dim);
	return ((get_index(facet) != -1) &&
	        (get_index(get_zero_pivot_cofacet(ripser, facet, dim - 1)) == get_index(simplex)))
	           ? facet
	           : index_diameter_t(-1, 0);
}

index_diameter_t get_zero_apparent_cofacet(ripser& ripser, index_diameter_t simplex, index_t dim) {
	index_diameter_t cofacet = get_zero_pivot_cofacet(ripser, simplex, dim);
	return ((get_index(cofacet) != -1) &&
	        (get_index(get_zero_pivot_facet(ripser, cofacet, dim + 1)) == get_index(simplex)))
	           ? cofacet
	           : index_diameter_t(-1, 0);
}

bool is_in_zero_apparent_pair(ripser& ripser, index_diameter_t simplex, index_t dim) {
	if(dim == 0) {
		return get_index(get_zero_apparent_cofacet(ripser, simplex, dim)) != -1;
	//} else if(dim == ripser.dim_max) {
	//	return (get_index(get_zero_apparent_facet(ripser, simplex, dim)) != -1);
	} else {
		return (get_index(get_zero_apparent_facet(ripser, simplex, dim)) != -1) ||
		       (get_index(get_zero_apparent_cofacet(ripser, simplex, dim)) != -1);
	}
}

// https://stackoverflow.com/questions/59817055/resident-set-size-remains-the-same-after-delete
long get_memory_usage() {
	long vmSize = 0;
	std::ifstream buffer("/proc/self/statm");
	buffer >> vmSize;
	buffer.close();
	return vmSize;
}

#endif
