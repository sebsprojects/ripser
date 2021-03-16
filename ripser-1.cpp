#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>


/* **************************************************************************
 * Types and Utility
 * *************************************************************************/

typedef float value_t;
typedef int64_t index_t;

const value_t INF = std::numeric_limits<value_t>::infinity();

void check_overflow(index_t i) {
	if (i < 0)
		throw std::overflow_error("simplex index " +
			std::to_string((uint64_t)i) +
			" in filtration is out of bounds");
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

bool greater_diameter_or_smaller_index(const index_diameter_t& a,
                                       const index_diameter_t& b) {
	return (get_diameter(a) > get_diameter(b)) ||
	       ((get_diameter(a) == get_diameter(b)) &&
	        (get_index(a) < get_index(b)));
}

class greater_diameter_or_smaller_index_comp {

public:
  bool operator()(const index_diameter_t& a, const index_diameter_t& b) {
	return greater_diameter_or_smaller_index(a, b);
  }
};

// A working column is represented by this type. The entries are ordered with
// respect to reverse filtration order
typedef std::priority_queue<index_diameter_t,
                            std::vector<index_diameter_t>,
                            greater_diameter_or_smaller_index_comp> Column;

// Binomial lookup table
class binomial_coeff_table {
	std::vector<std::vector<index_t>> B;

public:
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

// Type to store the input data as a distance matrix
class compressed_lower_distance_matrix {

public:
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

class compressed_sparse_matrix {

private:
	std::vector<size_t> bounds;
	std::vector<index_diameter_t> entries;

public:
	size_t size() const {
		return bounds.size();
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
 * Ripser
 * *************************************************************************/

class ripser {

public:
	const DistanceMatrix dist;
	const index_t n;
	const index_t dim_max;
	const value_t threshold;
	const float ratio;
	const binomial_coeff_table binomial_coeff;
	
	ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold, float _ratio)
		: dist(std::move(_dist)), n(dist.size()),
		dim_max(std::min(_dim_max, index_t(dist.size() - 2))),
		threshold(_threshold), ratio(_ratio),
		binomial_coeff(n, dim_max + 2)
	{ }
	
	// For a k-simplex idx, return the vertex (of idx) of maximal index
	// n is an upper bound for the search
	index_t get_max_vertex(const index_t idx, const index_t k, const index_t n) const {
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

};


/* **************************************************************************
 * Enumerating simplices
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

	bool has_next() {
		return k >= 0;
	}

	index_diameter_t next() {
		j = parent.get_max_vertex(idx_below, k + 1, j);
		std::cout << "j=" << j << " idx_above=" << idx_above
		          << " binom=" << parent.binomial_coeff(j, k + 1)
		          << " idx_below=" << idx_below << std::endl;
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


/* **************************************************************************
 * Matrix Reduction Operations
 * *************************************************************************/

// If the same row index appears twice, they sum up to 0 (in F2) and we continue
index_diameter_t pop_pivot(Column& column) {
	index_diameter_t pivot(-1, -1);
	while(!column.empty()) {
		pivot = column.top();
		column.pop();
		if(column.empty()) {
			return pivot;
		}
		if(get_index(column.top()) != get_index(pivot)) { // different index on top
			return pivot;
		} else { // same index on top
		  column.pop();
		}
	}
	return index_diameter_t(-1, -1);
}

// Note: May reduce to size of column by popping 'canceling' pivots but only
// replacing one
index_diameter_t get_pivot(Column& column) {
	index_diameter_t pivot = pop_pivot(column);
	if(get_index(pivot) != -1) {
		column.push(pivot); // push back the popped pivot
	}
	return pivot;
}

index_diameter_t init_coboundary_and_get_pivot(ripser &ripser,
                                               const index_diameter_t simplex,
                                               const index_t dim,
                                               Column working_coboundary) {
	simplex_coboundary_enumerator cofacets(ripser);
	std::vector<index_diameter_t> cofacet_entries;
	cofacets.set_simplex(simplex, dim);
	while(cofacets.has_next()) {
		index_diameter_t cofacet = cofacets.next();
		//TODO(seb): Check diam <= threshold
		//cofacet_entries.push_back(cofacet);
	//}
	//for(index_diameter_t cofacet : cofacet_entries) {
		working_coboundary.push(cofacet);
	}
	return get_pivot(working_coboundary);
}

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
		//TODO(seb): Check diam <= threshold
		working_coboundary.push(cofacet);
	}
}

// TODO(seb): param order should be uniform
void add_coboundary(ripser& ripser,
                    compressed_sparse_matrix& reduction_matrix,
                    Column& working_reduction_column,
                    Column& working_coboundary,
                    const std::vector<index_diameter_t>& columns_to_reduce,
                    const size_t index_column_to_add,
                    const size_t dim) {
	//TODO(seb): Do we need the correct diameter here?
	index_diameter_t column_to_add(columns_to_reduce.at(index_column_to_add), -1);
	add_simplex_coboundary(ripser, column_to_add, dim, working_reduction_column,
	                       working_coboundary);
	for(index_diameter_t simplex : reduction_matrix.subrange()) {
		add_simplex_coboundary(ripser, simplex, dim, working_reduction_column,
                               working_coboundary);
	}
}


/* **************************************************************************
 * Core Functionality
 * *************************************************************************/

// Takes a vector of dim-simplices as input and
// returns the ordered columns of the coboundary matrix
// Updates simplices to be TODO
void assemble_coboundaries_to_reduce(ripser &ripser,
                                     std::vector<index_diameter_t>& simplices,
                                     std::vector<index_diameter_t>& columns_to_reduce,
                                     const index_t dim) {
	columns_to_reduce.clear();
	std::vector<index_diameter_t> next_simplices;
	simplex_coboundary_enumerator cofacets(ripser);
	for(index_diameter_t& simplex : simplices) {
		cofacets.set_simplex(simplex, dim - 1);
		while(cofacets.has_next(false)) {
			//TODO(seb): Check diam <= threshold?
			index_diameter_t cofacet = cofacets.next();
			next_simplices.push_back(cofacet);
			columns_to_reduce.push_back(cofacet);
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              greater_diameter_or_smaller_index);
}

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   const index_t dim) {

}

void compute_barcodes() {
	std::vector<index_diameter_t> simplices;
	std::vector<index_diameter_t> columns_to_reduce;
}


/* **************************************************************************
 * Debugging and Printing
 * *************************************************************************/

void print_simplices(ripser& ripser, std::vector<index_diameter_t>& simplices, index_t d) {
  std::vector<index_t> vertices(ripser.n, -1);
  for(auto s : simplices) {
		ripser.get_simplex_vertices(get_index(s), d, ripser.n, vertices);
		std::cout << "idx=" << get_index(s) << " diam=" << get_diameter(s) << " vs=( ";
		for(auto i : vertices) {
			if(i >= 0) {
			  std::cout << i << " ";
			}
		}
		std::cout << ")" << std::endl;
		std::fill(vertices.begin(), vertices.end(), -1);
  }
}

void list_all_simplices(ripser& ripser) {
	std::vector<index_diameter_t> simpl_prev;
	std::vector<index_diameter_t> simpl_curr;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simpl_prev.push_back(index_diameter_t(i, diam));
	}
	print_simplices(ripser, simpl_prev, 0);
	std::cout << std::endl;
	
	simplex_coboundary_enumerator e(ripser);
	index_t dim = 1;
	for(auto i : simpl_prev) {
		e.set_simplex(i, dim - 1);
		while(e.has_next(false)) {
		  simpl_curr.push_back(e.next());
		}
	}
	print_simplices(ripser, simpl_curr, dim);
	std::cout << std::endl;
	
	simpl_prev.swap(simpl_curr);
	simpl_curr.clear();
	dim = 2;
	for(auto i : simpl_prev) {
		e.set_simplex(i, dim - 1);
		while(e.has_next(false)) {
		  simpl_curr.push_back(e.next());
		}
	}
	print_simplices(ripser, simpl_curr, dim);
}


/* **************************************************************************
 * Input Parsing and Main
 * *************************************************************************/

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}
	return compressed_lower_distance_matrix(std::move(distances));
}


int main(int argc, char** argv) {
	const char* filename = nullptr;
	if(argc == 2) {
		filename = argv[1];
	} else {
		std::cerr << "error: specify path to lower-distance matrix file as only arg"
				  << std::endl;
		exit(-1);
	}

	// Parameter setup
	// value_t threshold = std::numeric_limits<value_t>::max();
	index_t dim_max = 1;
	float ratio = 1;

	// Reading the distance matrix from file
	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "error: couldn't open file " << filename << std::endl;
		exit(-1);
	}
	DistanceMatrix dist = read_lower_distance_matrix(file_stream);

	// Compute enclosing radius and distance bounds
	value_t min = INF,
			max = -INF,
			max_finite = max;

	value_t enclosing_radius = INF;
	for (size_t i = 0; i < dist.size(); ++i) {
		value_t r_i = -INF;
		for (size_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
		enclosing_radius = std::min(enclosing_radius, r_i);
	}

	for (auto d : dist.distances) {
		min = std::min(min, d);
		max = std::max(max, d);
		if (d != INF) max_finite = std::max(max_finite, d);
	}
	std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

	std::cout << "distance matrix with " << dist.size()
			  << " points, using threshold at enclosing radius " << enclosing_radius
			  << std::endl;

	// Ripser entry point
	ripser r(std::move(dist), dim_max, enclosing_radius, ratio);
	list_all_simplices(r);
	exit(0);
}
