#ifndef RIPSER_CORE
#define RIPSER_CORE

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>


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

	std::vector<size_t> bounds;
	std::vector<index_diameter_t> entries;

	compressed_sparse_matrix() {}

	size_t column_start(const index_t index) const {
		return index == 0 ? 0 : bounds.at(index - 1);
	}

	size_t column_end(const index_t index) const {
		return bounds.at(index);
	}

	index_diameter_t get_entry(const size_t index) {
		return entries.at(index);
	}

	index_diameter_t get_entry(const size_t row_index, const size_t col_index) {
		return entries.at(column_start(col_index) + row_index);
	}

	bool has_entry_at(const size_t row_index, const size_t col_index) {
		if(bounds.size() <= col_index) {
			return false;
		}
		for(size_t i = column_start(col_index); i < column_end(col_index); ++i) {
			if(get_index(entries.at(i)) == row_index) {
				return true;
			}
		}
		return false;
	}

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

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}
	return compressed_lower_distance_matrix(std::move(distances));
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

// hash map with key = pivot index and value = (column index, basis element)
// the column index is relative to the current dimension since we only have
// a slice of V w.r.t to that dimension
// We also store the basis element corresponding to V_j to reconstruct to
// persistence pair in ripser_hom (needed in dim+1 where the column index is no
// longer meaningful)
typedef std::unordered_map<index_t,
                           size_t,
                           index_hash,
                           equal_index> entry_hash_map;


/* **************************************************************************
 * Ripser
 * *************************************************************************/

struct barcode {
	index_t dim;
	std::vector<std::pair<value_t, value_t>> persistence_intervals;
	size_t clearing_count;
	size_t emergent_count;
	size_t apparent_count;

	barcode(index_t _dim)
		: dim(_dim), persistence_intervals(),
		  clearing_count(0),
		  emergent_count(0),
		  apparent_count(0)
	{ }

	void add_interval(value_t birth, value_t death) {
		persistence_intervals.push_back(std::make_pair(birth, death));
	}
};

bool barcode_order(const barcode& a, const barcode& b) {
	return a.dim < b.dim;
}

bool persistence_interval_order(std::pair<value_t, value_t>& a,
                                std::pair<value_t, value_t>& b) {
	return (a.first < b.first) || (a.first == b.first && a.second < b.second);
}

struct ripser {

	const DistanceMatrix dist;
	const index_t n;
	const index_t dim_max;
	const value_t threshold;
	const float ratio;
	const binomial_coeff_table binomial_coeff;
	
	// Output
	mutable std::vector<barcode> barcodes;

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


/* **************************************************************************
 * Matrix Operations
 * *************************************************************************/

// Take a column (formal sum of simplices) and return the pivot element,
// i.e. the largest simplex that does not cancel out
// If and only if the column is zero, return (-1, -1)
template <typename Column>
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
		  // The same row index appears twice, they sum up to 0 (in F2), continue
		  column.pop();
		}
	}
	return index_diameter_t(-1, -1);
}

// Note: May reduce to size of column by popping 'canceling' pivots but only
// replacing one
template <typename Column>
index_diameter_t get_pivot(Column& column) {
	index_diameter_t pivot = pop_pivot(column);
	if(get_index(pivot) != -1) {
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
	add_simplex_coboundary(ripser,
	                       column_to_add,
	                       dim,
	                       working_reduction_column,
	                       working_coboundary);
	// Computation of R_j due to implicit reduction
	for(size_t i = reduction_matrix.column_start(index_column_to_add);
	    i < reduction_matrix.column_end(index_column_to_add);
	    ++i) {
		index_diameter_t simplex = reduction_matrix.get_entry(i);
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
		if(get_diameter(facet) <= ripser.threshold) {
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
	add_simplex_boundary(ripser,
	                     column_to_add,
	                     dim,
	                     working_reduction_column,
	                     working_boundary);
	// Computation of R_j due to implicit reduction
	for(size_t i = reduction_matrix.column_start(index_column_to_add);
	    i < reduction_matrix.column_end(index_column_to_add);
	    ++i) {
		index_diameter_t simplex = reduction_matrix.get_entry(i);
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
		for (index_t i = 0; i < n; ++i) parent[i] = i;
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
 * Other
 * *************************************************************************/

value_t compute_enclosing_radius(DistanceMatrix& dist) {
	// Compute enclosing radius and distance bounds
	value_t min = INF;
	value_t max = -INF;
	value_t max_finite = max;
	value_t enclosing_radius = INF;

	for(size_t i = 0; i < dist.size(); ++i) {
		value_t r_i = -INF;
		for (size_t j = 0; j < dist.size(); ++j) {
			r_i = std::max(r_i, dist(i, j));
		}
		enclosing_radius = std::min(enclosing_radius, r_i);
	}
	for(auto d : dist.distances) {
		min = std::min(min, d);
		max = std::max(max, d);
		if(d != INF) {
			max_finite = std::max(max_finite, d);
		}
	}
	std::cout << "info: value range: [" << min << "," << max_finite << "]" << std::endl;
	std::cout << "info: distance matrix with " << dist.size()
			  << " points, using threshold at enclosing radius " << enclosing_radius
			  << std::endl;
	return enclosing_radius;
}

#endif
