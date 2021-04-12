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

typedef int64_t index_t;
typedef float value_t;

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
public:
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

typedef std::unordered_map<index_t, size_t, index_hash, equal_index> entry_hash_map;


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

// Take a column (formal sum of simplices) and return the pivot element,
// i.e. the largest simplex that does not cancel out
// If and only if the column is zero, return (-1, -1)
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
index_diameter_t get_pivot(Column& column) {
	index_diameter_t pivot = pop_pivot(column);
	if(get_index(pivot) != -1) {
		column.push(pivot); // push back the popped pivot
	}
	return pivot;
}

// Take a dim-simplex as input and assemble to coboundary matrix column
// corresponding to that simplex. Return the pivot element of that column
index_diameter_t init_coboundary_and_get_pivot(ripser &ripser,
                                               const index_diameter_t simplex,
                                               const index_t dim,
                                               Column& working_coboundary) {
	simplex_coboundary_enumerator cofacets(ripser);
	// std::vector<index_diameter_t> cofacet_entries;
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

// Add the coboundary column of simplex to the coboundary column working_coboundary
// Add simplex to the working_reduction_column (matrix V)
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

// 
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


/* **************************************************************************
 * Core Functionality
 * *************************************************************************/

// Takes a vector of (dim-1)-simplices as input and sets columns_to_reduce
// to contain all cofacets of those simplices
// ordered column indices (indexed by dim-simplices) of the coboundary matrix
// Sets simplices to contain all dim-simplices
void assemble_columns_to_reduce(ripser &ripser,
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
			if(get_diameter(cofacet) <= ripser.threshold) {
				next_simplices.push_back(cofacet);
				columns_to_reduce.push_back(cofacet);
			}
		}
	}
	simplices.swap(next_simplices);
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              greater_diameter_or_smaller_index);
}

void print_column(ripser& ripser, Column &column, index_t dim);
void print_simplex(ripser& ripser, index_t simplex, index_t dim);
void print_simplices(ripser& ripser, std::vector<index_diameter_t>& simplices, index_t d);
void print_v(compressed_sparse_matrix& mat,
             std::vector<index_diameter_t> columns_to_reduce);
void print_mat(compressed_sparse_matrix& mat);

void compute_pairs(ripser &ripser,
                   const std::vector<index_diameter_t>& columns_to_reduce,
                   entry_hash_map& pivot_column_index,
                   const index_t dim) {
	std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
	compressed_sparse_matrix reduction_matrix; // V
	std::cout << "Num of col to red: " << columns_to_reduce.size() << std::endl;
	for(size_t j = 0; j < columns_to_reduce.size(); ++j) { // For j in J
		reduction_matrix.append_column();
		Column working_reduction_column; // V_j
		Column working_coboundary;       // R_j
		// Assemble the column j (corresponding to the simplex column_to_reduce)
		// and get the pivot of that column
		index_diameter_t column_to_reduce = columns_to_reduce.at(j);
		index_diameter_t pivot = init_coboundary_and_get_pivot(ripser,
		                                                       column_to_reduce,
		                                                       dim,
		                                                       working_coboundary);
		value_t birth = get_diameter(column_to_reduce);
		if(dim == 0 && birth == -INF) {
			birth = 0;
		}
		//print_column(ripser, working_coboundary, dim);
		// The reduction
		while(true) {
			// Check if the column is not (did not get reduced to) a 0-column
			if(get_index(pivot) != -1) {
				auto pair = pivot_column_index.find(get_index(pivot));
				// Check ?
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
				} else {
					value_t death = get_diameter(pivot);
					if(death > birth * ripser.ratio) {
						std::cout << " [" << birth << "," << death << ")" << std::endl;
					}
					pivot_column_index.insert({get_index(pivot), j});
					// Write V_j to V
					index_diameter_t e = pop_pivot(working_reduction_column);
					while(get_index(e) != -1) {
						reduction_matrix.push_back(e);
						e = pop_pivot(working_reduction_column);
					}
					break;
				}
			} else {
				// Zero column, persistent homology "pair"
				std::cout << " [" << birth << ", )" << std::endl;
				break;
			}
		}
		//print_mat(reduction_matrix);
		//print_v(reduction_matrix, columns_to_reduce);
		//std::cout << "-------------------------------------" << std::endl;
	}
}

void compute_barcodes(ripser& ripser) {
	std::vector<index_diameter_t> simplices;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simplices.push_back(index_diameter_t(i, diam));
	}
	for(index_t dim = 0; dim <= ripser.dim_max; ++dim) {
		std::vector<index_diameter_t> columns_to_reduce;
		entry_hash_map pivot_column_index;
		pivot_column_index.reserve(columns_to_reduce.size());
		if(dim == 0) {
			columns_to_reduce = std::vector<index_diameter_t>(simplices);
		} else {
			assemble_columns_to_reduce(ripser, simplices, columns_to_reduce, dim);
		}
		//std::vector<index_diameter_t> ctr;
		//for(auto c : columns_to_reduce) {
		//	if(get_index(c) == 1 || get_index(c) == 4 || get_index(c) == 0) {
		//		ctr.push_back(c);
		//	}
		//}
		compute_pairs(ripser, columns_to_reduce, pivot_column_index, dim);
	}
}


/* **************************************************************************
 * Debugging and Printing
 * *************************************************************************/

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
	std::cout << "info: list of all simplices (id :: diam :: vertices)" << std::endl;
	std::vector<index_diameter_t> simpl_prev;
	std::vector<index_diameter_t> simpl_curr;
	for(index_t i = 0; i < ripser.n; i++) {
		value_t diam = ripser.compute_diameter(i, 0);
		simpl_prev.push_back(index_diameter_t(i, diam));
	}
	std::cout << "  dim 0" << std::endl;
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
	offs += sprintf(offs + buf, "%i-", simplex);
	for(auto v : vertices) {
	  if(v == -1) continue;
	  offs += sprintf(offs + buf, "%i'", v);
	}
	return offs;
}

void print_simplex(ripser& ripser, index_t simplex, index_t dim) {
	char buf[1024]; buf[0] ='\0';
	sprint_simplex(buf, ripser, simplex, dim);
	printf("%s\n", buf);
}

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
	char buf[1024]; buf[0] = '\0';
	int offs = 0;
	int pad = 3;
	for(size_t row = 0; row < columns_to_reduce.size(); ++row) {
		index_t row_ele = get_index(columns_to_reduce.at(row));
		for(size_t col = 0; col < columns_to_reduce.size(); ++col) {
			if(col >= v.size()) {
				offs += sprint_pad(buf + offs, pad);
			} else {
				size_t count = 0;
				for(size_t i = v.column_start(col); i < v.column_end(col); ++i) {
					if(get_index(v.get_entry(i)) == row_ele) {
						count++;
					}
				}
				if(count % 2 == 1) {
					offs += sprint_element(buf + offs, 1, pad);
				} else {
					offs += sprint_pad(buf + offs, pad);
				}
			}
		}
		offs += sprintf(offs + buf, "\n");
	}
	printf("%s", buf);
}

void print_mat(compressed_sparse_matrix& mat) {
	char buf[1024]; buf[0] = '\0';
	int offs = 0;
	int pad = 3;
	size_t max_row_index = 0;
	for(size_t i = 0; i < mat.size(); ++i) {
		max_row_index = std::max(max_row_index,
		                         mat.column_end(i) - mat.column_start(i));
	}
//	for(auto e : mat.bounds) {
//		std::cout << e;
//	}
//	std::cout << std::endl;
//	for(auto e : mat.entries) {
//		std::cout << get_index(e) << " ";
//	}
	std::cout << std::endl;
	for(size_t row = 0; row < max_row_index; ++row) {
		for(size_t col = 0; col < mat.size(); ++col) {
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

void test_compressed_sparse_matrix_print() {
	compressed_sparse_matrix m;
	m.append_column();
	m.push_back(index_diameter_t(0, 0));
	m.push_back(index_diameter_t(1, 0));
	m.append_column();
	m.push_back(index_diameter_t(2, 0));
	m.push_back(index_diameter_t(3, 0));
	m.push_back(index_diameter_t(4, 0));
	m.append_column();
	m.push_back(index_diameter_t(5, 0));
	m.append_column();
	//print_matrix(m);
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
	std::cout << "info: value range: [" << min << "," << max_finite << "]" << std::endl;

	std::cout << "info: distance matrix with " << dist.size()
			  << " points, using threshold at enclosing radius " << enclosing_radius
			  << std::endl;
	
	ripser ripser(std::move(dist), dim_max, enclosing_radius, ratio);
	list_all_simplices(ripser);
	
	compute_barcodes(ripser);
	exit(0);
}
