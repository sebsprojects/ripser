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

void check_overflow(index_t i) {
	if (i < 0)
		throw std::overflow_error("simplex index " +
								  std::to_string((uint64_t)i) +
		                          " in filtration is out of bounds");
}

// The basic data-type for simplices, given by the index of the simplex
// w.r.t. binomial numbering system and its (cached) diameter
typedef std::pair<index_t, value_t> index_diameter_t;

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
		// TODO: Unsafe cases?
		assert((size_t) n < B.size() && (size_t) k < B[n].size() && n >= k - 1);
		return B[n][k];
	}
};

// Type to store the input data as a distance matrix
struct compressed_lower_distance_matrix {
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	compressed_lower_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
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
typedef struct compressed_lower_distance_matrix DistanceMatrix ;

/* **************************************************************************
 * Ripser
 * *************************************************************************/

class ripser {
	const DistanceMatrix dist;
	const index_t n;
	const index_t dim_max;
	const value_t threshold;
	const float ratio;
	const binomial_coeff_table binomial_coeff;

	mutable std::vector<index_t> vertices;

public:
	ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold, float _ratio)
		: dist(std::move(_dist)), n(dist.size()),
			   dim_max(std::min(_dim_max, index_t(dist.size() - 2))),
			   threshold(_threshold), ratio(_ratio),
			   binomial_coeff(n, dim_max + 2)
	{ }
};


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
	if(argc == 1) {
		filename = argv[1];
	} else {
		std::cout << "specify path to lower-distance matrix file as only arg"
				  << std::endl;
	}

	// Parameter setup
	// value_t threshold = std::numeric_limits<value_t>::max();
	index_t dim_max = 1;
	float ratio = 1;

	// Reading the distance matrix from file
	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}
	DistanceMatrix dist = read_lower_distance_matrix(file_stream);

	// Compute enclosing radius and distance bounds
	value_t min = std::numeric_limits<value_t>::infinity(),
			max = -std::numeric_limits<value_t>::infinity(),
			max_finite = max;

	value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
	for (size_t i = 0; i < dist.size(); ++i) {
		value_t r_i = -std::numeric_limits<value_t>::infinity();
		for (size_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
		enclosing_radius = std::min(enclosing_radius, r_i);
	}

	for (auto d : dist.distances) {
		min = std::min(min, d);
		max = std::max(max, d);
		if (d != std::numeric_limits<value_t>::infinity()) max_finite = std::max(max_finite, d);
	}
	std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

	std::cout << "distance matrix with " << dist.size()
			  << " points, using threshold at enclosing radius " << enclosing_radius
			  << std::endl;

	// Ripser entry point
	ripser r(std::move(dist), dim_max, enclosing_radius, ratio);
	exit(0);
}
