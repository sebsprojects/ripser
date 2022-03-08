This branch of (my fork of) Ripser contains all implemenetations, configs, scripts and most data set developed and used as part of my master thesis project.

## Compile and Run

All `ripser_*` source files can be compiled with `make`. See `Makefile` for all available builds. All versions are run as `./ripser_* path/to/config.ini`.

## Configuration: ripser_config.cpp

Instead of passing program parameters to Ripser directly, all versions expect a single parameter containing the path to `.ini` configuration file. The first line of this file is expected to contain `[ripser]`. Each subsequent line should comprise a single key-value pair formatted `key=value`. Empty lines, lines starting with a `#` and keys without value (e.g. `key=` are ignored. The available configuration options are the following:
- `input_path` (string) specifies the path to the data set relative to the current working directory.
- `input_type` (string) specifies whether the data set is a `point_cloud`, `lower_distance_matrix` or `full_distance_matrix`.
- `output_path` (string) specifies the output directory to which all output files are written.
- `dim_max` (non-negative integer) specifies the maximum considered dimension with respect to the (co)boundary matrix reduction. If a homology reduction is performed this means that the barcode decomposition up to dim_max - 1 is computed.
- `ratio` (non-negative float) specifies a cut-off ratio of the diameters of birth and death simplex. Non-essential pairs below this ratio are omitted from the output. If 0 is specified, zero-persistence pairs are part of the output (if not omitted due to the apparent pairs shortcut).
- `threshold` (float) specifies the diameter threshold at which Vietoris-Rips filtration is cutoff. If omitted or negative the full Vietoris-Rips filtration is considered.
- `use_enclosing_threshold` (bool) specifies whether the threshold is to be the enclosing threshold. If set to true (or 1) this overwrites the parameter `threshold`.
- `use_union_find` (bool) specifies whether the reduction in degree 0 should be replaced by a connected components computation. This parameter only has an effect if the version at hand implements this optimization.
- `print_progress` (bool) specifies whether the program prints progress updates to stdout.
- `absolute_subcomplex` (range) specifies which parts (points / submatrix) of the input data are considered for the computation. If omitted, the whole data set is used.
- `relative_subcomplex` (range) specifies which parts (point / submatrix) of the input data are used to induce the relative subcomplex. This parameter only has an effect if the version implements relative (co)homology. If omitted, the relative subcomplex is the empty complex (corresponding to absolute (co)homology).

Values of the type string must be unquoted. The type bool has values true, false (alternatively 1 and 0). Ranges are a comma-separated list of closed intervals specified by bounds (both bounds are included). For instance `0-99,200-299` specifies a 200 element range.

The parameters `input_path`, `inpute_type` and `dim_max` must be specified. All other parameters are optional.

## Ripser's Core Functionality: ripser_core.hpp

[ripser_core.hpp](./ripser_core.hpp) contains infrastructure code that is used for all versions of the reduction algorithm. It contains key types and functions implementing
- Data types,
- The combinatoric numbering system,
- Facet and cofacet enumeration,
- Matrix operations required by the reduction algorithm,
- Functionality for union-find, apparent pairs and relative (co)homology,
- Input processing,
- Output assembly and storage,
- Performance metric collection.

Most code concerning the functionality of Ripser is very close to that of Ripser v.1.2.1. However, many functions have been separated from the `ripser` class and templates have been reduced to a minimum.


## Reduction Algorithm: ripser_*.cpp

Each file contains an implementation of the reduction algorithm comprising the functions
- `init_(co)boundary_and_get_pivot`
- `assemble_columns_to_reduce`
- `compute_pairs` (the main reduction algorithm)
- `compute_barcodes`

All of the following versions implement a degreewise reduction, either in increasing or decreasing degree.

- `ripser_hom.cpp` implements the reduction in homology with increasing degree without any further optimizations
- `ripser_hom_optimized.cpp` implements the reduction in homology with increasing degree using the emergent and apparent pairs shortcut (no clearing)
- `ripser_hom_clearing.cpp` implements the in homology with decreasing degree using clearing
- `ripser_hom_clearing_optimized.cpp` adds the emergent and apparent pairs shortcut
- `ripser_cohom.cpp` implements the reduction in cohomology with decreasing degree without any further optimizations
- `ripser_cohom_clearing.cpp` implements the reduction in cohomolgoy with increasing degree using clearing
- `ripser_cohom_clearing_optimized.cpp` adds the emergent and apparent pairs shortcut
- `ripser_cohom_rephom.cpp` implements a secondary homology reduction to compute homology representatives
- `ripser_cohom_repinverse.cpp` implements matrix inversion to compute homology representatives
- `ripser_cohom_rel.cpp` implements relative cohomology with increasing degree using clearing and the emergent pairs shortcut (but not the apparent pairs shortcut)
- `ripser_cohom_rel_optimized.cpp` adds the apparent pairs shortcut
- `ripser_cohom_rel_rephom.cpp` implements a secondary reduction in (relative) homology to compute relative homology representatives


## Printing and Output: print_utils.hpp

[print_utils.hpp](./print_utils.hpp) contains all code concerning printing and output. The primary output functionality comprises the following functions:
- `output_barcode` prints the barcode in all (or a single specified) degrees, optionally with the associated (co)homology representative.
- `output_info` prints coarse information about the computation, for each degree separately. This includes the elapsed runtime for assembly, reduction and representative duration as well as various counts, such as the number of cleared columns, apparent pairs and more.
- `output_config` prints the configuration used for the computation.


## License

This fork of Ripser is licensed under the MIT license (`COPYING.txt`), with an extra clause (`CONTRIBUTING.txt`) clarifying the license for modifications released without an explicit written license agreement. Please contact the author if you want to use Ripser in your software under a different license.
