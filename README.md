# Ripser EDU

This is a modified version of [Ripser](https://github.com/Ripser/ripser) (by Ulrich Bauer) with the aim to educate. Refer to [Ripser's README](https://github.com/Ripser/ripser/blob/master/README.md) for general information.

The current version of Ripser includes a large number of optimizations and features in a single source file. This comes at the cost of source code readability and makes it more difficult to connect theory and implementation.

In this repository the original code is restructured, simplified and feature-reduced in the interest of education. More specifically:
- Restructured source code separate core functions and data types (in [ripser\_core.hpp](./ripser_core.hpp)) from specific implementations of the reduction algorithm
- Simplified types and only a bare minimum of templates
- Many added comments
- Implementation of persistent homology and cohomology **without** optimizations
- Optimizations such as clearing, emergent and apparent pairs are each moved to a separate file

The following features of Ripser have been **removed**:
- Support of coefficients and computation over prime fields
- Support of a multitude of input formats (only `lower_distance_matrix` remains)
- Progress indication
- Support for robinhoodhash

The following features have been **added**:
- The barcode (main output of Ripser) is not printed to stdout during the reduction but stored for further analysis and printing after completion

Despite the numerous structural changes, this version is in close correspondence with the original version. In particular data type, function and variable names have been (mostly) left unchanged. All implementations of the reduction use an implicit representation of the (co)boundary matrix. Most functions as well as the reduction algorithm are very close to the original if not unchanged.

# Outline

The original Ripser code is split into three types of files.

### 1. Core Functionality: ripser\_core

[ripser\_core.hpp](./ripser_core.hpp) contains infrastructure code that is used for all versions of the reduction algorithm. It contains
- The definition of most data types including the `ripser` class and `[co]boundary_enumerator`
- All functionality regarding the binomial numbering system
- All functionality regarding the enumeration of (co)faces
- The implementation of all basic matrix operations used in the reduction algorithm (column addition, pivot extraction)
- Functionality related to union-find

### 2. Implementation of the Reduction Algorithm: ripser\_*.cpp

The implementation of the reduction algorithm and computation of the persistence barcode is split into several files. Each implements the following functions separately depending on the particular choice of (co)homology and possibly a single optimization:
- `init_(co)boundary_and_get_pivot`
- `assemble_columns_to_reduce`
- `compute_pairs` (the main reduction algorithm)
- `compute_barcodes`
- `main`

The following implementations are currently available:
- Basic homology (in decreasing dimension)
- Basic cohomology (in increasing dimension)
- (Co)homology with clearing
- Cohomology with the emergent pairs optimization
- Cohomology with the apparent pairs optimization

Note, that all implementations use an implicit representation of the (co)boundary matrix (to keep memory-efficiency and a close correspondence with the original version).

### 3. Printing: print\_utils.hpp

[print\_utils.hpp](./print_utils.hpp) contains all code concerning printing to stdout such as printing simplices, columns, barcodes and auxiliary information


# License

Copyright © 2015–2021 [Ulrich Bauer]

Ripser is licensed under the [MIT] license (`COPYING.txt`), with an extra clause (`CONTRIBUTING.txt`) clarifying the license for modifications released without an explicit written license agreement.  Please contact the author if you want to use Ripser in your software under a different license.
