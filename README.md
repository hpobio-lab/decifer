# decifer

DeCiFer is a method for copy-aware clustering of somatic single nucleotide variants by the descendant cell fraction.

![Overview of DeCiFer](doc/overview.png)

## Contents

  1. [Compilation instructions](#compilation)
     * [Dependencies](#dep)
     * [Compilation](#comp)
  2. [Usage instructions](#usage)
     * [I/O format](#io)
     * [DeCiFer](#dcf)

<a name="compilation"></a>
## Compilation instructions

<a name="dep"></a>
### Dependencies

DeCiFer is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang). In addition, DeCiFer has the following dependencies.

* [CMake](http://www.cmake.org/) (>= 3.0)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* A recent ILP Solver:
  * [Gurobi](http://www.gurobi.com) (>= 6.0)
  * [CPLEX](http://www.gurobi.com) (>= 12.0)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting state tree DOT files, but is not required for compilation.
In case [doxygen](http://www.stack.nl/~dimitri/doxygen/) is available, extended source code documentation may be generated.

<a name="comp"></a>
### Compilation

To compile DeCiFer, execute the following commands from the root of the repository:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

By default DeCiFer will use the CPLEX ILP solver. To use Gurobi instead of CPLEX,  run the following command with adjusted paths:

    $ cmake -DLIBLEMON_ROOT=~/lemon -DCPLEX=OFF \
    -DGUROBI_INCLUDE_DIR=/usr/local/gurobi702/linux64/include \
    -DGUROBI_CPP_LIB=/usr/local/gurobi702/linux64/lib/libgurobi_c++.a \
    -DGUROBI_LIB=/usr/local/lib/libgurobi70.so ..

The compilation results in the following files in the `build` directory:

EXECUTABLE | DESCRIPTION
-----------|-------------
`decifer` | Infer evolutionary history and cellular prevalence of somatic single-nucleotide variants using a combinatorial algorithm based on the descendant cell fraction.
`summarize` | Extracts hard clustering from DeCiFer soft clustering.
`visualizepoststatetree` | Visualizes a state tree identified by DeCiFer.
