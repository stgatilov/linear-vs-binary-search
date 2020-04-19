# linear-vs-binary-search

These are the materials for the blog post:
  https://dirtyhandscoding.github.io/posts/performance-comparison-linear-search-vs-binary-search.html

It compares several implementations of linear or binary search intended to find position within a sorted array.
Here is the main plot with results (Broadwell 2 Ghz CPU), see the blog post for more information.

<img src="https://dirtyhandscoding.github.io/images/2_lin_bin_search/plot_search_65536_thr.png" width="500" height="500">
<img src="https://dirtyhandscoding.github.io/images/2_lin_bin_search/plot_search_65536_lat.png" width="500" height="500">

Contents
--------

All the C++ code is in `search.cpp` file, see more comments inside it.

Tiny script `c.bat` gives the command line for compiling `search.cpp` using MSVC compiler. Aside from the obvious `/O2` setting there, also `/arch:AVX2` and `/D NDEBUG` are important for proper performance measurements.
Also it contains commented command line for GCC compilation, which also needs one specific flag `-fno-strict-overflow`.

The python scripts work as follows (Python 3 is required).
First you run `collectdata.py`, which takes `search.cpp` source, copies it into `search_tmp.cpp` with some changes, compiles it using `c.bat` script and runs `search_tmp.exe` to generate text logs into `/res` subdirectory.
Then you run `makeplots.py`, which takes all logs in `/res` subdirectory and generates png images with plots (in the same subdirectory).
