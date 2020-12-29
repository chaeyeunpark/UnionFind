# UnionFind
C++ implementation of Union-Find decoder [arXiv:1709:06218](https://arxiv.org/abs/1709.06218). 
Python interface is also implemented using [pybind11](https://github.com/pybind/pybind11). 

Based on Python implemtation by Kai Meinerz.

Under LGPL lisence. 


# Implementation detail
Our implementation mostly follows the original paper but slightly differ in implementation of peeling decoder. 

The peeling decoder introduced in [arXiv:1703.01517](https://arxiv.org/abs/1703.01517) decodes erasure errors in linear-time by finding a minimum spanning tree from possible error configurations. This peeling decoder is utilized at the final stage of Union-Find decoder. 
There are lots of algorithms for finding a minimum spanning tree but here we implemented a greedy type of alogrithm as it is the most simple. 
Especially, our greedy algorithm does not introduce any additional hash set/map that may be required in other alrogrithms, saves some contant overheads for initializing such data structures.
In addition, even though the worst-time complexity of our greedy algorithm is not linear, we have not observed this part slows down the decoding from our profiling results. More than 90% of time is taken by the main Union-Find algorithm regardless of the lattice size and error rates.

# Usage
First, set up the source tree.
```bash
git clone https://github.com/chaeyeunpark/UnionFind.git
cd UnionFind
git submodule update --init --recursive
```

You can compile the code as below.
```bash
mkdir build && cd build
cmake ..
make union_find
```

After that, you can see `union_find.cpython-[some extra string].so` file in your `build` directory. Copy this file into your python code directory. Then you can use it as below
```python
from union_find import UnionFindDecoder
decoder = UnionFindDecoder(lattice_size) # L = lattice_size
decoder.decode(syndromes) # syndromes is a list of size L^2
decoder.clear() # should be called for reuse
```

# Notes
This repository does not contain an implementation of weighted Union-Find decoder. 
