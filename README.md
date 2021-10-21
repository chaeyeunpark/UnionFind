# UnionFind
[![Documentation Status](https://readthedocs.org/projects/unionfind/badge/?version=latest)](https://unionfind.readthedocs.io/en/latest/?badge=latest)

C++ implementation of the Union-Find decoder [arXiv:1709:06218](https://arxiv.org/abs/1709.06218). 
Python interface is also implemented using [pybind11](https://github.com/pybind/pybind11). 

Based on a Python implementation by Kai Meinerz.

Under the LGPL lisence. 

This repository includes codes for [arXiv:2101.07285](https://arxiv.org/abs/2101.07285) which explores a machine learning assisted preprocessing combined with conventional decoders such as minimum-weight perfect matching (MWPM) and the Union-Find decoder. 

## Implementation detail
Our implementation mostly follows the original paper but slightly differs in the implementation of peeling decoder. 

The peeling decoder introduced in [arXiv:1703.01517](https://arxiv.org/abs/1703.01517) decodes erasure errors in linear-time by finding a minimum spanning tree from possible error configurations. This peeling decoder is utilized at the final stage of the Union-Find decoder. 
There are lots of algorithms for finding a minimum spanning tree but here we implemented a greedy type of one as it is the most simple. 
Especially, our greedy algorithm does not introduce any additional hash set/map that may be required in other algorithms, saves some constant overheads for initializing such data structures.
In addition, even though the worst-time complexity of our greedy algorithm is not linear, we have not observed this part slows down the overall decoding time from our profiling results. More than 90% of time is taken by the main Union-Find algorithm regardless of lattice sizes and error rates.

## Usage
First, set up the source tree
```bash
git clone https://github.com/chaeyeunpark/UnionFind.git
cd UnionFind
git submodule update --init --recursive
```

and compile the code.
```bash
mkdir build && cd build
cmake ..
make union_find
```

After that, you can see `union_find.cpython-[some extra string].so` file in your `build` directory. Copy this file into your python code directory. Then you can use it.
```python
from union_find import UnionFindDecoder
decoder = UnionFindDecoder(lattice_size) # L = lattice_size
decoder.decode(syndromes) # syndromes is a list of size L^2
decoder.clear() # should be called for reuse
```

You can also compile binary executables:
```bash
cmake -DBUILD_EXECUTABLES=ON ..
make all
```

This builds the toric code examples with and without syndrome measurement errors for the bit-flip and depolarizing noise models. If MPI is found, the MPI version of code is also compiled.


## Notes
This repository does not contain an implementation of weighted Union-Find decoder. 

## Changes from v0.1
* Qubit ordering and interfaces have been changed. Now it is compatible with [PyMatching](https://github.com/oscarhiggott/PyMatching).
* 3D version is added to address noisy syndrome measurements.


## Reference
When you cite this repository, please use the following:
```
@misc{UnionFindCPP,
  author = {Chae-Yeun Park and Kai Meinerz},
  title = {Open-source C++ implementation of the Union-Find decoder},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/chaeyeunpark/UnionFind}}
}

