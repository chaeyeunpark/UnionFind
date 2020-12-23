# UnionFind
C++ implementation of Union-Find decoder [arXiv:1709:06218](https://arxiv.org/abs/1709.06218). 
Python interface is also available thanks to [pybind11](https://github.com/pybind/pybind11). 

Based on Python implemtation by Kai Meinerz.

Under LGPL lisence. 


# Implementation detail
Our implementation mostly follows the original paper but slightly differ in implementation of peeling decoder. 

The peeling decoder introduced in [arXiv:1703.01517](https://arxiv.org/abs/1703.01517) decodes erasure errors in linear-time by finding a minimum spanning tree from possible error configurations. This peeling decoder is utilized at the final stage of Union-Find decoder. 
There are lots of algorithms for finding a minimum spanning tree but here we implemented a greedy type of alogrithm as it is the most simple. 
Especially, our greedy algorithm does not introduce any additional hash set/map that may be required in other alrogrithms, saves some contant overheads for initializing such data structures.
In addition, even though the worst-time complexity of our greedy algorithm is not linear, we have not observed slowdone due to this part in our profiling results. More than 90% of time is taken by the main Union-Find algorithm regardless of the lattice size and error rates.

# Usage
to be added..
