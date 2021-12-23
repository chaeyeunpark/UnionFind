
[![Documentation Status](https://readthedocs.org/projects/unionfind/badge/?version=latest)](https://unionfind.readthedocs.io/en/latest/?badge=latest)
[![CodeFactor](https://www.codefactor.io/repository/github/chaeyeunpark/unionfind/badge)](https://www.codefactor.io/repository/github/chaeyeunpark/unionfind)


# UnionFind
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

and install the package.
```bash
pip install -r requirements.txt
python3 setup.py install
```

After that, you can see `union_find.cpython-[some extra string].so` file in your `build` directory. Copy this file into your python code directory. Then you can use it.
```python
from union_find import Decoder
decoder = Decoder(parity_matrix)
decoder.decode(syndromes) # syndromes is a list of measurment outcomes of each parity operator
```

For details, check the [document](https://unionfind.readthedocs.io/en/latest/?badge=latest).

## Notes
This repository does not contain an implementation of weighted Union-Find decoder. 

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
