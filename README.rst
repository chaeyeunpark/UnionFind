UnionFind
=========

|Documentation Status| |CodeFactor|

Introduction
------------

C++ implementation of the Union-Find decoder
`arXiv:1709:06218 <https://arxiv.org/abs/1709.06218>`__. Python
interface is also implemented using
`pybind11 <https://github.com/pybind/pybind11>`__.

Based on a Python implementation by Kai Meinerz.

Under the LGPL lisence.

This repository includes codes for
`arXiv:2101.07285 <https://arxiv.org/abs/2101.07285>`__ which explores a
machine learning assisted preprocessing combined with conventional
decoders such as minimum-weight perfect matching (MWPM) and the
Union-Find decoder.

Installation
------------

.. installation-start-inclusion-marker-do-not-remove

Currently, you can install this project by cloning the source code tree and compiling it.
To clone the source tree, use ``git clone`` as

.. code-block:: bash 

	git clone https://github.com/chaeyeunpark/UnionFind.git && cd UnionFind
	git submodule update --init --recursive

Then you can compile the code using 

.. code-block:: bash

    pip install -r requirements.txt
    python3 setup.py install

Note that a compiler with some C++20 support (e.g. GCC version => 10 or Clang++ version => 12) is required. For example, if you are using Ubuntu

.. code-block:: bash

    apt install -y g++-10
    CXX=g++-10 python3 setup.py install

will work.

PyPI support will be available soon. 

.. installation-end-inclusion-marker-do-not-remove

Basic Usage
-----------

.. code-block:: python

   from UnionFindPy import Decoder
   decoder = Decoder(parity_matrix)
   decoder.decode(syndromes) # syndromes is a list of measurment outcomes of each parity operator

For details, check the
`document <https://unionfind.readthedocs.io/en/latest/?badge=latest>`__.

Notes
-----

This repository does not contain an implementation of weighted
Union-Find decoder.

Reference
---------

When you cite this repository, please use the following:

.. code-block:: bibtex

    @misc{UnionFindCPP, 
        author = {Chae-Yeun Park and Kai Meinerz}, 
        title = {Open-source C++ implementation of the Union-Find decoder}, 
        year = {2020}, 
        publisher = {GitHub}, 
        journal = {GitHub repository},
        howpublished = {:raw-latex:`\url{https://github.com/chaeyeunpark/UnionFind}`}
    }



.. |Documentation Status| image:: https://readthedocs.org/projects/unionfind/badge/?version=latest
   :target: https://unionfind.readthedocs.io/en/latest/?badge=latest
.. |CodeFactor| image:: https://www.codefactor.io/repository/github/chaeyeunpark/unionfind/badge
   :target: https://www.codefactor.io/repository/github/chaeyeunpark/unionfind
