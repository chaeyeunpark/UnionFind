Installation
===============

.. highlight:: bash 

Currently, you can install this project by cloning the source code tree and compiling it.
To clone the source tree, use ``git clone`` as ::

	git clone https://github.com/chaeyeunpark/UnionFind.git && cd UnionFind
	git submodule update --init --recursive

Then you can compile the code using ::

    pip install -r requirements
    python3 setup.py install

Note that a compiler with some C++20 support (e.g. GCC version => 10 or Clang++ version => 12) is required. For example, if you are using Ubuntu ::

    apt install -y g++-10
    CXX=g++-10 python3 setup.py install

will work.

PyPI support will be available soon. 
