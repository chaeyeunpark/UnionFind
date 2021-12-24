Development
====================

This project consists of C++ implementation of `the Union-Find decoder <https://quantum-journal.org/papers/q-2021-12-02-595/>`_  and its python binding. For C++ code, one may look ``UnionFindPy/cpp`` directory. We note that this directory itself is a complete C++ project you can build with ``CMake``.

To compile C++ examples, you may run

.. code-block:: shell

    $ cd UnionFindPy/cpp
    $ mkdir build && cd build
    $ cmake -DBUILD_EXAMPLES=ON ..
    $ make

Other supported ``cmake`` options are  ``-DENABLE_AVX=ON``, ``-DBUILD_TESTS=ON``, ``-DCLANG_TIDY``.


For a contribution, I ask you install ``clang-tidy-12`` and ``clang-format-12``. You can format C++ source files with:

.. code-block:: shell

    $ cd UnionFindPy/cpp
    $ make format

where ``clang-tidy`` can be called with:

.. code-block:: shell

    $ cd UnionFindPy/cpp
    $ mkdir build && cd build
    $ cmake -DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON -DCLANG_TIDY=ON ..
    $ make
