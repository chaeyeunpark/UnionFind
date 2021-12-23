Development
====================

This project consists of C++ implementation of `the Union-Find decoder <https://quantum-journal.org/papers/q-2021-12-02-595/>`_  and its python binding. For C++ code, one may look ``union_find/cpp`` directory. We note that this directory itself is a complete C++ project you can build with ``CMake``.

To compile C++ examples, you may run

.. code-block:: shell

    cd union_find/cpp
    mkdir build && cd build
    cmake -DBUILD_EXAMPLES=ON ..
    make

Other supported ``cmake`` options are  ``-DENABLE_AVX=ON``, ``-DBUILD_TESTS=ON``, ``-DCLANG_TIDY``.
