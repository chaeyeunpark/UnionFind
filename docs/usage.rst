Usage
================

Python interface of this project works almost the same as PyMatching.
As in `the PyMatching document for decoding the toric code <https://pymatching.readthedocs.io/en/latest/toric-code-example.html>`_, we create stabilizers for the toric code

.. code-block:: python

    # This code is based on PyMatching document
    # https://pymatching.readthedocs.io/en/stable/toric-code-example.html
    # Modified for the C++ UnionFind Project <https://github.com/chaeyeunpark/UnionFind>

    import numpy as np
    from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag

    def repetition_code(n):
            ...
        return csr_matrix((data, (row_ind, col_ind)))

    def toric_code_x_stabilisers(L):
            ...
        return csr_matrix(H)

and the logical operator

.. code-block:: python

    def toric_code_x_logicals(L):
            ...
        return csr_matrix(x_logicals)

Then you can call the decoder with

.. code-block:: python

    decoder = Decoder(toric_code_x_stabilisers(L))
    correction = decoder.decode(syndrome) 

Noisy version also works almost exactly same as PyMatching except that a syndrome array saves a result of syndrome measurement of each time-slice in row (instead of column as in PyMatching example).

See code inside ``examples`` directory to see working examples.
