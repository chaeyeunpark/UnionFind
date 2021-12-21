import _union_find_py
import logging

class Decoder(_union_find_py.DecoderFromParity):
    def __init__(parity_matrix):
        """
        parity_matrix in CSR format
        """
        super().__init__(parity_matrix.shape[0], parity_matrix.shape[1], 
                parity_matrix.nnz, parity_matrix.indices, parity_matrix.indptr)

