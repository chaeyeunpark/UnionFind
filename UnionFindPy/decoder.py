from ._union_find_py import DecoderFromParity
import logging
from scipy.sparse import csr_matrix
import numpy as np

class Decoder:
    """Union-Find decoder class

    :param parity_matrix (scipy.sparse.csr_matrix): a parity matrix in CSR format
    """

    def __init__(self, parity_matrix, repetitions = None):
        """Create a decoder from a parity matrix"""

        if not isinstance(parity_matrix, csr_matrix):
            raise ValueError('Parameter parity_matrix must be a csr matrix.')

        if np.issubdtype(parity_matrix.dtype, int):
            raise ValueError('Data type of parity_matrix must be an integer, but {} is given'.format(parity_matrix.dtype))

        if not np.all(parity_matrix.data == 1):
            raise ValueError('Any non-zero value of the partiy matrix must be 1.')
        
        if repetitions is None:
            self._decoder = DecoderFromParity(parity_matrix.shape[0], 
                    parity_matrix.shape[1], parity_matrix.indices, parity_matrix.indptr)
        else:
            self._decoder = DecoderFromParity(parity_matrix.shape[0], 
                    parity_matrix.shape[1], parity_matrix.indices, parity_matrix.indptr,
                    repetitions)

    def decode(self, syndrome_arr):
        """Decode a given syndrome array.

        :param syndrome_arr: for a given parity index `i`, syndrome_arr[i] must be 0 or 1. 
        """
        res = self._decoder.decode(syndrome_arr)
        self._decoder.clear()
        return res
