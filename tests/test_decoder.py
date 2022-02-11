import pytest
from UnionFindPy import Decoder
import numpy as np
from scipy.sparse import csr_matrix

def test_toric33():
    parity_matrix = np.zeros((9, 18), dtype=np.int8)
    # P0
    parity_matrix[0,0] = 1
    parity_matrix[0,2] = 1
    parity_matrix[0,3] = 1
    parity_matrix[0,15] = 1
    # P1
    parity_matrix[1,0] = 1
    parity_matrix[1,1] = 1
    parity_matrix[1,4] = 1
    parity_matrix[1,16] = 1
    # P2
    parity_matrix[2,1] = 1
    parity_matrix[2,2] = 1
    parity_matrix[2,5] = 1
    parity_matrix[2,17] = 1
    # P3
    parity_matrix[3,3] = 1
    parity_matrix[3,6] = 1
    parity_matrix[3,8] = 1
    parity_matrix[3,9] = 1
    # P4
    parity_matrix[4,4] = 1
    parity_matrix[4,6] = 1
    parity_matrix[4,7] = 1
    parity_matrix[4,10] = 1
    # P5
    parity_matrix[5,5] = 1
    parity_matrix[5,7] = 1
    parity_matrix[5,8] = 1
    parity_matrix[5,11] = 1
    # P6
    parity_matrix[6,9] = 1
    parity_matrix[6,12] = 1
    parity_matrix[6,14] = 1
    parity_matrix[6,15] = 1
    # P7
    parity_matrix[7,10] = 1
    parity_matrix[7,12] = 1
    parity_matrix[7,13] = 1
    parity_matrix[7,16] = 1
    # P8
    parity_matrix[8,11] = 1
    parity_matrix[8,13] = 1
    parity_matrix[8,14] = 1
    parity_matrix[8,17] = 1

    parity_matrix = csr_matrix(parity_matrix)

    decoder = Decoder(parity_matrix)

    syndrom_arr = [0]*9

    # error = 3
    syndrom_arr[0] = 1
    syndrom_arr[3] = 1

    expected = [0] * 18
    expected[3] = 1

    assert np.all(decoder.decode(syndrom_arr) == expected)
