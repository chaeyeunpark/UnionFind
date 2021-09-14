# This code is based on PyMatching document 
# https://pymatching.readthedocs.io/en/stable/toric-code-example.html
# Modified for the C++ UnionFind Project <https://github.com/chaeyeunpark/UnionFind>

import numpy as np
import matplotlib.pyplot as plt

from utils import repetition_code, toric_code_x_logicals, toric_code_x_stabilisers
from union_find import UnionFindToric


def num_decoding_failures(L, H, logicals, p, num_trials):
    decoder = UnionFindToric(L)
    num_errors = 0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, 2*L*L)
        syndrome = H@noise % 2
        correction = decoder.decode(syndrome)
        decoder.clear()
        error = (noise + correction) % 2
        if np.any(error@logicals.T % 2):
            num_errors += 1
    return num_errors

if __name__ == '__main__':
    print(toric_code_x_stabilisers(4))

