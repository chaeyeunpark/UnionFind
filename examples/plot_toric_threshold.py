# This code is based on PyMatching document
# https://pymatching.readthedocs.io/en/stable/toric-code-example.html
# Modified for the C++ UnionFind Project <https://github.com/chaeyeunpark/UnionFind>

import numpy as np
import matplotlib.pyplot as plt

from utils import repetition_code, toric_code_x_logicals, toric_code_x_stabilisers
from UnionFindPy import Decoder

def num_decoding_failures(L, H, logicals, p, num_trials):
    decoder = Decoder(toric_code_x_stabilisers(L))
    num_errors = 0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, 2 * L * L)
        syndrome = H @ noise % 2
        correction = decoder.decode(syndrome)
        error = (noise + correction) % 2
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return num_errors


if __name__ == "__main__":
    num_trials = 20000
    Ls = range(3, 19, 2)
    ps = np.linspace(0.01, 0.15, 15)
    np.random.seed(2)
    log_errors_all_L = []
    for L in Ls:
        print("Simulating L={}...".format(L))
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        log_errors = []
        for p in ps:
            num_errors = num_decoding_failures(L, Hx, logX, p, num_trials)
            log_errors.append(num_errors / num_trials)
        log_errors_all_L.append(np.array(log_errors))

    print(log_errors_all_L)

    plt.figure()
    for L, logical_errors in zip(Ls, log_errors_all_L):
        std_err = (logical_errors * (1 - logical_errors) / num_trials) ** 0.5
        plt.errorbar(ps, logical_errors, yerr=std_err, label="L={}".format(L))
    plt.xlabel("Physical error rate")
    plt.ylabel("Logical error rate")
    plt.legend(loc=0)
    plt.show()
