# This code is based on PyMatching document
# https://pymatching.readthedocs.io/en/stable/toric-code-example.html
# Modified for the C++ UnionFind Project <https://github.com/chaeyeunpark/UnionFind>

import numpy as np
import matplotlib.pyplot as plt
from utils import repetition_code, toric_code_x_logicals, toric_code_x_stabilisers
from union_find import Decoder


def num_decoding_failures_noisy_syndromes(L, H, logicals, p, q, num_trials, repetitions):
    num_stabilisers, num_qubits = H.shape
    num_errors = 0
    decoder = UnionFind3D(L)
    for i in range(num_trials):
        noise_new = (np.random.rand(repetitions, num_qubits) < p).astype(np.uint8)
        noise_cumulative = (np.cumsum(noise_new, axis=0) % 2).astype(np.uint8)
        noise_total = noise_cumulative[-1, :]
        syndrome = noise_cumulative @ H.T % 2
        syndrome_error = (np.random.rand(repetitions, num_stabilisers) < q).astype(
            np.uint8
        )
        syndrome_error[-1,:] = 0 # Perfect measurements in last round to ensure even parity
        noisy_syndrome = (syndrome + syndrome_error) % 2
        # Convert to difference syndrome
        noisy_syndrome[1:,:] = (noisy_syndrome[1:, :] - noisy_syndrome[0:-1, :]) % 2
        correction = decoder.decode(noisy_syndrome.flatten())
        decoder.clear()

        correction_2d = np.zeros((2 * L * L,), dtype=np.int32)
        for i in range(repetitions):
            correction_2d += correction[i * 3 * L * L : i * 3 * L * L + 2 * L * L]
        error = (noise_total + correction_2d) % 2
        assert not np.any(H @ error % 2)
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return num_errors


if __name__ == "__main__":
    num_trials = 5000
    Ls = range(8, 13, 2)
    ps = np.linspace(0.02, 0.04, 7)
    log_errors_all_L = []
    for L in Ls:
        print("Simulating L={}...".format(L))
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        log_errors = []
        for p in ps:
            num_errors = num_decoding_failures_noisy_syndromes(
                L, Hx, logX, p, p, num_trials, L
            )
            log_errors.append(num_errors / num_trials)
        log_errors_all_L.append(np.array(log_errors))

    for L, logical_errors in zip(Ls, log_errors_all_L):
        std_err = (logical_errors * (1 - logical_errors) / num_trials) ** 0.5
        plt.errorbar(ps, logical_errors, yerr=std_err, label="L={}".format(L))
    plt.yscale("log")
    plt.xlabel("Physical error rate")
    plt.ylabel("Logical error rate")
    plt.legend(loc=0)
    plt.show()
