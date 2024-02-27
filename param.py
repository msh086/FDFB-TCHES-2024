import argparse
from math import log, ceil, log2, sqrt, gcd

encstd = 3.19
normbnd = 6.338
# Q54 = 18014398509404161
# P54 = 18014398509506561
# Q27 = 134215681

Q53 = 9007199254614017
P53 = 9007199254781953
Q26 = 67104769  # TODO: maybe larger?
P26 = 67127297

n35 = 1340
n25 = 955
n20 = 760


# def var_ACC(n, N, Q, B_g, std):
#     d_g = ceil(log(Q, B_g))
#     msd = ceil(Q / B_g**(d_g - 1))
#     return 4 * ((d_g - 1) * B_g**2 + msd**2) * n * N * std**2 / 6  # 1/12, MSD
def var_ACC(n, N, Q, B_g, std):
    d_g = ceil(log(Q, B_g))
    # msd = ceil(Q / B_g**(d_g - 1))
    return 4 * d_g * B_g**2 * n * N * std**2 / 6  # 1/12, MSD


def time_ACC(n, N, Q, B_g):
    d_g = ceil(log(Q, B_g))
    return 2 * n * 2 * (d_g + 1) * N * log(N)  # constant of NTT


def mem_ACC(n, N, Q, B_g):
    d_g = ceil(log(Q, B_g))
    return 2 * n * 2 * d_g * 2 * N * log2(Q)  # bits


# TODO: more fine-grained
# NOTE: Key switching is similar to FHEW/AP style bootstrapping, i.e. by using many keys and table lookup
def var_KS(N, q_ks, B_ks, std):
    d_ks = ceil(log(q_ks, B_ks))
    return N * d_ks * (1 - 1/B_ks) * std**2  # 1/B_ks


def time_KS(n, N, q_ks, B_ks):
    d_ks = ceil(log(q_ks, B_ks))
    return N * d_ks * (n+1)  # constant of add


def mem_KS(n, N, q_ks, B_ks):
    d_ks = ceil(log(q_ks, B_ks))
    return N * d_ks * (B_ks - 1) * log2(q_ks) * (n+1)  # bits


def var_MS(Qfrom, Qto, dim):
    denom = Qfrom / gcd(Qfrom, Qto)
    return (1-denom**-2)/12 * (1 + dim * 2 / 3)


# TODO: more find-grained
def var_KS_mult(N, q_ks, B_ks, std):
    d_ks = ceil(log(q_ks, B_ks))
    # 1/12 instead of 1/4 because LWE ctxt is nearly random
    return N * d_ks * B_ks**2 / 12 * std**2


def mem_KS_mult(n, N, q_ks, B_ks):
    d_ks = ceil(log(q_ks, B_ks))
    return N * d_ks * log2(q_ks) * (n+1)  # bits


def var_PK_scaled(N, qfrom, qto, B_pk, std):
    d_pk = ceil(log(qfrom, B_pk))
    # rounding error + ctxt error
    return N * d_pk * (1 - 1/B_pk) * (std**2 + 1/12)


"""
Qin: LWE modulus    # given
n: LWE dimension    # fixed 1305, TODO search?
q_ks: KS modulus    # search in 2^i, [Qin, 2^35]
B_ks: KS base       # search in 2^i, [2, sqrt()]
B_g: RGSW base      # search in 2^i
N: RLWE             # search in {2^10, 2^11}
"""


def search_params(Qin):
    Qin_bits = int(log2(Qin))
    n = 1305
    for N, Q in [
        (2**10, 2**27 - 1),
        # (2**11, 2**54 - 1),
    ]:  # set Q to odd number so var_MS can work
        Q_bits = ceil(log2(Q))
        print(f"N={N}, Q=2^{Q_bits}")
        ACC_map = dict()
        for B_g_bits in range(1, 1 + ((Q_bits + 1) >> 1)):
            B_g = 2**B_g_bits
            ACC_map[B_g_bits] = {
                "var": var_ACC(n, N, Q, B_g, encstd),
                "time": time_ACC(n, N, Q, B_g),
                "mem": mem_ACC(n, N, Q, B_g),
            }
        # KS is relatively lightweight compared to ACC
        for q_ks_bits in range(Qin_bits, 36):
            q_ks = 2**q_ks_bits
            for B_ks_bits in range(1, 1 + ((q_ks_bits + 1) >> 1)):
                B_ks = 2**B_ks_bits
                for B_g_bits, stat in ACC_map.items():
                    var_tmp = (
                        stat["var"] * (q_ks / Q) ** 2
                        + var_MS(Q, q_ks, N)
                        + var_KS(N, q_ks, B_ks, encstd)
                    )
                    std_boot = sqrt(var_tmp * (Qin / q_ks) **
                                    2 + var_MS(q_ks, Qin, n))
                    print(
                        f"B_g=2^{B_g_bits}, q_ks=2^{q_ks_bits}, B_ks=2^{B_ks_bits}, std_bt={std_boot}, bnd={std_boot*6.37}"
                    )
                print("----------")
            print("++++++++++")


# NOTE: upper bound for beta = 64 is var = 101.96 (beta = 128 -> var = 407.86)


def get_bt_std(Qin, n, Q, N, B_g, q_ks, B_ks):
    acc = var_ACC(n, N, Q, B_g, encstd)
    tmp = acc * (q_ks / Q) ** 2 + var_MS(Q, q_ks, N) + \
        var_KS(N, q_ks, B_ks, encstd)
    acc_contrib = acc * (Qin / Q)**2
    ks_contrib = (var_MS(Q, q_ks, N) + var_KS(N, q_ks,
                  B_ks, encstd)) * (Qin / q_ks) ** 2
    ms_contrib = var_MS(q_ks, Qin, n)
    print(f"contribs: ACC={acc_contrib}, ks={ks_contrib}, ms={ms_contrib}")
    return sqrt(tmp * (Qin / q_ks)**2 + var_MS(q_ks, Qin, n))


# OLD ALGORITHMS


# LMP22 params
# standard params:
#   p = 8
#   sqrt(2) * get_bt_std(2**11, n35, Q53, 2**11, 2**27, 2**20, 2**5) * normbnd = 77.56
# small params:
#   p = 16
#   sqrt(2) * get_bt_std(2**11, n20, Q53, 2**11, 2**27, 2**20, 2**5) * normbnd = 58.53

# FDFB-COMPRESS params
# standard params:
#   p = 16
#   get_bt_std(2**12, n35, Q53, 2**11, 2**27, 2**20, 2**5) * normbnd = 55.20
# small params:
#   p = 16
#   get_bt_std(2**12, n20, Q53, 2**11, 2**27, 2**20, 2**5) * normbnd = 41.85

# KS21 FDFB, decomp on Q
# NOTE: DONE
# NOTE: final KS & MS noise not included
# standard params:
# p = 32
#   (2**12, n35, Q53, 2**11, 2**18, 2**6, 2**30, 2**5)  |-> 11.27 + 75.77
#   (2**12, n35, Q53, 2**11, 2**14, 2**10, 2**35, 2**5) |-> 10.05 + 75.77
#   (2**12, n35, Q53, 2**11, 2**11, 2**13, 2**35, 2**5) |-> 11.55 + 75.77
# p = 16
#   (2**12, n35, Q53, 2**11, 2**18, 2**8, 2**30, 2**5)  |-> 140.22 + 75.77
#   (2**12, n35, Q53, 2**11, 2**14, 2**12, 2**35, 2**5) |-> 133.93 + 75.77
#   (2**12, n35, Q53, 2**11, 2**11, 2**15, 2**35, 2**5) |-> 147.82 + 75.77
# small params:
# p = 32
#   (2**12, n20, Q53, 2**11, 2**18, 2**7, 2**30, 2**5)  |-> 22.77 + 43.54
#   (2**12, n20, Q53, 2**11, 2**14, 2**11, 2**35, 2**5) |-> 19.00 + 43.54
#   (2**12, n20, Q53, 2**11, 2**11, 2**14, 2**35, 2**5) |-> 21.35 + 43.54
# p = 16
#   (2**12, n20, Q53, 2**11, 2**18, 2**9, 2**30, 2**5)  |-> 273.23 + 43.54
#   (2**12, n20, Q53, 2**11, 2**14, 2**13, 2**35, 2**5) |-> 303.96 + 43.54
#   (2**12, n20, Q53, 2**11, 2**11, 2**16, 2**35, 2**5) |-> 341.52 + 43.54
# N = 2^10 is impossible
def KS21_var(Qin, n, Q, N, B_g, B_g1, q_pk, B_pk):
    d_g1 = ceil(log(Q, B_g1))
    var = var_ACC(n, N, Q, B_g, encstd)  # blind rotation to get sgn
    var += var_MS(Q, q_pk, N) * (Q/q_pk)**2  # modulus switch to q_pk
    # LWE to RLWE packing (only constant term)
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    # NOTE: we use B_g1^2/4 instead of B_g1^2/12 because the LUT may not be random
    # implementation: product = RLWE'(sgn) * Decomp(m1 - m0) + m0
    # LUT x RLWE' product. TODO: use approximate decomp to achieve smaller noise?.
    var = B_g1**2 * N * var * d_g1 / 4
    print(f"var contrib of phase 1 = {var * (Qin / Q)**2}")
    var += var_ACC(n, N, Q, B_g, encstd)  # second blind rotation
    # var *= (q_ks / Q)**2  # modulus switch to q_ks
    # var += var_MS(Q, q_ks, N)
    # var += var_KS(N, q_ks, B_ks, encstd)  # key switch to LWE scheme
    # var *= (Qin / q_ks)**2  # modulus switch to Qin
    # var += var_MS(q_ks, Qin, n)
    return var * (Qin / Q) ** 2


# WoP-PBS npieces
# NOTE: done
# NOTE: final KS & MS noise not included
# standard params:
#   (32, 2**12, n35, Q53, P53, 2**11, 2**18, 2**25, 2**5, 2**27, 2) |-> 0.54 + 75.77 # ver 2
#   (16, 2**11, n35, Q53, P53, 2**11, 2**18, 2**25, 2**5, 2**27, 1) |-> 0.04 + 74.84 # ver 1
# small params:
#   (32, 2**12, n25, Q53, P53, 2**11, 2**18, 2**25, 2**5, 2**27, 2)  |-> 0.50 + 54.38  # ver 2
#   (16, 2**11, n20, Q53, P53, 2**11, 2**18, 2**20, 2**5, 2**27, 1)  |-> 25.36 + 42.62 # ver 1
# N=2^10 is impossible
def CLOT21_var(p, Qin, n, Q, P, N, B_g, q_pk, B_pk, B_rl, npieces=1):
    d_rl = ceil(log(Q, B_rl))
    # blind rotation + multi-value bootstrap
    # NOTE: the MSB needs to be computed using ptxt space=2p, so to use multi-value bts, we need to view the LUTs as Z_p -> Z_2p mappings
    var_sgn = var_ACC(n, N, Q, B_g, encstd)
    var_msg = var_ACC(n, N, Q, B_g, encstd)
    # XXX: debug
    print(
        f"log2(std sgn) = {log2(sqrt(var_sgn))}, log2(var msg) = {log2(sqrt(var_msg))}")

    # we only need to perform 1 bfv mult, regardless of npieces
    # i.e. MSB(m_- - m_+) + m_+
    var_sgn += (Q / q_pk)**2 * var_MS(Q, q_pk, N)
    var_msg += (Q / q_pk)**2 * var_MS(Q, q_pk, N)

    # background noise: LWE to RLWE packing
    var_pk = var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    var_ms = N / 18 + 1 / 12

    vars = [0] * 9
    vars[0] = var_msg + (p/2)**2 * var_sgn
    vars[1] = (Q/P)**2 * ((P/Q)**2 * var_pk + var_ms) + (p/2)**2 * var_pk
    vars[2] = (p/P)**2 * N * var_pk * ((P/Q)**2 * var_pk + var_ms)
    vars[3] = (p/Q)**2 * var_pk * var_msg + (p/P)**2 * var_sgn * ((P/Q)**2 * var_pk + var_ms)
    vars[4] = (p/P)**2 * var_sgn * var_msg
    vars[5] = 1/12 + N/18 + N/12 * 4/9 * N
    vars[6] = (N/18 + 1/12) * p**2 * (var_msg + var_sgn)
    vars[7] = (N/18 + 1/12) * p**2 * ((2 * var_pk + (Q/P)**2 * var_ms) * N)
    vars[8] = d_rl * B_rl**2 / 12 * N * encstd**2
    # Q/p*MSB + e_bt + e_bg
    # Q/p*diff + e_bt + e_bg -> P/p*diff + P/Q*e_bt + (P/Q*e_bg + e_ms)
    var = var_msg + (p/2)**2 * var_sgn \
        + (Q/P)**2 * ((P/Q)**2 * var_pk + var_ms) + (p/2)**2 * var_pk \
        + (p/P)**2 * N * var_pk * ((P/Q)**2 * var_pk + var_ms) \
        + (p/Q)**2 * var_pk * var_msg + (p/P)**2 * var_sgn * ((P/Q)**2 * var_pk + var_ms) \
        + (p/P)**2 * var_sgn * var_msg \
        + 1/12 + N/18 + N/12 * 4/9 * N \
        + (N/18 + 1/12) * p**2 * (var_msg + var_sgn +
                                  (2 * var_pk + (Q/P)**2 * var_ms) * N)
    totalvar = sum(vars)
    for i, ele in enumerate(vars):
        print(f"var[{i}] / total = {ele/totalvar}")
    # var = var + var_pk + npieces * (p**2 / 4 * var + p**2 / 4 * var_pk * N + (2 * var * var_pk * N) * (p/Q)**2) + var_MS(Q**2, Q, N)

    print(f"log2(std beforeKS) = {log2(sqrt(var))}")

    var += d_rl * B_rl**2 / 12 * N * encstd**2  # BFV relin
    if npieces > 1:
        var += var_msg
    return var * (Qin / Q)**2


# WoP-PBS npieces + multi-value bootstrap
# NOTE: done
# NOTE: final KS & MS noise not included
# standard params:
#   (32, 2**12, n35, Q53, P53, 2**11, 2**11, 2**25, 2**5, 2**27, 2) |-> 1.32 + 75.77   # ver2
#   (16, 2**12, n35, Q53, P53, 2**11, 2**18, 2**25, 2**5, 2**27, 2) |-> 250.06 + 75.77 # ver2
#   (16, 2**11, n35, Q53, P53, 2**11, 2**18, 2**25, 2**5, 2**27, 1) |-> 15.66 + 75.77  # ver1
# small params:
#   (32, 2**12, n20, Q53, P53, 2**11, 2**14, 2**25, 2**5, 2**27, 2) |-> 25.75 + 43.54 # ver2
#   (16, 2**12, n20, Q53, P53, 2**11, 2**18, 2**20, 2**5, 2**27, 2)  |-> 243.17 + 43.54 # ver2
#   (16, 2**11, n20, Q53, P53, 2**11, 2**18, 2**20, 2**5, 2**27, 1)  |-> 34.22 + 43.54 # ver1
# N=2^10 is impossible
def CLOT21_var_mv(p, Qin, n, Q, P, N, B_g, q_pk, B_pk, B_rl, npieces=1):
    d_rl = ceil(log(Q, B_rl))
    # blind rotation + multi-value bootstrap
    # NOTE: the MSB needs to be computed using ptxt space=2p, so to use multi-value bts, we need to view the LUTs as Z_p -> Z_2p mappings
    var_sgn = var_ACC(n, N, Q, B_g, encstd) * 4  # TV = 1....1,-1...-1
    # Z_p -> Z_2p if npieces == 2, otherwise Z_p -> Z_p
    var_msg = var_ACC(n, N, Q, B_g, encstd) * p * \
        (4 if npieces == 2 else 1) * (p - 1)**2
    # XXX: debug
    print(
        f"log2(std sgn) = {log2(sqrt(var_sgn))}, log2(var msg) = {log2(sqrt(var_msg))}")

    # we only need to perform 1 bfv mult, regardless of npieces
    # i.e. MSB(m_- - m_+) + m_+
    var_sgn += (Q / q_pk)**2 * var_MS(Q, q_pk, N)
    var_msg += (Q / q_pk)**2 * var_MS(Q, q_pk, N)

    # background noise: LWE to RLWE packing
    var_pk = var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    var_ms = N / 18 + 1 / 12
    # Q/p*MSB + e_bt + e_bg
    # Q/p*diff + e_bt + e_bg -> P/p*diff + P/Q*e_bt + (P/Q*e_bg + e_ms)
    vars = [0] * 9
    vars[0] = var_msg + (p/2)**2 * var_sgn
    vars[1] = (Q/P)**2 * ((P/Q)**2 * var_pk + var_ms) + (p/2)**2 * var_pk
    vars[2] = (p/P)**2 * N * var_pk * ((P/Q)**2 * var_pk + var_ms)
    vars[3] = (p/Q)**2 * var_pk * var_msg + (p/P)**2 * var_sgn * ((P/Q)**2 * var_pk + var_ms)
    vars[4] = (p/P)**2 * var_sgn * var_msg
    vars[5] = 1/12 + N/18 + N/12 * 4/9 * N
    vars[6] = (N/18 + 1/12) * p**2 * (var_msg + var_sgn)
    vars[7] = (N/18 + 1/12) * p**2 * ((2 * var_pk + (Q/P)**2 * var_ms) * N)
    vars[8] = d_rl * B_rl**2 / 12 * N * encstd**2
    totalvar = sum(vars)
    for i, ele in enumerate(vars):
        print(f"var[{i}] / total = {ele/totalvar}")

    var = var_msg + (p/2)**2 * var_sgn \
        + (Q/P)**2 * ((P/Q)**2 * var_pk + var_ms) + (p/2)**2 * var_pk \
        + (p/P)**2 * N * var_pk * ((P/Q)**2 * var_pk + var_ms) \
        + (p/Q)**2 * var_pk * var_msg + (p/P)**2 * var_sgn * ((P/Q)**2 * var_pk + var_ms) \
        + (p/P)**2 * var_sgn * var_msg \
        + 1/12 + N/18 + N/12 * 4/9 * N \
        + (N/18 + 1/12) * p**2 * (var_msg + var_sgn +
                                  (2 * var_pk + (Q/P)**2 * var_ms) * N)
    # var = var + var_pk + npieces * (p**2 / 4 * var + p**2 / 4 * var_pk * N + (2 * var * var_pk * N) * (p/Q)**2) + var_MS(Q**2, Q, N)

    print(f"log2(std beforeKS) = {log2(sqrt(var))}")

    var += d_rl * B_rl**2 / 12 * N * encstd**2  # BFV relin
    if npieces > 1:
        var += var_msg
    return var * (Qin / Q)**2


# WoPPBS-i with old noise analysis
# NOTE: OLD version
# NOTE: done
# NOTE: final KS & MS not included
# standard params
#   (32, 2**12, n35, Q53, P53, 2**11, 2**14, 2**30, 2**5, 2**27, 2) |-> 4.37 + 75.77  # ver 2
#   (16, 2**12, n35, Q53, P53, 2**11, 2**18, 2**30, 2**5, 2**27, 2) |-> 142.71 + 75.77 # ver 2
#   (16, 2**11, n35, Q53, P53, 2**11, 2**18, 2**30, 2**5, 2**27, 1) |-> 17.84 + 75.77 # ver 1
# small params
#   (32, 2**12, n20, Q53, P53, 2**11, 2**14, 2**30, 2**5, 2**27, 2) |-> 3.28 + 43.54 # ver 2
#   (16, 2**12, n20, Q53, P53, 2**11, 2**18, 2**30, 2**5, 2**27, 2) |-> 81.11 + 43.54 # ver 2
#   (16, 2**11, n20, Q53, P53, 2**11, 2**18, 2**30, 2**5, 2**27, 1) |-> 10.14 + 43.54 # ver 1
# N=2^10 is impossible
def CLOT21_var_old(p, Qin, n, Q, P, N, B_g, q_pk, B_pk, B_rl, npieces=1):
    d_rl = ceil(log(Q, B_rl))
    # blind rotation + multi-value bootstrap
    # NOTE: the MSB needs to be computed using ptxt space=2p, so to use multi-value bts, we need to view the LUTs as Z_p -> Z_2p mappings
    var_sgn = var_ACC(n, N, Q, B_g, encstd)  # TV = 1....1,-1...-1
    # Z_p -> Z_2p if npieces == 2, otherwise Z_p -> Z_p
    var_msg = var_ACC(n, N, Q, B_g, encstd)

    # we only need to perform 1 bfv mult, regardless of npieces
    # i.e. MSB(m_- - m_+) + m_+
    var_sgn += (Q / q_pk)**2 * var_MS(Q, q_pk, N)
    var_msg += (Q / q_pk)**2 * var_MS(Q, q_pk, N)

    # background noise: LWE to RLWE packing
    var_pk = var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    var_sgn += var_pk
    var_msg += var_pk
    var_ms = N / 18 + 1 / 12
    # Q/p*MSB + e_bt + e_bg
    # Q/p*diff + e_bt + e_bg -> P/p*diff + P/Q*e_bt + (P/Q*e_bg + e_ms)
    vars = [0] * 7
    vars[0] = N * (var_msg + (p/2)**2 * var_sgn)
    vars[1] = N * ((p/Q)**2 * var_msg * var_sgn)
    vars[2] = N * ((Q/P)**2 * var_ms + (p/P)**2 * var_sgn * var_ms)
    vars[3] = 1/12 + N/18 + N/12 * 4/9 * N
    vars[4] = (N/18 + 1/12) * p**2 * (var_msg + var_sgn) * N
    vars[5] = (N/18 + 1/12) * p**2 * ((Q/P)**2 * var_ms) * N
    vars[6] = d_rl * B_rl**2 / 12 * N * encstd**2
    totalvar = sum(vars)
    for i, ele in enumerate(vars):
        print(f"var[{i}] / total = {ele/totalvar}")

    var = N * (var_msg + (p/2)**2 * var_sgn \
        + (p/Q)**2 * var_msg * var_sgn \
        + (Q/P)**2 * var_ms + (p/P)**2 * var_sgn * var_ms) \
        + 1/12 + N/18 + N/12 * 4/9 * N \
        + (N/18 + 1/12) * p**2 * (var_msg + var_sgn +
                                  (Q/P)**2 * var_ms) * N
    # var = var + var_pk + npieces * (p**2 / 4 * var + p**2 / 4 * var_pk * N + (2 * var * var_pk * N) * (p/Q)**2) + var_MS(Q**2, Q, N)
    var *= npieces

    print(f"log2(std beforeKS) = {log2(sqrt(var))}")

    var += d_rl * B_rl**2 / 12 * N * encstd**2  # BFV relin
    return var * (Qin / Q)**2


# WoPPBS-i with old noise analysis, nulti-value
# NOTE: OLD version
# NOTE: done
# NOTE: final KS & MS not included
# standard params
#   (32, 2**12, n35, Q53, P53, 2**11, 2**6, 2**30, 2**5, 2**27, 2)  |-> 9.23 + 75.77
#   (16, 2**12, n35, Q53, P53, 2**11, 2**11, 2**30, 2**5, 2**27, 2) |-> 113.22 + 75.77
#   (16, 2**11, n35, Q53, P53, 2**11, 2**11, 2**30, 2**5, 2**27, 1) |-> 3.58 + 75.77
# small params
#   (32, 2**12, n20, Q53, P53, 2**11, 2**8, 2**30, 2**5, 2**27, 2)  |-> 52.16 + 43.54
#   (16, 2**12, n20, Q53, P53, 2**11, 2**11, 2**30, 2**5, 2**27, 2) |-> 64.39 + 43.54
#   (16, 2**11, n20, Q53, P53, 2**11, 2**11, 2**25, 2**5, 2**27, 1) |-> 52.70 + 43.54
# N=2^10 is impossible
def CLOT21_var_mv_old(p, Qin, n, Q, P, N, B_g, q_pk, B_pk, B_rl, npieces=1):
    d_rl = ceil(log(Q, B_rl))
    # blind rotation + multi-value bootstrap
    # NOTE: the MSB needs to be computed using ptxt space=2p, so to use multi-value bts, we need to view the LUTs as Z_p -> Z_2p mappings
    var_sgn = var_ACC(n, N, Q, B_g, encstd) * 4  # TV = 1....1,-1...-1
    # Z_p -> Z_2p if npieces == 2, otherwise Z_p -> Z_p
    var_msg = var_ACC(n, N, Q, B_g, encstd) * p * \
        (4 if npieces == 2 else 1) * (p - 1)**2

    # we only need to perform 1 bfv mult, regardless of npieces
    # i.e. MSB(m_- - m_+) + m_+
    var_sgn += (Q / q_pk)**2 * var_MS(Q, q_pk, N)
    var_msg += (Q / q_pk)**2 * var_MS(Q, q_pk, N)

    # background noise: LWE to RLWE packing
    var_pk = var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    var_sgn += var_pk
    var_msg += var_pk
    var_ms = N / 18 + 1 / 12

    vars = [0] * 7
    vars[0] = N * (var_msg + (p/2)**2 * var_sgn)
    vars[1] = N * ((p/Q)**2 * var_msg * var_sgn)
    vars[2] = N * ((Q/P)**2 * var_ms + (p/P)**2 * var_sgn * var_ms)
    vars[3] = 1/12 + N/18 + N/12 * 4/9 * N
    vars[4] = (N/18 + 1/12) * p**2 * (var_msg + var_sgn) * N
    vars[5] = (N/18 + 1/12) * p**2 * ((Q/P)**2 * var_ms) * N
    vars[6] = d_rl * B_rl**2 / 12 * N * encstd**2
    totalvar = sum(vars)
    for i, ele in enumerate(vars):
        print(f"var[{i}] / total = {ele/totalvar}")
    # Q/p*MSB + e_bt + e_bg
    # Q/p*diff + e_bt + e_bg -> P/p*diff + P/Q*e_bt + (P/Q*e_bg + e_ms)
    var = N * (var_msg + (p/2)**2 * var_sgn \
        + (p/Q)**2 * var_msg * var_sgn \
        + (Q/P)**2 * var_ms + (p/P)**2 * var_sgn * var_ms) \
        + 1/12 + N/18 + N/12 * 4/9 * N \
        + (N/18 + 1/12) * p**2 * (var_msg + var_sgn +
                                  (Q/P)**2 * var_ms) * N
    # var = var + var_pk + npieces * (p**2 / 4 * var + p**2 / 4 * var_pk * N + (2 * var * var_pk * N) * (p/Q)**2) + var_MS(Q**2, Q, N)
    var *= npieces
    print(f"log2(std beforeKS) = {log2(sqrt(var))}")

    var += d_rl * B_rl**2 / 12 * N * encstd**2  # BFV relin
    return var * (Qin / Q)**2


# NEW ALGORITHMS

# KS21 FDFB, decomp on p + CIM19 factorization
# NOTE: done
# NOTE: final KS & MS noise not included
# standard params
#   (32, 2**12, n35, Q53, 2**11, 2**27, 2, 2**20, 2**5) |-> 17.07 + 75.77
#   (16, 2**12, n35, Q53, 2**11, 2**27, 4, 2**20, 2**5) |-> 17.07 + 75.77
# small params
#   (32, 2**12, n20, Q53, 2**11, 2**27, 8, 2**20, 2**5) |-> 52.18 + 43.54
#   (16, 2**12, n20, Q53, 2**11, 2**27, 8, 2**20, 2**5) |-> 26.12 + 43.54
# N=2^10 is impossible
def FDFB_preselect(p, Qin, n, Q, N, B_g, B_mv, q_pk, B_pk):
    d_mv = ceil(log(2 * p, B_mv))
    # use sgn to batch select Q/2p*sgn*B^i. 1/4 comes from rand round
    var = var_ACC(n, N, Q, B_g, encstd) + 1/4
    # use LWE to RLWE packing to get Q/2p*sgn*B^i*TV0.
    var += (Q / q_pk)**2 * var_MS(Q, q_pk, N)  # mod switch to q_pk
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd)  # scaled packing
    # NOTE: we use B_g1^2/4 instead of B_g1^2/12 because the LUT may not be random
    # implementation: product = RLWE'(sgn) * Decomp(m1 - m0) + m0
    # LUT x RLWE' product. TODO: use approximate decomp to achieve smaller noise?.
    var *= B_mv**2 * p * d_mv / 4
    print(f"var contrib of phase 1 = {var * (Qin / Q)**2}")
    var += var_ACC(n, N, Q, B_g, encstd)  # second blind rotation
    # var *= (q_ks / Q)**2  # modulus switch to q_ks
    # var += var_MS(Q, q_ks, N)
    # var += var_KS(N, q_ks, B_ks, encstd)  # key switch to LWE scheme
    # var *= (Qin / q_ks)**2  # modulus switch to Qin
    # var += var_MS(q_ks, Qin, n)
    return var * (Qin / Q) ** 2


# NOTE: done
# NOTE: final KS & MS noise not included
# standard params
#   p = 16
#   (2**11, n35, Q53, 2**11, 2**27, 2**15, 2**5) |-> 0.47 + 75.77
# small params
#   p = 16
#   (2**11, n20, Q53, 2**11, 2**27, 2**15, 2**5) |-> 0.56 + 43.54
# N=2^10 is impossible
def FDFB_cancelsign_var(Qin, n, Q, N, B_g, q_pk, B_pk):
    var = var_ACC(n, N, Q, B_g, encstd)
    var += (Q/q_pk)**2 * var_MS(Q, q_pk, N)
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    return var * (Qin / Q)**2


# NOTE: done
# NOTE: final KS & MS noise not included
# NOTE: the final selection BTS do not need small B_g
# standard params
#   (32, 2**12, n35, Q53, 2**11, 2**18, 2**15, 2**5, 64) |-> 3.58 + 75.77
#   (16, 2**12, n35, Q53, 2**11, 2**18, 2**15, 2**5, 32) |-> 3.57 + 75.77
# small params
#   (32, 2**12, n20, Q53, 2**11, 2**18, 2**15, 2**5, 64) |-> 3.57 + 43.54
#   (16, 2**12, n20, Q53, 2**11, 2**27, 2**15, 2**5, 32) |-> 205.01 + 43.54
# N=2^10 is impossible
def FDFB_select_mv_var(p, Qin, n, Q, N, B_g, q_pk, B_pk, B_mv):
    d_mv = ceil(log(2*p, B_mv))
    # blind rotation + decomposed multi-value bootstrap
    var = var_ACC(n, N, Q, B_g, encstd)
    # inner product NOTE: this has to be small enough for 2-out-of-1 selection to work
    var *= d_mv * B_mv**2 / 4 * p
    var += (Q/q_pk)**2 * var_MS(Q, q_pk, N) * 2  # packing two messages
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd) * 2
    var += var_ACC(n, N, Q, B_g, encstd)  # final selection
    return var * (Qin / Q)**2


# NOTE: done
# NOTE: final KS & MS noise not included
# standard params
#   p = 32 & p = 16
#   (2**12, n35, Q53, 2**11, 2**27, 2**15, 2**5) |-> 3.74 + 75.77
# small params
#   p = 32 & 16
#   (2**12, n20, Q53, 2**11, 2**27, 2**15, 2**5) |-> 3.66 + 43.54
# N = 2^10:
#   p = 16
#   (2**11, n20, Q26, 2**10, 2**5, 2**15, 2**5) |-> 50.85 + 42.46
def FDFB_select_var(Qin, n, Q, N, B_g, q_pk, B_pk):
    # 4 bts
    var = var_ACC(n, N, Q, B_g, encstd)
    var += (Q/q_pk)**2 * var_MS(Q, q_pk, N) * 2  # packing two messages
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd) * 2
    var += var_ACC(n, N, Q, B_g, encstd)  # final selection
    return var * (Qin / Q)**2


# assume B_mv = 4p (i.e. original version multi-value bts)
# NOTE: done
# NOTE: final KS & MS noise not included
# NOTE: the final selection BTS do not need small B_g
# standard params
#   p = 32 & 16
#   (32, 2**12, n35, Q53, 2**11, 2**18, 2**15, 2**5) |-> 3.67 + 75.77
# small params
#   p = 32 & p = 16
#   (32, 2**12, n20, Q53, 2**11, 2**18, 2**15, 2**5) |-> 3.63 + 43.54
# N = 2^10 is impossible
def FDFB_select_alt_mv_var(p, Qin, n, Q, N, B_g, q_pk, B_pk):
    var = var_ACC(n, N, Q, B_g, encstd)
    var_diff = var * p * (p-1)**2
    var_add = var * p * 4 * (p-1)**2
    # packing
    var_diff += (Q/q_pk)**2 * var_MS(Q, q_pk, N) + var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    # bts 
    var_out = var_diff + var + var_diff + var_add
    return var_out * (Qin / Q)**2


# NOTE: done
# NOTE: final KS & MS noise not included
# standard params
#   p = 32 & 16
#   (2**12, n35, Q53, 2**11, 2**27, 2**15, 2**5) |-> 2.13 + 75.77
# small params
#   p = 32 & 16
#   (2**12, n20, Q53, 2**11, 2**27, 2**15, 2**5) |-> 1.98 + 43.54
# N = 2**10:
#   p = 16
#   (2**11, n20, Q26, 2**10, 2**4, 2**15, 2**5) |-> 30.76 + 42.46
def FDFB_select_alt_var(Qin, n, Q, N, B_g, q_pk, B_pk):
    # 3 bts
    var = var_ACC(n, N, Q, B_g, encstd)
    # packing
    var += (Q/q_pk)**2 * var_MS(Q, q_pk, N) + var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    # bts and sum
    var = var + 3 * var_ACC(n, N, Q, B_g, encstd) # var + (var + var_bt) + var_bt
    return var * (Qin / Q)**2


# NOTE: done
# NOTE: final KS & MS noise not included
# NOTE: this is the error of preprocessing
# standard params
#   (32, 2**12, n35, Q53, 2**11, 2**18) |-> 0.01 + 75.77
#   (16, 2**12, n35, Q53, 2**11, 2**27)  |-> 27.40 + 74.53 (NOTE: q_ks = 2**25)
# small params
#   (32, 2**12, n20, Q53, 2**11, 2**27) |-> 56.45 + 43.54
#   (16, 2**12, n20, Q53, 2**11, 2**27)  |-> 15.54 + 43.54
# N=2^10 is impossible
def FDFB_comp_mv_var(p, Qin, n, Q, N, B_g):
    var = var_ACC(n, N, Q, B_g, encstd)
    var *= (4 * (p - 1) + p**2)  # NOTE: the TVs of Comp are independent of the LUT
    return var * (Qin / Q)**2


# NOTE: done
# NOTE: final KS & MS noise not included
# standard params
# p = 32 & p = 16
#   (2**12, n35, Q53, 2**11, 2**27) |-> 0.018 + 75.77
# small params
# p = 32 & p = 16
#   (2**12, n20, Q53, 2**11, 2**27) |-> 0.10 + 43.54
# N = 2**10:
#   (2**11, n20, Q26, 2**10, 2**5) |-> 50.40 + 42.46
def FDFB_comp_var(Qin, n, Q, N, B_g):
    return var_ACC(n, N, Q, B_g, encstd) * 2 * (Qin / Q)**2


def common_part_var(Qin, n, Q, N, q_ks, B_ks):
    return (var_MS(Q, q_ks, N) + var_KS(N, q_ks, B_ks, encstd)) * (Qin / q_ks) ** 2 + var_MS(q_ks, Qin, n)


######################## CKKS part ################################


# Qin = 2**12, Lipschitz constant = 1, deltain = 2**8, range = (-8, 8) ~ (-1, 1)
#   (2**60, 2**16, 2**12, n35, 2**35, 2**12, True)|-> 0.0338
#   (2**60, 2**16, 2**12, n25, 2**25, 2**5, True) |-> 0.0296
#   (2**60, 2**16, 2**12, n20, 2**20, 2, False)   |-> 0.0469
# Qin = 2**11, Lipschitz constant = 1, deltain = 2**7, range = (-8, 8) ~ (-1, 1)
#   (2**60, 2**16, 2**11, n35, 2**35, 2**12, True) |-> 0.0675
#   (2**60, 2**16, 2**11, n25, 2**25, 2**5, True) |-> 0.0576
#   (2**60, 2**16, 2**11, n20, 2**20, 2, False)   |-> 0.0644
def CKKS_input_err_std(Qckks, Nckks, Qin, n, q_ks, B_ks, use_mult=False):
    var_from_KS = var_KS_mult(Nckks, q_ks, B_ks, encstd) if use_mult else var_KS(
        Nckks, q_ks, B_ks, encstd)
    print(
        f"Memory usage = {(mem_KS_mult(n, Nckks, q_ks, B_ks) if use_mult else mem_KS(n, Nckks, q_ks, B_ks)) / 2**33} GB")
    print(f"contrib of final MS = {var_MS(q_ks, Qin, n)}")
    return sqrt((Qin / q_ks)**2 * (var_MS(Qckks, q_ks, Nckks) + var_from_KS) + var_MS(q_ks, Qin, n))


def CKKS_input_err_std_predef(n, qin):
    if n == n35:
        return CKKS_input_err_std(2**60, 2**16, qin, n35, 2**35, 2**12, True)
    elif n == n25:
        return CKKS_input_err_std(2**60, 2**16, qin, n25, 2**25, 2**5, True)
    elif n == n20:
        return CKKS_input_err_std(2**60, 2**16, qin, n20, 2**20, 2, False)
    else:
        raise KeyError()


## TODO: maybe it is unnecessary to consider the standard params case?


# NOTE: done
# standard params for output range (-2, 2) and Lipschitz constant 1
#   (2**11, n35, Q53, 2**10, 2**27, 2**25, 2**5, 2**7, 2**23, 2**25, 1) |-> 0.0675
# small params for output range (-2, 2) and Lipschitz constant 1
#   (2**11, n25, Q53, 2**10, 2**27, 2**25, 2**5, 2**7, 2**23, 2**25, 1) |-> 0.0570
#   (2**11, n20, Q53, 2**10, 2**27, 2**20, 2**5, 2**7, 2**18, 2**20, 1) |-> 0.0513
#
# standard params
#   (2**11, n35, Q53, 2**10, 2**27, 2**25, 2**5, 2**7, 2**23, 2**25, 1/4) |-> 0.0239
# small params
#   (2**11, n25, Q53, 2**10, 2**27, 2**25, 2**5, 2**7, 2**23, 2**25, 1/4) |-> 0.0203
#   (2**11, n20, Q53, 2**10, 2**27, 2**20, 2**5, 2**7, 2**18, 2**20, 1/4) |-> 0.0206
def LMP22_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, deltain, deltaout, qout, lipschitz):
    bt_std = get_bt_std(2 * Qin, n, Q, N, B_g, q_ks, B_ks)
    lip_var = bt_std**2
    lip_var *= (lipschitz / deltain)**2
    bt_var = var_ACC(n, N, Q, B_g, encstd) * (qout / Q)**2
    bt_var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    print(
        f"lipschitz contrib = {lip_var}, acc contrib = {bt_var / deltaout**2}")
    input_var = CKKS_input_err_std_predef(n, Qin)**2  * deltain**-2 * lipschitz**2
    return sqrt(lip_var + bt_var / deltaout**2), sqrt(lip_var + bt_var * deltaout**-2 + input_var)



# NOTE: done
# standard params for output range (-2, 2) and Lipschitz constant 1
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1) |-> std = 0.0757
# small params for output range (-2, 2) and Lipschitz constant 1
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1) |-> std = 0.0629
#   (2**12, n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**8, 2**18, 2**20, 1) |-> std = 0.0564
# N = 2**11 case:
#   (2**11, n20, Q26, 2**10, 2**4, 2**20, 2**5, 2**8, 2**18, 2**20, 1)  |-> std = 0.0675
#
# standard params
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0208
# small params
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0174
#   (2**12, n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**8, 2**18, 2**20, 1/4) |-> 0.0184
# N = 2**11:
#   (2**11, n20, Q26, 2**10, 2**4, 2**20, 2**5, 2**8, 2**18, 2**20, 1)    |-> 0.0195
def FDFB_compress_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, deltain, deltaout, qout, lipschitz):
    bt_std = get_bt_std(Qin, n, Q, N, B_g, q_ks, B_ks)
    beta = ceil(normbnd * bt_std)
    # print(f"beta = {beta}")
    lip_var = bt_std ** 2
    lip_var += 1 / 4  # randomnized rounding error of index
    lip_var *= ((N - 1) / (N / 2 - 2*beta) / deltain *
                lipschitz) ** 2   # this is the main contributor; q = 2N
    bt_var = var_ACC(n, N, Q, B_g, encstd)
    bt_var += 1 / 4  # randomnized rounding error of fval
    bt_var *= (qout / Q)**2
    bt_var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    # bt_var += var_MS(Q, q_ks, N)
    # bt_var += var_KS(N, q_ks, B_ks, encstd)
    # bt_var *= (qout / q_ks)**2
    # bt_var += var_MS(q_ks, qout, n)
    print(
        f"lipschitz contrib = {lip_var}, acc contrib = {bt_var / deltaout**2}")
    input_var = CKKS_input_err_std_predef(n, Qin)**2 * deltain**-2 * lipschitz**2
    return sqrt(lip_var + bt_var / deltaout**2), sqrt(lip_var + bt_var * deltaout**-2 + input_var)


# NOTE: done
# standard params for output range (-2, 2)
#   (n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.000291
# small params for output range (-2, 2)
#   (n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25)  |-> 0.000246
#   (n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20)  |-> 0.00111
# N = 2**11 is not considered
#
# standard params
#   (2**11, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**7, 2**23, 2**25, 1/4) |-> 0.0169
# small params
#   (2**11, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**7, 2**23, 2**25, 1/4) |-> 0.0144
#   (2**11, n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**20, 2**5, 2**7, 2**18, 2**20, 1/4) |-> 0.0162
def FDFB_cancelsign_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, deltain, deltaout, qout, lipschitz):
    core_var = FDFB_cancelsign_var(qout, n, Q, N, B_g, q_pk, B_pk)
    core_var_old = core_var
    core_var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    print(
        f"core var contrib = {core_var_old * deltaout**-2}, common contrib = {(core_var - core_var_old) * deltaout**-2}")
    input_var = CKKS_input_err_std_predef(n, Qin)**2 * lipschitz**2 * deltain**-2
    return sqrt(core_var) / deltaout, sqrt(core_var * deltaout**-2 + input_var)


# NOTE: done
# NOTE: no multi-value BTS here
# standard params for output range (-2, 2)
#   (n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.000409
# small params for output range (-2, 2)
#   (n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.000346
#   (n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20) |-> 0.00113
# N = 2**10, q = 2**11 case:
#   (n20, Q26, 2**10, 2**5, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20)  |-> 0.0139
#   (n20, Q26, 2**10, 2**4, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20)  |-> 0.00767
#
# standard params
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0085
# small params
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0075
def FDFB_select_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, deltain, deltaout, qout, lipschitz):
    core_var = FDFB_select_var(qout, n, Q, N, B_g, q_pk, B_pk)
    core_var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    input_var = CKKS_input_err_std_predef(n, Qin)**2  * deltain**-2 * lipschitz**2
    return sqrt(core_var) / deltaout, sqrt(core_var * deltaout**-2 + input_var)


# NOTE: done
# NOTE: no multi-value BTS here
# standard params for output range (-2, 2)
#   (n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.000577
# small params for output range (-2, 2)
#   (n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.000487
#   (n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20) |-> 0.00117
# N = 2**10, q = 2**11 case:
#   (n20, Q26, 2**10, 2**4, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20)  |-> 0.0109
#
# standard params
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0085
# small params
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0075
def FDFB_select_alt_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, deltain, deltaout, qout, lipschitz):
    core_var = FDFB_select_alt_var(qout, n, Q, N, B_g, q_pk, B_pk)
    core_var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    input_var = CKKS_input_err_std_predef(n, Qin)**2 * deltain**-2 * lipschitz**2
    return sqrt(core_var) / deltaout, sqrt(core_var * deltaout**-2 + input_var)


# NOTE: done
# NOTE: need intermediate modulus p_mid
# standard params for output range (-2, 2)
#   (2**8, n35, Q53, 2**11, 2**27, 2, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25) |-> 0.0385
# small params for output range (-2, 2)
#   (2**8, n25, Q53, 2**11, 2**27, 2, 2**25, 2**5, 2**25, 2**5, 2**23, 2**25)  |-> 0.0330
#   (2**8, n20, Q53, 2**11, 2**27, 2, 2**20, 2**5, 2**20, 2**5, 2**18, 2**20)  |-> 0.0304
# N = 2**11: no suitable params
#
# standard params
#   (2**12, 2**8, n35, Q53, 2**11, 2**27, 2, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0394
# small params
#   (2**12, 2**8, n25, Q53, 2**11, 2**27, 2, 2**25, 2**5, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0338
def FDFB_preselect_CKKS_std(Qin, p_mid, n, Q, N, B_g, B_mv, q_ks, B_ks, q_pk, B_pk, deltain, deltaout, qout, lipschitz):
    d_mv = ceil(log(p_mid, B_mv))
    var = var_ACC(n, N, Q, B_g, encstd)  # use sgn to batch select Q/2p*sgn*B^i
    # use LWE to RLWE packing to get Q/2p*sgn*B^i*TV0. TODO: rounding error here
    # var += (Q / q_pk)**2 * (var_KS(N, q_pk, B_pk, encstd) +
    #                         var_MS(Q, q_pk, N)) + var_MS(q_pk, Q, N)
    var += (Q / q_pk)**2 * var_MS(Q, q_pk, N)
    var += var_PK_scaled(N, q_pk, Q, B_pk, encstd)
    var = B_mv**2 * N * var * d_mv / 4 + \
        (Q / p_mid)**2 * 1/4 * 2  # p_mid ptxt contains (randomnized) rounding error
    var += var_ACC(n, N, Q, B_g, encstd)
    var *= (qout / Q)**2
    # common part, KS + MS
    var += common_part_var(qout, n, Q, N, q_ks, B_ks)
    lip_var = CKKS_input_err_std_predef(n, Qin)**2 * deltain**-2 * lipschitz**2
    return sqrt(var) / deltaout, sqrt(var * deltaout**-2 + lip_var)


# PSO(first half) mapping: 0 -> beta, q/2-1 -> q/2 - beta
# y = (q/2-2beta)/(q/2-1)x+beta
# PSE(first half) mapping: 0 -> q/4+beta, q/2-1 -> 3q/4 - beta
# y = (q/2-2beta)/(q/2-1)x+q/4+beta
# NOTE: done
# standard params for input range (-8,8), output range (-2, 2) and Lipschitz constant 1
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1) |-> std = 0.0358
# small params for input range (-8,8), output range (-2, 2) and Lipschitz constant 1
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1) |-> std = 0.0300
#   (2**12, n20, Q53, 2**11, 2**27, 2**20, 2**5, 2**8, 2**18, 2**20, 1) |-> std = 0.0270
# N = 2^10 case:
#   (2**11, n20, Q26, 2**10, 2**5, 2**20, 2**5, 2**8, 2**18, 2**20, 1)  |-> std = 0.0385
#
# (special or not doesn't matter)
# standard params
#   (2**12, n35, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0123
# small params
#   (2**12, n25, Q53, 2**11, 2**27, 2**25, 2**5, 2**8, 2**23, 2**25, 1/4) |-> 0.0106
def FDFB_comp_CKKS_std(Qin, n, Q, N, B_g, q_ks, B_ks, deltain, deltaout, qout, lipschitz, special=False):
    bt_std = get_bt_std(Qin, n, Q, N, B_g, q_ks, B_ks)
    beta = ceil(normbnd * bt_std)
    lip_var = bt_std ** 2 + 1/4
    lip_var *= ((N-1)/(N-2*beta)*lipschitz/deltain)**2  # q = 2N
    bt_var = (2 if not special else 1) * (var_ACC(n, N, Q, B_g, encstd) + 1/4)
    bt_var *= (qout / Q)**2
    bt_var += common_part_var(qout, n, Q, N, q_ks, B_ks)  # lazy KS + MS
    print(
        f"lipschitz contrib = {lip_var}, acc contrib = {bt_var / deltaout**2}, max compensation err = {0.5 / deltain * lipschitz}")
    input_var = CKKS_input_err_std_predef(n, Qin)**2 * deltain**-2 * lipschitz**2
    return sqrt(lip_var + bt_var / deltaout**2), sqrt(lip_var + bt_var * deltaout**-2 + input_var)


# NOTE: 
# standard params
#   (n35, Q53, 2**11, 2**27, 2**35, 2**5, 2**35, 2**5) |-> log2(std*bnd)+1 = 24.91 ## 10
#   (n35, Q53, 2**11, 2**18, 2**35, 2**5, 2**35, 2**5) |-> log2(std*bnd)+1 = 16.33 ## 18
#   (n35, Q53, 2**11, 2**14, 2**35, 2**5, 2**35, 2**5) |-> log2(std*bnd)+1 = 12.90  ## 22
def ReLU_std(n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, qin):
    var = (Q / q_pk)**2 * var_MS(qin, q_pk, N)
    var += var_PK_scaled(n, q_pk, Q, B_pk, encstd)  # NOTE: n here
    var += var_ACC(n, N, Q, B_g, encstd)
    # note that qout = q_ks
    var *= (q_ks / Q)**2
    var += common_part_var(q_ks, n, Q, N, q_ks, B_ks)
    return sqrt(var)


# qout defaults to q_ks
# standard params
#   input range = (-8, 8)
#   (n35, Q53, 2**11, 2**27, 2**35, 2**5, 2**35, 2**5, 2**28, 2**24) |-> std = 7.32e-7
#   (n35, Q53, 2**11, 2**18, 2**35, 2**5, 2**35, 2**5, 2**28, 2**24) |-> std = 3.75e-8
#   (n35, Q53, 2**11, 2**14, 2**35, 2**5, 2**35, 2**5, 2**28, 2**24) |-> std = 1.15e-8
#
# standard params
#   (n35, Q53, 2**11, 2**27, 2**35, 2**5, 2**35, 2**5, 2**24, 2**20) |-> std = 0.000779
def ReLU_CKKS_std(n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, qin, delta):
    core_var = ReLU_std(n, Q, N, B_g, q_ks, B_ks, q_pk, B_pk, qin)
    lip_var = CKKS_input_err_std_predef(n, qin) ** 2 * delta**-2  # lipschitz = 1
    return sqrt(core_var) / (delta * q_ks / qin), sqrt(core_var * (delta * q_ks / qin)**-2 + lip_var)


if __name__ == "__main__":
    # the following code finds the proper value of Bg for different qin of HomDecomp
    bases = [2**27, 2**18, 2**14]
    # HomFloor
    print("HomFloor:")
    i = 13
    curbase_idx = 0
    while True:
        base = bases[curbase_idx]
        changed = False
        bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
        while (bt_std * sqrt(2) * normbnd >= 128 or \
            sqrt(bt_std ** 2 * (1 + 2**-8) + var_MS(2**30, 2**26, n35)) * normbnd >= 128) and curbase_idx < 3:
            curbase_idx += 1
            if curbase_idx < 3:
                base = bases[curbase_idx]
                bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
            changed = True
        if curbase_idx == 3:
            print(f"stop at log(q) = {i}")
            break
        if changed:
            print(f"log(q) = {i}, log(Bg) = {log2(bases[curbase_idx])}")
        i = i+1
    
    # HomFloorAlt
    print("\nHomFloorAlt:")
    i = 13
    curbase_idx = 0
    while True:
        base = bases[curbase_idx]
        changed = False
        bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
        while (bt_std * sqrt(2) * normbnd >= 512 or \
            sqrt(bt_std ** 2 * 2**-10 + var_MS(2**30, 2**25, n35)) * normbnd >= 64) and curbase_idx < 3:
            curbase_idx += 1
            if curbase_idx < 3:
                base = bases[curbase_idx]
                bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
            changed = True
        if curbase_idx == 3:
            print(f"stop at log(q) = {i}")
            break
        if changed:
            print(f"log(q) = {i}, log(Bg) = {log2(bases[curbase_idx])}")
        i = i+1
    
    # HomReduce
    print("\nHomDecomp-Reduce:")
    i = 13
    curbase_idx = 0
    while True:
        base = bases[curbase_idx]
        changed = False
        bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
        while sqrt(bt_std ** 2 * 2**-8 + var_MS(2**30, 2**25, n35)) * normbnd + 64 >= 128 and curbase_idx < 3:
            curbase_idx += 1
            if curbase_idx < 3:
                base = bases[curbase_idx]
                bt_std = get_bt_std(2**i, n35, Q53, 2**11, base, 2**35, 2**5)
            changed = True
        if curbase_idx == 3:
            print(f"stop at log(q) = {i}")
            break
        if changed:
            print(f"log(q) = {i}, log(Bg) = {log2(bases[curbase_idx])}")
        i = i+1
    
    # decomp usign fdfb
    print("\nHomDecomp-FDFB:")
    i = 13
    curbase_idx = 0
    while True:
        base = bases[curbase_idx]
        changed = False
        eval_std = FDFB_compress_CKKS_std(2**12, n35, Q53, 2**11, base, 2**35, 2**5, 1, 1, 2**i, 1)[0]
        while sqrt(eval_std ** 2 * 2**-10 + var_MS(2**30, 2**25, n35)) * normbnd >= 64 and curbase_idx < 3:
            curbase_idx += 1
            if curbase_idx < 3:
                base = bases[curbase_idx]
                eval_std = FDFB_compress_CKKS_std(2**12, n35, Q53, 2**11, base, 2**35, 2**5, 1, 1, 2**i, 1)[0]
            changed = True
        if curbase_idx == 3:
            print(f"stop at log(q) = {i}")
            break
        if changed:
            print(f"log(q) = {i}, log(Bg) = {log2(bases[curbase_idx])}")
        i = i+1