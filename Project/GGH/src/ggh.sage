from sage.modules.free_module_integer import IntegerLattice
import sage.matrix.matrix_integer_dense_hnf as hnf
from sage.matrix.matrix_space import MatrixSpace

def key_gen(security_parameter, gen_param = 500, hadamard_threshold = 0.7, delta=4):
    R = MatrixSpace(ZZ, security_parameter)

    U = identity_matrix(ZZ, security_parameter)

    V = None

    while hadamard_ratio(V) < hadamard_threshold:

        V = R.random_element(x = -security_parameter,
                             y = security_parameter + 1)

    # for i in xrange(round(security_parameter * 10)):
    #     c1 = ZZ.random_element(security_parameter)
    #     c2 = ZZ.random_element(security_parameter)

    #     if c1 == c2:
    #         continue

    #     sign = -1 if ZZ.random_element(2) == 1 else 1

    #     U.add_multiple_of_column(c1, c2, sign)

    # W = U * V

    W, _ = hnf.hnf(V)

    W_inv = W.inverse()

    print("Hadamard ratio of PubK: {0:.5}".format(float(hadamard_ratio(W))))
    print("Hadamard ratio of PrivK: {0:.5}".format(float(hadamard_ratio(V))))

    return (W, delta), (V, W_inv)

def hadamard_ratio(basis):
    if basis is None:
        return 0

    vol = abs(det(basis))
    block = product(map(norm, basis))

    return n((vol/block)^(1/basis.rank()))

def encrypt(pubkey, plaintext):
    W, delta = pubkey

    m = matrix(map(ord, plaintext))

    r = MatrixSpace(ZZ, 1, W.rank()).random_element(x = -delta, y = delta + 1)

    print("perturbation: {}".format(r))

    return m * W + r


def decrypt(privkey, ciphertext):
    V, W_inv = privkey

    appr_vec = ciphertext * V.inverse()

    appr_vec = appr_vec.apply_map(round)

    latt_vec = appr_vec * V * W_inv

    print(latt_vec)

    return ''.join(map(chr, list(latt_vec[0])))
