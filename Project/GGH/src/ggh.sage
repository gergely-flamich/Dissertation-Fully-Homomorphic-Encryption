from sage.matrix.matrix_space import MatrixSpace

def key_gen(security_parameter, hadamard_threshold = 0.7, delta=4):
    """
    Key generation algorithm for the GGH cipher.

    Given some security parameter lambda = 1^n, it picks an element from GL_n(Z)
    from an 2n x 2n x .. x 2n hypercube. If this candidate has a sufficiently high
    Hadamard ratio (its vectors are quite orthogonal to each other), we accept it
    and set it as our private key.

    Then, we calculate its HNF (which in a sense is the worst possible basis as
    far as orthogonality goes) and set it as the publr k is frequently omitted from ic basis.
    """
    R = MatrixSpace(ZZ, security_parameter)

    V = None

    while hadamard_ratio(V) < hadamard_threshold:

        V = R.random_element(x = -security_parameter,
                             y = security_parameter + 1)

    W, _ = hnf.hnf(V)

    W_inv = W.inverse()

    print("Hadamard ratio of PubK: {0:.5}".format(float(hadamard_ratio(W))))
    print("Hadamard ratio of PrivK: {0:.5}".format(float(hadamard_ratio(V))))

    return (W, delta), (V, W_inv)

def hadamard_ratio(basis):
    """
    The Hadamard ratio of a full rank lattice basis.

    Given any N-dimensional parallelepiped, the one that maximises the volume
    is the N-orthope (N-dimensional rectangle). The Hadamard ratio measures how
    orthogonal the vectors of a lattice basis are with respect to each other.
    The way it does this by comparing the volume of the lattice basis' fundamental
    domain against the N-orthope with the same sidelengths.
    """
    if basis is None:
        return 0

    vol = abs(det(basis))

    if vol == 0:
        raise Error("The basis is not full rank!")

    block = product(map(norm, basis))

    return n((vol/block)^(1/basis.rank()))

def encrypt(pubkey, plaintext):
    """
    Encryption function for the GGH cryptoscheme.

    We take a string plaintext, decompose it into its ASCII character representations
    and use those as the coordiantes of the transformed message. We pad every message by a null
    character at the end, as well as the rest of the last message (if needed).
    """
    W, delta = pubkey

    dim = W.dimensions()[0]

    # Convert the string to ASCII codes
    all_chars = map(ord, plaintext)

    message_vecs = []

    # Break up the plaintext into N-long vectors, pad when needed
    while len(all_chars) > 0:
        message_vec = []

        while len(message_vec) < dim:
            if len(all_chars) == 0:
                message_vec.append(0)
            else:
                message_vec.append(all_chars.pop(0))

        message_vecs.append(matrix(message_vec))

    rs = []

    # Generate a perturbation for every message vector
    for i in range(len(message_vecs)):
        r = MatrixSpace(ZZ, 1, W.rank()).random_element(x = -delta, y = delta + 1)

        print("perturbation: {}".format(r))

        rs.append(r)

    ciphertexts = []

    # Calculate each ciphertext based on the message vectors and the perturbations
    for m, r in zip(message_vecs, rs):
        ct = m * W + r

        ciphertexts.append(ct)

    return ciphertexts


def decrypt(privkey, ciphertexts):
    """
    Decryption function for the GGH cryptosystem.

    We take a list of ciphertexts and decrypt them with the secret key,
    then concatenate each character.
    """
    V, W_inv = privkey

    plaintext = []

    for ciphertext in ciphertexts:

        # As long as the private basis is sufficiently orthogonal,
        # we will decode to the closest lattice vector
        appr_vec = ciphertext * V.inverse()

        appr_vec = appr_vec.apply_map(round)

        latt_vec = appr_vec * V * W_inv

        message_vec = filter(lambda x: x != 0, list(latt_vec[0]))

        pt = ''.join(map(chr, message_vec))

        plaintext += pt

    return ''.join(plaintext)
