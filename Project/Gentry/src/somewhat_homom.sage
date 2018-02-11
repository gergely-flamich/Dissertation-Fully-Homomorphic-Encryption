def ideal_gen(security_param, ring, ideal_i):
    pass

def key_gen(security_param, l_samp = None):
    """
    Generates the key for Craig Gentry's somewhat homomorphic scheme.

    security_param - dimension of the ideal lattice
    l_samp - bound on the norm of the random element from Z[x]/f
    """

    # Just have a shorter name for the security parameter internally
    n = security_param

    # Create the quotient ring Z[x]/f with f(x) = x^n - 1.
    # TODO: choose different f, with f irreducible. In particular, attempt Thm 7.4.2
    qring = PolynomialRing(ZZ, 'x').quotient(x^n - 1)

    # s will be the generator for the ideal I.
    s = qring(2)

    print(s)

    return v

def samp(ideal_i, plaintext):
    pass

def encrypt(pk, plaintext):
    pass

def decrypt(sk, ciphertext):
    pass

def evaluate(pk, circuit, ct_vec):
    pass

def add(pk, ct_1, ct_2):
    return ct_1 + ct_2

def mult(pk, ct_1, ct_2):
    return ct_1 * ct_2
