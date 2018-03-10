def key_gen(security_param):
    """
    Key generation algorithm for the Paillier encryption scheme.

    @param security_param - Bit length of the two factors of N

    @return pk, sk - pk is simply N, sk is N, phi(N) and its inverse.
    """

    # Choose p and q, remembering with the constraint that |p| = |q|
    magnitude = security_param

    p = random_prime(2^magnitude - 1, True, 2^(magnitude - 1))

    q = None

    while q is None or q == p:
        q = random_prime(2^magnitude - 1, True, 2^(magnitude - 1))

    # Calculate N = pq

    N = p * q

    # Calculate Phi(N) = (p - 1)(q - 1)
    # Note: since we require a proof that p and q are prime, we can be
    # absolutely certain that the above formula works.

    phi = (p - 1) * (q - 1)

    # Calculate the inverse of phi mod N. Note that by a theorem from ITMC, we know
    # that gcd(phi(N), N) = 1, so phi has an inverse mod N.

    phi_inv = phi.inverse_mod(N)

    return (N, (N, phi, phi_inv))


def encrypt(public_key, plaintext):
    """
    Encryption function. Note that this function is non-deterministic.
    """

    N = public_key

    r = 0

    # Pick a random element of Z_N^*
    while True:
        r = mod(Zmod(N).random_element(), N ^ 2)

        if gcd(r, N) == 1:
           break

    m = plaintext

    # Calculate the group isomorphism from Z_N x Z_N^* to Z_{N^2}
    return (1 + N) ^ m * mod(r ^ N, N^2)


def decrypt(private_key, ciphertext):
    """
    Decryption function.
    """

    N, phi, phi_inv = private_key

    c = mod(ciphertext, N ^ 2)

    # Invert the group isomorphism based on phi(N) being the trapdoor.
    return mod(int(c ^ phi - 1) / N * phi_inv, N)
