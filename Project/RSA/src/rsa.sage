def key_gen(security_param):
    """
    Implements the key generation algorithm for textbook RSA.

    @param security_param - the bit-length of the primes we should choose.
    The modulus will therefore be of length 2 * security_param
    """
    m = security_param

    # Pick two different primes
    p = random_prime(2 ^ (m + 1) - 1, True, 2 ^ m)

    while True:
        q = random_prime(2 ^ (m + 1) - 1, True, 2 ^ m)

        if q is not p:
            break

    N = p * q

    phi = (p - 1) * (q - 1)

    i = 0

    # Pick a random element of Z_N* by checking that the
    # gcd(i, phi(N)) = 1.
    while True:

        i = ZZ.random_element(phi)

        d, j, _ = xgcd(i, phi)

        if d == 1:
            break

    return (N, i), (N, j)


def encrypt(pubkey, message):
    """
    Encryption function for textbook RSA.

    @param pubkey - Public key under which we are performing the encryption.
    Should be of the form (N = modulus, i = exponent)

    @param message - Plaintext to be encrypted. We assume that the message is smaller than the modulus.
    """
    N, i = pubkey

    Zn = Zmod(N)

    # For convenience, the output will be a member of Sage's Z_N ring instead of a standard integer.
    return Zn(message) ^ i


def decrypt(privkey, ciphertext):
    """
    Decryption function for textbook RSA.

    @param privkey - Private key under which we are performing the decryption.
    Should be of the form (N = modulus, j = inverse of the exponent i modulo N)

    @param ciphertext - Ciphertext to be decrypted. We assume that it is of the form m^i mod N.
    """

    N, j = privkey

    Zn = Zmod(N)

    return Zn(ciphertext) ^ j
