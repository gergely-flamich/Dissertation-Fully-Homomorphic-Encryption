
def key_gen(security_param):

    # Choose p and q, remembering with the constraint that |p| = |q|
    magnitude = security_param

    p = random_prime(2^magnitude - 1, True, 2^(magnitude - 1))
    q = random_prime(2^magnitude - 1, True, 2^(magnitude - 1))

    # Generate N = pq

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

    N = public_key

    r = 83

    while True:
        r = mod(Zmod(N).random_element(), N ^ 2)

        if gcd(r, N) == 1:
           break

    m = plaintext

    return (1 + N) ^ m * mod(r ^ N, N^2)


def decrypt(private_key, ciphertext):

    N, phi, phi_inv = private_key

    c = mod(ciphertext, N ^ 2)

    return mod(int(c ^ phi - 1) / N * phi_inv, N)
