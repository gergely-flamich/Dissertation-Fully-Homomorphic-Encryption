def key_gen(security_param):

    m = security_param

    p = random_prime(2 ^ (m + 1) - 1, True, 2 ^ m)

    while True:
        q = random_prime(2 ^ (m + 1) - 1, True, 2 ^ m)

        if q is not p:
            break

    N = p * q

    phi = (p - 1) * (q - 1)

    i = 0

    while True:

        i = ZZ.random_element(phi)

        d, j, _ = xgcd(i, phi)

        if d == 1:
            break

    return (N, i), (N, j)


def encrypt(pubkey, message):

    N, i = pubkey

    Zn = Zmod(N)

    return Zn(message) ^ i


def decrypt(privkey, ciphertext):

    N, j = privkey

    Zn = Zmod(N)

    return Zn(ciphertext) ^ j
