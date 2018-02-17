def key_gen(N = 5, t = 32, num_retry = 10):
    Zx = PolynomialRing(ZZ, 'x')

    n = 2^N

    f = Zx(x^n + 1)

    R = Zx.quotient(f)

    lower_bound = -2^(t-1)
    upper_bound = 2^(t-1)

    counter = 0

    for i in xrange(num_retry):

        counter = i + 1

        v = R.random_element(x=lower_bound, y=upper_bound)

        V = v.matrix()

        (d, w, _) = xgcd(v.lift(), f)

        w = R(w)

        Zd = Zmod(d)

        w_prime = list(w)

        try:
            r = Zd(w[0]) * Zd(w[1]).inverse_of_unit()
        except ArithmeticError:
            continue

        if r^n == -1:
            break

    print("HNF in correct form after {} tries.".format(counter))

    w_i = filter(lambda x: abs(x) % 2 == 1, w_prime)

    return (N, centered_mod(r, d), d), (d, w_i[0], N, w, v)

def hadamard_ratio(basis):
    if basis is None:
        return 0

    vol = abs(det(basis))
    block = product(map(norm, basis))

    return n((vol/block)^(1/basis.rank()))

def construct_hnf(n, d, r):

    first_row = [d] + [0]*(n - 1)

    r_prime = mod(r, d)

    first_col = [-centered_mod(r_prime^i, d) for i in xrange(1, n)]

    return block_matrix([[matrix(ZZ, first_row)],
                         [matrix(ZZ, first_col).transpose(), identity_matrix(n-1,n-1)]],
                        subdivide = False)

def encrypt(pk, plaintext, u=None):

    if plaintext != 0 and plaintext != 1:
        raise ArithmeticError("Plaintext must be a bit!")

    N, r, d = pk

    Zx = PolynomialRing(ZZ, 'x')

    n = 2^N

    f = Zx(x^n + 1)

    R = Zx.quotient(f)

    p = 10/n

    if u is None:
        u_coeff_vec = [0 if abs(candidate) > p else sgn(candidate) for candidate in random_vector(RR, n)]
        u = R(u_coeff_vec)

    print("Perturbation is {} dimensional and has {} non-zero entries.".format(n, sum(map(abs, u))))

    m_prime = R([plaintext])

    a = 2 * u + m_prime

    return centered_mod(a.lift()(r), d)

def inefficient_encrypt(pk, plaintext):
    if plaintext != 0 and plaintext != 1:
        raise ArithmeticError("Plaintext must be a bit!")

    N, r, d = pk

    Zx = PolynomialRing(ZZ, 'x')

    n = 2^N

    f = Zx(x^n + 1)

    R = Zx.quotient(f)

    p = 10/n

    u_coeff_vec = [0 if abs(candidate) > p else sgn(candidate) for candidate in random_vector(RR, n)]

    u = R(u_coeff_vec)

    a = 2 * u + plaintext

    B = construct_hnf(n, d, r)

    Binv = B.inverse()

    a_prime = vector(a)

    amodB = (vector(map(round, (a_prime * Binv))) * B)

    ciphertext = a_prime - amodB

    print("a(r) mod d: {}. Is within required bounds ({}): {}".format(ciphertext[0], d//2, abs(ciphertext[0])< d//2))

    return ciphertext, u, amodB


def decrypt(sk, ciphertext):
    w_i, d, _, _, _ = sk

    b = centered_mod(ciphertext * w_i, d).mod(2)

    return b

def inefficient_decrypt(sk, ciphertext):
    _, d, N, w, v = sk

    W = w.matrix()
    V = v.matrix()

    print(ciphertext)

    c_prime = (ciphertext * W)/d

    print(c_prime)

    c_prime = vector(map(round, c_prime))
    print(c_prime)

    a = ciphertext - c_prime * V
    print(a)

    return a[0].mod(2)

def centered_mod(x, modulus):
    x = ZZ(x)
    res = x.mod(modulus)
    half_mod = (modulus + 1)//2
    return  res if res < half_mod else res - modulus

