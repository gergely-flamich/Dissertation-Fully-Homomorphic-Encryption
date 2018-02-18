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

        try:
            r = Zd(w[0]) * Zd(w[1]).inverse_of_unit()
        except ArithmeticError:
            continue

        if r^n == -1:
            break

    print("HNF in correct form after {} tries.".format(counter))

    w_i = filter(lambda x: x % 2 == 1, w)

    d = ZZ(d)
    r = ZZ(r)

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

    p = 0.5 if p > 0.5 else p

    if u is None:
        u_coeff_vec = [0 if abs(candidate) > p else sgn(candidate) for candidate in random_vector(RR, n)]
        u = R(u_coeff_vec)

    print("Perturbation is {} dimensional and has {} non-zero entries.".format(n, sum(map(abs, u))))
    print("Namely, u(x) = {}".format(u))

    a = 2 * u

    return (a.lift()(r) + plaintext).mod(d)

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

    amodB = (vector(map(c_round, (a_prime * Binv))) * B)

    ciphertext = a_prime - amodB

    print("a(r) mod d: {}. Is within required bounds ({}): {}".format(ciphertext[0], d//2, abs(ciphertext[0])< d//2))

    return ciphertext, u


def decrypt(sk, ciphertext):
    d, w_i, _, _, _ = sk

    e = (ciphertext * w_i).mod(d)

    b = e.mod(2)

    if d % 2 == 1 and 2*e >= d:
        b = 1 - b

    return b

def inefficient_decrypt(sk, ciphertext):
    d, w_i, N, w, v = sk

    W = w.matrix()
    V = v.matrix()

    print("Received ciphertext: {}".format(ciphertext))

    c_prime = (ciphertext * W)/d

    print("c * V^-1 = {}".format(c_prime))

    c_prime = vector(map(c_round, c_prime))
    print("Rounded: {}".format(c_prime))

    a = ciphertext - c_prime * V
    print("Fractional part is therefore {}".format(a))

    return a[0].mod(2)

def centered_mod(x, modulus):
    x = ZZ(x)
    res = x.mod(modulus)
    half_mod = (modulus + 1)//2
    return  res if res < half_mod else res - modulus


def c_round(x):
    return floor(x + 0.5)

def test_ct(pk, sk, n = 100):

    zero_counter = 0
    one_counter = 0

    for i in xrange(n):
        if decrypt(sk, encrypt(pk, 0)) != 0:
            zero_counter += 1
        if decrypt(sk, encrypt(pk, 1)) != 1:
            one_counter += 1

    print("Encryptions of 0 failed {}/{}".format(zero_counter, n))
    print("Encryptions of 1 failed {}/{}".format(one_counter, n))
