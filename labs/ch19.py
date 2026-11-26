"""Chapter 19 laboratory: the primary decomposition theorem.

Exact arithmetic over Q (fractions.Fraction).  A polynomial is a list of
Fractions, constant term first, with no trailing zero.  A matrix is a
list of rows.  The model is the chapter's running example: the operator
M on Q^5 with minimal polynomial (x-2)^2 (x+1) (x^2+1).  The lab builds
the projections E_j = h_j(M) from a Bezout identity, and checks every
claim of the theorem and of the diagonalizable-plus-nilpotent corollary.
"""
from fractions import Fraction

Q = Fraction


# ---------------------------------------------------------------- polys
def trim(p):
    p = [Fraction(c) for c in p]
    while p and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [0] * (n - len(a))
    b = list(b) + [0] * (n - len(b))
    return trim([x + y for x, y in zip(a, b)])


def pscale(a, c):
    return trim([c * x for x in a])


def pmul(a, b):
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def ppow(a, r):
    out = [Fraction(1)]
    for _ in range(r):
        out = pmul(out, a)
    return out


def pdivmod(a, b):
    """Quotient and remainder of a by b over Q, with b nonzero."""
    a = trim(a)
    b = trim(b)
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        shift = len(a) - len(b)
        c = a[-1] / b[-1]
        q[shift] = c
        a = padd(a, pscale([0] * shift + list(b), -c))
    return trim(q), a


def pinverse_mod(f, m):
    """The polynomial g with deg g < deg m and f g = 1 mod m.

    Requires gcd(f, m) = 1.  Extended Euclid with a running coefficient.
    """
    r0, r1 = trim(m), pdivmod(f, m)[1]
    t0, t1 = [], [Fraction(1)]
    while r1:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, padd(t0, pscale(pmul(q, t1), -1))
    assert len(r0) == 1, 'f and m are not relatively prime'
    return pscale(pdivmod(t0, m)[1], 1 / r0[0])


def pstr(p):
    if not p:
        return '0'
    terms = []
    for k in range(len(p) - 1, -1, -1):
        c = p[k]
        if c == 0:
            continue
        mono = '' if k == 0 else ('x' if k == 1 else 'x^%d' % k)
        if k == 0:
            coef = str(c)
        elif c == 1:
            coef = ''
        elif c == -1:
            coef = '-'
        else:
            coef = '(%s)' % c
        terms.append(coef + mono)
    return ' + '.join(terms)


# ------------------------------------------------------------- matrices
def ident(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def zeros(n):
    return [[Fraction(0)] * n for _ in range(n)]


def mmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def madd(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def mscale(a, c):
    return [[c * x for x in r] for r in a]


def peval(p, m):
    """p(M) by Horner's rule."""
    n = len(m)
    out = zeros(n)
    for c in reversed(p):
        out = madd(mmul(out, m), mscale(ident(n), c))
    return out


def rref(rows):
    """Reduced row echelon form and pivot columns, exact."""
    a = [list(map(Fraction, r)) for r in rows]
    piv, r = [], 0
    for c in range(len(a[0]) if a else 0):
        s = next((i for i in range(r, len(a)) if a[i][c] != 0), None)
        if s is None:
            continue
        a[r], a[s] = a[s], a[r]
        a[r] = [x / a[r][c] for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c] != 0:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[r])]
        piv.append(c)
        r += 1
        if r == len(a):
            break
    return a, piv


def rank(rows):
    return len(rref(rows)[1]) if rows else 0


def nullspace(a):
    """Basis (list of vectors) of {x : A x = 0}."""
    n = len(a[0])
    red, piv = rref(a)
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for row, pc in zip(red, piv):
            v[pc] = -row[f]
        basis.append(v)
    return basis


def columns(a):
    return [list(c) for c in zip(*a)]


def sum_dim(*bases):
    """Dimension of the sum of the spans of the given bases."""
    return rank([v for b in bases for v in b])


def restrict(m, basis):
    """Matrix of x -> M x on span(basis), written in that basis."""
    k = len(basis)
    aug = [[basis[j][i] for j in range(k)] +
           [sum(m[i][t] * basis[c][t] for t in range(len(m)))
            for c in range(k)] for i in range(len(m))]
    red, piv = rref(aug)
    assert piv[:k] == list(range(k)), 'basis is not independent'
    return [[red[i][k + c] for c in range(k)] for i in range(k)]


def minimal_polynomial(m):
    """Monic polynomial of least degree annihilating M."""
    n = len(m)
    powers = [ident(n)]
    for _ in range(n):
        powers.append(mmul(powers[-1], m))
    flat = [[x for row in p for x in row] for p in powers]
    for d in range(1, n + 1):
        sol = nullspace(columns(flat[:d + 1]))
        if sol:
            c = sol[0]
            return trim([x / c[d] for x in c])
    raise AssertionError('no annihilating polynomial of degree <= n')


def inverse(a):
    n = len(a)
    red, piv = rref([list(r) + ident(n)[i] for i, r in enumerate(a)])
    assert piv[:n] == list(range(n)), 'matrix is singular'
    return [r[n:] for r in red]


# ------------------------------------------------- primary decomposition
def primary_data(factors):
    """Return m, the powers p_i^r_i, the cofactors f_i, and the h_i.

    Here h_i = g_i f_i reduced modulo m, where g_i f_i = 1 modulo
    p_i^r_i and m = prod p_i^r_i.  Then sum h_i = 1 exactly.
    """
    powers = [ppow(p, r) for p, r in factors]
    m = [Fraction(1)]
    for q in powers:
        m = pmul(m, q)
    fs, hs = [], []
    for q in powers:
        f = pdivmod(m, q)[0]
        g = pinverse_mod(f, q)
        fs.append(f)
        hs.append(pdivmod(pmul(g, f), m)[1])
    return m, powers, fs, hs


def main():
    p1, p2, p3 = trim([-2, 1]), trim([1, 1]), trim([1, 0, 1])
    factors = [(p1, 2), (p2, 1), (p3, 1)]
    A = [[2, 0, 0], [1, 2, 0], [0, 0, -1]]
    C = [[0, -1], [1, 0]]
    n = 5
    M = zeros(n)
    for i in range(3):
        for j in range(3):
            M[i][j] = Q(A[i][j])
    for i in range(2):
        for j in range(2):
            M[3 + i][3 + j] = Q(C[i][j])

    # 1. minimal polynomial of the running example
    m, powers, fs, hs = primary_data(factors)
    mp = minimal_polynomial(M)
    assert mp == m, (pstr(mp), pstr(m))
    assert peval(m, M) == zeros(n)
    print('minimal polynomial      :', pstr(mp))

    # 2. Bezout identity and the polynomials h_j
    total = []
    for i, (f, h) in enumerate(zip(fs, hs)):
        print('f_%d = %s' % (i + 1, pstr(f)))
        print('h_%d = %s' % (i + 1, pstr(h)))
        total = padd(total, h)
    assert total == [Q(1)], total

    # 3. projections E_j = h_j(M)
    E = [peval(h, M) for h in hs]
    for j in range(3):
        assert mmul(E[j], E[j]) == E[j]
        assert mmul(E[j], M) == mmul(M, E[j])
        for i in range(3):
            if i != j:
                assert mmul(E[i], E[j]) == zeros(n)
    s = zeros(n)
    for e in E:
        s = madd(s, e)
    assert s == ident(n)
    block = [[0, 1], [2], [3, 4]]
    for j in range(3):
        want = zeros(n)
        for a in block[j]:
            want[a][a] = Q(1)
        assert E[j] == want
    print('E_1, E_2, E_3 are the coordinate projections onto')
    print('  span(e1,e2), span(e3), span(e4,e5): checked')

    # 4. W_j = ker p_j(M)^r_j equals im E_j, dimensions 2, 1, 2
    W = [nullspace(peval(q, M)) for q in powers]
    dims = [len(w) for w in W]
    print('dim W_1, W_2, W_3       :', dims)
    assert dims == [2, 1, 2]
    for j in range(3):
        img = columns(E[j])
        assert rank(img) == dims[j] == sum_dim(img, W[j])
    assert sum_dim(*W) == n == sum(dims)

    # 5. minimal polynomial of each restriction T_j is p_j^r_j
    for j in range(3):
        assert minimal_polynomial(restrict(M, W[j])) == powers[j]
    print('min poly of T_1, T_2, T_3:', [pstr(q) for q in powers])

    # 6. semisimple part S and nilpotent part N
    S = madd(mscale(E[0], Q(2)), mscale(E[1], Q(-1)))
    S = madd(S, mmul(M, E[2]))
    N = mmul(madd(M, mscale(ident(n), Q(-2))), E[0])
    assert madd(S, N) == M
    assert mmul(S, N) == mmul(N, S)
    assert mmul(N, N) == zeros(n) and N != zeros(n)
    TE3 = mmul(M, E[2])
    assert madd(mmul(TE3, TE3), E[2]) == zeros(n)
    print('N^2 = 0, N != 0, SN = NS, S + N = M: checked')
    print('N =', [[str(v) for v in r] for r in N])

    # 7. the formula does not depend on the basis
    P = [[1, 1, 0, 0, 0], [0, 1, 1, 0, 0], [0, 0, 1, 1, 0],
         [0, 0, 0, 1, 1], [0, 0, 0, 0, 1]]
    P = [[Q(v) for v in r] for r in P]
    Pinv = inverse(P)
    M2 = mmul(mmul(P, M), Pinv)
    E2 = [peval(h, M2) for h in hs]
    for j in range(3):
        assert E2[j] == mmul(mmul(P, E[j]), Pinv)
    assert E2[0] != E[0]
    print('E_j(P M P^-1) = P E_j P^-1 and differs from E_j: checked')

    # 8. boundary: the exponent of (x-2) is too small
    k1 = nullspace(peval(p1, M))
    print('dim ker(M-2I)           :', len(k1))
    assert len(k1) == 1
    assert sum_dim(k1, W[1], W[2]) == 4 < n

    # 9. boundary: factors that are not relatively prime
    q1, q2 = ppow(p1, 2), pmul(p1, p2)
    a, b = nullspace(peval(q1, M)), nullspace(peval(q2, M))
    print('dim ker q1, ker q2, sum :', len(a), len(b), sum_dim(a, b))
    assert sum_dim(a, b) < len(a) + len(b)

    # 10. x^2+1 has no root in Q (x^2+1 >= 1 for every rational x), so
    # T_3 has no eigenvalue over Q; its restriction has x^2+1 as minimal
    # polynomial, which is irreducible over Q.
    assert all(t * t + 1 != 0 for t in (Q(a, 7) for a in range(-70, 71)))
    assert all(len(nullspace(madd(restrict(M, W[2]),
                                  mscale(ident(2), -Q(a, 3))))) == 0
               for a in range(-30, 31))
    # 11. numbers quoted in the chapter and the exercises
    assert m == trim([4, 0, 1, 1, -3, 1])
    D = madd(mscale(E[0], Q(2)), mscale(E[1], Q(-1)))
    assert S == madd(D, mmul(M, E[2]))

    # exercise 1: the block A
    A3 = [[Q(v) for v in r] for r in A]
    h1 = pscale(pmul(trim([5, -1]), trim([1, 1])), Q(1, 9))
    h2 = pscale(ppow(p1, 2), Q(1, 9))
    assert padd(h1, h2) == [1]
    assert peval(h1, A3) == [[1, 0, 0], [0, 1, 0], [0, 0, 0]]
    assert peval(h2, A3) == [[0, 0, 0], [0, 0, 0], [0, 0, 1]]

    # exercise 2: complex projections of the quarter turn
    half = [[0.5, 0.5j], [-0.5j, 0.5]]
    Cc = [[0, -1], [1, 0]]
    for i in range(2):
        for j in range(2):
            v = (Cc[i][j] + (1j if i == j else 0)) / 2j
            assert v == half[i][j]
    assert mmul(half, half) == half

    # exercise 6: the cyclic shift over Q, F_5, F_2
    Sh = [[Q(int((i - j) % 4 == 1)) for j in range(4)] for i in range(4)]
    assert minimal_polynomial(Sh) == trim([-1, 0, 0, 0, 1])
    for c, v in ((1, (1, 1, 1, 1)), (4, (1, 4, 1, 4)),
                 (2, (1, 3, 4, 2)), (3, (1, 2, 4, 3))):
        image = [v[(i - 1) % 4] for i in range(4)]
        assert all((image[i] - c * v[i]) % 5 == 0 for i in range(4))
    N1 = madd(Sh, mscale(ident(4), Q(-1)))
    cube = mmul(mmul(N1, N1), N1)
    assert any(int(x) % 2 for row in cube for x in row)
    assert all(int(x) % 2 == 0 for row in mmul(cube, N1) for x in row)
    print('exercises 1, 2, 6: checked')
    print('all assertions passed')


if __name__ == '__main__':
    main()
