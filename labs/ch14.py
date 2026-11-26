"""Chapter 14 laboratory: Cayley-Hamilton and the minimal polynomial.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows.  A polynomial is a list of coefficients, lowest degree first:
[a0, a1, a2] means a0 + a1*x + a2*x^2.  Nothing here proves a theorem;
each routine checks the chapter's statements on concrete matrices.
numpy is used only for the floating-point eigenvalues of section 6.
"""
import random
from fractions import Fraction as Fr

import numpy as np


def Q(rows):
    return [[Fr(x) for x in row] for row in rows]


def eye(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def mul(A, B):
    return [[sum(a*b for a, b in zip(r, c)) for c in zip(*B)] for r in A]


def add(A, B):
    return [[a + b for a, b in zip(r, s)] for r, s in zip(A, B)]


def scal(c, A):
    return [[c*a for a in r] for r in A]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def zeros(n):
    return [[Fr(0)]*n for _ in range(n)]


def mpow(A, k):
    R = eye(len(A))
    for _ in range(k):
        R = mul(R, A)
    return R


def poly_at_matrix(p, A):
    """p(A) by Horner's rule."""
    R = zeros(len(A))
    for c in reversed(p):
        R = add(mul(R, A), scal(Fr(c), eye(len(A))))
    return R


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(p, q):
    r = [Fr(0)]*(len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a*b
    return r


def pdivmod(p, q):
    p, q = trim(p), trim(q)
    quo = [Fr(0)]*max(len(p) - len(q) + 1, 1)
    while len(p) >= len(q) and p:
        c, d = p[-1]/q[-1], len(p) - len(q)
        quo[d] = c
        for i, b in enumerate(q):
            p[i + d] -= c*b
        p = trim(p)
    return trim(quo), p


def pgcd(p, q):
    p, q = trim(p), trim(q)
    while q:
        p, q = q, pdivmod(p, q)[1]
    return [c/p[-1] for c in p]


def pderiv(p):
    return [i*c for i, c in enumerate(p)][1:]


def det(A):
    """Determinant by elimination with row swaps."""
    A, n, d = [r[:] for r in A], len(A), Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c]/A[c][c]
            A[r] = [a - f*b for a, b in zip(A[r], A[c])]
    return d


def rank(A):
    A, rk = [r[:] for r in A], 0
    for c in range(len(A[0])):
        piv = next((r for r in range(rk, len(A)) if A[r][c] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for r in range(len(A)):
            if r != rk and A[r][c] != 0:
                f = A[r][c]/A[rk][c]
                A[r] = [a - f*b for a, b in zip(A[r], A[rk])]
        rk += 1
    return rk


def inverse(A):
    n = len(A)
    M = [r[:] + e for r, e in zip(A, eye(n))]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        M[c] = [a/M[c][c] for a in M[c]]
        for r in range(n):
            if r != c:
                M[r] = [a - M[r][c]*b for a, b in zip(M[r], M[c])]
    return [r[n:] for r in M]


def minor_matrix(A, i, j):
    return [r[:j] + r[j + 1:] for k, r in enumerate(A) if k != i]


def adjugate(A):
    """Transpose of the cofactor matrix."""
    n = len(A)
    return [[(-1)**(i + j)*det(minor_matrix(A, j, i))
             for j in range(n)] for i in range(n)]


def char_poly(A):
    """Faddeev-LeVerrier: needs division by 1..n, so characteristic 0."""
    n = len(A)
    c, M = [Fr(0)]*n + [Fr(1)], eye(n)
    for k in range(1, n + 1):
        AM = mul(A, M)
        c[n - k] = -trace(AM)/k
        M = add(AM, scal(c[n - k], eye(n)))
    return c


def det_char_poly_at(A, t):
    n = len(A)
    return det(add(scal(Fr(t), eye(n)), scal(Fr(-1), A)))


def poly_at(p, t):
    return sum(c*Fr(t)**i for i, c in enumerate(p))


def min_poly(A):
    """Least-degree monic p with p(A) = 0, by the first dependent power.

    Incremental elimination on the flattened powers I, A, A^2, ...;
    each stored row remembers its combination of powers.  Cayley-Hamilton
    guarantees the loop stops at degree <= n.
    """
    n, basis, P, k = len(A), [], eye(len(A)), 0
    while True:
        vec = [x for r in P for x in r]
        coef = [Fr(0)]*k + [Fr(1)]
        for bvec, bcoef, piv in basis:
            if vec[piv] != 0:
                f = vec[piv]/bvec[piv]
                vec = [a - f*b for a, b in zip(vec, bvec)]
                for i, b in enumerate(bcoef):
                    coef[i] -= f*b
        if all(v == 0 for v in vec):
            return coef
        piv = next(i for i, v in enumerate(vec) if v != 0)
        basis.append((vec, coef, piv))
        P, k = mul(P, A), k + 1


def nullity(A):
    return len(A) - rank(A)


def shift(A, lam):
    return add(A, scal(Fr(-lam), eye(len(A))))


def multiplicity(p, lam):
    """Algebraic multiplicity of lam as a root of p."""
    m, p = 0, trim(p)
    while True:
        q, r = pdivmod(p, [Fr(-lam), Fr(1)])
        if r:
            return m
        m, p = m + 1, q


def section1_two_by_two():
    A = Q([[1, 2], [3, 4]])
    f = char_poly(A)
    assert f == [-2, -5, 1]
    assert poly_at_matrix(f, A) == zeros(2)
    B1, B0 = eye(2), [[Fr(-4), Fr(2)], [Fr(3), Fr(-1)]]
    adj = adjugate(add(scal(Fr(7), eye(2)), scal(Fr(-1), A)))
    assert adj == add(scal(Fr(7), B1), B0)
    assert add(B0, scal(Fr(-1), mul(A, B1))) == scal(f[1], eye(2))
    assert scal(Fr(-1), mul(A, B0)) == scal(f[0], eye(2))
    Ainv = scal(Fr(1, 2), add(A, scal(Fr(-5), eye(2))))
    assert mul(A, Ainv) == eye(2)
    print('2x2: c =', [int(c) for c in f], ' A^-1 =', Ainv)


def adjugate_coefficients(A):
    """B_k with adj(xI - A) = sum x^k B_k, by the proof's recurrence."""
    n, f = len(A), char_poly(A)
    B = {n - 1: eye(n)}
    for k in range(n - 1, 0, -1):
        B[k - 1] = add(scal(f[k], eye(n)), mul(A, B[k]))
    assert scal(Fr(-1), mul(A, B[0])) == scal(f[0], eye(n))
    return B


def section2_cayley_hamilton_checks():
    A = Q([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]])
    f = char_poly(A)
    assert f == [0, 0, -4, 0, 1]
    for t in range(-4, 5):
        assert det_char_poly_at(A, t) == poly_at(f, t)
    B = adjugate_coefficients(A)
    for t in range(-3, 4):
        adj = adjugate(shift(scal(Fr(-1), A), -t))
        tot = zeros(4)
        for k, Bk in B.items():
            tot = add(tot, scal(Fr(t)**k, Bk))
        assert tot == adj
    assert poly_at_matrix(f, A) == zeros(4)
    A2, A3, A4 = mpow(A, 2), mpow(A, 3), mpow(A, 4)
    assert A3 == scal(Fr(4), A) and A4 == scal(Fr(4), A2)
    m = min_poly(A)
    assert m == [0, -4, 0, 1]
    assert pdivmod(f, m)[1] == [] and pdivmod(f, m)[0] == [0, 1]
    print('4x4: c =', [int(c) for c in f], ' m =', [int(c) for c in m])
    print('A^2 =', [[int(x) for x in r] for r in A2])
    return A, f, m


def section3_multiplicities(A, f):
    out = {}
    for lam in (0, 2, -2):
        d = multiplicity(f, lam)
        g = nullity(shift(A, lam))
        assert 1 <= g <= d
        out[lam] = (d, g)
    assert out == {0: (2, 2), 2: (1, 1), -2: (1, 1)}
    assert rank(A) == 2 and rank(shift(A, 2)) == 3
    P = Q([[-1, 0, 1, 1], [0, -1, 1, -1], [1, 0, 1, 1], [0, 1, 1, -1]])
    D = mul(mul(inverse(P), A), P)
    assert D == Q([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, -2]])
    print('(alg, geom) at 0, 2, -2:', out, ' det P =', det(P))
    # powers by remainder: x^10 mod m, then evaluate
    m = [Fr(0), Fr(-4), Fr(0), Fr(1)]
    xk = [Fr(0)]*10 + [Fr(1)]
    r = pdivmod(xk, m)[1]
    assert r == [0, 0, 256]
    assert poly_at_matrix(r, A) == mpow(A, 10)
    assert mpow(A, 10) == scal(Fr(256), mpow(A, 2))
    print('x^10 mod m =', [int(c) for c in r])


def section4_three_by_three():
    A = Q([[6, -3, -2], [4, -1, -2], [10, -5, -3]])
    f = char_poly(A)
    assert f == [-2, 1, -2, 1]
    assert pmul([Fr(-2), Fr(1)], [Fr(1), Fr(0), Fr(1)]) == f
    assert det(A) == 2
    assert min_poly(A) == f
    Ainv = scal(Fr(1, 2), add(add(mpow(A, 2), scal(Fr(-2), A)), eye(3)))
    assert mul(A, Ainv) == eye(3)
    assert Ainv == inverse(A)
    assert rank(shift(A, 2)) == 2
    v = [[Fr(1)], [Fr(0)], [Fr(2)]]
    assert mul(A, v) == scal(Fr(2), v)
    assert pgcd(f, pderiv(f)) == [1]
    ev = np.linalg.eigvals(np.array([[6, -3, -2], [4, -1, -2],
                                     [10, -5, -3]], dtype=float))
    ev = sorted(ev, key=lambda z: (z.imag, z.real))
    assert np.allclose(ev, [-1j, 2, 1j], atol=1e-9)
    print('3x3: c =', [int(c) for c in f], ' det =', det(A))
    print('3x3 inverse =', [[str(x) for x in r] for r in Ainv])
    print('numpy eigenvalues:', np.round(ev, 6))


def rand_matrix(rng, n, vals):
    return Q([[rng.choice(vals) for _ in range(n)] for _ in range(n)])


def section5_ab_ba():
    A, B = Q([[0, 0], [1, 0]]), Q([[1, 0], [0, 0]])
    AB, BA = mul(A, B), mul(B, A)
    assert AB == A and BA == zeros(2)
    assert min_poly(AB) == [0, 0, 1] and min_poly(BA) == [0, 1]
    assert char_poly(AB) == char_poly(BA) == [0, 0, 1]
    rng, differ, total = random.Random(14), 0, 0
    for n in (2, 3, 4):
        for _ in range(150):
            A = rand_matrix(rng, n, (-1, 0, 0, 1, 2))
            B = rand_matrix(rng, n, (-1, 0, 0, 1, 2))
            AB, BA = mul(A, B), mul(B, A)
            assert char_poly(AB) == char_poly(BA)
            a, b = min_poly(AB), min_poly(BA)
            xa, xb = [Fr(0)] + a, [Fr(0)] + b
            assert a == b or xa == b or xb == a
            assert (det(A) == 0 and det(B) == 0) or a == b
            total += 1
            differ += a != b
    print('AB vs BA trials:', total, ' minimal polynomials differ:', differ)
    assert 0 < differ < total


def section6_geometric_vs_algebraic():
    rng, strict, total = random.Random(28), 0, 0
    for _ in range(200):
        n = 4
        A = [[Fr(0)]*n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                A[i][j] = Fr(rng.choice((0, 0, 1, 2)) if j > i
                             else rng.choice((0, 1, 2)))
        f = char_poly(A)
        for lam in (0, 1, 2):
            d = multiplicity(f, lam)
            if d == 0:
                continue
            g = nullity(shift(A, lam))
            assert 1 <= g <= d
            total += 1
            strict += g < d
        assert pdivmod(f, min_poly(A))[1] == []
    print('geometric <= algebraic checks:', total, ' strict:', strict)
    assert strict > 0


if __name__ == '__main__':
    section1_two_by_two()
    A4, f4, m4 = section2_cayley_hamilton_checks()
    section3_multiplicities(A4, f4)
    section4_three_by_three()
    section5_ab_ba()
    section6_geometric_vs_algebraic()
    print('all chapter 14 checks passed')
