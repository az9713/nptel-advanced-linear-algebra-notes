"""Chapter 16 laboratory: triangulability and diagonalizability tests.

Exact arithmetic over Q (fractions.Fraction); one small computation over
F_2 written with plain integers.  A matrix is a list of rows; a vector is
a list; a polynomial is a list of coefficients, constant term first,
with no trailing zero (the zero polynomial is the empty list).

The laboratory checks, on matrices chosen by the chapter:
  * the minimal polynomial from Krylov dependence,
  * the criterion "triangulable iff the minimal polynomial splits",
    including an explicit triangulating basis built from conductors,
  * the criterion "diagonalizable iff the minimal polynomial is a
    product of distinct linear factors", against eigenspace dimensions,
  * the boundary cases: no rational root (rotation), characteristic 2,
    and the integral operator on polynomials (infinite dimension).
It does not prove the theorems; it tests them on finitely many cases.
"""
from fractions import Fraction
from math import gcd

import numpy as np


# ---------------------------------------------------------------- matrices
def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fraction(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, v):
    return [sum((A[i][k] * v[k] for k in range(len(v))), Fraction(0))
            for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


def mat_sub_scalar(A, c):
    n = len(A)
    return [[A[i][j] - (c if i == j else 0) for j in range(n)]
            for i in range(n)]


def rref(M):
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0])
    pivots, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        piv = M[r][c]
        M[r] = [x / piv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return M, pivots


def rank(M):
    return len(rref(M)[1])


def inverse(A):
    n = len(A)
    aug = [A[i][:] + identity(n)[i] for i in range(n)]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n)), 'matrix is singular'
    return [row[n:] for row in R]


def solve_in_span(cols, b):
    """Coefficients x with sum x_j cols[j] = b, or None if b is outside."""
    if not cols:
        return [] if all(t == 0 for t in b) else None
    M = [[cols[j][i] for j in range(len(cols))] + [b[i]]
         for i in range(len(b))]
    R, piv = rref(M)
    if len(cols) in piv:
        return None
    x = [Fraction(0)] * len(cols)
    for row, c in zip(R, piv):
        x[c] = row[len(cols)]
    return x


def flat(M):
    return [x for row in M for x in row]


# -------------------------------------------------------------- polynomials
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
                 for i in range(n)])


def pmul(p, q):
    if not p or not q:
        return []
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def pdivmod(p, d):
    p, d = trim(p), trim(d)
    q = [Fraction(0)] * max(len(p) - len(d) + 1, 0)
    r = p[:]
    while len(r) >= len(d) and r:
        k = len(r) - len(d)
        f = r[-1] / d[-1]
        q[k] = f
        r = trim([r[i] - (f * d[i - k] if 0 <= i - k < len(d) else 0)
                  for i in range(len(r))])
    return trim(q), r


def monic(p):
    p = trim(p)
    return [x / p[-1] for x in p]


def deriv(p):
    return trim([i * p[i] for i in range(1, len(p))])


def pgcd(p, q):
    p, q = trim(p), trim(q)
    while q:
        p, q = q, pdivmod(p, q)[1]
    return monic(p)


def peval(p, x):
    out = Fraction(0)
    for a in reversed(p):
        out = out * x + a
    return out


def poly_of_matrix(p, A):
    n = len(A)
    out = [[Fraction(0)] * n for _ in range(n)]
    for a in reversed(p):
        out = matmul(out, A)
        for i in range(n):
            out[i][i] += a
    return out


def poly_str(p):
    terms = []
    for k in range(len(p) - 1, -1, -1):
        if p[k] != 0:
            terms.append('%s*x^%d' % (p[k], k))
    return ' + '.join(terms) if terms else '0'


def divisors(n):
    n = abs(n)
    return [d for d in range(1, n + 1) if n % d == 0]


def rational_roots(p):
    """Distinct rational roots of p (rational root theorem)."""
    p = trim(p)
    roots = []
    if p and p[0] == 0:
        roots.append(Fraction(0))
    k = 0
    while k < len(p) and p[k] == 0:
        k += 1
    p = p[k:]
    if len(p) <= 1:
        return roots
    den = 1
    for a in p:
        den = den * a.denominator // gcd(den, a.denominator)
    ints = [int(a * den) for a in p]
    for a in divisors(ints[0]):
        for b in divisors(ints[-1]):
            for s in (1, -1):
                r = Fraction(s * a, b)
                if r not in roots and peval(p, r) == 0:
                    roots.append(r)
    return roots


def linear_factorization(p):
    """Return (list of (root, multiplicity), cofactor) over Q."""
    p = monic(p)
    found = []
    for r in rational_roots(p):
        m = 0
        while True:
            qt, rem = pdivmod(p, [-r, Fraction(1)])
            if rem:
                break
            p, m = qt, m + 1
        found.append((r, m))
    return found, p


def splits_over_Q(p):
    return len(linear_factorization(p)[1]) <= 1


# -------------------------------------------------- minimal polynomial etc.
def minimal_polynomial(A):
    n = len(A)
    powers = [identity(n)]
    for k in range(1, n * n + 1):
        powers.append(matmul(powers[-1], A))
        coef = solve_in_span([flat(P) for P in powers[:-1]],
                             flat(powers[-1]))
        if coef is not None:
            return [-c for c in coef] + [Fraction(1)]
    raise AssertionError('unreachable: Cayley-Hamilton bounds the degree')


def eigenspace_dim(A, c):
    n = len(A)
    return n - rank(mat_sub_scalar(A, c))


def triangulable(A):
    return splits_over_Q(minimal_polynomial(A))


def diagonalizable_by_polynomial(A):
    """Splits over Q with all multiplicities equal to 1."""
    found, rest = linear_factorization(minimal_polynomial(A))
    return len(rest) <= 1 and all(m == 1 for _, m in found)


def diagonalizable_by_eigenspaces(A):
    """Sum of eigenspace dimensions over the rational eigenvalues is n."""
    found, _ = linear_factorization(minimal_polynomial(A))
    return sum(eigenspace_dim(A, r) for r, _ in found) == len(A)


def conductor(A, beta, W):
    """Monic generator g of {g : g(A) beta in span(W)} and its degree.

    W is a list of independent vectors that span an A-invariant space.
    The smallest k with A^k beta in span(W, beta, ..., A^(k-1) beta) gives
    the generator x^k - sum c_j x^j.
    """
    vs = [beta]
    for k in range(1, len(beta) + 2):
        nxt = matvec(A, vs[-1])
        coef = solve_in_span(W + vs, nxt)
        if coef is not None:
            c = coef[len(W):]
            return [-x for x in c] + [Fraction(1)]
        vs.append(nxt)
    raise AssertionError('unreachable in finite dimension')


def triangulating_basis(A):
    """Basis (alpha_1..alpha_n) with A alpha_j in span(alpha_1..alpha_j).

    Follows the proof: pick beta outside W, take its conductor g into W,
    write g = (x - c) h, and adjoin alpha = h(A) beta.
    Returns (basis, list of the eigenvalues c used).
    """
    n = len(A)
    basis, used = [], []
    while len(basis) < n:
        beta = next(e for e in ([Fraction(int(i == j)) for j in range(n)]
                                for i in range(n))
                    if solve_in_span(basis, e) is None)
        g = conductor(A, beta, basis)
        found, rest = linear_factorization(g)
        if len(rest) > 1:
            raise ValueError('minimal polynomial does not split over Q')
        c = found[0][0]
        h, rem = pdivmod(g, [-c, Fraction(1)])
        assert not rem
        alpha = matvec(poly_of_matrix(h, A), beta)
        assert solve_in_span(basis, alpha) is None
        assert solve_in_span(basis + [alpha], matvec(
            mat_sub_scalar(A, c), alpha)) is not None
        basis.append(alpha)
        used.append(c)
    return basis, used


def is_upper_triangular(M):
    return all(M[i][j] == 0 for i in range(len(M)) for j in range(i))


def change_basis(A, basis):
    P = transpose(basis)
    return matmul(inverse(P), matmul(A, P))


def vec(*xs):
    return [Fraction(x) for x in xs]


# -------------------------------------------------------------------- main
def main():
    # ---- running matrix: triangulable, not diagonalizable -------------
    A = Q([[0, 1, 0], [2, -2, 2], [2, -3, 2]])
    mA = minimal_polynomial(A)
    print('A =', [[str(x) for x in r] for r in A])
    print('minimal polynomial of A:', poly_str(mA))
    assert mA == [0, 0, 0, 1]
    assert matmul(A, A) != Q([[0] * 3] * 3)
    assert poly_of_matrix(mA, A) == Q([[0] * 3] * 3)
    assert triangulable(A) and not diagonalizable_by_polynomial(A)
    assert eigenspace_dim(A, Fraction(0)) == 1
    print('dim ker A =', eigenspace_dim(A, Fraction(0)), '(n = 3)')

    # hand-made basis (a1, a2, a3), as in the chapter
    a1, a2, a3 = vec(1, 0, -1), vec(1, 1, 0), vec(0, 2, Fraction(5, 2))
    assert matvec(A, a1) == vec(0, 0, 0)
    assert matvec(A, a2) == a1
    assert matvec(A, a3) == [x + y for x, y in zip(a1, a2)]
    T1 = change_basis(A, [a1, a2, a3])
    print('triangular form in (a1, a2, a3):',
          [[str(x) for x in r] for r in T1])
    assert T1 == Q([[0, 1, 1], [0, 0, 1], [0, 0, 0]])

    # basis from the conductor construction
    basis, used = triangulating_basis(A)
    T2 = change_basis(A, basis)
    print('conductor basis:', [[str(x) for x in b] for b in basis])
    print('triangular form in conductor basis:',
          [[str(x) for x in r] for r in T2])
    assert is_upper_triangular(T2)
    assert [T2[i][i] for i in range(3)] == [0, 0, 0]
    assert T2 == Q([[0, 1, 0], [0, 0, 1], [0, 0, 0]])

    # the lemma's alpha for W = span(a1), beta = e1
    g = conductor(A, vec(1, 0, 0), [a1])
    print('conductor of e1 into span(a1):', poly_str(g))

    # ---- idempotent: diagonalizable ----------------------------------
    P = transpose([a1, a2, a3])
    D = Q([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
    E = matmul(P, matmul(D, inverse(P)))
    print('E =', [[str(x) for x in r] for r in E])
    assert matmul(E, E) == E
    mE = minimal_polynomial(E)
    print('minimal polynomial of E:', poly_str(mE))
    assert mE == [0, -1, 1]
    assert diagonalizable_by_polynomial(E)
    assert eigenspace_dim(E, Fraction(1)) + eigenspace_dim(
        E, Fraction(0)) == 3
    assert sum(E[i][i] for i in range(3)) == rank(E) == 2

    # ---- agreement of the two tests on a family of matrices -----------
    F = Q([[2, 0, 0], [0, 2, 0], [0, 0, 3]])
    J = Q([[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 5, 0], [0, 0, 0, 5]])
    S = Q([[0, 1], [1, 0]])
    Rot = Q([[0, -1], [1, 0]])
    cases = {'A': A, 'E': E, 'F': F, 'J': J, 'S': S, 'Rot': Rot,
             'zero': Q([[0, 0], [0, 0]]), 'scalar': Q([[7, 0], [0, 7]]),
             'Jordan2': Q([[1, 1], [0, 1]])}
    for name, M in cases.items():
        d1 = diagonalizable_by_polynomial(M)
        d2 = diagonalizable_by_eigenspaces(M)
        t = triangulable(M)
        print('%-8s m=%-28s triangulable=%-5s diagonalizable=%s'
              % (name, poly_str(minimal_polynomial(M)), t, d1))
        assert d1 == d2, name
        if d1:
            assert t
        if t:
            b, _ = triangulating_basis(M)
            assert is_upper_triangular(change_basis(M, b)), name
        else:
            try:
                triangulating_basis(M)
                raise AssertionError(name)
            except ValueError:
                pass
    assert not diagonalizable_by_polynomial(J) and triangulable(J)
    assert diagonalizable_by_polynomial(F)
    assert minimal_polynomial(F) == pmul([Fraction(-2), Fraction(1)],
                                         [Fraction(-3), Fraction(1)])
    assert minimal_polynomial(S) == pmul([Fraction(-1), Fraction(1)],
                                         [Fraction(1), Fraction(1)])
    # triangulating basis for J (two eigenvalues, one repeated)
    bJ, usedJ = triangulating_basis(J)
    TJ = change_basis(J, bJ)
    assert is_upper_triangular(TJ)
    print('J eigenvalues used in order:', [str(c) for c in usedJ])
    assert sorted(usedJ) == [2, 2, 5, 5]
    assert sorted(TJ[i][i] for i in range(4)) == [2, 2, 5, 5]

    # ---- boundary 1: no rational root ----------------------------------
    assert not triangulable(Rot)
    try:
        triangulating_basis(Rot)
        raise AssertionError('rotation must not be triangulable over Q')
    except ValueError:
        print('Rot: no triangulating basis over Q (m = x^2 + 1)')
    ev = np.linalg.eigvals(np.array([[0.0, -1.0], [1.0, 0.0]]))
    assert np.allclose(sorted(ev.imag), [-1.0, 1.0])
    assert np.allclose(ev.real, 0.0)
    print('Rot over C: eigenvalues', ev)

    # ---- boundary 2: characteristic 2 -----------------------------------
    def mul2(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(2)) % 2
                 for j in range(2)] for i in range(2)]

    M2 = [[1, 1], [0, 1]]
    I2 = [[1, 0], [0, 1]]
    N2 = [[(M2[i][j] - I2[i][j]) % 2 for j in range(2)] for i in range(2)]
    assert mul2(M2, M2) == I2 and M2 != I2
    assert N2 != [[0, 0], [0, 0]] and mul2(N2, N2) == [[0, 0], [0, 0]]
    print('over F_2: M^2 = I, M != I, (M - I)^2 = 0, so m = (x - 1)^2')
    # over Q the same relation M^2 = I forces distinct roots 1, -1:
    assert diagonalizable_by_polynomial(S)

    # ---- boundary 3: integral operator on polynomials -------------------
    def integ(f):
        return trim([Fraction(0)] + [f[k] / (k + 1) for k in range(len(f))])

    f = [Fraction(-1, 2), Fraction(1)]          # x - 1/2
    Tf = integ(f)
    print('T(x - 1/2) =', poly_str(Tf), ' value at 1/2:',
          peval(Tf, Fraction(1, 2)))
    assert peval(f, Fraction(1, 2)) == 0
    assert Tf == [0, Fraction(-1, 2), Fraction(1, 2)]
    assert peval(Tf, Fraction(1, 2)) == Fraction(-1, 8)
    for deg in range(0, 6):                    # T raises degree by one
        g = [Fraction(0)] * deg + [Fraction(1)]
        assert len(integ(g)) == len(g) + 1
    for c in (Fraction(1, 3), Fraction(2), Fraction(-5, 7)):
        h = [-c, Fraction(1)]
        assert peval(integ(h), c) == -c * c / 2   # never 0 for c != 0
    # polynomials vanishing at 0 are invariant: T f(0) = 0 always
    assert all(peval(integ(g), 0) == 0 for g in ([1], [1, 2], [0, 0, 3]))
    # ---- commuting families (Section 8) --------------------------------
    assert matmul(A, E) != matmul(E, A)
    for M in (A, E):
        assert is_upper_triangular(change_basis(M, [a1, a2, a3]))
    N1, N2t = Q([[0, 1], [0, 0]]), Q([[0, 0], [1, 0]])
    for N, e, other in ((N1, vec(1, 0), N2t), (N2t, vec(0, 1), N1)):
        assert eigenspace_dim(N, Fraction(0)) == 1
        assert matvec(N, e) == vec(0, 0)
        assert solve_in_span([e], matvec(other, e)) is None
    print('A E != E A; the two nilpotents share no eigenvector')

    # ---- exercise 5: the derivative on polynomials of degree <= n ------
    def nil_index(n, p=0):
        # smallest k with D^k = 0, in characteristic p (0 means over Q)
        k = 0
        coeffs = [1] * (n + 1)          # coefficient left on D^k x^j
        while any(coeffs):
            k += 1
            coeffs = [c * (j - k + 1) if j >= k else 0
                      for j, c in enumerate(coeffs)]
            if p:
                coeffs = [c % p for c in coeffs]
        return k

    for n in range(7):
        assert nil_index(n) == n + 1
        for p in (2, 3, 5):
            assert nil_index(n, p) == min(n + 1, p)
    print('derivative: m_D = x^(n+1) over Q, x^min(n+1,p) mod p')
    print('lab checks passed')


if __name__ == '__main__':
    main()
