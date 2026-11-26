"""Chapter 13 laboratory: diagonalizability and minimal polynomials.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows.  A polynomial is a list of coefficients, lowest degree first:
[a0, a1, a2] stands for a0 + a1 x + a2 x^2.  Every routine returns new
objects and never changes its inputs.

The laboratory checks, on the running matrix A and the defective matrix
B of the chapter:
  * f(T)v = f(c)v for an eigenvector v (Lemma 2.1);
  * the dimension test for diagonalizability (Theorem 4.1);
  * the minimal polynomial from dependence of powers (Theorems 5.1, 6.1);
  * invariance of the minimal polynomial under similarity (Section 6);
  * roots of the minimal polynomial are the eigenvalues (Theorem 7.1);
  * the Lagrange operators and the criterion (Theorem 8.1).
Apart from a floating-point check of the rotation over C it computes
only over Q.
"""
from fractions import Fraction
import numpy as np


def Q(rows):
    """Matrix over Q from integer or Fraction entries."""
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def mat_mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def mat_add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def mat_scale(c, a):
    return [[c*x for x in row] for row in a]


def mat_sub(a, b):
    return mat_add(a, mat_scale(Fraction(-1), b))


def mat_inv(a):
    """Inverse by Gauss-Jordan elimination; raises on a singular matrix."""
    n = len(a)
    m = [list(row) + list(idr) for row, idr in zip(a, identity(n))]
    for col in range(n):
        piv = next((r for r in range(col, n) if m[r][col] != 0), None)
        if piv is None:
            raise ZeroDivisionError('singular matrix')
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        m[col] = [x/pv for x in m[col]]
        for r in range(n):
            if r != col and m[r][col] != 0:
                f = m[r][col]
                m[r] = [x - f*y for x, y in zip(m[r], m[col])]
    return [row[n:] for row in m]


def rref_rank_nullspace(a):
    """Return (rank, basis of the null space) of a, exactly."""
    rows, cols = len(a), len(a[0])
    m = [list(r) for r in a]
    pivots = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        m[r] = [x/pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f*y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fcol in free:
        v = [Fraction(0)]*cols
        v[fcol] = Fraction(1)
        for row, pc in enumerate(pivots):
            v[pc] = -m[row][fcol]
        basis.append(v)
    return len(pivots), basis


# ---- polynomials ------------------------------------------------------

def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def poly_mul(p, q):
    out = [Fraction(0)]*(len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a*b
    return trim(out)


def poly_divmod(p, d):
    """Quotient and remainder of p by d (d nonzero)."""
    p, d = trim(p), trim(d)
    q = [Fraction(0)]*max(len(p) - len(d) + 1, 1)
    r = list(p)
    while len(r) >= len(d) and r:
        k = len(r) - len(d)
        c = r[-1]/d[-1]
        q[k] = c
        for i, x in enumerate(d):
            r[i + k] -= c*x
        r = trim(r)
    return trim(q), r


def poly_eval(p, x):
    out = Fraction(0)
    for c in reversed(p):
        out = out*x + c
    return out


def poly_at_matrix(p, a):
    """p(A) by Horner's rule; the constant term multiplies I."""
    n = len(a)
    out = [[Fraction(0)]*n for _ in range(n)]
    for c in reversed(p):
        out = mat_add(mat_mul(out, a), mat_scale(c, identity(n)))
    return out


def from_roots(roots):
    p = [Fraction(1)]
    for r in roots:
        p = poly_mul(p, [Fraction(-r), Fraction(1)])
    return p


def show(p):
    return ' + '.join('%s*x^%d' % (c, i) for i, c in enumerate(p) if c != 0)


# ---- the objects of the chapter ---------------------------------------

def char_poly(a):
    """Faddeev-LeVerrier: coefficients of det(xI - A), exact over Q."""
    n = len(a)
    coeffs = [Fraction(1)]
    mk = [[Fraction(0)]*n for _ in range(n)]
    for k in range(1, n + 1):
        mk = mat_add(mat_mul(a, mk), mat_scale(coeffs[-1], identity(n)))
        am = mat_mul(a, mk)
        ck = -sum(am[i][i] for i in range(n))/k
        coeffs.append(ck)
    return list(reversed(coeffs))


def vec(a):
    return [x for row in a for x in row]


def min_poly(a):
    """Least-degree monic p with p(A) = 0, from linear dependence of
    I, A, A^2, ...; also returns the number of powers inspected."""
    n = len(a)
    powers = [identity(n)]
    for d in range(1, n*n + 1):
        powers.append(mat_mul(powers[-1], a))
        cols = [vec(p) for p in powers]
        mat = [[cols[j][i] for j in range(len(cols))]
               for i in range(n*n)]
        rank, null = rref_rank_nullspace(mat)
        if rank < len(cols):
            rel = null[0]
            return trim([x/rel[-1] for x in rel]), len(cols)
    raise AssertionError('no dependence among n^2 + 1 powers')


def eig_space_dim(a, c):
    n = len(a)
    rank, _ = rref_rank_nullspace(mat_sub(a, mat_scale(Fraction(c),
                                                       identity(n))))
    return n - rank


def integer_roots(p, bound=30):
    return [r for r in range(-bound, bound + 1) if poly_eval(p, r) == 0]


def multiplicity(p, r):
    k = 0
    while True:
        q, rem = poly_divmod(p, [Fraction(-r), Fraction(1)])
        if rem:
            return k
        p, k = q, k + 1


def lagrange(nodes, i):
    """f_i with f_i(c_j) = delta_ij; the nodes must be distinct."""
    p = [Fraction(1)]
    for j, c in enumerate(nodes):
        if j != i:
            p = poly_mul(p, [Fraction(-c)/(nodes[i] - c),
                             Fraction(1)/(nodes[i] - c)])
    return p


def is_diagonalizable_over_Q(a):
    """Dimension test of Theorem 4.1 with integer eigenvalues."""
    n = len(a)
    roots = integer_roots(char_poly(a))
    return sum(eig_space_dim(a, c) for c in roots) == n


def apply(a, v):
    return [sum(a[i][j]*v[j] for j in range(len(v)))
            for i in range(len(a))]


def main():
    A = Q([[5, -6, -6], [-1, 4, 2], [3, -6, -4]])
    B = Q([[2, 1, 0], [0, 2, 0], [0, 0, 1]])
    n = 3
    I3 = identity(n)
    zero = [[Fraction(0)]*n for _ in range(n)]

    # characteristic polynomials: equal for A and B
    cA, cB = char_poly(A), char_poly(B)
    print('c_A =', show(cA))
    assert cA == cB == from_roots([1, 2, 2])
    assert cA == [Fraction(-4), Fraction(8), Fraction(-5), Fraction(1)]

    # eigenspace dimensions
    dA = {c: eig_space_dim(A, c) for c in (1, 2)}
    dB = {c: eig_space_dim(B, c) for c in (1, 2)}
    print('dim E_c for A:', dA, ' for B:', dB)
    assert dA == {1: 1, 2: 2} and dB == {1: 1, 2: 1}
    assert is_diagonalizable_over_Q(A)
    assert not is_diagonalizable_over_Q(B)

    # Lemma 2.1: f(T)v = f(c)v on an eigenvector
    v1 = [Fraction(3), Fraction(-1), Fraction(3)]
    v2a = [Fraction(2), Fraction(1), Fraction(0)]
    v2b = [Fraction(2), Fraction(0), Fraction(1)]
    assert apply(A, v1) == v1
    assert apply(A, v2a) == [2*x for x in v2a]
    assert apply(A, v2b) == [2*x for x in v2b]
    f = [Fraction(c) for c in (7, -2, 0, 5, 1)]
    for c, v in ((1, v1), (2, v2a), (2, v2b)):
        assert apply(poly_at_matrix(f, A), v) == \
            [poly_eval(f, c)*x for x in v]
    print('f(1) =', poly_eval(f, 1), ' f(2) =', poly_eval(f, 2))

    # Theorem 4.1: the eigenvectors form a basis; P^-1 A P = D
    P = [[v1[i], v2a[i], v2b[i]] for i in range(n)]
    D = mat_mul(mat_inv(P), mat_mul(A, P))
    print('P^-1 A P =', [[str(x) for x in r] for r in D])
    assert D == Q([[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    assert rref_rank_nullspace(P)[0] == 3

    # Section 5: n^2 + 1 powers are dependent; minimal polynomial
    mA, countA = min_poly(A)
    mB, countB = min_poly(B)
    print('m_A =', show(mA), ' powers inspected:', countA)
    print('m_B =', show(mB), ' powers inspected:', countB)
    assert mA == from_roots([1, 2]) and mB == from_roots([1, 2, 2])
    assert countA <= n*n + 1 and countB <= n*n + 1
    assert poly_at_matrix(mA, A) == zero
    assert poly_at_matrix(from_roots([1, 2]), B) != zero
    # same characteristic polynomial, different minimal polynomial
    assert cA == cB and mA != mB
    # m divides c in both examples (an instance of Chapter 14)
    for m, c in ((mA, cA), (mB, cB)):
        assert poly_divmod(c, m)[1] == []
    # every annihilating polynomial is a multiple of m_A
    g = poly_mul(mA, [Fraction(3), Fraction(1), Fraction(4)])
    assert poly_at_matrix(g, A) == zero
    assert poly_divmod(g, mA)[1] == []
    assert poly_at_matrix([Fraction(1), Fraction(1)], A) != zero

    # Section 6: similar matrices share the minimal polynomial
    S = Q([[1, 2, 0], [0, 1, 3], [1, 0, 1]])
    A2 = mat_mul(mat_inv(S), mat_mul(A, S))
    assert min_poly(A2)[0] == mA
    assert A2 != A

    # Section 7: roots of m are the eigenvalues, multiplicities differ
    assert integer_roots(mA) == integer_roots(cA) == [1, 2]
    assert integer_roots(mB) == integer_roots(cB) == [1, 2]
    print('mult of 2 in c_A, m_A:', multiplicity(cA, 2),
          multiplicity(mA, 2))
    assert (multiplicity(cA, 2), multiplicity(mA, 2)) == (2, 1)
    assert (multiplicity(cB, 2), multiplicity(mB, 2)) == (2, 2)

    # Section 8: converse by Lagrange projections
    nodes = [Fraction(1), Fraction(2)]
    ps = [lagrange(nodes, i) for i in range(2)]
    for i in range(2):
        for j in range(2):
            assert poly_eval(ps[i], nodes[j]) == int(i == j)
    E = [poly_at_matrix(p, A) for p in ps]
    total = mat_add(E[0], E[1])
    assert total == I3
    for i, c in enumerate(nodes):
        assert mat_mul(mat_sub(A, mat_scale(c, I3)), E[i]) == zero
        assert mat_mul(E[i], E[i]) == E[i]
    ranks = [rref_rank_nullspace(e)[0] for e in E]
    print('ranks of p_i(A):', ranks, ' p_1 =', show(ps[0]))
    assert ranks == [1, 2]
    assert mat_mul(E[0], E[1]) == zero

    # powers of A from the minimal polynomial
    Ak = I3
    for k in range(1, 13):
        Ak = mat_mul(Ak, A)
        formula = mat_add(mat_scale(Fraction(2**k - 1), A),
                          mat_scale(Fraction(2 - 2**k), I3))
        assert Ak == formula
    print('A^12 row 1:', [str(x) for x in Ak[0]])

    # field dependence: the rotation has m = x^2 + 1
    R = Q([[0, -1], [1, 0]])
    mR, _ = min_poly(R)
    print('m_R =', show(mR))
    assert mR == [Fraction(1), Fraction(0), Fraction(1)]
    assert integer_roots(mR) == []
    assert not is_diagonalizable_over_Q(R)
    Rc = np.array([[0, -1], [1, 0]], dtype=complex)
    Pc = np.array([[1, 1], [1j, -1j]])
    Dc = np.linalg.inv(Pc) @ Rc @ Pc
    assert np.allclose(Dc, np.diag([-1j, 1j]), atol=1e-12)

    # infinite dimension: d/dx on Q[x] has no annihilating polynomial
    def deriv(coeffs):
        return [i*c for i, c in enumerate(coeffs)][1:] or [Fraction(0)]

    def p_of_D(p, coeffs):
        out = [Fraction(0)]*(len(coeffs))
        cur = list(coeffs)
        for a in p:
            for i, c in enumerate(cur):
                out[i] += a*c
            cur = deriv(cur) + [Fraction(0)]*(len(coeffs) - len(cur) + 1)
            cur = cur[:len(coeffs)]
        return out

    for p in ([Fraction(0), Fraction(2), Fraction(5)],
              [Fraction(3), Fraction(1)],
              [Fraction(0), Fraction(0), Fraction(-4), Fraction(1)]):
        j = next(i for i, c in enumerate(p) if c != 0)
        xj = [Fraction(0)]*j + [Fraction(1)]
        res = p_of_D(p, xj)
        assert res[0] == p[j]*_fact(j) != 0
    print('all checks passed')


def _fact(j):
    out = 1
    for i in range(2, j + 1):
        out *= i
    return out


if __name__ == '__main__':
    main()
