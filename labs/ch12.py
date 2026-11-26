"""Chapter 12 laboratory: eigenvalues, eigenvectors, characteristic polynomial.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of rows.
Every routine returns a new object and never changes its input.  A
polynomial is a list of coefficients, highest degree first: [1, -5, 8, -4]
stands for x^3 - 5x^2 + 8x - 4.
"""
from fractions import Fraction


def Q(rows):
    """Matrix over Q from integer or Fraction entries."""
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def mat_add(a, b, s=1):
    """Return a + s*b."""
    return [[x + s * y for x, y in zip(r, t)] for r, t in zip(a, b)]


def scalar_mul(s, a):
    return [[s * x for x in row] for row in a]


def mat_pow(a, k):
    out = identity(len(a))
    for _ in range(k):
        out = mat_mul(out, a)
    return out


def mat_vec(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def rref(a):
    """Reduced row echelon form and the list of pivot columns."""
    m = [row[:] for row in a]
    rows, cols = len(m), len(m[0])
    pivots, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        lead = m[r][c]
        m[r] = [x / lead for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return m, pivots


def rank(a):
    return len(rref(a)[1])


def null_basis(a):
    """Basis of {x : a x = 0}, one vector per free column."""
    m, pivots = rref(a)
    n = len(a[0])
    basis = []
    for f in (c for c in range(n) if c not in pivots):
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for row, pc in zip(m, pivots):
            v[pc] = -row[f]
        basis.append(v)
    return basis


def det(a):
    """Determinant by fraction-exact elimination."""
    m = [row[:] for row in a]
    n, d = len(m), Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if m[i][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            m[c], m[p] = m[p], m[c]
            d = -d
        d *= m[c][c]
        for i in range(c + 1, n):
            f = m[i][c] / m[c][c]
            m[i] = [x - f * y for x, y in zip(m[i], m[c])]
    return d


def inverse(a):
    n = len(a)
    aug = [list(r) + list(e) for r, e in zip(a, identity(n))]
    m, pivots = rref(aug)
    assert pivots[:n] == list(range(n)), 'matrix is singular'
    return [row[n:] for row in m]


def charpoly(a):
    """det(xI - A) by the Faddeev-LeVerrier recurrence, exact over Q."""
    n = len(a)
    coeffs = [Fraction(1)]
    m = identity(n)
    for k in range(1, n + 1):
        am = mat_mul(a, m)
        ck = -trace(am) / k
        coeffs.append(ck)
        m = mat_add(am, identity(n), ck)
    return coeffs


def poly_mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def poly_eval(p, x):
    out = Fraction(0)
    for c in p:
        out = out * x + c
    return out


def poly_div_linear(p, c):
    """Quotient and remainder of p by (x - c), synthetic division."""
    q = [p[0]]
    for coef in p[1:]:
        q.append(coef + c * q[-1])
    return q[:-1], q[-1]


def multiplicity(p, c):
    """Largest m with (x - c)^m dividing p."""
    m = 0
    while len(p) > 1:
        q, r = poly_div_linear(p, c)
        if r != 0:
            break
        p, m = q, m + 1
    return m


def eigenspace(a, c):
    """Basis of E_c = ker(A - cI)."""
    return null_basis(mat_add(a, scalar_mul(c, identity(len(a))), -1))


def is_eigenpair(a, c, v):
    return any(x != 0 for x in v) and mat_vec(a, v) == [c * x for x in v]


def columns_to_matrix(vectors):
    return [[v[i] for v in vectors] for i in range(len(vectors[0]))]


def show(m):
    return [[str(x) for x in row] for row in m]


def main():
    # ---- A and B, the two running matrices --------------------------
    A = Q([[3, 1, -1], [2, 2, -1], [2, 2, 0]])
    B = Q([[5, -6, -6], [-1, 4, 2], [3, -6, -4]])
    target = poly_mul(poly_mul([1, -1], [1, -2]), [1, -2])
    target = [Fraction(x) for x in target]
    print('target (x-1)(x-2)^2 =', [str(x) for x in target])
    for name, M in (('A', A), ('B', B)):
        cp = charpoly(M)
        print(name, 'charpoly coefficients:', [str(x) for x in cp])
        assert cp == target
        assert cp[1] == -trace(M) == -5
        assert cp[3] == -det(M) == -4
        assert multiplicity(cp, 1) == 1 and multiplicity(cp, 2) == 2
        # roots are exactly the eigenvalues (Theorem 2.3)
        for c in range(-3, 5):
            has_vec = len(eigenspace(M, Fraction(c))) > 0
            assert has_vec == (poly_eval(cp, Fraction(c)) == 0)
        assert det(mat_add(M, scalar_mul(2, identity(3)), -1)) == 0
    print('det A =', det(A), ' tr A =', trace(A), ' det B =', det(B),
          ' tr B =', trace(B))

    # ---- eigenvectors of A (the defective matrix) -------------------
    EA1 = eigenspace(A, Fraction(1))
    EA2 = eigenspace(A, Fraction(2))
    print('A: E_1 basis', show(EA1), ' E_2 basis', show(EA2))
    assert [len(EA1), len(EA2)] == [1, 1]
    assert is_eigenpair(A, Fraction(1), [1, 0, 2])
    assert is_eigenpair(A, Fraction(2), [1, 1, 2])
    assert rank(mat_add(A, identity(3), -1)) == 2
    assert rank(mat_add(A, scalar_mul(2, identity(3)), -1)) == 2
    assert rank(columns_to_matrix([[1, 0, 2], [1, 1, 2]])) == 2
    assert len(EA1) + len(EA2) < 3   # A is not diagonalizable
    print('A: sum of eigenspace dimensions =', len(EA1) + len(EA2), '< 3')

    # ---- eigenvectors of B and diagonalization ----------------------
    EB1 = eigenspace(B, Fraction(1))
    EB2 = eigenspace(B, Fraction(2))
    print('B: E_1 basis', show(EB1), ' E_2 basis', show(EB2))
    assert [len(EB1), len(EB2)] == [1, 2]
    a1, a2, a3 = Q([[3, -1, 3], [2, 1, 0], [2, 0, 1]])
    for c, v in ((1, a1), (2, a2), (2, a3)):
        assert is_eigenpair(B, Fraction(c), v)
    P = columns_to_matrix([a1, a2, a3])
    Pinv = inverse(P)
    D = mat_mul(mat_mul(Pinv, B), P)
    print('P =', show(P), ' det P =', det(P))
    print('P^-1 =', show(Pinv))
    print('P^-1 B P =', show(D))
    assert D == Q([[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    assert det(P) == -1
    assert mat_mul(B, P) == mat_mul(P, D)

    # ---- same charpoly, different similarity class ------------------
    rA = rank(mat_add(A, scalar_mul(2, identity(3)), -1))
    rB = rank(mat_add(B, scalar_mul(2, identity(3)), -1))
    print('rank(A-2I) =', rA, ' rank(B-2I) =', rB)
    assert (rA, rB) == (2, 1)       # rank is a similarity invariant
    I3 = identity(3)

    def f(M):
        return mat_mul(mat_add(M, I3, -1),
                       mat_add(M, scalar_mul(2, I3), -1))

    zero = Q([[0] * 3] * 3)
    assert f(B) == zero and f(A) != zero
    print('(A-I)(A-2I) =', show(f(A)))

    # ---- opening problem: powers -----------------------------------
    for k in (10, 20):
        Ak, Bk = mat_pow(A, k), mat_pow(B, k)
        assert Bk == mat_mul(mat_mul(P, mat_pow(D, k)), Pinv)
        mA = max(abs(x) for r in Ak for x in r) / 2 ** k
        mB = max(abs(x) for r in Bk for x in r) / 2 ** k
        print('k =', k, ' max|A^k|/2^k =', float(mA),
              ' max|B^k|/2^k =', float(mB))
    A10, B10 = mat_pow(A, 10), mat_pow(B, 10)
    print('A^10 =', show(A10))
    print('B^10 =', show(B10))
    assert max(abs(x) for r in mat_pow(A, 20) for x in r) / 2 ** 20 > \
        1.9 * max(abs(x) for r in A10 for x in r) / 2 ** 10
    assert max(abs(x) for r in mat_pow(B, 20) for x in r) / 2 ** 20 < \
        max(abs(x) for r in B10 for x in r) / 2 ** 10 * 1.01

    # ---- field dependence: the rotation ----------------------------
    R = Q([[0, -1], [1, 0]])
    cpR = charpoly(R)
    print('rotation charpoly:', [str(x) for x in cpR])
    assert cpR == Q([[1, 0, 1]])[0]
    for c in range(-5, 6):
        assert poly_eval(cpR, Fraction(c)) >= 1   # x^2 + 1 >= 1 on Q
    assert null_basis(mat_add(R, identity(2), 1)) == []
    for c, v in ((1j, [1, -1j]), (-1j, [1, 1j])):
        Rv = [R[0][0] * v[0] + R[0][1] * v[1],
              R[1][0] * v[0] + R[1][1] * v[1]]
        assert Rv == [c * x for x in v]
    print('R(1,-i) = i(1,-i) and R(1,i) = -i(1,i) verified')

    # ---- similarity invariance (Theorem 3.1) ----------------------
    S = Q([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
    assert det(S) == 2
    A_s = mat_mul(mat_mul(inverse(S), A), S)
    assert charpoly(A_s) == charpoly(A)
    v = Q([[1, 0, 2]])[0]
    assert is_eigenpair(A_s, Fraction(1), mat_vec(inverse(S), v))
    print('charpoly(S^-1 A S) = charpoly(A);  S^-1 (1,0,2) is a 1-eigenvector')

    # ---- coefficients of charpoly (Theorem 2.2), a 4x4 test -------
    M4 = Q([[2, -1, 0, 3], [1, 4, 1, 0], [0, 2, -3, 1], [5, 0, 1, 1]])
    c4 = charpoly(M4)
    print('4x4 charpoly:', [str(x) for x in c4], ' tr =', trace(M4),
          ' det =', det(M4))
    assert c4[1] == -trace(M4) and c4[4] == det(M4)

    # ---- geometric <= algebraic (Theorem 6.2) --------------------
    J = Q([[2, 1], [0, 2]])
    assert multiplicity(charpoly(J), Fraction(2)) == 2
    assert len(eigenspace(J, Fraction(2))) == 1
    for M in (A, B, J, R):
        cp = charpoly(M)
        for c in range(-3, 4):
            assert len(eigenspace(M, Fraction(c))) <= multiplicity(
                cp, Fraction(c))
    print('Jordan block: algebraic multiplicity 2, geometric 1')

    # ---- distinct eigenvalues (Theorem 6.1, Exercise 2) -----------
    C = Q([[1, 1, 0], [0, 2, 1], [0, 0, 3]])
    assert charpoly(C) == [Fraction(x) for x in (1, -6, 11, -6)]
    vecs = []
    for c in (1, 2, 3):
        E = eigenspace(C, Fraction(c))
        assert len(E) == 1
        vecs.append(E[0])
    print('C eigenvectors:', [show([v])[0] for v in vecs])
    PC = columns_to_matrix(Q([[1, 0, 0], [1, 1, 0], [1, 2, 2]]))
    for c, col in zip((1, 2, 3), ([1, 0, 0], [1, 1, 0], [1, 2, 2])):
        assert is_eigenpair(C, Fraction(c), Q([col])[0])
    DC = mat_mul(mat_mul(inverse(PC), C), PC)
    assert DC == Q([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    assert rank(columns_to_matrix(vecs)) == 3
    print('P_C^-1 C P_C =', show(DC))

    # ---- Exercise 1: p(T) alpha = p(c) alpha -----------------------
    pB = mat_add(mat_add(mat_pow(B, 2), B, -3), I3)       # B^2-3B+I
    for c, v in ((1, a1), (2, a2), (2, a3)):
        p_c = c * c - 3 * c + 1
        assert mat_vec(pB, v) == [p_c * x for x in v]
    print('p(x) = x^2-3x+1: p(B)v = p(c)v on all three eigenvectors')

    # ---- Exercise 3: rotation squared ------------------------------
    assert mat_mul(R, R) == Q([[-1, 0], [0, -1]])
    assert len(eigenspace(mat_mul(R, R), Fraction(-1))) == 2

    # ---- Exercise 4: field dependence ------------------------------
    K = Q([[0, 2], [1, 0]])
    cpK = charpoly(K)
    assert cpK == Q([[1, 0, -2]])[0]
    assert all(poly_eval(cpK, Fraction(c)) != 0
               for c in (-2, -1, 1, 2, Fraction(1, 2), Fraction(-1, 2)))
    root = 2 ** 0.5
    for lam, vec in ((root, [root, 1]), (-root, [-root, 1])):
        Kv = [2 * vec[1], vec[0]]
        assert all(abs(x - lam * y) < 1e-12 for x, y in zip(Kv, vec))
    print('K charpoly x^2-2: no rational root; real eigenvalues +-sqrt 2')

    # ---- Exercise 6: same charpoly, both diagonalizable => similar --
    D2 = Q([[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    B2 = mat_mul(mat_mul(S, D2), inverse(S))
    assert charpoly(B2) == charpoly(B) and B2 != B
    W = mat_mul(S, Pinv)                 # B2 = W B W^-1
    assert mat_mul(mat_mul(W, B), inverse(W)) == B2
    print('B2 = S diag(1,2,2) S^-1 =', show(B2), ' is similar to B')
    print('all assertions passed')


if __name__ == '__main__':
    main()
