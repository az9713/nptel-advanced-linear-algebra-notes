"""Chapter 18 laboratory: invariant direct sums and block-diagonal operators.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of rows,
a vector is a list, a polynomial is a list of coefficients from degree 0
upward.  Every routine returns new objects and never changes its input.

What the checks cover: the chapter's running operator on Q^3, the
commutation criterion, the spectral form T = c1 E1 + c2 E2, the Lagrange
formula E_j = p_j(T), the converse theorem on a second example, and three
boundary cases.  They verify finite examples only; they prove no theorem.
"""
from fractions import Fraction


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def zeros(m, n):
    return [[Fraction(0)] * n for _ in range(m)]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def matadd(A, B):
    return [[a + b for a, b in zip(r, s)] for r, s in zip(A, B)]


def scale(c, A):
    return [[Fraction(c) * a for a in row] for row in A]


def matsub(A, B):
    return matadd(A, scale(-1, B))


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def column(A, j):
    return [row[j] for row in A]


def matvec(A, v):
    return [sum(a * x for a, x in zip(row, v)) for row in A]


def rref(A):
    """Row-reduced echelon form and the list of pivot columns."""
    R = [row[:] for row in A]
    pivots = []
    r = 0
    for c in range(len(R[0]) if R else 0):
        p = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [x - f * y for x, y in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
        if r == len(R):
            break
    return R, pivots


def rank(A):
    return len(rref(A)[1])


def nullspace(A):
    """A basis (list of vectors) of {x : A x = 0}."""
    R, pivots = rref(A)
    n = len(A[0])
    basis = []
    for f in (c for c in range(n) if c not in pivots):
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for r, p in enumerate(pivots):
            v[p] = -R[r][f]
        basis.append(v)
    return basis


def inverse(A):
    n = len(A)
    aug = [row + e for row, e in zip(A, eye(n))]
    R, pivots = rref(aug)
    assert pivots[:n] == list(range(n)), 'matrix is singular'
    return [row[n:] for row in R]


def is_invariant(A, vectors):
    """True when A maps span(vectors) into span(vectors)."""
    base = [list(v) for v in vectors]
    image = base + [matvec(A, v) for v in vectors]
    return rank(image) == rank(base)


# Polynomials: coefficient lists, lowest degree first.

def poly_mul(f, g):
    out = [Fraction(0)] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i + j] += a * b
    return out


def poly_at_matrix(f, A):
    """f(A) by Horner's rule."""
    n = len(A)
    out = zeros(n, n)
    for c in reversed(f):
        out = matadd(matmul(out, A), scale(c, eye(n)))
    return out


def lagrange(cs, j):
    """p_j(x) = product over i != j of (x - c_i)/(c_j - c_i)."""
    p = [Fraction(1)]
    for i, ci in enumerate(cs):
        if i != j:
            d = Fraction(cs[j]) - ci
            p = poly_mul(p, [-Fraction(ci) / d, 1 / d])
    return p


def charpoly(A):
    """Faddeev-LeVerrier: monic det(xI - A), lowest degree first."""
    n = len(A)
    coeffs = [Fraction(0)] * (n + 1)
    coeffs[n] = Fraction(1)
    M = zeros(n, n)
    for k in range(1, n + 1):
        M = matadd(matmul(A, M), scale(coeffs[n - k + 1], eye(n)))
        coeffs[n - k] = -trace(matmul(A, M)) / k
    return coeffs


def min_poly(A):
    """Monic polynomial of least degree d with A^d in span(I..A^(d-1))."""
    n = len(A)
    powers = [eye(n)]
    for d in range(1, n * n + 1):
        powers.append(matmul(powers[-1], A))
        flat = [[x for row in P for x in row] for P in powers]
        cols = [[flat[k][i] for k in range(len(flat))]
                for i in range(n * n)]
        null = nullspace(cols)
        if null:
            v = null[0]
            return [x / v[-1] for x in v]
    raise AssertionError('no dependency found')


def poly_from_roots(roots):
    p = [Fraction(1)]
    for r in roots:
        p = poly_mul(p, [-Fraction(r), Fraction(1)])
    return p


def poly_str(f):
    terms = []
    for k in range(len(f) - 1, -1, -1):
        if f[k] != 0:
            terms.append('%s*x^%d' % (f[k], k))
    return ' + '.join(terms)


def show(M):
    """Compact text for a matrix with integer or Fraction entries."""
    return '[' + '; '.join(' '.join(str(x) for x in row)
                           for row in M) + ']'


def diag_selector(n, idx):
    D = zeros(n, n)
    for i in idx:
        D[i][i] = Fraction(1)
    return D


def is_zero(A):
    return all(x == 0 for row in A for x in row)


def commutes(A, B):
    return matmul(A, B) == matmul(B, A)


# ---------------------------------------------------------------------
# 1. The running operator T on Q^3, in the standard basis.
A = Q([[3, 0, -2], [0, 3, -2], [0, 0, 1]])
b1, b2, b3 = [1, 1, 1], [1, 0, 0], [0, 1, 0]
P = Q([[1, 1, 0], [1, 0, 1], [1, 0, 0]])  # columns b1, b2, b3
Pinv = inverse(P)
assert [column(P, j) for j in range(3)] == Q([b1, b2, b3])

# W1 = span(b1) and W2 = span(b2, b3) are invariant; matrix is block-diag.
assert is_invariant(A, Q([b1])) and is_invariant(A, Q([b2, b3]))
D = matmul(Pinv, matmul(A, P))
assert D == Q([[1, 0, 0], [0, 3, 0], [0, 0, 3]])
print('[T]_B =', show(D))

# 2. The projections of V = W1 + W2 along the decomposition.
E1 = matmul(P, matmul(diag_selector(3, [0]), Pinv))
E2 = matmul(P, matmul(diag_selector(3, [1, 2]), Pinv))
assert E1 == Q([[0, 0, 1], [0, 0, 1], [0, 0, 1]])
assert E2 == Q([[1, 0, -1], [0, 1, -1], [0, 0, 0]])
assert matmul(E1, E1) == E1 and matmul(E2, E2) == E2
assert is_zero(matmul(E1, E2)) and is_zero(matmul(E2, E1))
assert matadd(E1, E2) == eye(3)
assert rank(E1) == 1 and rank(E2) == 2
print('E1 =', show(E1))
print('E2 =', show(E2))

# 3. Commutation criterion: A commutes with E1, E2.
assert commutes(A, E1) and commutes(A, E2)
print('A E1 = E1 A and A E2 = E2 A: True')

# 4. Spectral form and powers.
assert matadd(E1, scale(3, E2)) == A
for n in (2, 5, 10):
    An = eye(3)
    for _ in range(n):
        An = matmul(An, A)
    assert An == matadd(E1, scale(3 ** n, E2))
print('A^10 =', show(matadd(E1, scale(3 ** 10, E2))))

# g(x) = x^3 - 2x gives g(A) = g(1) E1 + g(3) E2 = -E1 + 21 E2.
g = Q([[0, -2, 0, 1]])[0]
assert poly_at_matrix(g, A) == matadd(scale(-1, E1), scale(21, E2))
print('g(1), g(3) =', 1 - 2, 27 - 6)

# 5. Lagrange polynomials give the projections.
cs = [1, 3]
p1, p2 = lagrange(cs, 0), lagrange(cs, 1)
assert p1 == [Fraction(3, 2), Fraction(-1, 2)]
assert p2 == [Fraction(-1, 2), Fraction(1, 2)]
assert poly_at_matrix(p1, A) == E1 and poly_at_matrix(p2, A) == E2
print('p1(x) =', poly_str(p1), '; p2(x) =', poly_str(p2))

# Minimal and characteristic polynomials of A.
assert min_poly(A) == poly_from_roots([1, 3])
assert charpoly(A) == poly_from_roots([1, 3, 3])
print('m_A =', poly_str(min_poly(A)), '; c_A =', poly_str(charpoly(A)))

# 6. Converse theorem on a second example with three summands.
cs3 = [2, -1, 7]
Es = [matmul(P, matmul(diag_selector(3, [i]), Pinv)) for i in range(3)]
T = zeros(3, 3)
for c, E in zip(cs3, Es):
    T = matadd(T, scale(c, E))
total = zeros(3, 3)
for E in Es:
    total = matadd(total, E)
assert total == eye(3)
for i in range(3):
    for j in range(3):
        if i != j:
            assert is_zero(matmul(Es[i], Es[j]))
    assert matmul(Es[i], Es[i]) == Es[i] and rank(Es[i]) == 1
    kernel = nullspace(matsub(T, scale(cs3[i], eye(3))))
    assert len(kernel) == rank(Es[i])
    assert is_zero(matmul(matsub(T, scale(cs3[i], eye(3))), Es[i]))
assert charpoly(T) == poly_from_roots(cs3)
assert min_poly(T) == poly_from_roots(cs3)
print('T =', show(T), 'has eigenvalues', cs3)
print('m_T =', poly_str(min_poly(T)))

# 7. Boundary: W1 = span(e1) is invariant, W2 = span(e2, e3) is not.
e1, e2, e3 = [1, 0, 0], [0, 1, 0], [0, 0, 1]
assert is_invariant(A, Q([e1])) and not is_invariant(A, Q([e2, e3]))
F1 = diag_selector(3, [0])
assert not commutes(A, F1)
print('invariant W1 with non-invariant complement: A F1 != F1 A')

# 8. Boundary: invariant decomposition, block-diagonal, not diagonalizable.
N = Q([[2, 1, 0], [0, 2, 0], [0, 0, 3]])
assert is_invariant(N, Q([e1, e2])) and is_invariant(N, Q([e3]))
G1, G2 = diag_selector(3, [0, 1]), diag_selector(3, [2])
assert commutes(N, G1) and commutes(N, G2)
mN = min_poly(N)
assert mN == poly_from_roots([2, 2, 3])
assert charpoly(N) == poly_from_roots([2, 2, 3])
# m_N is the lcm of the minimal polynomials of the blocks (x-2)^2, (x-3).
assert min_poly(N) == poly_mul(poly_from_roots([2, 2]), [-3, 1])
# T = a G1 + b G2 is impossible: it is diagonal, N is not.
assert N[0][1] != 0
print('m_N =', poly_str(mN), '= (x-2)^2 (x-3): not diagonalizable')

# 9. Boundary: converse needs nonzero E_i and distinct c_i.
T_id = matadd(matadd(E1, E2), scale(5, zeros(3, 3)))  # 5 * (zero E)
assert T_id == eye(3)
assert charpoly(T_id) == poly_from_roots([1, 1, 1])  # 5 is not a root
T_eq = matadd(scale(4, E1), scale(4, E2))
assert T_eq == scale(4, eye(3))
assert len(nullspace(matsub(T_eq, scale(4, eye(3))))) == 3
assert rank(E1) == 1  # ker(T - 4I) has dimension 3, range E1 has 1
print('zero E_2 with c_2 = 5: 5 is not an eigenvalue of I')
print('c_1 = c_2 = 4: ker(T - 4I) has dimension 3, range E1 has dim 1')
# Section 6 matrix S and the exercises.
S = matadd(matadd(scale(2, Es[0]), scale(-1, Es[1])), scale(7, Es[2]))
assert S == Q([[-1, 0, 3], [0, 7, -5], [0, 0, 2]]) and S == T
Ainv = matadd(E1, scale(Fraction(1, 3), E2))
assert matmul(A, Ainv) == eye(3)
assert Ainv == [[Fraction(1, 3), 0, Fraction(2, 3)],
                [0, Fraction(1, 3), Fraction(2, 3)], [0, 0, 1]]
q1 = [x / -15 for x in poly_from_roots([-1, 7])]
assert poly_at_matrix(q1, S) == Es[0]
M2 = matsub(N, scale(2, eye(3)))
assert matmul(M2, M2) == diag_selector(3, [2])
print('S =', show(S), '; A^-1 =', show(Ainv))
print('all checks passed')
