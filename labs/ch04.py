"""Chapter 4 laboratory: independence, bases, dimension, coordinates.

Exact arithmetic over Q (fractions.Fraction).  A vector is a list of
Fractions; a list of vectors is a list of such lists.  Matrices are
lists of rows.  Vectors are written as columns in the chapter.
"""
from fractions import Fraction
from itertools import product, combinations
import random


def Q(vec):
    return [Fraction(x) for x in vec]


def rref(A):
    """Row-reduced echelon form of A (rows), with the pivot columns."""
    R = [row[:] for row in A]
    pivots, r = [], 0
    for c in range(len(R[0]) if R else 0):
        piv = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if piv is None:
            continue
        R[r], R[piv] = R[piv], R[r]
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


def cols_to_matrix(vectors):
    """Matrix whose columns are the given vectors."""
    return [list(row) for row in zip(*vectors)]


def rank(vectors):
    if not vectors:
        return 0
    return len(rref(cols_to_matrix(vectors))[1])


def is_independent(vectors):
    return rank(vectors) == len(vectors)


def null_space(A):
    """Basis of {x : A x = 0}, one vector per free column."""
    R, piv = rref(A)
    n = len(A[0])
    basis = []
    for f in (c for c in range(n) if c not in piv):
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, p in zip(R, piv):
            x[p] = -row[f]
        basis.append(x)
    return basis


def solve(vectors, v):
    """Coefficients c with sum c_j vectors[j] = v, or None."""
    A = cols_to_matrix(vectors)
    aug = [row + [y] for row, y in zip(A, v)]
    R, piv = rref(aug)
    n = len(vectors)
    if n in piv:
        return None
    c = [Fraction(0)] * n
    for row, p in zip(R, piv):
        c[p] = row[n]
    return c


def in_span(vectors, v):
    return solve(vectors, v) is not None


def basis_from_spanning(vectors):
    """Indices of a basis of the span: the pivot columns."""
    return rref(cols_to_matrix(vectors))[1]


def extend_to_basis(indep, n):
    """Extend an independent list of Q^n by standard vectors."""
    basis = [v[:] for v in indep]
    for i in range(n):
        e = [Fraction(int(i == j)) for j in range(n)]
        if not in_span(basis, e):
            basis.append(e)
    return basis


def matmul(A, B):
    return [[sum((a * b for a, b in zip(row, col)), Fraction(0))
             for col in zip(*B)] for row in A]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def inverse(P):
    n = len(P)
    aug = [row + idr for row, idr in zip(P, identity(n))]
    R, piv = rref(aug)
    if piv[:n] != list(range(n)):
        return None
    return [row[n:] for row in R]


def intersection_basis(U, W):
    """Basis of span(U) cap span(W) from the null space of [U | -W]."""
    M = cols_to_matrix(U + [[-x for x in w] for w in W])
    out = []
    for z in null_space(M):
        v = [sum((c * u[i] for c, u in zip(z[:len(U)], U)), Fraction(0))
             for i in range(len(U[0]))]
        out.append(v)
    return [out[i] for i in basis_from_spanning(out)] if out else []


def fmt(v):
    return '(' + ', '.join(str(x) for x in v) + ')'


# ------------------------------------------------------------------
# Running example (original): the list v1..v4 in Q^4.
v1, v2, v3 = Q([1, 1, 0, 0]), Q([0, 1, 1, 0]), Q([0, 0, 1, 1])
v4 = [a + 2 * b - c for a, b, c in zip(v1, v2, v3)]
assert v4 == Q([1, 3, 1, -1])
L = [v1, v2, v3, v4]

# 1. Dependence, redundancy, dimension of the span.
assert not is_independent(L) and rank(L) == 3
assert null_space(cols_to_matrix(L)) == [Q([-1, -2, 1, 1])]
assert is_independent([v1, v2, v3])
assert basis_from_spanning(L) == [0, 1, 2]
assert in_span([v1, v2, v3], v4)
print('rank of (v1, v2, v3, v4) =', rank(L))
print('relation: v1 + 2 v2 - v3 - v4 = 0; null space =',
      fmt(null_space(cols_to_matrix(L))[0]))

# 2. Lemma: S independent and beta outside span(S) => S + [beta]
# independent; and beta inside span(S) => dependent.
S = [v1, v2]
beta_out, beta_in = v3, [a + b for a, b in zip(v1, v2)]
assert not in_span(S, beta_out) and is_independent(S + [beta_out])
assert in_span(S, beta_in) and not is_independent(S + [beta_in])

# 3. Extension of an independent list to a basis of Q^4.
ext = extend_to_basis([v1, v2, v3], 4)
assert len(ext) == 4 and is_independent(ext)
print('extend (v1, v2, v3) to a basis of Q^4 by', fmt(ext[3]))
ext1 = extend_to_basis([v1], 4)
assert len(ext1) == 4 and is_independent(ext1)

# 4. Coordinates of v4 relative to B = (v1, v2, v3).
B = [v1, v2, v3]
x = solve(B, v4)
assert x == Q([1, 2, -1])
print('[v4]_B =', fmt(x))
assert solve([v1, v2, v3, v4], v4) is not None
# uniqueness of coordinates fails for the dependent list L
c1, c2 = Q([0, 0, 0, 1]), Q([1, 2, -1, 0])
assert [sum(c * v[i] for c, v in zip(c1, L)) for i in range(4)] == v4
assert [sum(c * v[i] for c, v in zip(c2, L)) for i in range(4)] == v4

# 5. Change of ordered basis.  B' = (v1, v1 + v2, v2 + v3).
Bp = [v1, [a + b for a, b in zip(v1, v2)], [a + b for a, b in zip(v2, v3)]]
assert is_independent(Bp)
P = cols_to_matrix([solve(B, w) for w in Bp])
assert P == [Q([1, 1, 0]), Q([0, 1, 1]), Q([0, 0, 1])]
xp = solve(Bp, v4)
assert xp == Q([-2, 3, -1])
Pxp = [sum(P[i][j] * xp[j] for j in range(3)) for i in range(3)]
assert Pxp == x
Pinv = inverse(P)
assert matmul(P, Pinv) == identity(3) and matmul(Pinv, P) == identity(3)
print('[v4]_B\' =', fmt(xp), ' P x\' =', fmt(Pxp))

# 6. Dimension formula for the running pair U, W in Q^4.
U = [v1, v2, v3]
w1, w2 = v1, Q([0, 1, 0, 1])
Wl = [w1, w2]
dU, dW, dS = rank(U), rank(Wl), rank(U + Wl)
inter = intersection_basis(U, Wl)
dI = rank(inter)
assert (dU, dW, dS, dI) == (3, 2, 4, 1) and dU + dW == dI + dS
assert inter and rank(inter + [v1]) == 1
print('dim U, dim W, dim(U+W), dim(U cap W) =', dU, dW, dS, dI)

# 7. Randomized check of the dimension formula over Q.
rng = random.Random(4)
for _ in range(200):
    n = rng.randint(2, 5)
    a = [Q([rng.randint(-2, 2) for _ in range(n)])
         for _ in range(rng.randint(1, 4))]
    b = [Q([rng.randint(-2, 2) for _ in range(n)])
         for _ in range(rng.randint(1, 4))]
    ia = intersection_basis(a, b)
    assert rank(a) + rank(b) == rank(ia) + rank(a + b)
print('dimension formula: 200 random pairs in Q^2..Q^5 agree')

# 8. Boundary: inclusion-exclusion fails for three subspaces.
l1, l2, l3 = Q([1, 0]), Q([0, 1]), Q([1, 1])
assert rank([l1, l2, l3]) == 2
assert 1 + 1 + 1 - 0 - 0 - 0 + 0 == 3 != 2
print('three lines in a plane: dim of sum = 2, inclusion-exclusion = 3')

# 9. Invertible matrix <=> columns form a basis (n = 3).
Pm = [Q([1, 0, 0]), Q([1, 1, 0]), Q([0, 1, 1])]     # columns
assert is_independent(Pm) and rank(Pm) == 3
assert inverse(cols_to_matrix(Pm)) is not None
sing = [Q([1, 2, 3]), Q([2, 4, 6]), Q([0, 0, 1])]
assert not is_independent(sing) and inverse(cols_to_matrix(sing)) is None

# 10. Polynomials of degree <= 2 over Q: coefficient vectors.
mono = [Q([1, 0, 0]), Q([0, 1, 0]), Q([0, 0, 1])]
assert is_independent(mono)
assert all(rank(mono + [Q(v)]) == 3 for v in [[1, 2, 3], [0, 0, 5]])
# A polynomial of degree 2 has at most 2 roots: a function that
# vanishes at three Q-points has all coefficients zero.
pts = [Fraction(0), Fraction(1), Fraction(2)]
vander = [[t ** k for k in range(3)] for t in pts]
assert inverse(vander) is not None


# 11. Boundary: over F_2 the functions 1, x, x^2 are dependent.
def fn(coeffs, p):
    return tuple(sum(c * t ** k for k, c in enumerate(coeffs)) % p
                 for t in range(p))


assert fn([0, 1, 1], 2) == (0, 0) == fn([0, 0, 0], 2)
print('over F_2: the polynomial x + x^2 is the zero function on F_2')

# 12. Exhaustive check of invariance of dimension over F_2^3.
vecs = [v for v in product((0, 1), repeat=3)]


def span2(sel):
    out = {(0, 0, 0)}
    for v in sel:
        out |= {tuple((a + b) % 2 for a, b in zip(o, v)) for o in out}
    return out


def indep2(sel):
    return len(span2(sel)) == 2 ** len(sel)


full = set(vecs)
nonzero = [v for v in vecs if any(v)]
sizes_basis, sizes_min_span, sizes_max_indep = set(), set(), set()
n_ordered = 0
for k in range(0, 8):
    for sel in combinations(nonzero, k):
        ind, sp = indep2(sel), span2(sel) == full
        if ind and sp:
            sizes_basis.add(k)
            n_ordered += 1
        if sp and all(span2(sel[:i] + sel[i + 1:]) != full
                      for i in range(k)):
            sizes_min_span.add(k)
        if ind and all(not indep2(sel + (v,)) for v in nonzero
                       if v not in sel):
            sizes_max_indep.add(k)
assert sizes_basis == sizes_min_span == sizes_max_indep == {3}
assert n_ordered * 6 == 7 * 6 * 4
print('F_2^3: unordered bases', n_ordered, '(ordered', n_ordered * 6, ')',
      '; every minimal spanning and maximal independent set has size 3')
print('ch04 laboratory: all assertions passed')
