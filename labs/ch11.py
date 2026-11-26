"""Chapter 11 laboratory: annihilators, the double dual, the transpose.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows.  A vector of F^n is a list; a linear functional on F^n is a ROW
(a list), so the functional f acts on the vector v as sum_i f_i v_i.
The map T(x) = A x has matrix A; its transpose T^t sends the functional
g to g o T, whose row is g A.
"""
from fractions import Fraction


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(col) for col in zip(*A)]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fraction(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def apply(A, x):
    """Column vector A x."""
    return [sum((a * b for a, b in zip(row, x)), Fraction(0)) for row in A]


def pair(f, v):
    """Value f(v) of the functional f (a row) at the vector v."""
    return sum((a * b for a, b in zip(f, v)), Fraction(0))


def rref(A):
    """Row-reduced echelon form and the list of pivot columns."""
    R = [row[:] for row in A]
    pivots, r = [], 0
    for c in range(len(R[0]) if R else 0):
        p = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                R[i] = [x - R[i][c] * y for x, y in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
    return R, pivots


def rank(A):
    return len(rref(A)[1])


def nullspace(A):
    """Basis of {x : A x = 0}, one vector per free column."""
    R, piv = rref(A)
    n = len(A[0])
    basis = []
    for free in (c for c in range(n) if c not in piv):
        x = [Fraction(0)] * n
        x[free] = Fraction(1)
        for row, c in zip(R, piv):
            x[c] = -row[free]
        basis.append(x)
    return basis


def inverse(M):
    n = len(M)
    aug = [row + e for row, e in zip(M, identity(n))]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n)), 'matrix is singular'
    return [row[n:] for row in R]


def annihilator(vectors, n):
    """Basis (rows) of S^0 for S = span(vectors) inside F^n."""
    if not vectors:
        return identity(n)
    return nullspace(vectors)


def dual_basis(cols):
    """Dual basis (rows) of the basis given as a list of vectors."""
    M = transpose(cols)
    return inverse(M)


def span_contains(vectors, y):
    """Is y in the span of vectors?  Compare ranks."""
    return rank(vectors) == rank(vectors + [y])


def matrix_in_bases(A, B, C):
    """[T]_{C<-B}: column j holds the C-coordinates of T(b_j)."""
    Cmat = transpose(C)
    Bmat = transpose(B)
    return matmul(inverse(Cmat), matmul(A, Bmat))


def rank_mod_p(A, p):
    """Rank over the prime field F_p of an integer matrix."""
    R = [[x % p for x in row] for row in A]
    r = 0
    for c in range(len(R[0])):
        piv = next((i for i in range(r, len(R)) if R[i][c]), None)
        if piv is None:
            continue
        R[r], R[piv] = R[piv], R[r]
        inv = pow(R[r][c], p - 2, p)
        R[r] = [x * inv % p for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c]
                R[i] = [(x - f * y) % p for x, y in zip(R[i], R[r])]
        r += 1
    return r


def show(M):
    def cell(x):
        return str(x)
    if M and isinstance(M[0], list):
        return '[' + '; '.join(' '.join(cell(x) for x in r)
                               for r in M) + ']'
    return '(' + ', '.join(cell(x) for x in M) + ')'


# Running example (an original construction): T(x) = A x, Q^4 -> Q^3.
A = Q([[1, 2, 0, 1], [0, 1, 1, -1], [1, 3, 1, 0]])
m, n = 3, 4

# 1. Rank, kernel, range, and the two annihilators.
r = rank(A)
ker = nullspace(A)
assert r == 2 and len(ker) == n - r
assert ker == Q([[2, -1, 1, 0], [-3, 1, 0, 1]])
cols = transpose(A)
ann_range = annihilator(cols, m)          # (im T)^0 as rows
assert len(ann_range) == m - r
ker_t = nullspace(transpose(A))           # ker T^t: g with A^T g = 0
assert ker_t == Q([[-1, -1, 1]]) == ann_range
ann_ker = annihilator(ker, n)             # (ker T)^0 as rows
assert len(ann_ker) == n - len(ker) == r
for f in ann_ker:
    assert all(pair(f, v) == 0 for v in ker)
assert rank(A + ann_ker) == r             # (ker T)^0 = row space of A
print('rank T =', r, '; ker T basis', show(ker))
print('ker T^t = (im T)^0 spanned by', show(ker_t[0]))
print('(ker T)^0 basis (rows)', show(ann_ker))

# 2. The consistency test: A x = y is solvable iff g(y) = 0.
g = ker_t[0]
for y, ok in ((Q([[1, 2, 3]])[0], True), (Q([[1, 2, 4]])[0], False)):
    assert span_contains(cols, y) == ok
    assert (pair(g, y) == 0) == ok
print('y=(1,2,3): g(y) =', pair(g, Q([[1, 2, 3]])[0]),
      '(solvable); y=(1,2,4): g(y) =', pair(g, Q([[1, 2, 4]])[0]))

# 3. Dimension formula dim S + dim S^0 = dim V for several subspaces.
for S in (Q([[1, 1, 0, 0], [0, 0, 1, 1]]), Q([[1, 2, 3, 4]]),
          Q([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), []):
    assert rank(S if S else Q([[0, 0, 0, 0]])) + \
        len(annihilator(S, 4)) == 4
print('dim S + dim S^0 = 4 checked for subspaces of dimension 2, 1, 4, 0')

# 4. Every basis of V* is dual to a basis of V (section 3 example).
F = Q([[1, 2, 0, 1], [0, 1, 1, -1], [0, 0, 1, 0], [0, 0, 0, 1]])
vecs = transpose(inverse(F))              # columns of F^{-1}
assert vecs == Q([[1, 0, 0, 0], [-2, 1, 0, 0],
                  [2, -1, 1, 0], [-3, 1, 0, 1]])
for i, f in enumerate(F):
    for j, v in enumerate(vecs):
        assert pair(f, v) == int(i == j)
print('dual basis of f1..f4 in V:', ' '.join(show(v) for v in vecs))

# 5. Matrix of T in adapted bases, and the transpose in dual bases.
Bb = Q([[1, 0, 0, 0], [0, 1, 0, 0], [2, -1, 1, 0], [-3, 1, 0, 1]])
Cb = Q([[1, 0, 1], [0, 1, 1], [0, 0, 1]])
phi = dual_basis(Bb)
gam = dual_basis(Cb)
assert phi == Q([[1, 0, -2, 3], [0, 1, 1, -1], [0, 0, 1, 0], [0, 0, 0, 1]])
assert gam == Q([[1, 0, 0], [0, 1, 0], [-1, -1, 1]])
MT = matrix_in_bases(A, Bb, Cb)
assert MT == Q([[1, 2, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]])
# Matrix of T^t built independently: entry (k, j) is (gam_j o T)(b_k).
Mt = [[pair(gam[j], apply(A, Bb[k])) for j in range(m)] for k in range(n)]
assert Mt == transpose(MT)
# and gam_j o T = sum_k Mt[k][j] phi_k as rows.
for j in range(m):
    row = [pair(gam[j], c) for c in transpose(A)]
    comb = [sum((Mt[k][j] * phi[k][i] for k in range(n)), Fraction(0))
            for i in range(n)]
    assert row == comb
print('[T] in adapted bases =', show(MT))
print('[T^t] in dual bases  =', show(Mt))
print('T^t(gamma_3) =', show([sum((Mt[k][2] * phi[k][i] for k in range(n)),
                                  Fraction(0)) for i in range(n)]))

# 6. Rank of T^t equals rank of T; ranks over other fields.
assert rank(transpose(A)) == rank(A) == 2
for M in (Q([[1, 2, 3], [2, 4, 6]]), Q([[0, 1], [0, 0], [0, 0]]),
          Q([[1, 0, 2, 1], [3, 1, 0, 0], [4, 1, 2, 1]])):
    assert rank(M) == rank(transpose(M))
N2 = [[1, 1], [1, 1]]
N3 = [[1, 2, 0], [2, 1, 0], [0, 0, 1]]
assert rank_mod_p(N2, 2) == rank_mod_p(transpose(N2), 2) == 1
assert rank_mod_p(N3, 3) == rank_mod_p(transpose(N3), 3) == 2
assert rank(Q(N3)) == 3     # the rank depends on the field, not on T^t
print('rank A = rank A^T over Q, F_2, F_3 (rank of N3: Q 3, F_3 2)')

# 7. T injective but T^t not injective; T^t is surjective.
Ainj = Q([[1], [0]])                       # T: Q -> Q^2, x -> (x, 0)
assert rank(Ainj) == 1 and nullspace(Ainj) == []
assert len(nullspace(transpose(Ainj))) == 1
assert rank(transpose(Ainj)) == 1          # onto Q^* (dimension 1)
print('T(x)=(x,0): ker T = 0 but ker T^t =',
      show(nullspace(transpose(Ainj))[0]))

# 8. Section 4 example: T(x1, x2) = (x1, 0); T^t(a, b) = (a, 0).
A22 = Q([[1, 0], [0, 0]])
for a, b in ((2, 5), (-1, 3)):
    gg = Q([[a, b]])[0]
    gt = [sum((gg[i] * A22[i][j] for i in range(2)), Fraction(0))
          for j in range(2)]
    assert gt == Q([[a, 0]])[0]
print('T(x1,x2)=(x1,0): T^t g has row (a, 0) for g = (a, b)')

# 9. (S T)^t = T^t S^t on matrices; I^t = I.
S = Q([[1, 1, 0], [0, 2, 1]])
assert transpose(matmul(S, A)) == matmul(transpose(A), transpose(S))
assert transpose(identity(4)) == identity(4)
print('(S T)^t = T^t S^t: checked for S 2x3, T 3x4')

# 10. Separation lemma: for v != 0 construct f with f(v) = 1.


def separating_functional(v):
    basis = [v]
    for e in identity(len(v)):
        if rank(basis + [e]) > len(basis):
            basis.append(e)
    return dual_basis(basis)[0]


for v in ker + [Q([[0, 0, 0, 1]])[0], Q([[3, -2, 5, 7]])[0]]:
    assert pair(separating_functional(v), v) == 1
print('separating functionals built for 4 nonzero vectors')

# 11. Exercises 2, 4, 5 of the chapter.
F2 = Q([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
half = Fraction(1, 2)
assert inverse(F2) == [[half, -half, half], [half, half, -half],
                       [-half, half, half]]
A2 = Q([[1, 2, 3], [2, 4, 6]])
assert rank(A2) == 1
assert nullspace(transpose(A2)) == Q([[-2, 1]])
assert annihilator(transpose(A2), 2) == Q([[-2, 1]])
assert nullspace(A2) == Q([[-2, 1, 0], [-3, 0, 1]])
third = Fraction(1, 3)
assert annihilator(nullspace(A2), 3) == [[third, 2 * third, Fraction(1)]]
Nsub = Q([[1, 1, 0, 0], [0, 0, 1, 1]])
assert annihilator(Nsub, 4) == Q([[-1, 1, 0, 0], [0, 0, -1, 1]])
print('exercise data: ker A2^t spanned by (-2, 1); (ker A2)^0 by (1,2,3)')
print('all checks passed')
