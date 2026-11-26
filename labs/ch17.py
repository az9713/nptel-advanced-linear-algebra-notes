"""Chapter 17 laboratory: direct sums and projections.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows; a vector is a list.  The routines implement the chapter model:
independence of subspaces by rank and by the criterion
W_j meets W_1 + ... + W_{j-1} only in 0; the projection E onto R along
N built from two bases; the projections E_j of a direct sum; and the
block-diagonal matrix of an operator with invariant summands.
Subspaces are given by lists of spanning vectors.
"""
from fractions import Fraction
from itertools import product


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fraction(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(col) for col in zip(*A)]


def madd(A, B):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(A, B)]


def msub(A, B):
    return [[x - y for x, y in zip(r, s)] for r, s in zip(A, B)]


def mscale(c, A):
    return [[c * x for x in row] for row in A]


def trace(A):
    return sum((A[i][i] for i in range(len(A))), Fraction(0))


def rref(A):
    """Reduced row echelon form and pivot columns, exact."""
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
    return len(rref(A)[1]) if A and A[0] else 0


def nullspace(A):
    """Basis of {x : A x = 0}."""
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


def inverse(A):
    n = len(A)
    aug = [row + e for row, e in zip(A, identity(n))]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n)), 'matrix is not invertible'
    return [row[n:] for row in R]


def span_dim(vectors):
    """Dimension of the span of a list of vectors."""
    return rank(vectors) if vectors else 0


def independent_by_rank(bases):
    """Criterion (c): the concatenated bases are linearly independent."""
    flat = [v for B in bases for v in B]
    return span_dim(flat) == len(flat)


def independent_by_intersections(spans):
    """Criterion (b): dim(W_j meets W_1+...+W_{j-1}) = 0 for every j.

    Uses dim(U meets S) = dim U + dim S - dim(U + S).
    """
    for j in range(1, len(spans)):
        S = [v for W in spans[:j] for v in W]
        meet = span_dim(spans[j]) + span_dim(S) - span_dim(spans[j] + S)
        if meet != 0:
            return False
    return True


def projection(range_basis, null_basis):
    """E with E = identity on span(range_basis), E = 0 on span(null_basis).

    E = P diag(I_r, 0) P^{-1}, with the vectors as columns of P.
    """
    cols = range_basis + null_basis
    n = len(cols)
    P = transpose(cols)
    D = [[Fraction(int(i == j and i < len(range_basis)))
          for j in range(n)] for i in range(n)]
    return matmul(matmul(P, D), inverse(P))


def apply(A, v):
    return [sum((a * x for a, x in zip(row, v)), Fraction(0)) for row in A]


def column_space_basis(A):
    R, piv = rref(transpose(A))
    return R[:len(piv)]


def show(M):
    return '[' + '; '.join(' '.join(str(x) for x in row)
                           for row in M) + ']'


I3 = identity(3)
Z3 = [[Fraction(0)] * 3 for _ in range(3)]

# 1. Pairwise zero intersection is not independence (three lines in Q^2).
l1, l2, l3 = Q([[1, 0]])[0], Q([[0, 1]])[0], Q([[1, 1]])[0]
for a, b in ((l1, l2), (l1, l3), (l2, l3)):
    assert span_dim([a, b]) == 2          # each pair meets only in 0
assert not independent_by_rank([[l1], [l2], [l3]])
assert not independent_by_intersections([[l1], [l2], [l3]])
assert l1[0] + l2[0] - l3[0] == 0 and l1[1] + l2[1] - l3[1] == 0
print('three lines in Q^2: pairwise zero intersection, not independent')

# 2. The running example: V = Q^3, W1 = span(p1, p2), W2 = span(p3).
p1, p2, p3 = Q([[1, 0, 1]])[0], Q([[0, 1, 1]])[0], Q([[1, 1, 1]])[0]
P = transpose([p1, p2, p3])
Pinv = inverse(P)
assert Pinv == Q([[0, -1, 1], [-1, 0, 1], [1, 1, -1]])
assert independent_by_rank([[p1, p2], [p3]])
assert independent_by_intersections([[p1, p2], [p3]])
E1 = projection([p1, p2], [p3])
E2 = projection([p3], [p1, p2])
assert E2 == Q([[1, 1, -1], [1, 1, -1], [1, 1, -1]])
assert E1 == Q([[0, -1, 1], [-1, 0, 1], [-1, -1, 2]])
assert matmul(E1, E1) == E1 and matmul(E2, E2) == E2
assert matmul(E1, E2) == Z3 and matmul(E2, E1) == Z3
assert madd(E1, E2) == I3
assert rank(E1) == 2 and rank(E2) == 1
assert trace(E1) == 2 and trace(E2) == 1      # trace = rank
assert apply(E2, p3) == p3 and apply(E2, p1) == [0, 0, 0]
print('E1 =', show(E1))
print('E2 =', show(E2))
print('E1 + E2 = I, E1 E2 = E2 E1 = 0, ranks 2 and 1, traces 2 and 1')

# E2 is not an orthogonal projection: it is not symmetric.
assert transpose(E2) != E2
# Null space of E2 is the plane x + y - z = 0 and equals W1.
assert [r[:] for r in nullspace(E2)] != [] and span_dim(
    nullspace(E2) + [p1, p2]) == 2
# Range of E1 = null space of E2 (the complement rule).
assert span_dim(column_space_basis(E1) + [p1, p2]) == 2

# 3. The operator A with invariant summands W1 and W2.
A = Q([[2, 1, 0], [1, 3, -1], [0, 1, 2]])
assert apply(A, p1) == [2 * x for x in p1]
assert apply(A, p2) == [x + 2 * y for x, y in zip(p1, p2)]
assert apply(A, p3) == [3 * x for x in p3]
B = matmul(matmul(Pinv, A), P)
assert B == Q([[2, 1, 0], [0, 2, 0], [0, 0, 3]])
assert matmul(A, E1) == matmul(E1, A) and matmul(A, E2) == matmul(E2, A)
assert madd(matmul(A, E1), matmul(A, E2)) == A     # T = T E1 + T E2
print('P^{-1} A P =', show(B))

# A is not diagonalizable: eigenspaces E_2, E_3 have dimension 1 each.
dim2 = len(nullspace(msub(A, mscale(Fraction(2), I3))))
dim3 = len(nullspace(msub(A, mscale(Fraction(3), I3))))
assert (dim2, dim3) == (1, 1)
# The sum of the eigenspaces is direct but is a proper subspace.
e2 = nullspace(msub(A, mscale(Fraction(2), I3)))
e3 = nullspace(msub(A, mscale(Fraction(3), I3)))
assert independent_by_rank([e2, e3]) and dim2 + dim3 < 3
# The invariant plane W1 is the kernel of (A - 2I)^2.
M2 = msub(A, mscale(Fraction(2), I3))
K = nullspace(matmul(M2, M2))
assert len(K) == 2 and span_dim(K + [p1, p2]) == 2
print('dim E_2 = 1, dim E_3 = 1; ker (A - 2I)^2 = W1 has dimension 2')

# 4. A complement that is not invariant: W3 = span(e1).
e1 = Q([[1, 0, 0]])[0]
assert independent_by_rank([[p1, p2], [e1]])
E3 = projection([e1], [p1, p2])
assert E3 == Q([[1, 1, -1], [0, 0, 0], [0, 0, 0]])
assert matmul(A, E3) != matmul(E3, A)
assert apply(A, e1) == [2, 1, 0]
print('span(e1) is not A-invariant: A e1 = (2, 1, 0); A E != E A')

# 5. The two-dimensional example: range (1,-1), null space (1,2).
E = projection([Q([[1, -1]])[0]], [Q([[1, 2]])[0]])
assert E == Q([[Fraction(2, 3), Fraction(-1, 3)],
               [Fraction(-2, 3), Fraction(1, 3)]])
assert matmul(E, E) == E and rank(E) == 1
assert apply(E, Q([[1, 2]])[0]) == [0, 0]
assert apply(E, Q([[1, -1]])[0]) == Q([[1, -1]])[0]
print('lecture example E =', show(E))

# 6. Matrix of a projection in an adapted basis: diag(I_r, 0).
assert matmul(matmul(Pinv, E1), P) == Q([[1, 0, 0], [0, 1, 0], [0, 0, 0]])

# 7. Symmetric + skew-symmetric in M_3(Q): dimensions 6 and 3.
n = 3


def unit(i, j):
    v = [Fraction(0)] * (n * n)
    v[i * n + j] = Fraction(1)
    return v


sym = [[x + y for x, y in zip(unit(i, j), unit(j, i))]
       for i in range(n) for j in range(i, n)]
sym = [[x / 2 if i == j else x for x in v] for v, (i, j)
       in zip(sym, [(i, j) for i in range(n) for j in range(i, n)])]
skew = [[x - y for x, y in zip(unit(i, j), unit(j, i))]
        for i in range(n) for j in range(i + 1, n)]
assert (len(sym), len(skew)) == (6, 3)
assert independent_by_rank([sym, skew]) and span_dim(sym + skew) == 9
Mx = Q([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
S = mscale(Fraction(1, 2), madd(Mx, transpose(Mx)))
K_ = mscale(Fraction(1, 2), msub(Mx, transpose(Mx)))
assert S == Q([[1, 3, 5], [3, 5, 7], [5, 7, 10]])
assert K_ == Q([[0, -1, -2], [1, 0, -1], [2, 1, 0]])
assert madd(S, K_) == Mx and transpose(S) == S
assert transpose(K_) == mscale(Fraction(-1), K_)
print('M = S + K with S =', show(S), 'and K =', show(K_))

# 8. Exercise 6: projection M -> (M + M^T)/2 on M_2(Q), basis E11,E12,E21,E22.
m = 2
cols = []
for idx in range(m * m):
    Eij = [[Fraction(int(r * m + c == idx)) for c in range(m)]
           for r in range(m)]
    img = mscale(Fraction(1, 2), madd(Eij, transpose(Eij)))
    cols.append([x for row in img for x in row])
TS = transpose(cols)
assert TS == Q([[1, 0, 0, 0],
                [0, Fraction(1, 2), Fraction(1, 2), 0],
                [0, Fraction(1, 2), Fraction(1, 2), 0],
                [0, 0, 0, 1]])
assert matmul(TS, TS) == TS and trace(TS) == 3 and rank(TS) == 3
print('symmetrizer on M_2: trace 3 = rank 3')

# 9. Exercise 2: R^3, range span(1,2,0), null plane x + y + z = 0.
r_, n1, n2 = (Q([[1, 2, 0]])[0], Q([[1, -1, 0]])[0], Q([[0, 1, -1]])[0])
Ex = projection([r_], [n1, n2])
third = Fraction(1, 3)
assert Ex == Q([[third, third, third], [2 * third] * 3, [0, 0, 0]])
assert matmul(Ex, Ex) == Ex and trace(Ex) == 1
print('exercise 2: E =', show(Ex))

# 10. Characteristic 2: three copies of I on F_2^2.  E_j^2 = E_j and
# E1 + E2 + E3 = 3 I = I, yet E_i E_j = I != 0 and V is not their direct
# sum (3 * 2 = 6 > 2 = dim V).
I2 = [[1, 0], [0, 1]]
tot = [[(3 * I2[i][j]) % 2 for j in range(2)] for i in range(2)]
assert tot == I2 and I2 != [[0, 0], [0, 0]]
assert 3 * 2 > 2
# Over Q the same three matrices give 3 I != I.
assert mscale(Fraction(3), identity(2)) != identity(2)
print('over F_2: I + I + I = I with I*I != 0; trace argument needs char 0')

# 11. Exercise 4: every projection on F_2^2; trace = rank mod 2.
vecs = list(product(range(2), repeat=2))
found = []
for e in product(range(2), repeat=4):
    M = [[e[0], e[1]], [e[2], e[3]]]
    M2_ = [[(M[i][0] * M[0][j] + M[i][1] * M[1][j]) % 2 for j in range(2)]
           for i in range(2)]
    if M2_ == M:
        image = {tuple((r[0] * v[0] + r[1] * v[1]) % 2 for r in M)
                 for v in vecs}
        rk = {1: 0, 2: 1, 4: 2}[len(image)]
        tr = (M[0][0] + M[1][1]) % 2
        assert tr == rk % 2
        found.append(rk)
assert sorted(found) == [0] + [1] * 6 + [2]
print('projections on F_2^2: 8 (I has trace 0 but rank 2)')

# Exercise 1: components of (1, 2, 5), and A^10 on the W2 part.
v = Q([[1, 2, 5]])[0]
v2 = apply(E2, v)
assert v2 == [-2, -2, -2] and apply(E1, v) == [3, 4, 7]
for _ in range(10):
    v2 = apply(A, v2)
assert v2 == [-118098] * 3
assert apply(E2, Q([[3, 1, 2]])[0]) == [2, 2, 2]
assert apply(E1, Q([[3, 1, 2]])[0]) == [1, -1, 0]
print('exercise 1: components (3,4,7) and (-2,-2,-2); A^10 -> -118098')

# 12. Exercise 5: p(A) acts block by block; m_A = lcm of block minimals.
# On W1 the block is [[2,1],[0,2]] with minimal polynomial (x-2)^2; on
# W2 it is [3] with minimal polynomial x-3.  So m_A = (x-2)^2 (x-3).
M2x = msub(A, mscale(Fraction(2), I3))
M3x = msub(A, mscale(Fraction(3), I3))
assert matmul(matmul(M2x, M2x), M3x) == Z3
assert matmul(M2x, M3x) != Z3 and matmul(M2x, M2x) != Z3
print('m_A = (x-2)^2 (x-3): annihilates A, no proper divisor does')
print('all checks passed')
