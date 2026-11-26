"""Chapter 9 laboratory: matrices of linear maps, change of basis, similarity.

Exact arithmetic over Q (fractions.Fraction).  A vector of F^k is a tuple,
a matrix is a list of rows, an ordered basis is a list of vectors, and a
linear map is a Python function from tuples to tuples.  Polynomials of
degree at most 3 are stored as coefficient tuples (p0, p1, p2, p3).

  coords(B, v)          the coordinate vector [v]_B
  matrix_of(T, B, C)    the matrix [T]_{C<-B}, column j = [T(B_j)]_C
  map_from_matrix(M, B, C)   the map U with [U]_{C<-B} = M
  change_matrix(B, Bp)  P_{B<-B'}, column j = [B'_j]_B

The checks verify the chapter's theorems on its examples.  They are
evidence for those examples and for random small cases, not proofs.
"""
import random
from fractions import Fraction as Fr


def mat(rows):
    return [[Fr(x) for x in row] for row in rows]


def vec(*xs):
    return tuple(Fr(x) for x in xs)


def identity(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(col) for col in zip(*A)]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, x):
    return tuple(sum(a * b for a, b in zip(row, x)) for row in A)


def matpow(A, k):
    R = identity(len(A))
    for _ in range(k):
        R = matmul(R, A)
    return R


def inverse(A):
    """Gauss-Jordan inverse over Q; ValueError when A is singular."""
    n = len(A)
    M = [list(A[i]) + identity(n)[i] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            raise ValueError('singular matrix')
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def rank(A):
    M = [list(r) for r in A]
    rk, rows, cols = 0, len(M), len(M[0])
    for c in range(cols):
        p = next((r for r in range(rk, rows) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(rk + 1, rows):
            f = M[r][c] / M[rk][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def det(A):
    """Determinant by Laplace expansion along row 0 (small matrices)."""
    if len(A) == 1:
        return A[0][0]
    return sum((-1) ** j * A[0][j] *
               det([r[:j] + r[j + 1:] for r in A[1:]])
               for j in range(len(A)))


def charpoly(A):
    """Coefficients of det(xI - A), highest degree first (Faddeev-LeVerrier).
    """
    n = len(A)
    coeffs, M = [Fr(1)], identity(n)
    for k in range(1, n + 1):
        AM = matmul(A, M)
        c = -trace(AM) / k
        coeffs.append(c)
        M = [[AM[i][j] + (c if i == j else 0) for j in range(n)]
             for i in range(n)]
    return coeffs


def basis_matrix(B):
    """Matrix whose columns are the vectors of the ordered basis B."""
    return transpose([list(b) for b in B])


def coords(B, v):
    """[v]_B: solve sum_j c_j B_j = v (B must be a basis of F^k)."""
    return matvec(inverse(basis_matrix(B)), v)


def combine(B, c):
    """The vector sum_j c_j B_j."""
    return tuple(sum(cj * b[i] for cj, b in zip(c, B))
                 for i in range(len(B[0])))


def matrix_of(T, B, C):
    """[T]_{C<-B}: column j is the coordinate vector [T(B_j)]_C."""
    cols = [coords(C, T(b)) for b in B]
    return [[cols[j][i] for j in range(len(B))] for i in range(len(C))]


def map_from_matrix(M, B, C):
    """The unique U with U(B_j) = sum_i M[i][j] C_i."""
    def U(v):
        return combine(C, matvec(M, coords(B, v)))
    return U


def change_matrix(B, Bp):
    """P_{B<-B'}: column j is [B'_j]_B, so [v]_B = P [v]_B'."""
    return transpose([list(coords(B, b)) for b in Bp])


def grid(k, lo=-2, hi=2):
    """All integer vectors of F^k with entries in lo..hi."""
    out = [()]
    for _ in range(k):
        out = [t + (Fr(x),) for t in out for x in range(lo, hi + 1)]
    return out


def show(name, A):
    print(name + ' =', [[str(x) for x in row] for row in A])


# ---------------------------------------------------------------- Part 1
# Section 1 and 3: T(x) = (x1 - 2 x2, x2 + x3) from F^3 to F^2.
def T1(x):
    return vec(x[0] - 2 * x[1], x[1] + x[2])


E3 = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)]
E2 = [vec(1, 0), vec(0, 1)]
Bs = [vec(1, 1, 1), vec(1, 1, 0), vec(1, 0, 0)]
Cp = [vec(1, 0), vec(1, 1)]

A1 = matrix_of(T1, E3, E2)
assert A1 == mat([[1, -2, 0], [0, 1, 1]])
show('Section 3  [T]_std', A1)
for x in grid(3):
    # Theorem 1.1: [T x]_C = [T]_{C<-B} [x]_B
    assert coords(E2, T1(x)) == matvec(A1, coords(E3, x))

A1s = matrix_of(T1, Bs, E2)
assert A1s == mat([[-1, -1, 1], [2, 1, 0]])
show('Section 3  [T]_{std<-B*}', A1s)
A1ss = matrix_of(T1, Bs, Cp)
assert A1ss == mat([[-3, -2, 1], [2, 1, 0]])
show('Section 3  [T]_{C\'<-B*}', A1ss)

# Section 7: equivalence formula A' = Q^-1 A P.
P1 = change_matrix(E3, Bs)
Q1 = change_matrix(E2, Cp)
assert P1 == mat([[1, 1, 1], [1, 1, 0], [1, 0, 0]])
assert Q1 == mat([[1, 1], [0, 1]])
assert matmul(matmul(inverse(Q1), A1), P1) == A1ss
assert matmul(A1, P1) == A1s

# ---------------------------------------------------------------- Part 2
# Section 2: every matrix is the matrix of exactly one map.
M2 = mat([[1, 0, 2], [3, 1, -1]])
U = map_from_matrix(M2, Bs, Cp)
assert U(Bs[0]) == vec(4, 3)
assert U(Bs[1]) == vec(1, 1)
assert U(Bs[2]) == vec(1, -1)
for x in grid(3):
    assert U(x) == vec(x[0] + 3 * x[2], -x[0] + 2 * x[1] + 2 * x[2])
assert matrix_of(U, Bs, Cp) == M2
print('Section 2  U(x) = (x1 + 3 x3, -x1 + 2 x2 + 2 x3); U(B*_1) =',
      [str(t) for t in U(Bs[0])])

random.seed(9)
for _ in range(30):
    m, n = random.randint(1, 3), random.randint(1, 3)
    M = mat([[random.randint(-3, 3) for _ in range(n)] for _ in range(m)])
    Ei = [tuple(Fr(int(i == j)) for j in range(n)) for i in range(n)]
    Fi = [tuple(Fr(int(i == j)) for j in range(m)) for i in range(m)]
    assert matrix_of(map_from_matrix(M, Ei, Fi), Ei, Fi) == M
print('Section 2  30 random matrices: matrix -> map -> matrix is identity')

# ---------------------------------------------------------------- Part 3
# Section 3: differentiation on polynomials of degree <= 3.


def D(p):
    return vec(p[1], 2 * p[2], 3 * p[3], 0)


Mono = [vec(1, 0, 0, 0), vec(0, 1, 0, 0), vec(0, 0, 1, 0), vec(0, 0, 0, 1)]
Dm = matrix_of(D, Mono, Mono)
assert Dm == mat([[0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3], [0, 0, 0, 0]])
Fact = [vec(1, 0, 0, 0), vec(0, 1, 0, 0), vec(0, 0, Fr(1, 2), 0),
        vec(0, 0, 0, Fr(1, 6))]
Nm = matrix_of(D, Fact, Fact)
assert Nm == mat([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]])
Pd = change_matrix(Mono, Fact)
assert matmul(matmul(inverse(Pd), Dm), Pd) == Nm
assert matpow(Dm, 3) != mat([[0] * 4] * 4)
assert matpow(Dm, 4) == mat([[0] * 4] * 4)
show('Section 3  [D]_(1,x,x^2,x^3)', Dm)
show('Section 3  [D]_(1,x,x^2/2,x^3/6)', Nm)

# Section 3: the matrix depends on the order of the basis.
A9 = mat([[1, 1, 1], [0, 2, 1], [0, 0, 3]])


def T9(x):
    return matvec(A9, x)


Bsw = [E3[1], E3[0], E3[2]]
Asw = matrix_of(T9, Bsw, Bsw)
assert Asw == mat([[2, 0, 1], [1, 1, 1], [0, 0, 3]])
show('Section 3  [T]_(e2,e1,e3)', Asw)
# The identity operator has matrix I only when both bases are equal.
assert matrix_of(lambda v: v, Bsw, E3) == mat([[0, 1, 0], [1, 0, 0],
                                              [0, 0, 1]])

# ---------------------------------------------------------------- Part 4
# Section 4: Psi: L(V,W) -> F^{m x n} is linear and bijective.


def T2(x):
    return vec(x[0] + x[1] + x[2], x[1] + x[2])


A2 = matrix_of(T2, E3, E2)
assert A2 == mat([[1, 1, 1], [0, 1, 1]])
c = Fr(5, 3)
Tc = (lambda x: tuple(c * a + b for a, b in zip(T1(x), T2(x))))
assert matrix_of(Tc, E3, E2) == [[c * a + b for a, b in zip(r1, r2)]
                                 for r1, r2 in zip(A1, A2)]


def unit_map(p, q, m, n):
    """The map with matrix E^{pq}: B_q -> C_p and B_j -> 0 for j != q."""
    def Eq(x):
        return tuple(x[q] if i == p else Fr(0) for i in range(m))
    return Eq


units = [(p, q) for p in range(2) for q in range(3)]
mats = [matrix_of(unit_map(p, q, 2, 3), E3, E2) for p, q in units]
flat = [[x for row in M for x in row] for M in mats]
assert rank(flat) == 6          # six independent maps: dim L(F^3,F^2) = 6
assert [flat[k][k] for k in range(6)] == [1] * 6
print('Section 4  dim L(F^3, F^2) = mn =', len(units))

# ---------------------------------------------------------------- Part 5
# Section 5: composition corresponds to the matrix product.


def S1(y):
    return vec(y[0] + y[1], y[0] - y[1])


AS = matrix_of(S1, E2, E2)
assert AS == mat([[1, 1], [1, -1]])
ST = matrix_of(lambda x: S1(T1(x)), E3, E2)
assert ST == matmul(AS, A1) == mat([[1, -1, 1], [1, -3, -1]])
show('Section 5  [S o T]', ST)
assert matmul(Dm, Dm) == matrix_of(lambda p: D(D(p)), Mono, Mono)

# ---------------------------------------------------------------- Part 6
# Section 6: row-vector convention gives T(x) = x A^T.
Arow = mat([[1, 1, 1], [0, 1, 1]])
for x in grid(3):
    assert T2(x) == tuple(matmul([list(x)], transpose(Arow))[0])
print('Section 6  T(x) = x A^T verified on', len(grid(3)), 'row vectors')

# ---------------------------------------------------------------- Part 7
# Sections 7-9: change of basis and similarity on the running example.
Bp = [vec(1, 0, 0), vec(1, 1, 0), vec(1, 1, 1)]
P9 = change_matrix(E3, Bp)
assert P9 == mat([[1, 1, 1], [0, 1, 1], [0, 0, 1]])
P9i = inverse(P9)
assert P9i == mat([[1, -1, 0], [0, 1, -1], [0, 0, 1]])
D9 = matrix_of(T9, Bp, Bp)
assert D9 == mat([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
assert D9 == matmul(matmul(P9i, A9), P9)            # Theorem 7.1
assert A9 == matmul(matmul(P9, D9), P9i)
assert matrix_of(lambda v: v, Bp, E3) == P9          # [I]_{std<-B'}
show('Section 9  P', P9)
show('Section 9  P^-1 A P', D9)
for n in range(0, 12):
    An = matpow(A9, n)
    assert An == [[1, 2 ** n - 1, 3 ** n - 2 ** n],
                  [0, 2 ** n, 3 ** n - 2 ** n],
                  [0, 0, 3 ** n]]
    assert An == matmul(matmul(P9, matpow(D9, n)), P9i)
show('Section 9  A^10', matpow(A9, 10))

# The converse (Section 8): a basis built from an invertible P.
for _ in range(25):
    while True:
        Pr = mat([[random.randint(-2, 2) for _ in range(3)]
                  for _ in range(3)])
        if det(Pr) != 0:
            break
    Bnew = [combine(E3, [Pr[i][j] for i in range(3)]) for j in range(3)]
    assert rank([list(b) for b in Bnew]) == 3          # a basis
    assert change_matrix(E3, Bnew) == Pr
    assert matrix_of(T9, Bnew, Bnew) == \
        matmul(matmul(inverse(Pr), A9), Pr)
print('Section 8  25 random invertible P: B\' basis, [T]_B\' = P^-1 A P')

# The exercise left open in the lecture: a second worked change of basis.
G = mat([[1, 2, -1], [1, 0, 1], [2, 1, 0]])
Bl = [vec(1, 1, 1), vec(1, 1, 0), vec(1, 0, 0)]
H = matrix_of(lambda x: matvec(G, x), Bl, Bl)
assert H == mat([[3, 3, 2], [-1, -2, -1], [0, 2, 0]])
Pl = change_matrix(E3, Bl)
assert H == matmul(matmul(inverse(Pl), G), Pl)
show('Section 9  H', H)

# ---------------------------------------------------------------- Part 8
# Section 10: invariants, and what similarity does not preserve.
for X, Y in ((A9, D9), (G, H)):
    assert rank(X) == rank(Y) and trace(X) == trace(Y)
    assert det(X) == det(Y) and charpoly(X) == charpoly(Y)
print('Section 10  A vs D: trace', trace(A9), 'det', det(A9),
      'charpoly', [str(x) for x in charpoly(A9)])
print('Section 10  G vs H: trace', trace(G), 'det', det(G),
      'charpoly', [str(x) for x in charpoly(G)])

# Same characteristic polynomial and rank, yet not similar.
I2 = identity(2)
J2 = mat([[1, 1], [0, 1]])
assert charpoly(I2) == charpoly(J2) and rank(I2) == rank(J2)
for _ in range(50):
    while True:
        Pr = mat([[random.randint(-3, 3) for _ in range(2)]
                  for _ in range(2)])
        if det(Pr) != 0:
            break
    assert matmul(matmul(inverse(Pr), I2), Pr) == I2     # conjugates of I
    assert I2 != J2

# Equivalent but not similar: I and diag(1, 2) have the same rank 2.
K2 = mat([[1, 0], [0, 2]])
assert rank(I2) == rank(K2) and trace(I2) != trace(K2)

# Rank normal form: bases with matrix [[I_r, 0], [0, 0]] for T1.
a1, a2, a3 = vec(1, 0, 0), vec(0, 1, 0), vec(2, 1, -1)
assert T1(a3) == vec(0, 0)                      # a3 spans ker T1
Bn = [a1, a2, a3]
Cn = [T1(a1), T1(a2)]
Nrm = matrix_of(T1, Bn, Cn)
assert Nrm == mat([[1, 0, 0], [0, 1, 0]])
assert matmul(matmul(inverse(change_matrix(E2, Cn)), A1),
              change_matrix(E3, Bn)) == Nrm
show('Section 10  rank normal form of [T]', Nrm)

# Field dependence: rotation by 90 degrees.
R = [[0, -1], [1, 0]]
assert charpoly(mat(R)) == mat([[1, 0, 1]])[0]    # x^2 + 1: no root in Q
Pc = [[1, 1], [-1j, 1j]]
Rc = [[complex(x) for x in r] for r in R]


def cmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


Pci = [[0.5, 0.5j], [0.5, -0.5j]]
assert cmul(Pc, Pci) == [[1, 0], [0, 1]]
assert cmul(cmul(Pci, Rc), Pc) == [[1j, 0], [0, -1j]]
print('Section 10  rotation R is similar to diag(i, -i) over C only')

print('ch09 laboratory: all assertions passed')
