"""Chapter 8 laboratory: the algebra of linear maps.

Exact arithmetic over Q (fractions.Fraction).  A linear map F^n -> F^m is
stored as its m x n matrix (a list of rows) acting on column vectors.
Everything is checked on the chapter's running example:

    A = [[1, 1, 0], [0, 1, 1]]     T : Q^3 -> Q^2
    B = [[2, 0, 1], [1, -1, 0]]    U : Q^3 -> Q^2
    C = [[1, 2], [1, 3]]           R : Q^2 -> Q^2 (invertible)
    N = [[0, 1], [0, 0]]           nilpotent operator on Q^2
"""
from fractions import Fraction
from itertools import product


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def show(M):
    """Compact text for a matrix or vector of Fractions."""
    if M and isinstance(M[0], list):
        return '[' + '; '.join(show(r) for r in M) + ']'
    return '(' + ', '.join(str(x) for x in M) + ')'


def shape(M):
    return len(M), len(M[0])


def add(M, N):
    assert shape(M) == shape(N), 'sum needs maps with one domain and codomain'
    return [[a + b for a, b in zip(r, s)] for r, s in zip(M, N)]


def scale(c, M):
    return [[c * a for a in row] for row in M]


def matmul(M, N):
    """Matrix of the composite M after N (apply N first)."""
    assert shape(N)[0] == shape(M)[1], 'composite is not defined'
    return [[sum(M[i][k] * N[k][j] for k in range(len(N)))
             for j in range(len(N[0]))] for i in range(len(M))]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def apply(M, x):
    return [sum(a * b for a, b in zip(row, x)) for row in M]


def rref(M):
    """Reduced row echelon form and pivot columns."""
    R = [row[:] for row in M]
    pivots, r = [], 0
    for c in range(len(R[0])):
        p = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = [a / R[r][c] for a in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                R[i] = [a - R[i][c] * b for a, b in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
        if r == len(R):
            break
    return R, pivots


def rank(M):
    return len(rref(M)[1])


def nullspace(M):
    """A basis of {x : Mx = 0} from the reduced form."""
    R, piv = rref(M)
    n = len(M[0])
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for i, c in enumerate(piv):
            v[c] = -R[i][f]
        basis.append(v)
    return basis


def det(M):
    """Determinant by Laplace expansion (small matrices only)."""
    n = len(M)
    if n == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * det([r[:j] + r[j + 1:] for r in M[1:]])
               for j in range(n))


def inverse(M):
    n = len(M)
    aug = [row + identity(n)[i] for i, row in enumerate(M)]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n)), 'matrix is not invertible'
    return [row[n:] for row in R]


def unit(p, q, m, n):
    """The matrix e_pq (1 at row p, column q), 1-based as in the text."""
    return [[Fraction(int(i == p - 1 and j == q - 1)) for j in range(n)]
            for i in range(m)]


def flat(M):
    return [x for row in M for x in row]


A = Q([[1, 1, 0], [0, 1, 1]])
B = Q([[2, 0, 1], [1, -1, 0]])
C = Q([[1, 2], [1, 3]])
N = Q([[0, 1], [0, 0]])
I2, I3 = identity(2), identity(3)

# 1. Sums and scalar multiples act pointwise.
S = add(A, B)
assert S == Q([[3, 1, 1], [1, 0, 1]])
assert scale(2, A) == Q([[2, 2, 0], [0, 2, 2]])
for x in product(range(-2, 3), repeat=3):
    assert apply(S, x) == [a + b for a, b in zip(apply(A, x), apply(B, x))]
print('T+U has matrix', show(S), ' 2T has matrix', show(scale(2, A)))

# 2. The set of invertible operators is not a subspace; I + (-I) = 0.
assert add(I2, scale(-1, I2)) == [[0, 0], [0, 0]]
assert rank(add(I2, scale(-1, I2))) == 0

# 3. dim L(V, W) = mn: the six units e_pq are independent and span.
units = [unit(p, q, 2, 3) for p in (1, 2) for q in (1, 2, 3)]
assert rank([flat(E) for E in units]) == 6
for M in (A, B, S):
    rebuilt = [[Fraction(0)] * 3 for _ in range(2)]
    for p in (1, 2):
        for q in (1, 2, 3):
            rebuilt = add(rebuilt, scale(M[p - 1][q - 1], unit(p, q, 2, 3)))
    assert rebuilt == M
# E^{pq}(alpha_j) = delta_{jq} beta_p for the standard bases.
for p, q in product((1, 2), (1, 2, 3)):
    for j in (1, 2, 3):
        col = apply(unit(p, q, 2, 3), I3[j - 1])
        assert col == [Fraction(int(j == q and i == p)) for i in (1, 2)]
print('dim L(Q^3, Q^2) = 6; coefficients of T:', show(flat(A)),
      'coefficients of U:', show(flat(B)))
n = 3
assert rank([flat(unit(p, q, n, n)) for p in range(1, n + 1)
             for q in range(1, n + 1)]) == n * n == 9

# 4. Other bases change the E^{pq}: beta' = ((1,0), (1,1)) in W.
# E'^{21} sends alpha_1 to beta'_2 = (1,1), not to (0,1).
P = Q([[1, 1], [0, 1]])
E21p = matmul(P, unit(2, 1, 2, 3))
assert apply(E21p, I3[0]) == Q([[1, 1]])[0]
assert apply(E21p, I3[1]) == Q([[0, 0]])[0]

# 5. Composition: UT is defined for T: V -> W, U: W -> Z only.
CA = matmul(C, A)
assert CA == Q([[1, 3, 2], [1, 4, 3]])
for x in product(range(-2, 3), repeat=3):
    assert apply(CA, x) == apply(C, apply(A, x))
try:
    matmul(A, C)
    raise SystemExit('A after C must be undefined')
except AssertionError:
    pass
assert rank(CA) == 2 <= min(rank(C), rank(A))
print('R after T has matrix', show(CA), ' rank', rank(CA))
CB = matmul(C, B)
assert matmul(C, add(A, B)) == add(CA, CB) == Q([[5, 1, 3], [6, 1, 4]])
assert matmul(add(A, B), I3) == add(A, B)
# Right distributivity (U1+U2)T = U1T + U2T with a second map into Z.
D2 = Q([[1, 0], [2, 5]])
assert matmul(add(C, D2), A) == add(matmul(C, A), matmul(D2, A))

# 6. Operators: non-commutativity and zero divisors in L(V).
CN, NC = matmul(C, N), matmul(N, C)
assert CN == Q([[0, 1], [0, 1]]) and NC == Q([[1, 3], [0, 0]]) and CN != NC
assert matmul(N, N) == [[0, 0], [0, 0]] and N != [[0, 0], [0, 0]]
e11, e12 = unit(1, 1, 2, 2), unit(1, 2, 2, 2)
assert matmul(e11, e12) == e12
assert matmul(e12, e11) == [[0, 0], [0, 0]]
print('CN =', show(CN), ' NC =', show(NC), ' N^2 = 0')

# 7. L(V) is an algebra with identity (spot check with fixed matrices).
X = Q([[1, -1], [2, 3]])
Y = Q([[0, 4], [1, 1]])
Z = Q([[2, 0], [-1, 5]])
c = Fraction(7, 3)
assert matmul(X, matmul(Y, Z)) == matmul(matmul(X, Y), Z)
assert matmul(X, add(Y, Z)) == add(matmul(X, Y), matmul(X, Z))
assert matmul(add(Y, Z), X) == add(matmul(Y, X), matmul(Z, X))
assert scale(c, matmul(X, Y)) == matmul(scale(c, X), Y)
assert scale(c, matmul(X, Y)) == matmul(X, scale(c, Y))
assert matmul(X, I2) == X == matmul(I2, X)

# 8. n^2 + 1 powers of an operator are dependent; for C one finds
# C^2 = 4C - I, hence C^{-1} = 4I - C.
C2 = matmul(C, C)
assert C2 == add(scale(4, C), scale(-1, I2)) == Q([[3, 8], [4, 11]])
Cinv = add(scale(4, I2), scale(-1, C))
assert Cinv == Q([[3, -2], [-1, 1]]) == inverse(C)
assert matmul(C, Cinv) == I2 == matmul(Cinv, C)
assert det(C) == 1
print('C^2 = 4C - I;  C^{-1} =', show(Cinv))
powers = [flat(identity(2))]
M = identity(2)
for _ in range(4):
    M = matmul(M, C)
    powers.append(flat(M))
assert rank(powers) == 2 < 5

# 9. One-sided inverses for T = L_A: right inverse exists, left does not.
G = scale(Fraction(1, 3), Q([[2, -1], [1, 1], [-1, 2]]))
assert matmul(A, G) == I2
GT = matmul(G, A)
assert GT == scale(Fraction(1, 3), Q([[2, 1, -1], [1, 2, 1], [-1, 1, 2]]))
assert GT != I3 and matmul(GT, GT) == GT and rank(GT) == 2
ker = nullspace(A)
assert ker == [Q([[1, -1, 1]])[0]] and apply(A, ker[0]) == [0, 0]
print('T G = I_2 but G T has rank', rank(GT), '; ker T spanned by',
      show(ker[0]))
print('3 x 3 G T =', show(GT))

# 10. Nonsingular maps preserve independence; singular ones may not.
assert nullspace(C) == []
assert rank([apply(C, I2[0]), apply(C, I2[1])]) == 2
assert rank([apply(A, ker[0])]) == 0          # {(1,-1,1)} -> {0}
assert rank(Q([[1, -1, 1]])) == 1             # the source set is independent

# 11. Coordinates give an isomorphism P_2(Q) -> Q^3.
# p = a + b x + c x^2 is the triple (a, b, c).
S3 = Q([[1, 1, 0], [0, 1, 1], [1, 0, 1]])     # 1+x, x+x^2, 1+x^2
assert det(S3) == 2 and rank(S3) == 3
for p in (2, 3, 5):
    assert (int(det(S3)) % p == 0) == (p == 2)
rows_mod2 = [[int(x) % 2 for x in row] for row in S3]
assert [sum(col) % 2 for col in zip(*rows_mod2)] == [0, 0, 0]
print('{1+x, x+x^2, 1+x^2}: det over Q =', det(S3),
      '; over F_2 the sum of the three is 0')

# 12. Infinite dimension: Q[x] with D (derivative) and J (integral).
# A polynomial is a list of coefficients [a0, a1, ...].


def D(p):
    return [k * a for k, a in enumerate(p)][1:] or [Fraction(0)]


def J(p):
    return [Fraction(0)] + [a / (k + 1) for k, a in enumerate(p)]


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


p = Q([[3, 2, 1]])[0]
assert trim(D(J(p))) == p                      # D J = I
assert trim(J(D(p))) == Q([[0, 2, 1]])[0]      # J D p = p - p(0)
assert trim(J(D(p))) != p                      # J D != I
assert trim(J(Q([[1]])[0])) == Q([[0, 1]])[0]  # J is injective on samples
print('D J = I but J D p = p - p(0); J injective, not onto Q[x]')

# 13. Counting maps and isomorphisms over F_q: |L| = q^(mn) and
# |GL_2(F_q)| = (q^2 - 1)(q^2 - q).
for q in (2, 3):
    mats = list(product(range(q), repeat=4))
    assert len(mats) == q ** 4
    inv = [m for m in mats if (m[0] * m[3] - m[1] * m[2]) % q != 0]
    assert len(inv) == (q ** 2 - 1) * (q ** 2 - q)
    print(f'F_{q}: {len(mats)} maps F^2 -> F^2, {len(inv)} isomorphisms')
assert len(list(product(range(2), repeat=6))) == 2 ** 6   # L(F_2^3, F_2^2)

# 14. Exercise 2: matrix of 2T - U.
assert add(scale(2, A), scale(-1, B)) == Q([[0, 2, -1], [-1, 3, 2]])
print('all checks passed')
