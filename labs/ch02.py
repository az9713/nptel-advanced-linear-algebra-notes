"""Chapter 2 laboratory: the row-reduced echelon form (RREF).

Exact arithmetic over Q (fractions.Fraction) and over the prime fields
F_p (class GF).  A matrix is a list of rows; no routine changes its
input.  Elementary row operations are tuples:
    ('swap', r, s)       rows r and s are exchanged
    ('scale', r, c)      row r becomes c * row r, with c != 0
    ('add', r, s, c)     row r becomes row r + c * row s, with r != s
Rows and columns are numbered from 0 in code and from 1 in the chapter.
Checks: the running example, existence and uniqueness of the RREF, the
solution sets of AX = 0 and AX = Y, the left-kernel certificate of
inconsistency, and an exhaustive test of the uniqueness proof over F_2
and F_3.
"""
from fractions import Fraction
from itertools import permutations, product
import random


class GF:
    """An element of the prime field F_p = {0, 1, ..., p-1}."""

    def __init__(self, value, p):
        self.p = p
        self.v = value % p

    def _lift(self, other):
        return other if isinstance(other, GF) else GF(other, self.p)

    def __add__(self, other):
        return GF(self.v + self._lift(other).v, self.p)

    def __sub__(self, other):
        return GF(self.v - self._lift(other).v, self.p)

    def __mul__(self, other):
        return GF(self.v * self._lift(other).v, self.p)

    def __truediv__(self, other):
        other = self._lift(other)
        assert other.v != 0, 'division by zero in F_p'
        return GF(self.v * pow(other.v, self.p - 2, self.p), self.p)

    def __rtruediv__(self, other):
        return self._lift(other) / self

    def __neg__(self):
        return GF(-self.v, self.p)

    def __eq__(self, other):
        return self.v == self._lift(other).v

    def __hash__(self):
        return hash((self.v, self.p))

    def __repr__(self):
        return str(self.v)


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def Fp(rows, p):
    return [[GF(x, p) for x in row] for row in rows]


def one_of(A):
    return A[0][0] * 0 + 1


def identity(m, one):
    zero = one - one
    return [[one if i == j else zero for j in range(m)] for i in range(m)]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(1, len(B))),
                 A[i][0] * B[0][j])
             for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, X):
    return [sum((a * x for a, x in zip(row[1:], X[1:])), row[0] * X[0])
            for row in A]


def apply_op(op, A):
    B = [row[:] for row in A]
    if op[0] == 'swap':
        _, r, s = op
        B[r], B[s] = A[s][:], A[r][:]
    elif op[0] == 'scale':
        _, r, c = op
        assert c != 0
        B[r] = [c * x for x in A[r]]
    else:
        _, r, s, c = op
        assert r != s
        B[r] = [x + c * y for x, y in zip(A[r], A[s])]
    return B


def fmt_op(op):
    if op[0] == 'swap':
        return f'R{op[1] + 1} <-> R{op[2] + 1}'
    if op[0] == 'scale':
        return f'R{op[1] + 1} <- ({op[2]}) R{op[1] + 1}'
    return f'R{op[1] + 1} <- R{op[1] + 1} + ({op[3]}) R{op[2] + 1}'


def show(M):
    return '[' + '; '.join(' '.join(str(x) for x in row)
                           for row in M) + ']'


def rref(A):
    """Gauss-Jordan elimination.  Returns (R, pivots, ops).

    Columns are scanned left to right.  In column j, a nonzero entry at
    or below the current row t is swapped up, scaled to 1, and cleared
    from EVERY other row.  Then t grows by one.
    """
    R = [row[:] for row in A]
    m, n = len(R), len(R[0])
    pivots, ops, t = [], [], 0
    for j in range(n):
        if t == m:
            break
        s = next((i for i in range(t, m) if R[i][j] != 0), None)
        if s is None:
            continue
        if s != t:
            ops.append(('swap', t, s))
            R = apply_op(ops[-1], R)
        if R[t][j] != 1:
            ops.append(('scale', t, 1 / R[t][j]))
            R = apply_op(ops[-1], R)
        for i in range(m):
            if i != t and R[i][j] != 0:
                ops.append(('add', i, t, -R[i][j]))
                R = apply_op(ops[-1], R)
        pivots.append(j)
        t += 1
    return R, pivots, ops


def left_factor(ops, m, one):
    """P = E_k ... E_1 with PA = R: apply the operations to I_m."""
    P = identity(m, one)
    for op in ops:
        P = apply_op(op, P)
    return P


def lead(row):
    return next((k for k, x in enumerate(row) if x != 0), None)


def is_row_reduced(A):
    for t, row in enumerate(A):
        k = lead(row)
        if k is None:
            continue
        if row[k] != 1:
            return False
        if any(A[i][k] != 0 for i in range(len(A)) if i != t):
            return False
    return True


def is_rref(A):
    """Row-reduced, zero rows last, leading columns strictly increase."""
    if not is_row_reduced(A):
        return False
    leads = [lead(row) for row in A]
    nz = [k for k in leads if k is not None]
    return (leads[:len(nz)] == nz and sorted(set(nz)) == nz)


def kernel_basis(R, pivots):
    """One vector per free column j: x_j = 1, other free x = 0,
    x_{k_i} = -R[i][j].  This is the parametrization of Section 5."""
    n = len(R[0])
    zero = R[0][0] - R[0][0]
    free = [j for j in range(n) if j not in pivots]
    basis = []
    for j in free:
        X = [zero] * n
        X[j] = one_of(R)
        for i, k in enumerate(pivots):
            X[k] = -R[i][j]
        basis.append(X)
    return basis


def solve(A, Y):
    """Return None if AX = Y is inconsistent, else (X0, kernel basis)."""
    R, pivots, ops = rref(A)
    P = left_factor(ops, len(A), one_of(A))
    Z = matvec(P, Y)
    r = len(pivots)
    if any(z != 0 for z in Z[r:]):
        return None
    X0 = [R[0][0] - R[0][0]] * len(A[0])
    for i, k in enumerate(pivots):
        X0[k] = Z[i]
    return X0, kernel_basis(R, pivots)


def witness(A, Y):
    """A row vector w with wA = 0 and wY != 0, or None if consistent."""
    R, pivots, ops = rref(A)
    P = left_factor(ops, len(A), one_of(A))
    for w in P[len(pivots):]:
        if sum((a * y for a, y in zip(w[1:], Y[1:])), w[0] * Y[0]) != 0:
            return w
    return None


def one_based(cols):
    return tuple(k + 1 for k in cols)


# ---------------------------------------------------------------- 1
# Running example (an original construction), 4 x 5 over Q.
A = Q([[0, 0, 0, 0, 1],
       [1, 2, 1, 7, 0],
       [2, 4, 1, 10, 1],
       [3, 6, 2, 17, 1]])
R, pivots, ops = rref(A)
P = left_factor(ops, 4, Fraction(1))
R_expected = Q([[1, 2, 0, 3, 0], [0, 0, 1, 4, 0],
                [0, 0, 0, 0, 1], [0, 0, 0, 0, 0]])
assert R == R_expected and is_rref(R)
assert matmul(P, A) == R
assert pivots == [0, 2, 4]
print('operations:', '; '.join(fmt_op(op) for op in ops))
print('R =', show(R))
print('P =', show(P))
print('pivot columns k =', one_based(pivots), ' r =', len(pivots))

# ---------------------------------------------------------------- 2
# The row-reduced matrix reached without any interchange is NOT an
# echelon matrix: its leading columns are (5, 1, 3).
R0 = Q([[0, 0, 0, 0, 1], [1, 2, 0, 3, 0],
        [0, 0, 1, 4, 0], [0, 0, 0, 0, 0]])
assert is_row_reduced(R0) and not is_rref(R0)
assert [lead(row) for row in R0] == [4, 0, 2, None]
assert rref(R0)[0] == R and rref(A)[0] == rref(R0)[0]
R1 = apply_op(('swap', 0, 1), R0)
R2 = apply_op(('swap', 1, 2), R1)
assert not is_rref(R1) and R2 == R
# Non-examples of the definition.
assert not is_row_reduced(Q([[1, 2], [0, 1]]))
assert is_row_reduced(Q([[1, 2], [0, 0]]))
assert not is_rref(Q([[0, 0, 0], [1, 0, 0]]))
assert is_rref(Q([[0, 0], [0, 0]])) and is_rref(identity(3, Fraction(1)))
print('R0 leading columns (5, 1, 3): row-reduced, not echelon;',
      'two swaps give R')

# ---------------------------------------------------------------- 3
# Homogeneous solutions and counting (Section 5).
basis = kernel_basis(R, pivots)
assert basis == Q([[-2, 1, 0, 0, 0], [-3, 0, -4, 1, 0]])
Z4 = [Fraction(0)] * 4
for a, b in product(range(-3, 4), repeat=2):
    X = [a * u + b * v for u, v in zip(*basis)]
    assert matvec(A, X) == Z4 and matvec(R0, X) == Z4
print('ker A = x2*(-2,1,0,0,0) + x4*(-3,0,-4,1,0); free columns (2, 4)')
assert len(A[0]) - len(pivots) == 2

# ---------------------------------------------------------------- 4
# Non-homogeneous systems (Section 7).
Y = Q([[1, 2, 5, 7]])[0]
Yb = Q([[1, 2, 5, 8]])[0]
Z = matvec(P, Y)
assert Z == Q([[2, 0, 1, 0]])[0]
sol = solve(A, Y)
assert sol is not None and sol[0] == Q([[2, 0, 0, 0, 1]])[0]
assert matvec(A, sol[0]) == Y
assert solve(A, Yb) is None
w = witness(A, Yb)
assert w == Q([[0, -1, -1, 1]])[0] and witness(A, Y) is None
wA = [sum((w[i] * A[i][j] for i in range(1, 4)), w[0] * A[0][j])
      for j in range(5)]
assert wA == [0] * 5
assert sum(wi * yi for wi, yi in zip(w, Yb)) == 1
print('Y = (1,2,5,7): Z = PY =', [str(z) for z in Z],
      '; X0 =', [str(x) for x in sol[0]])
print('Y = (1,2,5,8): inconsistent; w =', [str(x) for x in w],
      ' wA = 0, wY =', sum(wi * yi for wi, yi in zip(w, Yb)))
# Last r..m rows of P: every w with wA = 0 is a combination of them.
assert all(sum((a * b for a, b in zip(P[3], col)), 0) == 0
           for col in zip(*A))

# Three lines in the plane, no two parallel, no common point.
L = Q([[1, 0], [0, 1], [1, 1]])
YL = Q([[0, 0, 1]])[0]
assert solve(L, YL) is None
wL = witness(L, YL)
assert wL is not None
assert [sum(wL[i] * L[i][j] for i in range(3)) for j in range(2)] == [0, 0]
print('lines x1=0, x2=0, x1+x2=1: witness w =', [str(x) for x in wL],
      ', wY =', sum(a * b for a, b in zip(wL, YL)))
# The two parallel lines of the text: x1+2x2 = 1 and x1+2x2 = 4.
assert solve(Q([[1, 2], [1, 2]]), Q([[1, 4]])[0]) is None

# ---------------------------------------------------------------- 5
# Uniqueness: all row-reduced matrices row-equivalent to A.
rows = [tuple(row) for row in R]
perms = {tuple(p) for p in permutations(rows)}
assert len(perms) == 24
for p in perms:
    M = [list(row) for row in p]
    assert is_row_reduced(M) and rref(M)[0] == R
assert sum(1 for p in perms if is_rref([list(r) for r in p])) == 1
print('row-reduced matrices that are row-equivalent to A (row',
      'orders of R):', len(perms), '; of these in echelon form: 1')

# Random elementary operations never change the RREF.
rng = random.Random(2)
for _ in range(60):
    M = [row[:] for row in A]
    for _ in range(25):
        kind = rng.choice(['swap', 'scale', 'add'])
        r, s = rng.sample(range(4), 2)
        c = Fraction(rng.choice([-3, -2, -1, 1, 2, 3]), rng.choice([1, 2]))
        M = apply_op((kind, r, s) if kind == 'swap' else
                     ('scale', r, c) if kind == 'scale' else
                     ('add', r, s, c), M)
    assert rref(M)[0] == R
print('60 random chains of 25 elementary operations: same RREF')

# A staircase matrix without the "clear above" rule is not unique.
S1, S2 = Q([[1, 2], [0, 1]]), identity(2, Fraction(1))
assert rref(S1)[0] == S2 and rref(S2)[0] == S2 and S1 != S2


# ---------------------------------------------------------------- 6
# Exhaustive test of the uniqueness proof over F_p.
def kernel_bruteforce(A, p):
    n = len(A[0])
    return [x for x in product(range(p), repeat=n)
            if all(sum(a * b for a, b in zip(row, x)) % p == 0
                   for row in A)]


def rref_from_kernel(kernel, m, n, p):
    """Rebuild the RREF from the solution set ALONE (proof of Theorem 8.1).

    Column j is free iff some solution has its last nonzero coordinate
    in position j.  For a free column j, X(j) is the solution with
    x_j = 1 and every other free coordinate 0; then R[i][j] = -X(j)[k_i].
    """
    free = sorted({max(k for k, c in enumerate(x) if c)
                   for x in kernel if any(x)})
    piv = [k for k in range(n) if k not in free]
    R = [[0] * n for _ in range(m)]
    for i, k in enumerate(piv):
        R[i][k] = 1
    for j in free:
        xs = [x for x in kernel if x[j] == 1 and
              all(x[q] == 0 for q in free if q != j)]
        assert len(xs) == 1
        for i, k in enumerate(piv):
            R[i][j] = (-xs[0][k]) % p
    return R


def rref_ints(A, p):
    Rp = rref(Fp(A, p))[0]
    return [[x.v for x in row] for row in Rp]


for p, m, n, expect in ((2, 2, 2, 5), (2, 2, 3, 15), (3, 2, 2, 6),
                        (3, 2, 3, 27)):
    classes = {}
    for flat in product(range(p), repeat=m * n):
        M = [list(flat[i * n:(i + 1) * n]) for i in range(m)]
        ker = kernel_bruteforce(M, p)
        Rm = rref_ints(M, p)
        assert rref_from_kernel(ker, m, n, p) == Rm
        classes.setdefault(frozenset(ker), set()).add(
            tuple(map(tuple, Rm)))
    assert all(len(v) == 1 for v in classes.values())
    assert len(classes) == expect
    print(f'{m}x{n} matrices over F_{p}: {p ** (m * n)} matrices,',
          f'{len(classes)} distinct RREFs, same kernel => same RREF')

# The running example reduced modulo 3: 3^(5-3) = 9 solutions of AX = 0.
A3 = [[int(x) for x in row] for row in A]
ker3 = kernel_bruteforce(A3, 3)
assert len(ker3) == 9
assert rref_ints(A3, 3) == [[int(x) % 3 for x in row] for row in R]
assert rref_from_kernel(ker3, 4, 5, 3) == rref_ints(A3, 3)
print('A over F_3: |ker A| =', len(ker3), '= 3^(n-r)')

# ---------------------------------------------------------------- 7
# Field dependence of the RREF.
M2 = [[1, 1], [1, -1]]
assert rref(Q(M2))[0] == identity(2, Fraction(1))
assert rref(Fp(M2, 3))[0] == identity(2, GF(1, 3))
R2f = rref(Fp(M2, 2))[0]
assert [[x.v for x in row] for row in R2f] == [[1, 1], [0, 0]]
M3 = [[1, 2], [2, 1]]
assert rref(Q(M3))[0] == identity(2, Fraction(1))
R3f = rref(Fp(M3, 3))[0]
assert [[x.v for x in row] for row in R3f] == [[1, 2], [0, 0]]
print('RREF of (1 1; 1 -1): I over Q and F_3, (1 1; 0 0) over F_2')
print('RREF of (1 2; 2 1): I over Q, (1 2; 0 0) over F_3')
# Q -> C: the same elimination steps are valid, so the RREF is the same.
Cm = [[complex(x) for x in row] for row in A]
Rc = rref([[x for x in row] for row in Cm])[0]
assert [[complex(x) for x in row] for row in R] == Rc

# Counting with m < n (needs infinitely many scalars to be infinite).
assert len(kernel_bruteforce([[1, 1]], 2)) == 2
assert solve(Q([[1, 1, 1], [1, 1, 1]]), Q([[1, 2]])[0]) is None

# ---------------------------------------------------------------- 8
# Exercise checks.
B = Q([[1, 2, -1, 3], [2, 4, 1, 0], [3, 6, 0, 3]])
RB, pivB, opsB = rref(B)
PB = left_factor(opsB, 3, Fraction(1))
assert RB == Q([[1, 2, 0, 1], [0, 0, 1, -2], [0, 0, 0, 0]])
assert matmul(PB, B) == RB and pivB == [0, 2]
YB = Q([[1, 5, 6]])[0]
solB = solve(B, YB)
assert solB[0] == Q([[2, 0, 1, 0]])[0]
assert solve(B, Q([[1, 5, 7]])[0]) is None
print('exercise 2: R =', show(RB), ' P =', show(PB), ' X0 =',
      [str(x) for x in solB[0]])
# Exercise 6: column j of A equals sum_i R[i][j] * (column k_i of A).
for j in range(5):
    col = [A[i][j] for i in range(4)]
    comb = [sum(R[i][j] * A[t][pivots[i]] for i in range(3))
            for t in range(4)]
    assert col == comb
print('column relation: a4 = 3 a1 + 4 a3 holds in A')
# Exercise 1: RREF(A) = I for square A iff A is invertible.
Mi = Q([[2, 1], [5, 3]])
assert rref(Mi)[0] == identity(2, Fraction(1))
Ms = Q([[2, 1], [4, 2]])
assert rref(Ms)[0] != identity(2, Fraction(1))
print('all checks passed')
