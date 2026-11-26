"""Chapter 1 laboratory: linear systems and elementary row operations.

Exact arithmetic over Q (fractions.Fraction) and over the prime fields
F_p (class GF).  A matrix is a list of rows; every routine returns a new
matrix and never changes its input.  An elementary row operation is a
tuple:
    ('scale', r, c)      row r becomes c * row r, with c != 0
    ('add', r, s, c)     row r becomes row r + c * row s, with r != s
    ('swap', r, s)       rows r and s are exchanged
Rows are numbered from 0 in code and from 1 in the chapter.
"""
from fractions import Fraction
from itertools import product


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
        if other.v == 0:
            raise ZeroDivisionError('division by zero in F_p')
        return GF(self.v * pow(other.v, self.p - 2, self.p), self.p)

    def __neg__(self):
        return GF(-self.v, self.p)

    def __eq__(self, other):
        return self.v == self._lift(other).v

    def __hash__(self):
        return hash((self.v, self.p))

    def __repr__(self):
        return str(self.v)


def Q(rows):
    """Matrix over Q from integer or Fraction entries."""
    return [[Fraction(x) for x in row] for row in rows]


def Fp(rows, p):
    """Matrix over F_p from integer entries."""
    return [[GF(x, p) for x in row] for row in rows]


def identity(m, one=Fraction(1)):
    zero = one - one
    return [[one if i == j else zero for j in range(m)] for i in range(m)]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(1, len(B))),
                 A[i][0] * B[0][j])
             for j in range(len(B[0]))] for i in range(len(A))]


def apply_op(op, A):
    """Return e(A) for the elementary row operation op."""
    B = [row[:] for row in A]
    if op[0] == 'scale':
        _, r, c = op
        assert c != 0, 'scaling by zero is not an elementary operation'
        B[r] = [c * x for x in A[r]]
    elif op[0] == 'add':
        _, r, s, c = op
        assert r != s, 'type-2 operation needs two different rows'
        B[r] = [x + c * y for x, y in zip(A[r], A[s])]
    else:
        _, r, s = op
        B[r], B[s] = A[s][:], A[r][:]
    return B


def inverse_op(op):
    """The elementary operation of the same type that undoes op."""
    if op[0] == 'scale':
        return ('scale', op[1], 1 / op[2])
    if op[0] == 'add':
        return ('add', op[1], op[2], -op[3])
    return op


def elementary_matrix(op, m, one=Fraction(1)):
    """E = e(I) for the m x m identity."""
    return apply_op(op, identity(m, one))


def row_reduce(A):
    """Reduce A to a row-reduced matrix without interchanging rows.

    Rows are processed from top to bottom.  A nonzero row is scaled so
    that its first nonzero entry becomes 1; that entry's column is then
    cleared in EVERY other row, above and below.  Returns (R, ops) with
    R = e_k(...e_1(A)) for ops = [e_1, ..., e_k].
    """
    R = [row[:] for row in A]
    ops = []
    for t in range(len(R)):
        lead = next((k for k, x in enumerate(R[t]) if x != 0), None)
        if lead is None:
            continue
        if R[t][lead] != 1:
            op = ('scale', t, 1 / R[t][lead])
            R = apply_op(op, R)
            ops.append(op)
        for i in range(len(R)):
            if i != t and R[i][lead] != 0:
                op = ('add', i, t, -R[i][lead])
                R = apply_op(op, R)
                ops.append(op)
    return R, ops


def is_row_reduced(A):
    """Condition (a): leading entry of each nonzero row is 1.
    Condition (b): a column holding a leading entry is otherwise zero."""
    for t, row in enumerate(A):
        lead = next((k for k, x in enumerate(row) if x != 0), None)
        if lead is None:
            continue
        if row[lead] != 1:
            return False
        if any(A[i][lead] != 0 for i in range(len(A)) if i != t):
            return False
    return True


def is_solution(A, X, Y):
    return all(sum((a * x for a, x in zip(row, X)), y - y) == y
               for row, y in zip(A, Y))


def augment(A, Y):
    return [row + [y] for row, y in zip(A, Y)]


def fmt_op(op):
    """Write an operation in the chapter's 1-based notation."""
    if op[0] == 'scale':
        return f'R{op[1] + 1} <- ({op[2]}) R{op[1] + 1}'
    if op[0] == 'add':
        return f'R{op[1] + 1} <- R{op[1] + 1} + ({op[3]}) R{op[2] + 1}'
    return f'R{op[1] + 1} <-> R{op[2] + 1}'


def show(M):
    return '[' + '; '.join(' '.join(str(x) for x in row)
                           for row in M) + ']'


# Running example of the chapter (an original construction).
A = Q([[0, 0, 2, 2], [2, 4, 1, 7], [2, 4, 3, 9]])
I3 = identity(3)

# 1. The elementary-matrix identity e(A) = e(I) A for every type.
sample_ops = [('scale', 1, Fraction(-3)), ('add', 0, 2, Fraction(5, 2)),
              ('add', 2, 0, Fraction(-1)), ('swap', 0, 2)]
for op in sample_ops:
    assert apply_op(op, A) == matmul(elementary_matrix(op, 3), A)
    assert apply_op(inverse_op(op), apply_op(op, A)) == A
    E, F = elementary_matrix(op, 3), elementary_matrix(inverse_op(op), 3)
    assert matmul(E, F) == I3 and matmul(F, E) == I3
print('e(A) = e(I)A and inverse operations: checked for',
      len(sample_ops), 'operations')

# Boundary: "row r + c * row r" with c = -1 kills a row; not invertible.
bad = [row[:] for row in A]
bad[0] = [x + (-1) * x for x in A[0]]
assert all(x == 0 for x in bad[0])

# 2. Reduce the running example and record P with PA = R.
R, ops = row_reduce(A)
P = I3
for op in ops:
    P = matmul(elementary_matrix(op, 3), P)
assert matmul(P, A) == R and is_row_reduced(R)
assert R == Q([[0, 0, 1, 1], [1, 2, 0, 3], [0, 0, 0, 0]])
print('operations:', '; '.join(fmt_op(op) for op in ops))
print('R =', show(R))
print('P =', show(P))

# R is row-reduced but not echelon: the leading 1 of row 1 lies to the
# right of the leading 1 of row 2.  Swapping gives a second row-reduced
# matrix that is row-equivalent to A (exercise 3).
R_swapped = apply_op(('swap', 0, 1), R)
assert is_row_reduced(R_swapped) and R_swapped != R

# 3. Homogeneous solutions: x3 = -x4, x1 = -2 x2 - 3 x4.
Z3 = [Fraction(0)] * 3
h1, h2 = Q([[-2, 1, 0, 0]])[0], Q([[-3, 0, -1, 1]])[0]
for a, b in product(range(-2, 3), repeat=2):
    X = [a * u + b * v for u, v in zip(h1, h2)]
    assert is_solution(A, X, Z3) and is_solution(R, X, Z3)
print('AX = 0 solved by x2*(-2,1,0,0) + x4*(-3,0,-1,1)')

# 4. Consistency: after the same operations Y becomes PY, and the
# zero row of R demands y3 - y1 - y2 = 0.
for Y in (Q([[2, 9, 11]])[0], Q([[2, 9, 12]])[0]):
    RY, _ = row_reduce(augment(A, Y))
    print('Y =', [str(y) for y in Y], '-> reduced augmented', show(RY))
consistent = row_reduce(augment(A, Q([[2, 9, 11]])[0]))[0]
assert consistent[2] == Q([[0, 0, 0, 0, 0]])[0]
assert row_reduce(augment(A, Q([[2, 9, 12]])[0]))[0][2][4] == 1
X0 = Q([[4, 0, 1, 0]])[0]
assert is_solution(A, X0, Q([[2, 9, 11]])[0])
print('particular solution for Y=(2,9,11): X0 = (4,0,1,0)')

# 5. The opening failure: E1 <- E1+E2 and E2 <- E1+E2 simultaneously.
S = [[x + y for x, y in zip(A[0], A[1])]] * 2 + [A[2][:]]
assert S[0] == A[2]
Xbad = Q([[0, 0, 3, -1]])[0]
assert is_solution(S, Xbad, Z3) and not is_solution(A, Xbad, Z3)
print('simultaneous replacement admits X = (0,0,3,-1); A*X =',
      [str(sum(a * x for a, x in zip(row, Xbad))) for row in A])

# 6. The two-line example: 2x1 + x2 = 4, x1 - 2x2 = -3.
L = Q([[2, 1], [1, -2]])
RL, _ = row_reduce(augment(L, Q([[4, -3]])[0]))
assert RL == Q([[1, 0, 1], [0, 1, 2]])
print('two-line example: unique solution (1, 2)')
Lswap = apply_op(('swap', 0, 1), L)
Rbad, _ = row_reduce(augment(Lswap, Q([[4, -3]])[0]))
assert [r[2] for r in Rbad] == [Fraction(-2, 5), Fraction(-11, 5)]
print('swapping coefficient rows only gives (-2/5, -11/5)')

# 7. Three pairwise non-parallel planes with no common point.
planes = Q([[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 1]])
RP, _ = row_reduce(planes)
assert RP[2] == Q([[0, 0, 0, 1]])[0]
print('planes x1=0, x2=0, x1+x2=1: reduced last row', show([RP[2]]))

# 8. Field dependence: x1 + x2 = 1, x1 - x2 = 0 over Q, F_3, F_2.
for p in (3, 2):
    M = Fp([[1, 1, 1], [1, -1, 0]], p)
    sols = [(a, b) for a, b in product(range(p), repeat=2)
            if is_solution([r[:2] for r in M], [GF(a, p), GF(b, p)],
                           [r[2] for r in M])]
    print(f'over F_{p}: solutions {sols}')
    assert sols == ([(2, 2)] if p == 3 else [])
RQ, _ = row_reduce(augment(Q([[1, 1], [1, -1]]), Q([[1, 0]])[0]))
assert RQ == Q([[1, 0, Fraction(1, 2)], [0, 1, Fraction(1, 2)]])
R2, _ = row_reduce(Fp([[1, 1, 1], [1, -1, 0]], 2))
assert [x.v for x in R2[1]] == [0, 0, 1]
print('over Q: (1/2, 1/2); over F_2 the reduced row is 0 = 1')

# 9. Two distinct solutions force |F| solutions; over F_2 that is two.
two = [(a, b) for a, b in product(range(2), repeat=2)
       if (a + b) % 2 == 1]
assert two == [(0, 1), (1, 0)]
print('x1 + x2 = 1 over F_2 has exactly', len(two), 'solutions')

# 10. Exercise 2 matrix: the "clear below only" rule fails here.
M = Q([[1, 2, 1, 0], [2, 4, 3, 1], [0, 1, 1, 2]])
RM, opsM = row_reduce(M)
PM = I3
for op in opsM:
    PM = matmul(elementary_matrix(op, 3), PM)
assert matmul(PM, M) == RM == Q([[1, 0, 0, -3], [0, 0, 1, 1],
                                 [0, 1, 0, 1]])
assert is_solution(M, Q([[3, -1, -1, 1]])[0], Z3)
print('exercise 2 operations:', '; '.join(fmt_op(op) for op in opsM))
print('exercise 2 R =', show(RM), ' P =', show(PM))


def reduce_below_only(A):
    """The faulty variant: clear each leading column only below."""
    B = [row[:] for row in A]
    for t in range(len(B)):
        lead = next((k for k, x in enumerate(B[t]) if x != 0), None)
        if lead is None:
            continue
        B = apply_op(('scale', t, 1 / B[t][lead]), B)
        for i in range(t + 1, len(B)):
            B = apply_op(('add', i, t, -B[i][lead]), B)
    return B


below_only = reduce_below_only(M)
assert below_only == Q([[1, 2, 1, 0], [0, 0, 1, 1], [0, 1, 0, 1]])
assert not is_row_reduced(below_only)
print('clearing below only gives', show(below_only), '(not row-reduced)')

# 11. Column dependencies are solutions of AX = 0 and survive reduction.
for c in (Q([[2, -1, 0, 0]])[0], Q([[3, 0, 1, -1]])[0]):
    assert is_solution(A, c, Z3) and is_solution(R, c, Z3)
print('column relations a2 = 2 a1 and a4 = 3 a1 + a3 hold in A and R')

# 12. Exercise 3: the six row orders of R are distinct, row-reduced,
# and row-equivalent to A.
orders = {tuple(map(tuple, (R[i] for i in perm)))
          for perm in [(0, 1, 2), (0, 2, 1), (1, 0, 2),
                       (1, 2, 0), (2, 0, 1), (2, 1, 0)]}
assert len(orders) == 6
assert all(is_row_reduced([list(r) for r in m]) for m in orders)
print('row-reduced matrices row-equivalent to A found:', len(orders))

# 13. Exercise 5: a column operation A -> AQ changes the solutions by
# X = Q X'.  Swap the columns of the two-line example.
Qswap = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
LQ = matmul(L, Qswap)
RLQ = row_reduce(augment(LQ, Q([[4, -3]])[0]))[0]
assert RLQ == Q([[1, 0, 2], [0, 1, 1]])
print('after swapping columns the solution is (2, 1); Q(2, 1) = (1, 2)')
print('all checks passed')
