"""Chapter 3 laboratory: vector-space axioms, subspaces, span, sums.

Part A works over the rational field Q with fractions.Fraction (exact
row reduction, null spaces, span membership).  Part B enumerates small
vector spaces over the prime fields F_2 and F_3 and tests statements
exhaustively.  A vector is a tuple; a matrix is a list of rows.
"""
from fractions import Fraction
from itertools import product, combinations


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def rref(M):
    """Return (R, pivots): reduced row echelon form of M over Q."""
    R = [row[:] for row in M]
    pivots, r = [], 0
    for c in range(len(R[0])):
        p = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
        if r == len(R):
            break
    return R, pivots


def rank(M):
    return len(rref(M)[1])


def matvec(A, x):
    return tuple(sum((a * b for a, b in zip(row, x)), Fraction(0))
                 for row in A)


def nullspace(A):
    """Basis of {x : Ax = 0}, one vector per free column."""
    R, piv = rref(A)
    n = len(A[0])
    basis = []
    for f in (c for c in range(n) if c not in piv):
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, c in zip(R, piv):
            x[c] = -row[f]
        basis.append(tuple(x))
    return basis


def solve(A, y):
    """One solution of Ax = y, or None when the system is inconsistent."""
    aug = [list(row) + [Fraction(b)] for row, b in zip(A, y)]
    R, piv = rref(aug)
    n = len(A[0])
    if n in piv:
        return None
    x = [Fraction(0)] * n
    for row, c in zip(R, piv):
        x[c] = row[n]
    return tuple(x)


def columns(A):
    return [tuple(row[j] for row in A) for j in range(len(A[0]))]


def in_span(vectors, y):
    """Is y a linear combination of the vectors?  A linear system."""
    A = [[v[i] for v in vectors] for i in range(len(y))]
    return solve(A, y) is not None


def vec(*xs):
    return tuple(Fraction(x) for x in xs)


def parallel(u, v):
    return rank([list(u), list(v)]) <= 1


# ---- Part A: the running example over Q --------------------------------
A = Q([[1, 1, 0, -1], [0, 0, 1, 2], [1, 1, 1, 1]])
B = Q([[0, 1, 0, -1]])
Y_ok, Y_bad = vec(2, 3, 5), vec(2, 3, 4)


def part_a():
    R, piv = rref(A)
    assert piv == [0, 2] and rank(A) == 2
    print('rank A =', rank(A), '; pivot columns (from 0):', piv)

    N = nullspace(A)
    assert N == [vec(-1, 1, 0, 0), vec(1, 0, -2, 1)]
    assert all(matvec(A, z) == vec(0, 0, 0) for z in N)
    print('basis of N(A):', [tuple(map(str, z)) for z in N])

    x0 = solve(A, Y_ok)
    assert matvec(A, x0) == Y_ok
    assert matvec(A, vec(2, 0, 3, 0)) == Y_ok
    assert solve(A, Y_bad) is None
    assert matvec(A, vec(0, 0, 0, 0)) != Y_ok
    for s, t in product(range(-2, 3), repeat=2):
        z = tuple(s * a + t * b for a, b in zip(*N))
        x = tuple(p + q for p, q in zip(x0, z))
        assert matvec(A, x) == Y_ok
    print('AX=(2,3,5): solvable, x0 =', tuple(map(str, x0)))
    print('AX=(2,3,4): inconsistent;  0 is not a solution of AX=Y')

    cols = columns(A)
    assert len(set(cols)) == 3
    for y in product(range(-2, 4), repeat=3):
        y = vec(*y)
        assert in_span(cols, y) == (y[2] == y[0] + y[1])
    gens = [vec(1, 0, 1), vec(0, 1, 1)]
    assert all(in_span(gens, c) for c in cols)
    assert all(in_span(cols, g) for g in gens)
    print('column space = {y : y3 = y1 + y2}, checked on 216 integer y')

    stacked = A + B
    assert rank(stacked) == 3
    (m,) = nullspace(stacked)
    assert parallel(m, vec(0, 1, -2, 1))
    print('W1 meet W2 = span of', tuple(map(str, m)))

    e = [vec(*(int(i == j) for j in range(4))) for i in range(4)]
    W1b = N
    W2b = [vec(1, 0, 0, 0), vec(0, 1, 0, 1), vec(0, 0, 1, 0)]
    assert all(x[1] == x[3] for x in W2b)
    total = rank([list(v) for v in W1b + W2b])
    assert total == 4
    assert all(in_span(W1b + W2b, ei) for ei in e)
    print('rank of W1 and W2 generators together =', total, '(sum = Q^4)')

    u, w = vec(-1, 1, 0, 0), vec(1, 0, 0, 0)
    s = tuple(a + b for a, b in zip(u, w))

    def in_w1(x):
        return matvec(A, x) == vec(0, 0, 0)

    def in_w2(x):
        return x[1] == x[3]

    assert in_w1(u) and in_w2(w)
    assert not (in_w1(s) or in_w2(s))
    print('union fails: u in W1, w in W2, u+w =', tuple(map(str, s)))

    M = Q([[1, 2, 3], [0, 4, 5], [6, 7, 8]])
    n = 3
    Sy = [[(M[i][j] + M[j][i]) / 2 for j in range(n)] for i in range(n)]
    Sk = [[(M[i][j] - M[j][i]) / 2 for j in range(n)] for i in range(n)]
    assert all(Sy[i][j] == Sy[j][i] for i in range(n) for j in range(n))
    assert all(Sk[i][j] == -Sk[j][i] for i in range(n) for j in range(n))
    assert all(Sy[i][j] + Sk[i][j] == M[i][j]
               for i in range(n) for j in range(n))
    print('Sym + Skew split of a 3x3 matrix over Q: checked')

    U = [vec(1, 2, 3), vec(0, 1, 1)]
    for a, b in product(range(-3, 4), repeat=2):
        x = tuple(a * p + b * q for p, q in zip(*U))
        assert (x[0] - x[1] + x[2] == 0) == (a == 0)
    assert all(in_span(U, vec(0, t, t)) for t in range(-2, 3))
    print('exercise 2: U meet W = span{(0,1,1)}, checked on a grid')


# ---- Part B: finite spaces and exhaustive tests ------------------------
def space(p, n):
    return list(product(range(p), repeat=n))


def make_ops(p):
    def add(u, v):
        return tuple((a + b) % p for a, b in zip(u, v))

    def smul(c, u):
        return tuple((c * a) % p for a in u)

    return add, smul


def failed_axioms(elems, add, smul, p):
    """Names of axioms that fail; closure is tested by membership."""
    E = set(elems)
    bad = set()
    if any(add(u, v) not in E for u in elems for v in elems):
        bad.add('A0')
    if any(smul(c, u) not in E for c in range(p) for u in elems):
        bad.add('S0')
    if bad:
        return bad
    if any(add(add(u, v), w) != add(u, add(v, w))
           for u in elems for v in elems for w in elems):
        bad.add('A1')
    zeros = [e for e in elems if all(add(u, e) == u for u in elems)]
    if not zeros:
        bad.add('A2')
    elif any(not any(add(u, v) == zeros[0] for v in elems)
             for u in elems):
        bad.add('A3')
    if any(add(u, v) != add(v, u) for u in elems for v in elems):
        bad.add('A4')
    sc = range(p)
    if any(smul((a + b) % p, u) != add(smul(a, u), smul(b, u))
           for a in sc for b in sc for u in elems):
        bad.add('S1')
    if any(smul(c, add(u, v)) != add(smul(c, u), smul(c, v))
           for c in sc for u in elems for v in elems):
        bad.add('S2')
    if any(smul((a * b) % p, u) != smul(a, smul(b, u))
           for a in sc for b in sc for u in elems):
        bad.add('S3')
    if any(smul(1, u) != u for u in elems):
        bad.add('S4')
    return bad


def one_step(W, add, smul, p):
    return bool(W) and all(add(smul(c, u), v) in W
                           for c in range(p) for u in W for v in W)


def subspaces(elems, add, smul, p):
    out = []
    for r in range(1, len(elems) + 1):
        for sub in combinations(elems, r):
            W = set(sub)
            if one_step(W, add, smul, p):
                out.append(frozenset(W))
    return out


def part_b():
    for p, n in ((3, 2), (2, 3)):
        add, smul = make_ops(p)
        V = space(p, n)
        assert failed_axioms(V, add, smul, p) == set()
        count, agree = 0, 0
        for r in range(1, len(V) + 1):
            for sub in combinations(V, r):
                W = list(sub)
                by_def = failed_axioms(W, add, smul, p) == set()
                agree += by_def == one_step(set(W), add, smul, p)
                count += 1
        assert agree == count
        subs = subspaces(V, add, smul, p)
        print('F_%d^%d: all axioms hold; one-step test agrees with the'
              ' definition on %d subsets; %d subspaces'
              % (p, n, count, len(subs)))
        assert len(subs) == {(3, 2): 6, (2, 3): 16}[(p, n)]
        pairs = nested = 0
        for W1, W2 in combinations(subs, 2):
            is_sub = one_step(set(W1 | W2), add, smul, p)
            assert is_sub == (W1 <= W2 or W2 <= W1)
            pairs += 1
            nested += is_sub
        print('   union of two subspaces is a subspace iff nested:'
              ' %d pairs, %d nested' % (pairs, nested))

    add, smul = make_ops(2)
    V = space(2, 2)
    proper = [W for W in subspaces(V, add, smul, 2)
              if 1 < len(W) < len(V)]
    assert len(proper) == 3
    assert set().union(*proper) == set(V)
    assert not any(set().union(*pair) == set(V)
                   for pair in combinations(proper, 2))
    print('F_2^2 is the union of 3 proper subspaces, not of 2')

    add3, smul3 = make_ops(3)
    V3 = space(3, 2)
    lines = [W for W in subspaces(V3, add3, smul3, 3) if len(W) == 3]
    assert len(lines) == 4
    assert not any(set().union(*t) == set(V3)
                   for t in combinations(lines, 3))
    assert set().union(*lines) == set(V3)
    print('F_3^2: 4 lines cover, no 3 of them do')

    def zero_smul(c, u):
        return tuple(0 for _ in u)

    add3, _ = make_ops(3)
    bad = failed_axioms(space(3, 2), add3, zero_smul, 3)
    assert bad == {'S4'}
    print('scalar multiplication c*u = 0 fails exactly:', sorted(bad))

    def poly_function(coeffs, p):
        return tuple(sum(c * x ** k for k, c in enumerate(coeffs)) % p
                     for x in range(p))

    assert poly_function((0, 1, 1), 2) == (0, 0)
    V = Q([[x ** k for k in range(3)] for x in (0, 1, 2)])
    assert rank(V) == 3
    print('x^2 + x over F_2 is the zero function; over Q the evaluation'
          ' map on P_2 is injective (Vandermonde rank 3)')

    mats = list(product(range(2), repeat=4))

    def transpose(m):
        return (m[0], m[2], m[1], m[3])

    sym = {m for m in mats if transpose(m) == m}
    skew = {m for m in mats
            if transpose(m) == tuple((-x) % 2 for x in m)}
    total = {tuple((a + b) % 2 for a, b in zip(s, k))
             for s in sym for k in skew}
    assert sym == skew and len(sym) == 8 and len(total) == 8
    assert len(total) < len(mats)
    print('F_2^(2x2): Sym = Skew has 8 elements; Sym + Skew has 8 of 16')


if __name__ == '__main__':
    part_a()
    part_b()
    print('ch03 laboratory: all assertions passed')
