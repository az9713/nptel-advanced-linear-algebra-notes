"""Chapter 15 laboratory: transposes, invariant subspaces, conductors.

Exact arithmetic over Q (fractions.Fraction).  A vector is a list of
Fractions, a matrix is a list of rows, a polynomial is a list of
coefficients with the constant term first (so [c0, c1, c2] is
c0 + c1 x + c2 x^2).  Nothing is changed in place.

Sections of the laboratory:
  A. the transpose of a map: ker T^t = (im T)^0, im T^t = (ker T)^0,
     rank T^t = rank T, and [T^t] = A^T in dual bases;
  B. invariant subspaces and annihilators: W is T-invariant exactly
     when W^0 is T^t-invariant;
  C. the T-conductor S(v; W), computed by linear algebra;
  D. the lemma (a conductor into an invariant proper subspace is a
     linear polynomial when the minimal polynomial splits) and the
     construction of a triangularizing basis from it;
  E. the boundary: a rotation of Q^2 has no invariant line.
"""
from fractions import Fraction


def Q(rows):
    """Matrix over Q from integer or Fraction entries."""
    return [[Fraction(x) for x in row] for row in rows]


def vec(xs):
    return [Fraction(x) for x in xs]


def matmul(A, B):
    return [[sum(a * b for a, b in zip(row, col)) for col in zip(*B)]
            for row in A]


def matvec(A, v):
    return [sum(a * x for a, x in zip(row, v)) for row in A]


def transpose(A):
    return [list(col) for col in zip(*A)]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def sub(A, B):
    return [[a - b for a, b in zip(r, s)] for r, s in zip(A, B)]


def scal(c, A):
    return [[c * a for a in row] for row in A]


def rref(A):
    """Reduced row echelon form and the list of pivot columns."""
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
    return len(rref(A)[1]) if A else 0


def nullspace(A):
    """Basis of {x : A x = 0}, as a list of vectors."""
    n = len(A[0])
    R, piv = rref(A)
    basis = []
    for f in (c for c in range(n) if c not in piv):
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, p in zip(R, piv):
            x[p] = -row[f]
        basis.append(x)
    return basis


def colspace(A):
    """Basis of the column space: the pivot columns of A."""
    _, piv = rref(A)
    return [[row[c] for row in A] for c in piv]


def in_span(vs, v):
    """Is v in span(vs)?  (vs may be empty.)"""
    if not vs:
        return all(x == 0 for x in v)
    M = transpose(vs)
    return rank(M) == rank([r + [x] for r, x in zip(M, v)])


def indep(vs):
    return rank(vs) == len(vs) if vs else True


def invariant(A, W):
    """Is span(W) invariant under x -> A x?"""
    return all(in_span(W, matvec(A, w)) for w in W)


def annihilator(W, n):
    """Basis of W^0 as row vectors g (g . w = 0 for all w in W)."""
    return nullspace(W) if W else [vec(e) for e in identity(n)]


def row_act(A, g):
    """T^t g in the dual basis: the row vector g A."""
    return [sum(g[i] * A[i][j] for i in range(len(g)))
            for j in range(len(A[0]))]


# ---- polynomials (constant term first) ---------------------------

def ptrim(p):
    p = p[:]
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def pdivmod(p, d):
    p = ptrim(p)
    d = ptrim(d)
    q = [Fraction(0)] * max(len(p) - len(d) + 1, 1)
    while len(p) >= len(d) and p:
        c = p[-1] / d[-1]
        k = len(p) - len(d)
        q[k] = c
        for i, x in enumerate(d):
            p[i + k] -= c * x
        p = ptrim(p)
    return ptrim(q), p


def peval(p, x):
    out = Fraction(0)
    for c in reversed(p):
        out = out * x + c
    return out


def pmat(p, A):
    """The matrix p(A) by Horner's rule."""
    n = len(A)
    out = [[Fraction(0)] * n for _ in range(n)]
    for c in reversed(p):
        out = matmul(out, A)
        out = [[out[i][j] + (c if i == j else 0) for j in range(n)]
               for i in range(n)]
    return out


def pshow(p):
    p = ptrim(p)
    terms = []
    for k in reversed(range(len(p))):
        c = p[k]
        if c == 0:
            continue
        mono = 'x' + ('^%d' % k if k > 1 else '') if k else ''
        coef = '' if (c == 1 and k) else str(c)
        if c == -1 and k:
            coef = '-'
        terms.append((coef + mono))
    s = ' + '.join(terms) if terms else '0'
    return s.replace('+ -', '- ')


def linear_factors(g, candidates):
    """Split g (monic) into linear factors (x - c); None if impossible.

    Returns the list of roots with multiplicity, using the candidate
    roots only.  A None answer means g does not split over the
    candidates, which for rational roots of an integer-coefficient
    polynomial is a proof once the candidates cover the divisors.
    """
    roots = []
    g = ptrim(g)
    while len(g) > 1:
        c = next((c for c in candidates if peval(g, c) == 0), None)
        if c is None:
            return None
        g, rem = pdivmod(g, [-c, Fraction(1)])
        assert not rem
        roots.append(c)
    return roots


# ---- conductors ----------------------------------------------------

def krylov(A, v, k):
    """[v, Av, ..., A^k v]."""
    out = [v]
    for _ in range(k):
        out.append(matvec(A, out[-1]))
    return out


def conductor(A, v, W):
    """Monic generator of S(v; W) = {g : g(A) v in span W}.

    W is a list of vectors spanning a T-invariant subspace (maybe empty).
    Search for the least d such that v, Av, ..., A^d v are dependent
    modulo W.  Then the dependence is unique up to a scalar and its
    coefficient on A^d v is not zero.
    """
    n = len(A)
    for d in range(0, n + 2):
        K = krylov(A, v, d)
        # relation: sum c_i K_i + sum w_j W_j = 0; unknowns c then w.
        cols = K + list(W)
        M = transpose(cols) if cols else [[]]
        ns = nullspace(M)
        rel = [x for x in ns if any(x[:d + 1])]
        if rel:
            assert len(rel) == 1
            c = rel[0][:d + 1]
            assert c[d] != 0
            return [x / c[d] for x in c]
    raise AssertionError('no conductor found (W not invariant?)')


def triangularize(A, roots):
    """Basis alpha_1..alpha_n with T alpha_j in span(alpha_1..alpha_j).

    The construction of the chapter: with W = span(alpha_1..alpha_{j-1})
    pick beta not in W, factor its conductor g = (x - c) h, and put
    alpha_j = h(T) beta.  Then (T - c) alpha_j lies in W.
    """
    n = len(A)
    basis = []
    steps = []
    while len(basis) < n:
        beta = next(vec(e) for e in identity(n)
                    if not in_span(basis, vec(e)))
        g = conductor(A, beta, basis)
        fac = linear_factors(g, roots)
        assert fac is not None
        c = fac[0]
        h, rem = pdivmod(g, [-c, Fraction(1)])
        assert not rem
        alpha = matvec(pmat(h, A), beta)
        assert not in_span(basis, alpha)
        shifted = matvec(sub(A, scal(c, identity(n))), alpha)
        assert in_span(basis, shifted)
        steps.append((beta, g, c, alpha))
        basis.append(alpha)
    return basis, steps


def fmt(v):
    return '(' + ', '.join(str(x) for x in v) + ')'


def fmtM(A):
    return '[' + '; '.join(' '.join(str(x) for x in r) for r in A) + ']'


def lines_mod(m, p):
    """Invariant lines of the 2x2 matrix m = (a, b, c, d) over F_p."""
    a, b, c, d = m
    found = []
    for x, y in [(0, 1)] + [(1, t) for t in range(p)]:
        tx, ty = a * x + b * y, c * x + d * y
        if any((tx - s * x) % p == 0 and (ty - s * y) % p == 0
               for s in range(p)):
            found.append((x, y))
    return found


def main():
    A = Q([[2, 0, 1, 0], [-2, 4, 0, 1], [-3, 3, 0, 2], [-1, 1, -1, 3]])
    one = identity(4)
    N, N3 = sub(A, scal(2, one)), sub(A, scal(3, one))
    w1, w2 = vec([1, 1, 0, 0]), vec([0, 1, 1, 0])
    W, e4, ev3 = [w1, w2], vec([0, 0, 0, 1]), vec([1, 2, 1, 0])
    assert matvec(A, w1) == vec([2, 2, 0, 0])
    assert matvec(A, w2) == vec([1, 4, 3, 0])
    assert invariant(A, W) and not invariant(A, [vec([1, 0, 0, 0])])

    # block form in the adapted basis w1, w2, w3, w4
    P = Q([[1, 0, 0, 0], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
    M = Q([[2, 1, 1, 0], [0, 3, 0, 1], [0, 0, 2, 1], [0, 0, 0, 2]])
    assert matmul(A, P) == matmul(P, M) and rank(P) == 4
    print('block form P^-1 A P =', fmtM(M))

    N2 = matmul(N, N)
    low = matmul(N2, N3)
    assert low[0] == vec([1, -1, 1, -1])
    assert matmul(N, low) == scal(0, one)
    print('(A-2I)^2(A-3I) first row', fmt(low[0]), '; m_T=(x-2)^3(x-3)')
    assert [len(nullspace(X)) for X in (N, N2, matmul(N2, N))] == [1, 2, 3]
    assert in_span(nullspace(N), w1) and in_span(nullspace(N3), ev3)
    assert in_span(nullspace(N2), vec([0, 0, 1, 1]))

    # conductors of e4 along 0 < span(w1) < W < V
    chain = [[], [w1], W, [vec(e) for e in one]]
    got = [conductor(A, e4, U) for U in chain]
    assert got == [[24, -44, 30, -9, 1], [-12, 16, -7, 1], [4, -4, 1], [1]]
    for g in got:
        print('conductor of e4:', pshow(g))
    for big, small in zip(got, got[1:]):
        assert not pdivmod(big, small)[1]
    powers = krylov(N, e4, 4)
    assert powers[1:] == [vec([0, 1, 2, 1]), vec([2, 3, 1, 0]), ev3, ev3]

    # the lemma step and the triangular basis
    alpha = powers[1]
    assert not in_span(W, alpha) and in_span(W, matvec(N, alpha))
    basis, _ = triangularize(A, vec([2, 3]))
    assert basis == [w1, vec([-3, -3, 1, 1]), vec([-1, -2, -3, -1]),
                     vec([1, 0, 0, 0])]
    tri = Q([[2, 1, 0, 0], [0, 2, 1, 0], [0, 0, 2, 1], [0, 0, 0, 3]])
    S = transpose(basis)
    assert matmul(A, S) == matmul(S, tri)
    print('triangular basis', ' '.join(fmt(v) for v in basis))
    print('triangular matrix', fmtM(tri))

    # annihilator of W under the transpose
    f1, f2 = vec([1, -1, 1, 0]), vec([0, 0, 0, 1])
    assert all(sum(a * b for a, b in zip(f, w)) == 0
               for f in (f1, f2) for w in W)
    assert row_act(A, f1) == [a + b for a, b in zip(f1, f2)]
    assert row_act(A, f2) == [-a + 3 * b for a, b in zip(f1, f2)]
    print('T^t on W^0: [1 -1; 1 3], char poly (x-2)^2')

    # rotation: no rational root; lines over F_5 and F_2
    assert linear_factors(vec([1, 0, 1]), vec([1, -1])) is None
    for p, expect in ((5, [(1, 2), (1, 3)]), (2, [(1, 1)])):
        lines = lines_mod((0, -1, 1, 0), p)
        assert lines == expect
        print('invariant lines of the rotation over F_%d:' % p, lines)

    # exercise 1
    g = conductor(A, e4, [ev3])
    assert g == [-8, 12, -6, 1]
    print('conductor of e4 into E_3:', pshow(g))
    # exercise 3, checked over F_5 by exhaustion
    assert lines_mod((2, 0, 0, 3), 5) == [(0, 1), (1, 0)]
    assert lines_mod((2, 1, 0, 2), 5) == [(1, 0)]
    assert len(lines_mod((2, 0, 0, 2), 5)) == 6
    # exercise 5: derivative, coefficient lists with constant term first
    q = [0, 0, 0, 0, 1]
    degrees = []
    for _ in range(3):
        q = [k * q[k] for k in range(1, len(q))]
        degrees.append(len(q) - 1)
    assert degrees == [3, 2, 1] and q == [0, 24]
    # exercise 6: only W invariant, m_T = (x-2)^2 is not lcm = x-2
    J = sub(Q([[2, 1], [0, 2]]), scal(2, identity(2)))
    assert rank(J) == 1 and rank(matmul(J, J)) == 0
    print('all checks passed')


if __name__ == '__main__':
    main()
