"""Chapter 20 laboratory: diagonalizable plus nilpotent, canonical forms.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows, a vector a list.  The laboratory checks the running example (the
3 x 3 matrix with minimal polynomial (x-1)(x-2)^2), the chain relations
behind the solutions of X' = AX, the nilpotent block counts, the chain
basis built by induction on random nilpotent matrices, companion
matrices, and the computations of the exercises.
"""
from fractions import Fraction
import random


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def vec(*xs):
    return [Fraction(x) for x in xs]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def zeros(n):
    return [[Fraction(0)] * n for _ in range(n)]


def mul(A, B):
    return [[sum((a * b for a, b in zip(r, c)), Fraction(0))
             for c in zip(*B)] for r in A]


def add(A, B, c=1):
    return [[x + c * y for x, y in zip(r, s)] for r, s in zip(A, B)]


def apply(A, v):
    return [sum((a * x for a, x in zip(r, v)), Fraction(0)) for r in A]


def power(A, k):
    out = eye(len(A))
    for _ in range(k):
        out = mul(out, A)
    return out


def cols(vectors):
    return [list(r) for r in zip(*vectors)]


def rref(A):
    R = [r[:] for r in A]
    piv, r = [], 0
    for c in range(len(R[0]) if R else 0):
        s = next((i for i in range(r, len(R)) if R[i][c] != 0), None)
        if s is None:
            continue
        R[r], R[s] = R[s], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [x - f * y for x, y in zip(R[i], R[r])]
        piv.append(c)
        r += 1
        if r == len(R):
            break
    return R, piv


def rank(vectors):
    return len(rref(cols(vectors))[1]) if vectors else 0


def nullspace(A):
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


def solve(vectors, b):
    """Coefficients x with sum x_j vectors[j] = b (vectors independent)."""
    M = [list(row) + [y] for row, y in zip(cols(vectors), b)]
    R, piv = rref(M)
    k = len(vectors)
    assert k not in piv, 'b is not in the span'
    x = [Fraction(0)] * k
    for row, p in zip(R, piv):
        x[p] = row[k]
    return x


def inverse(A):
    n = len(A)
    R, piv = rref([r + e for r, e in zip(A, eye(n))])
    assert piv[:n] == list(range(n))
    return [r[n:] for r in R]


def nullity(A):
    return len(nullspace(A))


def jordan(blocks):
    """Block diagonal matrix of Jordan blocks given as (c, size)."""
    n = sum(s for _, s in blocks)
    J = zeros(n)
    at = 0
    for c, s in blocks:
        for i in range(s):
            J[at + i][at + i] = Fraction(c)
            if i + 1 < s:
                J[at + i][at + i + 1] = Fraction(1)
        at += s
    return J


def block_counts(N):
    """Number of nilpotent blocks of each size, from the nullities."""
    n = len(N)
    m = [nullity(power(N, i)) for i in range(n + 2)]
    return {i: 2 * m[i] - m[i - 1] - m[i + 1]
            for i in range(1, n + 1) if 2 * m[i] - m[i - 1] - m[i + 1]}


def chain_basis(N, space):
    """Chains [v, Nv, ...] forming a basis of the invariant span(space).

    This is the induction of the chapter: build chains for N(space),
    lift each head one step, then complete the chain ends to a basis of
    the kernel of N inside the space.
    """
    if not space:
        return []
    images = [apply(N, v) for v in space]
    piv = rref(cols(images))[1]
    image = [images[j] for j in piv]
    assert len(image) < len(space), 'N is not nilpotent on the space'
    chains = []
    for chain in chain_basis(N, image):
        head = chain[0]
        coef = solve(image, head)
        lift = [sum(c * space[piv[i]][t] for i, c in enumerate(coef))
                for t in range(len(head))]
        assert apply(N, lift) == head
        chains.append([lift] + chain)
    ends = [c[-1] for c in chains]
    # kernel of N restricted to span(space)
    kernel = [[sum(c * space[i][t] for i, c in enumerate(z))
               for t in range(len(space[0]))]
              for z in nullspace(cols(images))]
    for w in kernel:
        if rank(ends + [w]) > len(ends):
            ends.append(w)
            chains.append([w])
    return chains


def companion(coeffs):
    """Companion matrix of x^k + a_{k-1} x^{k-1} + ... + a_0.

    coeffs lists a_0, ..., a_{k-1}.
    """
    k = len(coeffs)
    C = zeros(k)
    for i in range(1, k):
        C[i][i - 1] = Fraction(1)
    for i in range(k):
        C[i][k - 1] = Fraction(-coeffs[i])
    return C


def block_diag(*blocks):
    n = sum(len(b) for b in blocks)
    M = zeros(n)
    at = 0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)):
                M[at + i][at + j] = b[i][j]
        at += len(b)
    return M


def poly_at(coeffs, A):
    """p(A) for p given by coefficients, constant term first."""
    out = zeros(len(A))
    for c in reversed(coeffs):
        out = add(mul(out, A), eye(len(A)), Fraction(c))
    return out


def main():
    # running example
    A = Q([[3, 1, -1], [2, 2, -1], [2, 2, 0]])
    one = eye(3)
    M2 = add(A, one, -2)
    E1 = mul(M2, M2)
    E2 = add(one, E1, -1)
    D = add(E1, E2, 2)
    N = add(A, D, -1)
    assert mul(E1, E1) == E1 and mul(E1, E2) == zeros(3)
    assert mul(D, N) == mul(N, D) and mul(N, N) == zeros(3) != N
    assert E1 == Q([[1, -1, 0], [0, 0, 0], [2, -2, 0]])
    assert D == Q([[1, 1, 0], [0, 2, 0], [-2, 2, 2]])
    assert N == Q([[2, 0, -1], [2, 0, -1], [4, 0, -2]])
    print('E1 =', [[str(x) for x in r] for r in E1])
    print('D =', [[str(x) for x in r] for r in D])
    print('N =', [[str(x) for x in r] for r in N])
    u1, u, w = vec(1, 0, 2), vec(1, 1, 2), vec(1, 1, 1)
    assert apply(A, u1) == u1
    assert apply(M2, u) == vec(0, 0, 0) and apply(M2, w) == u
    P = cols([u1, u, w])
    J = jordan([(1, 1), (2, 2)])
    assert mul(A, P) == mul(P, J)
    for k in range(1, 12):
        lhs = apply(power(A, k), w)
        assert lhs == [2**k * b + k * 2**(k - 1) * a for a, b in zip(u, w)]
    assert solve([u1, u, w], vec(0, 0, 1)) == [0, 1, -1]
    print('Jordan basis (1,0,2), (1,1,2), (1,1,1): checked')

    # first recording: the 2 x 2 system
    B = Q([[-2, 1], [0, -2]])
    v1, v2 = vec(1, 0), vec(0, 1)
    assert apply(add(B, eye(2), 2), v2) == v1
    assert apply(add(B, eye(2), 2), v1) == vec(0, 0)

    # nilpotent block counts: sizes (3,1,1) and (3,2)
    random.seed(20)
    for sizes in ((3, 1, 1), (3, 2), (3, 2, 2), (3, 3, 1), (2, 2)):
        n = sum(sizes)
        J = jordan([(0, s) for s in sizes])
        while True:
            S = Q([[random.randint(-2, 2) for _ in range(n)]
                   for _ in range(n)])
            if rank(cols(S)) == n:
                break
        Nm = mul(mul(S, J), inverse(S))
        expected = {s: sizes.count(s) for s in set(sizes)}
        assert block_counts(Nm) == expected
        chains = chain_basis(Nm, eye(n))
        assert sorted(map(len, chains)) == sorted(sizes)
        assert rank([v for c in chains for v in c]) == n
        for c in chains:
            assert [apply(Nm, v) for v in c] == c[1:] + [[0] * n]
        print('sizes', sizes, 'nullity', nullity(Nm), 'm2',
              nullity(mul(Nm, Nm)), 'recovered', expected)

    # rational forms: p1 = (t-1)^4 - 8, p2 = t - 3 on an 8-dim space
    p1 = [-7, -4, 6, -4, 1]
    C1 = companion(p1[:-1])
    C2 = companion([9, -6])
    c3 = companion([-3])
    assert poly_at(p1, C1) == zeros(4)
    for M in (block_diag(C1, C2, C2), block_diag(C1, C2, c3, c3)):
        a = poly_at(p1, M)
        b = add(M, eye(8), -3)
        assert mul(a, mul(b, b)) == zeros(8)
        assert mul(a, b) != zeros(8) and mul(b, b) != zeros(8)
    print('companion of p1, last column:',
          [str(r[-1]) for r in C1], '; two rational forms checked')
    # (t-1)^4 = 8 has no rational root
    assert all((Fraction(a, 7) - 1)**4 != 8 for a in range(-60, 61))

    # exercise 4: non-commuting split of [[1,1],[0,2]]
    T = Q([[1, 1], [0, 2]])
    Dg, Ng = Q([[1, 0], [0, 2]]), Q([[0, 1], [0, 0]])
    assert add(Dg, Ng) == T and mul(Dg, Ng) != mul(Ng, Dg)
    # exercise 5: two copies of x^2 + 1
    R = block_diag(companion([1, 0]), companion([1, 0]))
    assert mul(R, R) == add(zeros(4), eye(4), -1)
    print('all checks passed')


if __name__ == '__main__':
    main()
