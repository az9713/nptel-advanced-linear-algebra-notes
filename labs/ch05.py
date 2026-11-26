"""Chapter 5 laboratory: coordinates and change of ordered basis.

Exact arithmetic over Q (fractions.Fraction) and brute force over F_p.
A matrix is a list of rows.  A basis of F^n is given by the matrix M
whose columns are its vectors, written in standard coordinates.  For the
space of polynomials of degree at most 2 the monomial basis (1, x, x^2)
plays the role of the standard basis.

Convention (as in the chapter): P(B <- B') is the matrix with
[v]_B = P(B <- B') [v]_B', and its column j is [alpha'_j]_B.
"""
from fractions import Fraction
from itertools import product


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(r) for r in zip(*A)]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def rref(A):
    """Return (R, pivot columns) for A over Q, by Gauss-Jordan."""
    R = [row[:] for row in A]
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
                R[i] = [x - f * y for x, y in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
        if r == len(R):
            break
    return R, pivots


def rank(A):
    return len(rref(A)[1])


def inverse(A):
    n = len(A)
    aug = [row + identity(n)[i] for i, row in enumerate(A)]
    R, piv = rref(aug)
    if piv[:n] != list(range(n)) or len(piv) < n:
        raise ValueError('matrix is singular')
    return [row[n:] for row in R]


def coords(M, v):
    """Coordinates of the standard-coordinate vector v in basis M."""
    return matvec(inverse(M), v)


def change_matrix(MB, MBp):
    """P(B <- B') = MB^(-1) MBp (columns are [alpha'_j]_B)."""
    return matmul(inverse(MB), MBp)


def change_matrix_rref(MB, MBp):
    """The same matrix by row reducing [MB | MBp] to [I | P]."""
    n = len(MB)
    aug = [a + b for a, b in zip(MB, MBp)]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n))
    return [row[n:] for row in R]


def basis_from_P(MB, P):
    """New basis alpha'_j = sum_i P[i][j] alpha_i, as a matrix of columns."""
    return matmul(MB, P)


def det_mod(A, p):
    """Determinant of a small integer matrix modulo p (Laplace)."""
    if len(A) == 1:
        return A[0][0] % p
    total = 0
    for j in range(len(A)):
        minor = [r[:j] + r[j + 1:] for r in A[1:]]
        total += (-1) ** j * A[0][j] * det_mod(minor, p)
    return total % p


def count_ordered_bases(n, p):
    """Number of ordered bases of F_p^n, by brute force.

    A tuple of n vectors is a basis iff the matrix with those columns has
    nonzero determinant mod p.
    """
    vecs = list(product(range(p), repeat=n))
    bases = 0
    for tup in product(vecs, repeat=n):
        cols = [list(v) for v in tup]
        if det_mod(transpose(cols), p) != 0:
            bases += 1
    return bases


def gl_order(n, q):
    out = 1
    for i in range(n):
        out *= q ** n - q ** i
    return out


def show(M):
    return [[str(x) for x in r] for r in M]


def main():
    # Running example: P_2(Q), B = (1, x, x^2), B' = (1, x-1, (x-1)^2).
    MB = identity(3)
    MBp = Q([[1, -1, 1], [0, 1, -2], [0, 0, 1]])
    P = change_matrix(MB, MBp)
    Pinv = inverse(P)
    assert P == MBp
    assert Pinv == Q([[1, 1, 1], [0, 1, 2], [0, 0, 1]])
    X = [Fraction(2), Fraction(3), Fraction(5)]      # 2 + 3x + 5x^2
    Xp = matvec(Pinv, X)
    assert Xp == [10, 13, 5]                  # p(1), p'(1), p''(1)/2
    assert matvec(P, Xp) == X
    wrong = matvec(P, X)                      # the direction error
    assert wrong == [4, -7, 5] and wrong != Xp
    print('P =', show(P))
    print('P^-1 =', show(Pinv))
    print('X =', [str(x) for x in X], ' Xp =', [str(x) for x in Xp])
    print('P X (wrong direction) =', [str(x) for x in wrong])

    # Taylor check: coefficients of p in powers of (x-1).
    def p(t):
        return 2 + 3 * t + 5 * t * t
    for t in [Fraction(k, 3) for k in range(-4, 6)]:
        s = t - 1
        assert p(t) == Xp[0] + Xp[1] * s + Xp[2] * s * s

    # Composition: B'' = (1, x-2, (x-2)^2).
    MBpp = Q([[1, -2, 4], [0, 1, -4], [0, 0, 1]])
    P2 = change_matrix(MB, MBpp)
    P_bp_bpp = change_matrix(MBp, MBpp)
    assert P_bp_bpp == P                      # shift by 1 again
    assert matmul(P, P_bp_bpp) == P2 == matmul(P, P)
    assert change_matrix(MBp, MB) == Pinv     # P(B'<-B) = P(B<-B')^-1
    print('P^2 =', show(P2))

    # Linear functional f(p) = p(1): coordinate row a' = a P.
    a = Q([[1, 1, 1]])
    ap = matmul(a, P)
    assert ap == Q([[1, 0, 0]])
    print('f = p(1): a =', [str(x) for x in a[0]], ' ap =',
          [str(x) for x in ap[0]])

    # Second example in Q^3 with two non-standard bases.
    B1 = Q([[1, 0, 1], [1, 1, 0], [0, 1, 1]])    # columns (1,1,0),(0,1,1),..
    B2 = Q([[1, 0, 0], [1, 1, 0], [1, 1, 1]])    # (1,1,1),(0,1,1),(0,0,1)
    P12 = change_matrix(B1, B2)
    assert change_matrix_rref(B1, B2) == P12
    v = [Fraction(4), Fraction(-2), Fraction(6)]
    x1, x2 = coords(B1, v), coords(B2, v)
    assert matvec(P12, x2) == x1 and matvec(inverse(P12), x1) == x2
    assert matvec(B1, x1) == v and matvec(B2, x2) == v
    print('B1 columns:', show(transpose(B1)))
    print('B2 columns:', show(transpose(B2)))
    print('P(B1<-B2) =', show(P12))
    print('v =', [str(x) for x in v], ' [v]_B1 =', [str(x) for x in x1],
          ' [v]_B2 =', [str(x) for x in x2])
    # In closed form, [(x, y, z)]_B2 = (x, y - x, z - y).
    assert x2 == [4, -6, 8]
    # Swapping the order of a basis permutes the coordinates.
    swapped = [r[:] for r in B2]
    for r in swapped:
        r[0], r[1] = r[1], r[0]
    assert coords(swapped, v) == [x2[1], x2[0], x2[2]]

    # Converse: build B' from an invertible P (B = monomials).
    Pc = Q([[1, 2, 0], [0, 1, 1], [1, 0, 1]])
    Mnew = basis_from_P(MB, Pc)
    assert transpose(Mnew) == Q([[1, 0, 1], [2, 1, 0], [0, 1, 1]])
    assert rank(Mnew) == 3
    Xc = matvec(inverse(Pc), X)
    assert matvec(Pc, Xc) == X
    print('converse: P^-1 =', show(inverse(Pc)), ' [p]_B2 =',
          [str(x) for x in Xc])

    # Boundary: a singular P gives dependent vectors, no basis.
    Ps = Q([[1, 1], [1, 1]])
    assert rank(basis_from_P(identity(2), Ps)) == 1
    try:
        inverse(Ps)
        raise AssertionError('singular matrix inverted')
    except ValueError:
        pass
    assert rank(Q([[1, 1, 0], [0, 0, 1], [0, 0, 1]])) == 2

    # Ordered bases correspond to invertible matrices (counts).
    for n, q in [(2, 2), (3, 2), (2, 3)]:
        nb = count_ordered_bases(n, q)
        assert nb == gl_order(n, q)
        print('ordered bases of F_%d^%d: %d' % (q, n, nb))
    assert gl_order(2, 2) == 6 and gl_order(3, 2) == 168
    assert gl_order(2, 3) == 48
    print('all checks passed')


if __name__ == '__main__':
    main()
