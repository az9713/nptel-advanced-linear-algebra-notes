"""Chapter 7 laboratory: kernel, image, rank-nullity, prescribed values.

Exact arithmetic over Q (fractions.Fraction).  A matrix is a list of
rows; a vector is a list.  Routines return new objects and do not change
their input.  The laboratory checks:

* the running example A: image, kernel, rank 2, nullity 1;
* rank + nullity = number of columns on 300 integer matrices, and
  row rank = column rank on the same matrices;
* the consistency criterion for prescribed values (the criterion of
  Theorem 7.4), on one failing and one succeeding prescription;
* the infinite-dimensional boundary: integration and differentiation of
  polynomials.
"""
from fractions import Fraction
import random


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def rref(A):
    """Return (R, pivots): reduced row echelon form and pivot columns."""
    R = [row[:] for row in A]
    m, n = len(R), len(R[0])
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i in range(m):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
        if r == m:
            break
    return R, pivots


def transpose(A):
    return [list(col) for col in zip(*A)]


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def rank(A):
    return len(rref(A)[1])


def kernel_basis(A):
    """Basis of {x : A x = 0}, one vector per free column."""
    R, pivots = rref(A)
    n = len(A[0])
    basis = []
    for f in (c for c in range(n) if c not in pivots):
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, p in zip(R, pivots):
            x[p] = -row[f]
        basis.append(x)
    return basis


def image_basis(A):
    """Pivot columns of A: a basis of the column space = image."""
    _, pivots = rref(A)
    return [[row[c] for row in A] for c in pivots]


def row_space_basis(A):
    R, pivots = rref(A)
    return R[:len(pivots)]


def image_conditions(A):
    """Basis of the left null space: y with y^T A = 0.

    v lies in the image iff y . v = 0 for every such y."""
    return kernel_basis(transpose(A))


def prescription_exists(alphas, betas):
    """Theorem 7.4: T with T(alpha_i) = beta_i exists iff every linear
    relation among the alphas holds among the betas."""
    cols = transpose(alphas)
    for c in kernel_basis(cols):
        total = [sum(ci * b[j] for ci, b in zip(c, betas))
                 for j in range(len(betas[0]))]
        if any(t != 0 for t in total):
            return False
    return True


def poly_integrate(p):
    """J p: coefficients [a0, a1, ...] -> integral from 0 to x."""
    return [Fraction(0)] + [a / (k + 1) for k, a in enumerate(p)]


def poly_differentiate(p):
    return [k * a for k, a in enumerate(p)][1:] or [Fraction(0)]


def strip(p):
    p = p[:]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def main():
    # --- running example (Sections 2 and 7) ---------------------------
    A = Q([[1, -1, 2], [2, 1, 0], [-1, -2, 2]])
    R, piv = rref(A)
    r = len(piv)
    ker = kernel_basis(A)
    assert r == 2 and len(ker) == 1
    assert ker[0] == [Fraction(-2, 3), Fraction(4, 3), 1]
    v = [Fraction(-2), Fraction(4), Fraction(3)]
    assert matvec(A, v) == [0, 0, 0]
    img = image_basis(A)
    assert img == [[1, 2, -1], [-1, 1, -2]]
    ys = image_conditions(A)
    assert len(ys) == 1
    y = [Fraction(1), Fraction(-1), Fraction(-1)]
    assert ys[0] == [-1, 1, 1] and [-t for t in ys[0]] == y
    for col in transpose(A):
        assert sum(a * b for a, b in zip(y, col)) == 0
    assert rank(transpose(A)) == r
    print('A rank', r, 'nullity', len(ker), 'kernel', [str(x) for x in v])
    print('R(A) pivots', piv, 'image basis',
          [[str(x) for x in c] for c in img])
    print('image condition y =', [str(x) for x in y], '(a - b - c = 0)')
    print('row space basis', [[str(x) for x in w]
                              for w in row_space_basis(A)])

    # --- rank-nullity and row rank = column rank, random matrices -----
    rng = random.Random(7)
    seen = set()
    for _ in range(300):
        m, n = rng.randint(1, 5), rng.randint(1, 6)
        B = Q([[rng.choice([-2, -1, 0, 0, 1, 2]) for _ in range(n)]
               for _ in range(m)])
        rb = rank(B)
        assert rb + len(kernel_basis(B)) == n
        assert rb == rank(transpose(B))
        assert rb == len(image_basis(B)) == len(row_space_basis(B))
        for x in kernel_basis(B):
            assert all(t == 0 for t in matvec(B, x))
        seen.add((n, rb))
    print('rank-nullity and row rank = column rank verified on 300',
          'matrices;', len(seen), 'distinct (n, rank) pairs')

    # --- prescribed values (Section 4) --------------------------------
    alphas = Q([[1, 1], [0, 1], [1, 0]])
    bad = Q([[1, 1, 0], [0, 1, 1], [1, 2, 1]])
    good = Q([[1, 1, 0], [0, 1, 1], [1, 0, -1]])
    assert not prescription_exists(alphas, bad)
    assert prescription_exists(alphas, good)
    print('relation space of (1,1),(0,1),(1,0):',
          [[str(x) for x in c] for c in kernel_basis(transpose(alphas))])
    t_e1 = [bad[0][j] - bad[1][j] for j in range(3)]
    assert t_e1 == [1, 0, -1] and t_e1 != bad[2]
    print('forced T(1,0) =', [str(x) for x in t_e1],
          'prescribed', [str(x) for x in bad[2]])
    # unique map from the two basis vectors, as a matrix in standard basis
    M = Q([[1, 0], [0, 1], [-1, 1]])
    for a, b in zip(alphas, good):
        assert matvec(M, a) == b
    # prescribed images need not be independent: kernel appears
    assert prescription_exists(Q([[1, 0], [0, 1]]), Q([[1, 2], [2, 4]]))
    # affine map x -> A x + b fails T(0) = 0, hence is not linear
    assert matvec(A, [0, 0, 0])[0] + 1 != 0

    # --- infinite-dimensional boundary (Section 5) --------------------
    for p in ([Fraction(1)], Q([[3, 0, 2, 5]])[0], Q([[0, 1]])[0]):
        assert strip(poly_differentiate(poly_integrate(p))) == strip(p)
    one = [Fraction(1)]
    assert strip(poly_integrate(poly_differentiate(one))) == [0]
    assert strip(poly_integrate(poly_differentiate(one))) != one
    print('D J = identity on polynomials; J D (1) = 0 != 1:',
          'J is injective but not onto, D is onto but not injective')
    print('ch07 laboratory: all assertions passed')


if __name__ == '__main__':
    main()
