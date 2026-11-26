"""Chapter 22 laboratory: best approximation and orthogonal projection.

Exact arithmetic over Q (fractions.Fraction) for the real inner product
<u, v> = sum u_i v_i on Q^n.  A vector is a list of Fractions.  The
floating-point parts (the infinite-dimensional example and one complex
check) use explicit tolerances.

What is checked: Gram-Schmidt, the coefficient formula, the normal
equations, the characterisation "beta - alpha is orthogonal to W", the
Pythagorean identity behind minimality, idempotence, kernel and direct
sum of the orthogonal projection, a non-orthogonal projection that is not
a best approximation, and Bessel's inequality.  What is not checked: any
statement about a general inner product space; the proofs are in the
chapter.
"""
import math
from fractions import Fraction as Fr

import numpy as np


def vec(*xs):
    return [Fr(x) for x in xs]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def add(u, v):
    return [a + b for a, b in zip(u, v)]


def sub(u, v):
    return [a - b for a, b in zip(u, v)]


def scale(c, u):
    return [c * a for a in u]


def norm2(u):
    return dot(u, u)


def matvec(m, v):
    return [dot(row, v) for row in m]


def gram_schmidt(vectors):
    """Orthogonal (not normalised) basis of span(vectors)."""
    out = []
    for v in vectors:
        w = list(v)
        for a in out:
            w = sub(w, scale(dot(v, a) / norm2(a), a))
        if all(x == 0 for x in w):
            raise ValueError('vectors are linearly dependent')
        out.append(w)
    return out


def coefficients(beta, ortho):
    return [dot(beta, a) / norm2(a) for a in ortho]


def project_orthogonal(beta, ortho):
    """Best approximation of beta in span(ortho); ortho is orthogonal."""
    alpha = [Fr(0)] * len(beta)
    for c, a in zip(coefficients(beta, ortho), ortho):
        alpha = add(alpha, scale(c, a))
    return alpha


def solve(m, b):
    """Solve m x = b exactly for a square invertible m."""
    n = len(m)
    aug = [list(row) + [b[i]] for i, row in enumerate(m)]
    for col in range(n):
        piv = next(r for r in range(col, n) if aug[r][col] != 0)
        aug[col], aug[piv] = aug[piv], aug[col]
        p = aug[col][col]
        aug[col] = [x / p for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    return [aug[i][n] for i in range(n)]


def project_normal(beta, cols):
    """Best approximation via the normal equations A^T A x = A^T beta."""
    g = [[dot(a, b) for b in cols] for a in cols]
    rhs = [dot(a, beta) for a in cols]
    x = solve(g, rhs)
    alpha = [Fr(0)] * len(beta)
    for c, a in zip(x, cols):
        alpha = add(alpha, scale(c, a))
    return x, alpha


def nullspace(rows, n):
    """Basis of {x in Q^n : row . x = 0 for every row} by RREF."""
    m = [list(r) for r in rows]
    pivots = []
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        p = m[r][c]
        m[r] = [x / p for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    basis = []
    for free in (c for c in range(n) if c not in pivots):
        x = [Fr(0)] * n
        x[free] = Fr(1)
        for i, pc in enumerate(pivots):
            x[pc] = -m[i][free]
        basis.append(x)
    return basis


def projection_matrix(ortho, n):
    cols = []
    for j in range(n):
        e = [Fr(1) if i == j else Fr(0) for i in range(n)]
        cols.append(project_orthogonal(e, ortho))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def show(v):
    return '(' + ', '.join(str(x) for x in v) + ')'


def running_example():
    u1, u2 = vec(1, 1, 0), vec(1, 0, 1)
    beta = vec(1, 2, 4)
    a1, a2 = gram_schmidt([u1, u2])
    assert a1 == vec(1, 1, 0) and a2 == vec(Fr(1, 2), Fr(-1, 2), 1)
    assert dot(a1, a2) == 0
    assert norm2(a1) == 2 and norm2(a2) == Fr(3, 2)
    c = coefficients(beta, [a1, a2])
    assert c == [Fr(3, 2), Fr(7, 3)]
    alpha = project_orthogonal(beta, [a1, a2])
    assert alpha == vec(Fr(8, 3), Fr(1, 3), Fr(7, 3))
    resid = sub(beta, alpha)
    assert resid == vec(Fr(-5, 3), Fr(5, 3), Fr(5, 3))
    # characterisation: beta - alpha is orthogonal to W
    assert dot(resid, u1) == 0 and dot(resid, u2) == 0
    print('alpha =', show(alpha), ' beta - alpha =', show(resid))
    print('|beta|^2 =', norm2(beta), ' |alpha|^2 =', norm2(alpha),
          ' |beta - alpha|^2 =', norm2(resid))
    assert norm2(beta) == norm2(alpha) + norm2(resid)
    print('distance^2 =', norm2(resid), ' distance =',
          round(math.sqrt(norm2(resid)), 6))

    # normal equations in the non-orthogonal basis (u1, u2)
    x, alpha2 = project_normal(beta, [u1, u2])
    assert x == [Fr(1, 3), Fr(7, 3)] and alpha2 == alpha
    print('normal-equation coordinates in (u1, u2):', show(x))

    # W-perp by null space, and the matrix of E
    perp = nullspace([u1, u2], 3)
    assert len(perp) == 1
    n = vec(1, -1, -1)
    assert perp[0] == scale(Fr(-1), n)   # a multiple of n
    e = projection_matrix([a1, a2], 3)
    third = Fr(1, 3)
    assert e == [vec(2 * third, third, third), vec(third, 2 * third, -third),
                 vec(third, -third, 2 * third)]
    assert matvec(e, beta) == alpha
    ee = [[dot(row, [r[j] for r in e]) for j in range(3)] for row in e]
    assert ee == e                                   # idempotent
    assert all(e[i][j] == e[j][i] for i in range(3) for j in range(3))
    assert matvec(e, n) == vec(0, 0, 0)              # kernel contains W-perp
    assert sum(e[i][i] for i in range(3)) == 2       # trace = dim W
    print('E = (1/3) *', [[int(3 * x) for x in row] for row in e])

    # Pythagorean identity: |beta - g|^2 = |beta - alpha|^2 + |alpha - g|^2
    count = 0
    for a in range(-4, 5):
        for b in range(-4, 5):
            g = add(scale(Fr(a), u1), scale(Fr(b), u2))
            lhs = norm2(sub(beta, g))
            rhs = norm2(resid) + norm2(sub(alpha, g))
            assert lhs == rhs and lhs >= norm2(resid)
            count += 1
    print('Pythagorean identity verified at', count, 'points of W')

    # Bessel: partial and full orthogonal sets
    s1 = dot(beta, a1) ** 2 / norm2(a1)
    s2 = s1 + dot(beta, a2) ** 2 / norm2(a2)
    s3 = s2 + dot(beta, n) ** 2 / norm2(n)
    assert (s1, s2, s3) == (Fr(9, 2), Fr(38, 3), Fr(21))
    assert s1 <= s2 <= norm2(beta) and s3 == norm2(beta)
    print('Bessel sums:', s1, s2, s3, ' |beta|^2 =', norm2(beta))
    # direct sum: beta = E beta + (beta - E beta), unique components
    assert add(alpha, resid) == beta
    assert dot(alpha, resid) == 0


def one_dimensional_example():
    beta, w = vec(-10, 2, 8), vec(3, 12, -1)
    assert norm2(w) == 154 and dot(beta, w) == -14
    alpha = project_orthogonal(beta, [w])
    assert alpha == scale(Fr(-1, 11), w)
    assert dot(sub(beta, alpha), w) == 0
    print('line example: coefficient', dot(beta, w) / norm2(w),
          ' alpha =', show(alpha))


def non_orthogonal_projection():
    p = [vec(1, 1), vec(0, 0)]
    pp = [[dot(row, [r[j] for r in p]) for j in range(2)] for row in p]
    assert pp == p                                   # idempotent
    assert p[0][1] != p[1][0]                        # not symmetric
    beta = vec(0, 1)
    pb = matvec(p, beta)
    best = project_orthogonal(beta, [vec(1, 0)])
    assert pb == vec(1, 0) and best == vec(0, 0)
    assert norm2(sub(beta, pb)) == 2 and norm2(sub(beta, best)) == 1
    assert dot(sub(beta, pb), vec(1, 0)) != 0
    print('non-orthogonal projection: |beta - P beta|^2 = 2 against',
          'best-approximation distance^2 = 1')


def gram_schmidt_four():
    vs = [vec(1, 1, 1, 1), vec(1, 2, 3, 4), vec(2, 0, 1, 5)]
    g = gram_schmidt(vs)
    for i in range(3):
        for j in range(i):
            assert dot(g[i], g[j]) == 0
    # each span(v_1..v_k) equals span(g_1..g_k): v_k - g_k in earlier span
    for k in range(3):
        diff = sub(vs[k], g[k])
        assert project_orthogonal(diff, g[:k]) == diff if k else diff == \
            vec(0, 0, 0, 0)
    try:
        gram_schmidt([vec(1, 2), vec(2, 4)])
        raise AssertionError('dependent input must be rejected')
    except ValueError:
        pass
    print('Gram-Schmidt in Q^4 gave', [show(x) for x in g])


def infinite_dimensional_failure():
    # V = finitely supported real sequences, W = kernel of
    # phi(x) = sum x_k / k.  Check j e_j - k e_k lies in W (exactly).
    for j in range(1, 7):
        for k in range(1, 7):
            assert j * Fr(1, j) - k * Fr(1, k) == 0
    s = 0.0
    seen = {}
    prev = None
    limit = math.sqrt(6) / math.pi
    for k in range(1, 10 ** 6 + 1):
        s += 1.0 / (k * k)
        if k in (1, 10, 100, 10 ** 4, 10 ** 6):
            dist = 1.0 / math.sqrt(s)
            seen[k] = dist
            assert dist > limit                      # never attained
            assert prev is None or dist < prev       # decreasing
            prev = dist
    assert abs(seen[10 ** 6] - limit) < 1e-5
    print('c00 example: sqrt(6)/pi =', round(limit, 6))
    for k, d in seen.items():
        print('   n =', k, ' |e1 - w_n| =', round(d, 6))


def complex_check():
    rng = np.random.default_rng(22)
    a = rng.standard_normal((5, 2)) + 1j * rng.standard_normal((5, 2))
    q, _ = np.linalg.qr(a)                  # orthonormal basis of W
    beta = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    # <x, y> = y^* x, linear in the first slot
    coef = q.conj().T @ beta
    alpha = q @ coef
    assert np.allclose(q.conj().T @ (beta - alpha), 0, atol=1e-12)
    assert abs(np.vdot(beta, beta).real
               - np.vdot(alpha, alpha).real
               - np.vdot(beta - alpha, beta - alpha).real) < 1e-12
    assert np.sum(np.abs(coef) ** 2) <= np.vdot(beta, beta).real + 1e-12
    for _ in range(200):
        gam = q @ (rng.standard_normal(2) + 1j * rng.standard_normal(2))
        assert np.linalg.norm(beta - gam) >= np.linalg.norm(beta - alpha)
    print('complex check passed: residual orthogonal, Pythagoras holds')


def exercises():
    x, alpha = project_normal(vec(1, 3, 4), [vec(1, 1, 1), vec(0, 1, 2)])
    assert x == [Fr(7, 6), Fr(3, 2)]
    resid = sub(vec(1, 3, 4), alpha)
    assert resid == [Fr(-1, 6), Fr(1, 3), Fr(-1, 6)]
    w, b = [1, 1j], [1, 0]
    c = sum(a * z.conjugate() for a, z in zip(b, w)) / 2
    d = [a - c * z for a, z in zip(b, w)]
    assert c == 0.5 and sum(a * z.conjugate() for a, z in zip(d, w)) == 0
    print('exercise 2: x = 7/6, 3/2; exercise 5: coefficient 1/2')


def main():
    exercises()
    running_example()
    one_dimensional_example()
    non_orthogonal_projection()
    gram_schmidt_four()
    infinite_dimensional_failure()
    complex_check()
    print('all chapter 22 checks passed')


if __name__ == '__main__':
    main()
