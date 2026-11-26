"""Chapter 21 laboratory: inner products, norms, Cauchy-Schwarz.

Exact arithmetic with fractions.Fraction for the polynomial space P_2
(inner product: integral over [0, 1]) and for R^2, R^3.  Complex checks
use Python complex numbers with small integer entries, which are exact
in binary floating point.  Floating-point (numpy) is used only for the
trigonometric orthonormality check, with an explicit tolerance.

Conventions: the inner product is linear in the first argument and
conjugate-linear in the second.  The Gram matrix of an ordered basis
(b_1, ..., b_n) has entries G[j][k] = <b_k, b_j>, so that
<x, y> = Y* G X for coordinate columns X and Y.
"""
from fractions import Fraction as Fr
import math

import numpy as np


# ---------------------------------------------------------------- spaces
def poly_ip(p, q):
    """<p, q> = integral over [0, 1] of p(x) q(x) dx.

    A polynomial is a list of coefficients, constant term first.
    """
    return sum(Fr(a) * Fr(b) / (i + j + 1)
               for i, a in enumerate(p) for j, b in enumerate(q))


def rdot(x, y):
    """Real standard inner product, exact on Fractions."""
    return sum(Fr(a) * Fr(b) for a, b in zip(x, y))


# ------------------------------------------------------ exact linear algebra
def det(m):
    """Determinant by fraction-exact elimination."""
    a = [[Fr(x) for x in row] for row in m]
    n, d = len(a), Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            d = -d
        d *= a[c][c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return d


def leading_minors(m):
    return [det([row[:k] for row in m[:k]]) for k in range(1, len(m) + 1)]


def gram(basis, ip):
    """G[j][k] = <b_k, b_j> (the convention of Section 5)."""
    n = len(basis)
    return [[ip(basis[k], basis[j]) for k in range(n)] for j in range(n)]


def gram_schmidt(betas, ip):
    """Orthogonalize linearly independent vectors (lists of Fractions).

    Returns (alphas, C) with C[j][k] = <beta_k, alpha_j>/||alpha_j||^2
    for j < k, so that beta_k = alpha_k + sum_j C[j][k] alpha_j.
    """
    alphas, n = [], len(betas)
    coef = [[Fr(0)] * n for _ in range(n)]
    for k, b in enumerate(betas):
        a = [Fr(x) for x in b]
        for j, aj in enumerate(alphas):
            c = ip(b, aj) / ip(aj, aj)
            coef[j][k] = c
            a = [x - c * y for x, y in zip(a, aj)]
        assert any(x != 0 for x in a), 'dependent input'
        alphas.append(a)
    for k in range(n):
        coef[k][k] = Fr(1)
    return alphas, coef


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def transpose(a):
    return [list(r) for r in zip(*a)]


# --------------------------------------------------------------- checks
def pad_add(u, v):
    n = max(len(u), len(v))
    u, v = list(u) + [0] * (n - len(u)), list(v) + [0] * (n - len(v))
    return [Fr(a) + Fr(b) for a, b in zip(u, v)]


def lp_norm(x, p):
    if p == 'inf':
        return max(abs(t) for t in x)
    return sum(abs(t) ** p for t in x) ** (1 / p)


def parallelogram_gap(norm, a, b):
    s = [x + y for x, y in zip(a, b)]
    d = [x - y for x, y in zip(a, b)]
    return norm(s) ** 2 + norm(d) ** 2 - 2 * (norm(a) ** 2 + norm(b) ** 2)


def nsq(z):
    """|z|^2 for a complex number, exact for small integer parts."""
    return z.real * z.real + z.imag * z.imag


def cdot(x, y):
    return sum(a * b.conjugate() for a, b in zip(x, y))


def cnorm_sq(x):
    return sum(nsq(a) for a in x)


def bilinear_c2(x, y):
    """Bilinear (not conjugate-linear) form on C^2."""
    return x[0] * y[0] + x[1] * y[1]


def endpoint_form(p, q):
    """p(0)q(0) + p(1)q(1) on coefficient lists of length 3."""
    return p[0] * q[0] + sum(p) * sum(q)


def pullback_singular(a, b):
    """<a, b>_T = (Ta).(Tb) for the singular map T(x1, x2) = x1."""
    return a[0] * b[0]


def apply(m, v):
    return [sum(r[i] * v[i] for i in range(len(v))) for r in m]


def pullback_invertible(a, b):
    """(Ma).(Mb) for M = [[1, 1], [0, 1]]."""
    m = [[1, 1], [0, 1]]
    return rdot(apply(m, a), apply(m, b))


def twisted(a, b):
    """x1y1 - x1y2 - x2y1 + 4 x2y2 on R^2."""
    return a[0] * b[0] - a[0] * b[1] - a[1] * b[0] + 4 * a[1] * b[1]


def n1(v):
    return lp_norm(v, 1)


def ninf(v):
    return lp_norm(v, 'inf')


def n2(v):
    return lp_norm(v, 2)


def main():
    # ---- Section 1: why the conjugate is needed -----------------------
    z = [1, 1j]
    # bilinear fails
    assert bilinear_c2(z, z) == 0 and z != [0, 0]
    assert cdot(z, z) == 2                                # sesquilinear ok
    r = [0, -1, 1]                                        # x^2 - x
    # degenerate
    assert endpoint_form(r, r) == 0 and r != [0, 0, 0]

    # ---- running example P_2, integral inner product -------------------
    mono = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    H = gram(mono, poly_ip)
    assert H == [[Fr(1, j + k + 1) for k in range(3)] for j in range(3)]
    assert H[1][1] == Fr(1, 3) and H[2][2] == Fr(1, 5)
    minors = leading_minors(H)
    assert minors == [Fr(1), Fr(1, 12), Fr(1, 2160)]
    assert det(H) == Fr(1, 2160)
    print('Hilbert H3 =', H)
    print('leading minors =', minors)

    # norm, cosine, angle between 1 and x
    one, xx = [1, 0, 0], [0, 1, 0]
    assert poly_ip(one, xx) == Fr(1, 2)
    assert poly_ip(xx, xx) == Fr(1, 3)
    cos = float(poly_ip(one, xx)) / math.sqrt(float(poly_ip(xx, xx)))
    assert abs(cos - math.sqrt(3) / 2) < 1e-15
    angle_deg = math.degrees(math.acos(cos))
    assert abs(angle_deg - 30.0) < 1e-9
    print('cos(1,x) = %.6f, angle = %.4f degrees' % (cos, angle_deg))
    # Cauchy-Schwarz exactly, and the equality case
    assert poly_ip(one, xx) ** 2 <= poly_ip(one, one) * poly_ip(xx, xx)
    assert Fr(1, 4) < Fr(1, 3)
    three_one = [3, 0, 0]
    assert (poly_ip(one, three_one) ** 2
            == poly_ip(one, one) * poly_ip(three_one, three_one))
    # triangle inequality for 1 and x
    s = pad_add(one, xx)
    n_s = math.sqrt(float(poly_ip(s, s)))
    assert poly_ip(s, s) == Fr(7, 3)
    assert n_s <= 1 + math.sqrt(1 / 3)
    print('||1+x|| = %.6f <= %.6f' % (n_s, 1 + math.sqrt(1 / 3)))
    # parallelogram law for the induced norm (exact)
    d = [Fr(1), Fr(-1), Fr(0)]
    lhs = poly_ip(s, s) + poly_ip(d, d)
    rhs = 2 * (poly_ip(one, one) + poly_ip(xx, xx))
    assert lhs == rhs == Fr(8, 3)

    # ---- Gram-Schmidt on (1, x, x^2) -----------------------------------
    alphas, C = gram_schmidt(mono, poly_ip)
    assert alphas[0] == [1, 0, 0]
    assert alphas[1] == [Fr(-1, 2), 1, 0]
    assert alphas[2] == [Fr(1, 6), -1, 1]
    dd = [poly_ip(a, a) for a in alphas]
    assert dd == [Fr(1), Fr(1, 12), Fr(1, 180)]
    for i in range(3):
        for j in range(i):
            assert poly_ip(alphas[i], alphas[j]) == 0
    assert dd[0] * dd[1] * dd[2] == det(H)
    D = [[dd[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    assert matmul(matmul(transpose(C), D), C) == H       # H = C^T D C
    assert C == [[1, Fr(1, 2), Fr(1, 3)], [0, 1, 1], [0, 0, 1]]
    print('alphas =', alphas)
    print('||alpha_k||^2 =', dd, ' det H3 =', det(H))

    # ---- the same space, another inner product (coordinate pullback) ----
    assert rdot(one, xx) == 0 and poly_ip(one, xx) != 0
    Gc = gram(mono, rdot)
    assert Gc == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    # ---- Section 3 pullback: singular T gives only semi-definite -------
    assert pullback_singular([0, 1], [0, 1]) == 0  # T singular
    M = [[1, 1], [0, 1]]                                  # T invertible
    GM = gram([[1, 0], [0, 1]], pullback_invertible)
    assert GM == matmul(transpose(M), M) == [[1, 1], [1, 2]]
    assert leading_minors(GM) == [1, 1]

    # ---- Section 5: Gram matrix of R^2 twisted inner product ------------
    Gt = [[Fr(1), Fr(-1)], [Fr(-1), Fr(4)]]
    assert gram([[1, 0], [0, 1]], twisted) == Gt
    assert leading_minors(Gt) == [1, 3]
    for X, Y in (([2, 3], [5, -1]), ([1, 1], [1, 1])):
        yGx = sum(Y[j] * Gt[j][k] * X[k] for j in range(2) for k in range(2))
        assert yGx == twisted(X, Y)
    # invertible and Hermitian but not positive definite
    Gbad = [[1, 0], [0, -1]]
    assert det(Gbad) != 0 and leading_minors(Gbad) == [1, -1]
    # orthogonality is relative to the inner product
    t = (-3 + math.sqrt(13)) / 2
    alpha, beta = [1.0, t], [-t, 1.0]
    val = (alpha[0] * beta[0] - alpha[0] * beta[1]
           - alpha[1] * beta[0] + 4 * alpha[1] * beta[1])
    assert abs(val) < 1e-12
    print('twisted orthogonal slope y/x = %.6f' % t)

    # ---- Section 4: polarization over C ----------------------------------
    x, y = [1, 1j], [2, 1 + 1j]
    ip_xy = cdot(x, y)
    assert ip_xy == 3 + 1j

    def sq(ws):
        return cnorm_sq(ws)

    pol = 0
    for n in range(1, 5):
        w = [a + (1j ** n) * b for a, b in zip(x, y)]
        pol += (1j ** n) * sq(w)
    assert pol / 4 == ip_xy
    xp = [a + b for a, b in zip(x, y)]
    xm = [a - b for a, b in zip(x, y)]
    assert (sq(xp) - sq(xm)) / 4 == ip_xy.real == 3
    xip = [a + 1j * b for a, b in zip(x, y)]
    xim = [a - 1j * b for a, b in zip(x, y)]
    assert sq(xip) == 10 and sq(xim) == 6
    assert (sq(xip) - sq(xim)) / 4 == ip_xy.imag == 1
    print('<x,y> =', ip_xy, ' by polarization =', pol / 4)

    # ---- parallelogram law fails for l1 and sup norms ----------------------
    e1, e2 = [1, 0], [0, 1]
    assert parallelogram_gap(n1, e1, e2) == 4
    assert parallelogram_gap(ninf, e1, e2) == -2
    assert abs(parallelogram_gap(n2, e1, e2)) < 1e-12
    print('parallelogram gaps (l1, sup, l2):',
          parallelogram_gap(n1, e1, e2), parallelogram_gap(ninf, e1, e2),
          parallelogram_gap(n2, e1, e2))

    # ---- Pythagoras and its converse over C
    u, v = [1], [1j]
    assert cdot(u, v) != 0                                # not orthogonal
    assert cnorm_sq([a + b for a, b in zip(u, v)]) == 2 == (
        cnorm_sq(u) + cnorm_sq(v))                        # equation holds

    # ---- Section 10: trigonometric orthonormal system on [0, 1]
    N = 64
    xs = np.arange(N) / N
    fam = [np.ones(N)]
    for k in (1, 2, 3):
        fam.append(np.sqrt(2) * np.cos(2 * np.pi * k * xs))
        fam.append(np.sqrt(2) * np.sin(2 * np.pi * k * xs))
    Gtrig = np.array([[np.dot(a, b) / N for b in fam] for a in fam])
    assert np.allclose(Gtrig, np.eye(len(fam)), atol=1e-12)
    print('max |Gram - I| for 7 trigonometric functions: %.2e'
          % np.abs(Gtrig - np.eye(len(fam))).max())

    # ---- Section 12: the Gram-Schmidt example in R^3
    b = [[3, 0, 4], [-1, 0, 7], [2, 9, 11]]
    assert det(b) != 0
    a, Cr = gram_schmidt(b, rdot)
    assert a == [[3, 0, 4], [-4, 0, 3], [0, 9, 0]]
    assert [rdot(t, t) for t in a] == [25, 25, 81]
    assert Cr[0][1] == 1 and Cr[0][2] == 2 and Cr[1][2] == 1
    print('R^3 Gram-Schmidt:', [[int(t) for t in r] for r in a])

    # ---- Section 11: coefficients of an orthogonal expansion
    beta = [Fr(2) * p + Fr(-3) * q for p, q in zip(a[0], a[1])]
    coeffs = [rdot(beta, t) / rdot(t, t) for t in a]
    assert coeffs == [2, -3, 0]
    # dependent set is rejected by Gram-Schmidt
    try:
        gram_schmidt([[1, 0], [2, 0]], rdot)
        raise SystemExit('dependent input was accepted')
    except AssertionError:
        pass
    # ---- exercises ------------------------------------------------------
    sq_ = [0, 0, 1]
    c2 = poly_ip(xx, sq_) ** 2 / (poly_ip(xx, xx) * poly_ip(sq_, sq_))
    assert poly_ip(xx, sq_) == Fr(1, 4) and c2 == Fr(15, 16)
    deg = math.degrees(math.acos(math.sqrt(15) / 4))
    assert abs(deg - 14.4775) < 1e-3
    print('cos^2(x, x^2) = 15/16, angle = %.4f degrees' % deg)
    # x - x^2/... : distance data used in the opening
    best = [Fr(1, 2), 0, 0]
    gap = pad_add(xx, [-c for c in best])
    assert poly_ip(gap, one) == 0 and poly_ip(gap, gap) == Fr(1, 12)
    print('all chapter 21 checks passed')


if __name__ == '__main__':
    main()
