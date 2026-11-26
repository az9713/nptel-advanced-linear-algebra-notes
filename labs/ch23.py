"""Chapter 23 laboratory: representing vectors and adjoints.

Exact arithmetic over Q and Q(i).  A polynomial c0 + c1 x + c2 x^2 is the
list [c0, c1, c2], with <f, g> equal to the integral of f g over [0, 1].
Complex numbers use Python complex with small integer entries, which are
exact in floating point for the operations used here.  The statements about
the space of all polynomials are proved in the chapter, not checked here.
"""
from fractions import Fraction as Fr

import numpy as np

H = [[Fr(1, i + j + 1) for j in range(3)] for i in range(3)]


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(r) for r in zip(*a)]


def solve(m, rhs):
    n = len(m)
    aug = [list(m[i]) + list(rhs[i]) for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if aug[r][c] != 0)
        aug[c], aug[p] = aug[p], aug[c]
        aug[c] = [x / aug[c][c] for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    return [row[n:] for row in aug]


def ip(f, g):
    return sum(f[i] * g[j] * H[i][j] for i in range(3) for j in range(3))


def representer(values):
    """g with <x^k, g> = values[k]; solves the Gram system."""
    return [r[0] for r in solve(H, [[v] for v in values])]


def functionals():
    g0 = representer([Fr(1), Fr(0), Fr(0)])
    assert g0 == [9, -36, 30]
    ortho = [([1, 0, 0], 1), ([Fr(-1, 2), 1, 0], 12),
             ([Fr(1, 6), -1, 1], 180)]
    for a, w in ortho:
        assert ip(a, a) == Fr(1, w)
    for t in (Fr(0), Fr(1, 3), Fr(2)):
        g = [Fr(0)] * 3
        for a, w in ortho:
            value = sum(c * t ** k for k, c in enumerate(a))
            g = [x + w * value * c for x, c in zip(g, a)]
        assert g == representer([t ** k for k in range(3)])
    g1 = representer([Fr(0), Fr(1), Fr(0)])
    assert g1 == [-36, 192, -180]
    print('representer of f(0):', [str(c) for c in g0])
    print("representer of f'(0):", [str(c) for c in g1])


def derivative_adjoint():
    d = [[Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(2)], [Fr(0)] * 3]
    b = solve(H, mul(transpose(d), H))
    assert b == [[-6, 2, 3], [12, -24, -26], [0, 30, 30]]
    basis = [[Fr(int(i == j)) for i in range(3)] for j in range(3)]
    for j in range(3):
        for k in range(3):
            df = [d[i][j] for i in range(3)]
            dsg = [b[i][k] for i in range(3)]
            assert ip(df, basis[k]) == ip(basis[j], dsg)
    assert [b[i][0] for i in range(3)] == [-6, 12, 0]
    e3 = [Fr(1, 6), Fr(-1), Fr(1)]
    assert all(sum(b[i][k] * e3[k] for k in range(3)) == 0
               for i in range(3))
    print('[D*] =', [[int(x) for x in row] for row in b])


def matrices():
    w = [3, 12, -1]
    assert sum(x * x for x in w) == 154
    e = [[Fr(a * b, 154) for b in w] for a in w]
    assert mul(e, e) == e and transpose(e) == e
    assert [[int(154 * x) for x in row] for row in e] == [
        [9, 36, -3], [36, 144, -12], [-3, -12, 1]]
    g = [[Fr(1), Fr(0)], [Fr(0), Fr(2)]]
    a = [[Fr(0), Fr(1)], [Fr(0), Fr(0)]]
    assert solve(g, mul(transpose(a), g)) == [[0, 0], [Fr(1, 2), 0]]
    print('projection and weighted-adjoint examples verified')


def complex_checks():
    # beta = conj(f(gamma)) gamma / |gamma|^2 for f(a) = a1 + i a2
    f = np.array([1, 1j])
    gamma = (2 - 1j) * f.conj()
    assert abs(f @ np.array([1j, -1])) == 0          # (i, -1) in ker f
    assert np.vdot(gamma, np.array([1j, -1])) == 0
    beta = np.conj(f @ gamma) * gamma / np.vdot(gamma, gamma)
    assert np.allclose(beta, f.conj())
    m = np.array([[1, 1j], [0, 2]])
    a = np.array([[1, 2], [3j, 1 - 1j]])
    b = np.array([[2j, 0], [1, 1 + 1j]])

    def tr(x, y):
        return np.trace(y.conj().T @ x)

    assert tr(m @ a, b) == tr(a, m.conj().T @ b)
    assert tr(a @ m, b) == tr(a, b @ m.conj().T)
    u1 = (m + m.conj().T) / 2
    u2 = (m - m.conj().T) / 2j
    assert np.array_equal(u1, u1.conj().T)
    assert np.array_equal(u2, u2.conj().T)
    assert np.array_equal(u1 + 1j * u2, m)
    p = np.array([[1, 0], [0, 0]]) @ np.array([[0, 1], [1, 0]])
    assert not np.array_equal(p, p.T)
    turn = np.array([[0, -1], [1, 0]])
    v = np.array([3, 5])
    assert (turn @ v) @ v == 0 and np.array_equal(turn.T, -turn)
    print('complex representer, trace adjoints, decomposition verified')


def main():
    functionals()
    derivative_adjoint()
    matrices()
    complex_checks()
    print('all chapter 23 checks passed')


if __name__ == '__main__':
    main()
