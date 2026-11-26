"""Chapter 25 laboratory: spectral resolutions.

Exact arithmetic.  Rational matrices use fractions.Fraction.  Complex
matrices with integer or half-integer entries use Python complex, which is
exact for these values.  The general theorems are proved in the chapter;
this file checks the printed instances.
"""
from fractions import Fraction as Fr

import numpy as np


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def comb(c, a, d, b):
    return [[c * x + d * y for x, y in zip(r, s)] for r, s in zip(a, b)]


def running_example():
    n = range(3)
    eye = [[Fr(int(i == j)) for j in n] for i in n]
    s = [[Fr(2 if i == j else 1) for j in n] for i in n]
    e1 = [[Fr(1, 3) for j in n] for i in n]
    e2 = comb(1, eye, -1, e1)
    zero = [[0] * 3 for i in n]
    assert mul(e1, e1) == e1 and mul(e2, e2) == e2
    assert mul(e1, e2) == zero
    assert comb(4, e1, 1, e2) == s
    assert comb(Fr(1, 3), s, Fr(-1, 3), eye) == e1
    assert comb(Fr(-1, 3), s, Fr(4, 3), eye) == e2
    power = eye
    for _ in range(10):
        power = mul(power, s)
    assert power == comb(4 ** 10, e1, 1, e2)
    for a in (2, -2):
        for b in (1, -1):
            r = comb(a, e1, b, e2)
            assert mul(r, r) == s
    assert comb(2, e1, 1, e2) == comb(1, eye, 1, e1)
    u, w = [1, -1, 0], [1, 0, -1]
    dot = sum(x * y for x, y in zip(u, w))
    v = [y - Fr(dot, 2) * x for x, y in zip(u, w)]
    assert dot == 1 and [2 * x for x in v] == [1, 1, -2]
    for x in ([1, 1, 1], u, [1, 1, -2]):
        lam = 4 if x == [1, 1, 1] else 1
        assert [sum(p * q for p, q in zip(row, x)) for row in s] == [
            lam * t for t in x]
    print('S = 4 E1 + E2; S^10 entry (1,1) =', power[0][0])


def complex_examples():
    c = np.array([[1, -1], [1, 1]])
    assert np.array_equal(c @ c.T, c.T @ c)
    p = np.array([[1, 1j], [-1j, 1]]) / 2
    q = p.conj()
    assert np.array_equal((1 + 1j) * p + (1 - 1j) * q, c)
    assert np.array_equal(p @ p, p) and np.array_equal(p @ q, 0 * p)
    for sign in (1, -1):
        v = np.array([sign * 1j, 1])
        assert np.array_equal(c @ v, (1 + sign * 1j) * v)
    a = np.array([[2, 1 + 1j], [1 - 1j, 3]])
    eye = np.eye(2)
    e1, e4 = (4 * eye - a) / 3, (a - eye) / 3
    assert np.allclose(e1 + 4 * e4, a)
    assert np.allclose(e1 @ e1, e1) and np.allclose(e1 @ e4, 0)
    assert np.allclose(3 * e1, [[2, -1 - 1j], [-1 + 1j, 1]])
    m = np.array([[1, 1], [0, 2]])
    assert not np.array_equal(m @ m.T, m.T @ m)
    assert np.array_equal(m @ [1, 1], [2, 2])
    x = np.array([0, 1])
    assert np.vdot(m @ x, m @ x) != np.vdot(m.T @ x, m.T @ x)
    assert np.vdot(c @ x, c @ x) == np.vdot(c.T @ x, c.T @ x)
    print('complex resolutions and normality tests verified')


def main():
    running_example()
    complex_examples()
    print('all chapter 25 checks passed')


if __name__ == '__main__':
    main()
