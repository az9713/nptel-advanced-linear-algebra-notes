"""Chapter 24 laboratory: unitary matrices and triangular factorization.

Exact arithmetic over the Gaussian rationals: a number is a pair of
Fractions (class GQ).  The general theorems are proved in the chapter;
this file checks the printed instances.
"""
from fractions import Fraction as Fr
from math import isqrt


class GQ:
    """a + b i with rational a and b."""

    def __init__(self, re=0, im=0):
        self.re, self.im = Fr(re), Fr(im)

    @staticmethod
    def lift(x):
        return x if isinstance(x, GQ) else GQ(x)

    def __add__(self, o):
        o = GQ.lift(o)
        return GQ(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __sub__(self, o):
        o = GQ.lift(o)
        return GQ(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        o = GQ.lift(o)
        return GQ(self.re * o.re - self.im * o.im,
                  self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = GQ.lift(o)
        d = o.re ** 2 + o.im ** 2
        n = self * o.conj()
        return GQ(n.re / d, n.im / d)

    def __eq__(self, o):
        o = GQ.lift(o)
        return (self.re, self.im) == (o.re, o.im)

    def conj(self):
        return GQ(self.re, -self.im)


def mat(rows):
    return [[GQ.lift(x) for x in r] for r in rows]


def mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), GQ())
             for j in range(len(b[0]))] for i in range(len(a))]


def tr(a):
    return [list(r) for r in zip(*a)]


def adj(a):
    return [[x.conj() for x in r] for r in tr(a)]


def eye(n):
    return mat([[int(i == j) for j in range(n)] for i in range(n)])


def unitary(a):
    return mul(adj(a), a) == eye(len(a[0]))


def ip(x, y):
    return sum((a * b.conj() for a, b in zip(x, y)), GQ())


def lq(a):
    """Lower triangular M, positive diagonal, with M a unitary.  Raises
    if some Gram-Schmidt norm is irrational."""
    orth, comb, m = [], [], []
    n = len(a)
    for k in range(n):
        v, t = list(a[k]), [GQ(int(i == k)) for i in range(n)]
        for u, s in zip(orth, comb):
            c = ip(a[k], u) / ip(u, u)
            v = [x - c * y for x, y in zip(v, u)]
            t = [x - c * y for x, y in zip(t, s)]
        q = ip(v, v).re
        norm = Fr(isqrt(q.numerator), isqrt(q.denominator))
        assert norm ** 2 == q and norm > 0
        orth.append(v)
        comb.append(t)
        m.append([x / norm for x in t])
    return m


U = mat([[Fr(3, 5), GQ(0, Fr(4, 5))], [GQ(0, Fr(4, 5)), Fr(3, 5)]])


def running_example():
    assert unitary(U) and mul(U, adj(U)) == eye(2)
    for sign in (1, -1):
        lam = GQ(Fr(3, 5), sign * Fr(4, 5))
        v = mat([[1], [sign]])
        assert mul(U, v) == [[lam * v[0][0]], [lam * v[1][0]]]
        assert lam.re ** 2 + lam.im ** 2 == 1
    a = mat([[3, GQ(0, 4)], [0, 5]])
    m = lq(a)
    assert m == mat([[Fr(1, 5), 0], [GQ(0, Fr(4, 15)), Fr(1, 3)]])
    assert mul(m, a) == U
    assert mul(mat([[5, 0], [GQ(0, -4), 3]]), U) == a
    other = mul(mat([[1, 0], [0, -1]]), m)
    assert unitary(mul(other, a)) and other[1][1] == GQ(Fr(-1, 3))
    print('running example: unitary, eigenpairs, A = L U verified')


def other_matrices():
    f35, f45 = Fr(3, 5), Fr(4, 5)
    r = mat([[f35, -f45], [f45, f35]])
    p, pinv = mat([[1, 1], [0, 1]]), mat([[1, -1], [0, 1]])
    rb = mul(pinv, mul(r, p))
    assert rb == mat([[Fr(-1, 5), Fr(-8, 5)], [Fr(4, 5), Fr(7, 5)]])
    assert not unitary(rb)
    assert ip(tr(rb)[0], tr(rb)[0]) == Fr(17, 25)
    h = Fr(3, 4)
    c = mat([[Fr(5, 4), GQ(0, -h)], [GQ(0, h), Fr(5, 4)]])
    assert mul(tr(c), c) == eye(2) and not unitary(c)
    assert ip(tr(c)[0], tr(c)[0]) == Fr(34, 16)
    for w in (mat([[f35, f45], [-f45, f35]]),
              mat([[f35, f45], [f45, -f35]])):
        assert mul(tr(w), w) == eye(2)
    a, cc = GQ(f35), GQ(0, f45)
    for lam in (GQ(1), GQ(0, 1), GQ(f35, f45)):
        w = [[a, GQ() - lam * cc.conj()], [cc, lam * a.conj()]]
        assert unitary(w)
        assert w[0][0] * w[1][1] - w[0][1] * w[1][0] == lam
    a3 = mat([[f35, f45, 0], [0, 1, 0], [0, 0, 1]])
    m3 = lq(a3)
    assert m3 == mat([[1, 0, 0], [Fr(-4, 3), Fr(5, 3), 0], [0, 0, 1]])
    assert mul(m3, a3) == mat([[f35, f45, 0], [-f45, f35, 0], [0, 0, 1]])
    print('basis change, orthogonal forms, unitary family, 3x3 verified')


def shift_and_exercises():
    x = [2, 1, -3]
    assert ([0] + x)[1:] == x
    e1 = [1, 0, 0]
    assert [0] + e1[1:] != e1
    a = mat([[4, 3], [0, 5]])
    m = lq(a)
    assert m == mat([[Fr(1, 5), 0], [Fr(-3, 20), Fr(1, 4)]])
    assert mul(mat([[5, 0], [3, 4]]), mul(m, a)) == a
    h = mat([[0, GQ(0, 1)], [GQ(0, -1), 0]])
    assert unitary(h) and adj(h) == h
    assert Fr(1, 4) + Fr(1, 12) == Fr(1, 3)
    nil = mat([[0, 1], [0, 0]])
    assert mul(nil, adj(nil)) != mul(adj(nil), nil)
    print('shift identities and exercise data verified')


def main():
    running_example()
    other_matrices()
    shift_and_exercises()
    print('all chapter 24 checks passed')


if __name__ == '__main__':
    main()
