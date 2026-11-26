"""Chapter 26 laboratory: singular values and forms.

Integer and Gaussian-integer identities are checked exactly (Python
complex is exact for these small values).  Where square roots occur the
check uses floating point with numpy.allclose.  The general theorems are
proved in the chapter.
"""
from fractions import Fraction as Fr

import numpy as np


def running_example():
    a = np.array([[3, 0], [4, 5]])
    assert np.array_equal(a.T @ a, [[25, 20], [20, 25]])
    for v, lam in (([1, 1], 45), ([1, -1], 5)):
        assert np.array_equal(a.T @ a @ v, lam * np.array(v))
    v = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    u = np.array([[1, 3], [3, -1]]) / np.sqrt(10)
    d = np.diag([3 * np.sqrt(5), np.sqrt(5)])
    assert np.allclose(u.T @ u, np.eye(2))
    assert np.allclose(u @ d @ v.T, a)
    assert np.allclose(np.linalg.svd(a)[1], np.diag(d))
    assert np.isclose(d[0, 0] * d[1, 1], 15)
    for x1, x2 in ((1, 0), (2, -3), (5, 7)):
        s2, t2 = Fr((x1 + x2) ** 2, 2), Fr((x1 - x2) ** 2, 2)
        y = a @ [x1, x2]
        q = 25 * x1 ** 2 + 40 * x1 * x2 + 25 * x2 ** 2
        assert int(y @ y) == q == 45 * s2 + 5 * t2
    print('running example: sigma^2 = 45, 5; U D V^T = A')


def rank_one_example():
    b = np.array([[1, 1 - 1j], [1 + 1j, 2], [1, 1 - 1j]])
    h = b.conj().T @ b
    assert np.array_equal(h, 4 * np.array([[1, 1 - 1j], [1 + 1j, 2]]))
    v1 = np.array([1, 1 + 1j])
    v0 = np.array([1 - 1j, -1])
    assert np.array_equal(h @ v1, 12 * v1)
    assert np.array_equal(h @ v0, 0 * v0)
    assert np.vdot(v1, v0) == 0
    assert np.array_equal(b @ v1, 3 * np.array([1, 1 + 1j, 1]))
    u1 = np.array([1, 1 + 1j, 1]) / 2
    x2 = np.array([-1 + 1j, 2, -1 + 1j]) / 4
    x3 = np.array([-1, 0, 1]) / 2
    e2, e3 = np.array([0, 1, 0]), np.array([0, 0, 1])
    assert np.allclose(x2, e2 - np.vdot(u1, e2) * u1)
    w = e3 - np.vdot(u1, e3) * u1
    w = w - np.vdot(x2, e3) / np.vdot(x2, x2) * x2
    assert np.allclose(x3, w)
    u = np.column_stack([u1, np.sqrt(2) * x2, np.sqrt(2) * x3])
    v = np.column_stack([v1, v0]) / np.sqrt(3)
    d = np.zeros((3, 2))
    d[0, 0] = 2 * np.sqrt(3)
    assert np.allclose(u.conj().T @ u, np.eye(3))
    assert np.allclose(u.conj().T @ b @ v, d)
    print('rank-one example: sigma = 2 sqrt 3, U completed')


def forms():
    # f(a, b) = <T a, b> in the basis (1,0), (1,1)
    t = np.array([[1, 2], [0, 3]])
    p = np.array([[1, 1], [0, 1]])
    pinv = np.array([[1, -1], [0, 1]])

    def f(x, y):
        return np.vdot(y, t @ x)

    m = np.array([[f(p[:, k], p[:, j]) for k in range(2)]
                  for j in range(2)])
    gram = p.T @ p
    op = pinv @ t @ p
    assert np.array_equal(m, gram @ op)
    assert np.array_equal(m, p.T @ t @ p)
    assert not np.array_equal(m, op)
    s = np.array([[0, 1j], [2, 1]])
    c, sn = Fr(3, 5), Fr(4, 5)
    basis = np.array([[c, -sn], [sn, c]], dtype=float)
    total = sum(np.vdot(basis[:, j], s @ basis[:, k])
                * np.conj(np.vdot(basis[:, j], t @ basis[:, k]))
                for j in range(2) for k in range(2))
    assert np.isclose(total, np.trace(s @ t.conj().T))
    print('form matrices and basis independence verified')


def exercises():
    n = np.array([[1, 1], [0, 0]])
    v = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    assert np.allclose(np.diag([np.sqrt(2), 0]) @ v.T, n)
    q = np.array([[5, -3], [-3, 5]])
    assert np.array_equal(q @ [1, 1], [2, 2])
    assert np.array_equal(q @ [1, -1], [8, -8])
    m = np.array([[0, -1], [1, 0]])
    assert np.array_equal(m @ [1, -1j], 1j * np.array([1, -1j]))
    assert np.array_equal(1j * m, (1j * m).conj().T)
    print('exercise data verified')


def main():
    running_example()
    rank_one_example()
    forms()
    exercises()
    print('all chapter 26 checks passed')


if __name__ == '__main__':
    main()
