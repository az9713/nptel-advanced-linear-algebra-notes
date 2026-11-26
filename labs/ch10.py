"""Chapter 10 laboratory: linear functionals, dual bases, annihilators.

Exact arithmetic over Q with fractions.Fraction.  A vector in Q^n is a
list of n Fractions.  A linear functional on Q^n is stored by its row of
coefficients, so f(x) = sum_i c_i x_i.  A matrix is a list of rows.  Every
routine returns new objects and never changes its input.

The laboratory checks, in exact arithmetic and on one fixed example each:
  1. the dual basis of a basis of Q^3 is given by the rows of P^{-1};
  2. coordinates and functional expansions (Theorem 10.2);
  3. the trace functional and its dual basis on 2x2 matrices;
  4. evaluation functionals on polynomials of degree at most 2, and the
     dual basis they induce (the Lagrange polynomials);
  5. the boundary cases: repeated nodes, F_2, infinite dimension;
  6. dim W + dim W^0 = dim V and W as an intersection of hyperspaces;
  7. the converse problem: the common null space of given functionals.
It does not prove anything for general fields or dimensions.
"""
from fractions import Fraction


def Q(rows):
    """Matrix over Q from integer or Fraction entries."""
    return [[Fraction(x) for x in row] for row in rows]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fraction(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(col) for col in zip(*A)]


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def rref(A):
    """Return (R, pivots): reduced row echelon form and pivot columns."""
    R = [row[:] for row in A]
    m, n = len(R), len(R[0])
    pivots = []
    r = 0
    for c in range(n):
        if r == m:
            break
        p = next((i for i in range(r, m) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        piv = R[r][c]
        R[r] = [x / piv for x in R[r]]
        for i in range(m):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        pivots.append(c)
        r += 1
    return R, pivots


def rank(A):
    return len(rref(A)[1])


def nullspace(A):
    """Basis of {x : A x = 0} as a list of vectors."""
    R, pivots = rref(A)
    n = len(A[0])
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for f in free:
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, p in zip(R, pivots):
            x[p] = -row[f]
        basis.append(x)
    return basis


def inverse(P):
    """Inverse of a square matrix; raises ValueError if singular."""
    n = len(P)
    aug = [P[i][:] + identity(n)[i] for i in range(n)]
    R, pivots = rref(aug)
    if pivots[:n] != list(range(n)):
        raise ValueError('matrix is singular')
    return [row[n:] for row in R]


def apply(f, x):
    """Value of the functional with coefficient row f at the vector x."""
    return sum((a * b for a, b in zip(f, x)), Fraction(0))


def dual_basis(basis):
    """Coefficient rows (f_1, ..., f_n) with f_i(basis_j) = delta_ij.

    The basis vectors are the columns of P; the rows of P^{-1} are the
    functionals, because P^{-1} P = I.
    """
    P = transpose(basis)
    return inverse(P)


def coords(basis, x):
    """Coordinates of x relative to the ordered basis, via the dual basis."""
    return [apply(f, x) for f in dual_basis(basis)]


def poly_eval(p, t):
    """Value at t of the polynomial with coefficient list p = (p0, p1, ...)."""
    return sum((c * t ** k for k, c in enumerate(p)), Fraction(0))


def check_dual_basis_q3():
    alpha = Q([[1, 0, 1], [1, 1, 0], [0, 1, 1]])
    duals = dual_basis(alpha)
    print('dual basis of B in Q^3 (rows):')
    for row in duals:
        print('  ', [str(c) for c in row])
    for i, f in enumerate(duals):
        for j, a in enumerate(alpha):
            assert apply(f, a) == (1 if i == j else 0)
    v = Q([[3, 1, 2]])[0]
    c = coords(alpha, v)
    print('coordinates of (3,1,2):', [str(x) for x in c])
    assert c == [2, 1, 0]
    rebuilt = [sum(c[j] * alpha[j][k] for j in range(3)) for k in range(3)]
    assert rebuilt == v
    # Theorem 10.2: f = sum f(alpha_j) f_j for f(x, y, z) = 2x - y + 3z
    f = Q([[2, -1, 3]])[0]
    comp = [apply(f, a) for a in alpha]
    print('f = 2x - y + 3z has coordinates in B*:', [str(x) for x in comp])
    assert comp == [5, 1, 2]
    combo = [sum(comp[i] * duals[i][k] for i in range(3)) for k in range(3)]
    assert combo == f
    # The standard basis is its own dual basis: c_i = f(e_i).
    e = identity(3)
    assert dual_basis(e) == e
    assert [apply(f, ei) for ei in e] == f
    # Annihilator of W = span(alpha_1, alpha_2) is spanned by f_3.
    W = [alpha[0], alpha[1]]
    ann = nullspace(W)
    print('annihilator of span(alpha1, alpha2):',
          [[str(x) for x in r] for r in ann])
    assert len(ann) == 1
    assert ann[0] == [c * ann[0][0] / duals[2][0] for c in duals[2]]
    # A hyperplane is the kernel of one nonzero functional (rank 1).
    assert rank([duals[2]]) == 1 and len(nullspace([duals[2]])) == 2


def check_trace():
    # tr on 2x2 matrices, coordinates (a11, a12, a21, a22).
    tr = Q([[1, 0, 0, 1]])[0]
    A = Q([[1, 2, 3, 4]])[0]
    B = Q([[5, -1, 0, 7]])[0]
    c = Fraction(3)
    lin = [c * a + b for a, b in zip(A, B)]
    assert apply(tr, lin) == c * apply(tr, A) + apply(tr, B)
    # AB and BA need not be equal, but tr(AB) = tr(BA): not a functional
    # identity of one variable, so only the linear check is a claim here.
    ker = nullspace([tr])
    print('trace on 2x2 matrices: rank 1, nullity', len(ker))
    assert len(ker) == 3


def check_evaluation(nodes=(0, 1, 3)):
    t = [Fraction(x) for x in nodes]
    # Row i of the Vandermonde matrix is L_i applied to 1, x, x^2.
    V = [[ti ** k for k in range(3)] for ti in t]
    Vinv = inverse(V)
    print('Vandermonde inverse for nodes', nodes)
    for row in Vinv:
        print('  ', [str(x) for x in row])
    # Column j of Vinv holds the coefficients of p_j (L_i(p_j) = delta_ij).
    polys = [[Vinv[k][j] for k in range(3)] for j in range(3)]
    for i in range(3):
        for j in range(3):
            assert poly_eval(polys[j], t[i]) == (1 if i == j else 0)
    print('dual basis polynomials p_j (coefficients p0, p1, p2):')
    for p in polys:
        print('  ', [str(x) for x in p])
    # They are the Lagrange polynomials prod_{k != j} (x - t_k)/(t_j - t_k).
    for j in range(3):
        for x in [Fraction(-2), Fraction(1, 2), Fraction(5)]:
            lag = Fraction(1)
            for k in range(3):
                if k != j:
                    lag *= (x - t[k]) / (t[j] - t[k])
            assert poly_eval(polys[j], x) == lag
    # Expansion p = sum p(t_j) p_j for p = 2 - x + 4x^2.
    p = Q([[2, -1, 4]])[0]
    rebuilt = [sum(poly_eval(p, t[j]) * polys[j][k] for j in range(3))
               for k in range(3)]
    print('p(t_j) for p = 2 - x + 4x^2:',
          [str(poly_eval(p, tj)) for tj in t])
    assert rebuilt == p
    # det of the Vandermonde matrix is prod_{i<j} (t_j - t_i).
    det = Fraction(1)
    for i in range(3):
        for j in range(i + 1, 3):
            det *= t[j] - t[i]
    print('prod (t_j - t_i) =', det)
    assert det == 6 and rank(V) == 3
    # Boundary: a repeated node makes the evaluation functionals dependent.
    bad = [[x ** k for k in range(3)] for x in Q([[0, 1, 1]])[0]]
    assert rank(bad) == 2
    try:
        inverse(bad)
        raise AssertionError('singular matrix was inverted')
    except ValueError:
        pass
    # Boundary: over F_2 the functions 1, x, x^2 on {0, 1} have rank 2.
    rows = [[1, 1], [0, 1], [0, 1]]   # (p(0), p(1)) for 1, x, x^2 over F_2
    r = 0
    for c in range(2):
        p = next((i for i in range(r, 3) if rows[i][c] % 2), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(3):
            if i != r and rows[i][c] % 2:
                rows[i] = [(a + b) % 2 for a, b in zip(rows[i], rows[r])]
        r += 1
    print('rank of (1, x, x^2) as functions on F_2:', r)
    assert r == 2


def check_infinite():
    # On R[x], f_j(p) = coefficient of x^j are independent, but ev_1 is
    # not a finite combination: on x^(N+1) every f_j (j <= N) is zero
    # while ev_1(x^(N+1)) = 1.
    ok = 0
    for N in range(0, 21):
        mono = [Fraction(0)] * (N + 2)
        mono[N + 1] = Fraction(1)
        assert poly_eval(mono, Fraction(1)) == 1
        assert all(mono[j] == 0 for j in range(N + 1))
        ok += 1
    print('ev_1 separated from span(f_0..f_N) for N = 0..20:', ok, 'cases')


def check_annihilator():
    U = Q([[1, 2, 0, 1, 1],
           [0, 1, 1, -1, 2],
           [0, 0, 1, 1, -1],
           [1, 4, 1, -2, 6]])
    n = 5
    R, piv = rref(U)
    k = len(piv)
    print('RREF of the four spanning vectors of W (nonzero rows):')
    for row in R[:k]:
        print('  ', [str(x) for x in row])
    print('dim W =', k)
    ann = nullspace(U)
    print('dim W0 =', len(ann))
    assert k == 3 and len(ann) == n - k
    assert k + len(ann) == n
    print('basis of W0 (one functional per free column):')
    for row in ann:
        print('  ', [str(x) for x in row])
    # Every functional in ann kills every spanning vector.
    for f in ann:
        for u in U:
            assert apply(f, u) == 0
    # W equals the intersection of the two kernels: the common null space
    # of ann has dimension n - 2 = 3 = dim W and contains W.
    common = nullspace(ann)
    assert len(common) == k
    assert rank(common + U) == k
    # Fewer than n - k hyperspaces do not cut out W.
    assert len(nullspace([ann[0]])) == n - 1 != k
    # S^0 for a non-subspace S equals (span S)^0.
    S = [U[0], U[1]]
    assert len(nullspace(S)) == n - 2
    comb = [a + b for a, b in zip(U[0], U[1])]
    assert all(apply(f, comb) == 0 for f in nullspace(S))


def check_converse():
    # m = 3 functionals on Q^4 given by the rows of A; the third is the
    # sum of the first two, so r = 2 and dim W = 4 - 2 = 2.
    A = Q([[1, 2, 0, 1],
           [0, 1, 1, -1],
           [1, 3, 1, 0]])
    r = rank(A)
    W = nullspace(A)
    print('converse problem: r =', r, ', dim W =', len(W))
    for row in W:
        print('  ', [str(x) for x in row])
    assert r == 2 and len(W) == 4 - r
    for w in W:
        assert all(apply(f, w) == 0 for f in A)
    # Dropping the dependent functional changes nothing.
    assert len(nullspace(A[:2])) == len(W)
    assert rank(W + nullspace(A[:2])) == len(W)


def main():
    check_dual_basis_q3()
    check_trace()
    check_evaluation()
    check_infinite()
    check_annihilator()
    check_converse()
    print('all chapter 10 checks passed')


if __name__ == '__main__':
    main()
