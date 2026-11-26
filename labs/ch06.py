"""Chapter 6 laboratory: row space and computations with subspaces.

Exact arithmetic over Q (fractions.Fraction), plus one brute-force
count over the field F_2.  A matrix is a list of rows.  Rows and columns
are numbered from 0 in code and from 1 in the chapter.
"""
from fractions import Fraction
from itertools import product


def Q(rows):
    return [[Fraction(x) for x in row] for row in rows]


def matmul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))), Fraction(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(col) for col in zip(*A)]


def identity(m):
    return [[Fraction(int(i == j)) for j in range(m)] for i in range(m)]


def rref(A):
    """Return (R, P, pivots) with P*A = R, R row reduced echelon.

    Row operations are applied to [A | I]; the right block ends as P.
    pivots lists the leading-entry columns k_1 < ... < k_r.
    """
    m, n = len(A), len(A[0])
    M = [list(A[i]) + identity(m)[i] for i in range(m)]
    row, pivots = 0, []
    for col in range(n):
        piv = next((i for i in range(row, m) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[row], M[piv] = M[piv], M[row]
        c = M[row][col]
        M[row] = [x/c for x in M[row]]
        for i in range(m):
            if i != row and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f*y for x, y in zip(M[i], M[row])]
        pivots.append(col)
        row += 1
        if row == m:
            break
    R = [r[:n] for r in M]
    P = [r[n:] for r in M]
    return R, P, pivots


def is_rref(R):
    last = -1
    seen_zero = False
    pivots = []
    for row in R:
        lead = next((k for k, x in enumerate(row) if x != 0), None)
        if lead is None:
            seen_zero = True
            continue
        if seen_zero or row[lead] != 1 or lead <= last:
            return False
        last = lead
        pivots.append(lead)
    for k in pivots:
        if sum(1 for row in R if row[k] != 0) != 1:
            return False
    return True


def rank(A):
    return len(rref(A)[2])


def row_basis(A):
    R, _, piv = rref(A)
    return R[:len(piv)]


def coordinates(beta, R, pivots):
    """Coordinates of beta in the basis of nonzero rows of an RREF R."""
    c = [beta[k] for k in pivots]
    combo = [sum((c[i]*R[i][j] for i in range(len(c))), Fraction(0))
             for j in range(len(beta))]
    return c if combo == list(beta) else None


def null_basis(A):
    """Basis of {x : Ax = 0}: one vector per free column."""
    R, _, piv = rref(A)
    n = len(A[0])
    basis = []
    for f in (j for j in range(n) if j not in piv):
        v = [Fraction(0)]*n
        v[f] = Fraction(1)
        for i, k in enumerate(piv):
            v[k] = -R[i][f]
        basis.append(v)
    return basis


def in_row_space(A, beta):
    return rank(A + [list(beta)]) == rank(A)


def same_row_space(A, B):
    return row_basis(A) == row_basis(B)


def mat_vec(A, x):
    return [sum((a*b for a, b in zip(row, x)), Fraction(0)) for row in A]


def fp_matmul(A, B, p=2):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(len(B))) % p
                       for j in range(len(B[0]))) for i in range(len(A)))


def count_P_over_F2(A, R):
    """Count invertible P over F_2 with P*A = R by brute force."""
    m = len(A)
    cnt = 0
    for entries in product((0, 1), repeat=m*m):
        P = [entries[i*m:(i+1)*m] for i in range(m)]
        if fp_matmul(P, A) != tuple(tuple(r) for r in R):
            continue
        if invertible_f2(P):
            cnt += 1
    return cnt


def invertible_f2(P):
    m = len(P)
    M = [list(r) + [int(i == j) for j in range(m)]
         for i, r in enumerate(P)]
    for col in range(m):
        piv = next((i for i in range(col, m) if M[i][col]), None)
        if piv is None:
            return False
        M[col], M[piv] = M[piv], M[col]
        for i in range(m):
            if i != col and M[i][col]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[col])]
    return True


def main():
    # ---- running example: the 5 x 5 matrix of Section 5 -------------
    A = Q([[1, 2, 0, 3, 0],
           [1, 2, -1, -1, 0],
           [0, 0, 1, 4, 0],
           [2, 4, 1, 10, 1],
           [0, 0, 0, 0, 1]])
    R, P, piv = rref(A)
    assert is_rref(R) and matmul(P, A) == R
    assert piv == [0, 2, 4]
    assert R == Q([[1, 2, 0, 3, 0], [0, 0, 1, 4, 0], [0, 0, 0, 0, 1],
                   [0]*5, [0]*5])
    print('pivot columns (1-based):', [k + 1 for k in piv])
    print('P from this lab:', [[str(x) for x in r] for r in P])

    # the matrix P of Section 5 (one of many valid choices)
    P5 = Q([[1, 0, 0, 0, 0],
            [1, -1, 0, 0, 0],
            [0, 0, 0, 0, 1],
            [-3, 1, 0, 1, -1],
            [-1, 1, 1, 0, 0]])
    assert matmul(P5, A) == R
    assert rank(P5) == 5
    print('P of Section 5 reproduces R; rank P =', rank(P5))

    # ---- row space, basis, coordinates -----------------------------
    assert rank(A) == 3
    basis = R[:3]
    # beta = (b1, 2 b1, b3, 3 b1 + 4 b3, b5) with b1=2, b3=3, b5=5
    beta = [Fraction(2), Fraction(4), Fraction(3), Fraction(18),
            Fraction(5)]
    assert in_row_space(A, beta)
    c = coordinates(beta, R, piv)
    assert c == [2, 3, 5]
    print('coordinates of', [int(x) for x in beta], 'are',
          [int(x) for x in c])
    bad = [Fraction(x) for x in (1, 0, 0, 0, 0)]
    assert not in_row_space(A, bad) and coordinates(bad, R, piv) is None
    # every row of A lies in the row space of R and conversely
    assert all(in_row_space(R[:3], row) for row in A)
    assert all(in_row_space(A, row) for row in R[:3])

    # ---- null space --------------------------------------------------
    N = null_basis(A)
    assert N == [Q([[-2, 1, 0, 0, 0]])[0], Q([[-3, 0, -4, 1, 0]])[0]]
    assert all(mat_vec(A, v) == [0]*5 for v in N)
    assert len(N) == 5 - rank(A)
    print('null space basis:', [[int(x) for x in v] for v in N])

    # ---- solvability conditions for Ax = y ---------------------------
    z4, z5 = P5[3], P5[4]
    print('conditions: z4 =', [int(x) for x in z4], ' z5 =',
          [int(x) for x in z5])
    # column space = {y : z4.y = 0 and z5.y = 0}
    cols = transpose(A)
    for col in cols:
        assert sum(a*b for a, b in zip(z4, col)) == 0
        assert sum(a*b for a, b in zip(z5, col)) == 0
    y_good = mat_vec(A, [Fraction(x) for x in (1, 1, 1, 1, 1)])
    assert sum(a*b for a, b in zip(z4, y_good)) == 0
    y_bad = [Fraction(x) for x in (1, 0, 0, 0, 0)]
    assert sum(a*b for a, b in zip(z5, y_bad)) != 0
    # row rank = column rank, and the left null space has dimension 2
    assert rank(transpose(A)) == rank(A) == 3
    assert len(null_basis(transpose(A))) == 5 - 3
    for z in null_basis(transpose(A)):
        assert all(x == 0 for x in
                   [sum(a*b for a, b in zip(z, col)) for col in cols])
    print('row rank = column rank =', rank(A))

    # ---- P is not unique: P' = Q P with Q R = R ---------------------
    Qm = identity(5)
    Qm[0][4] = Fraction(7)      # add 7 * (zero row 5 of R) to row 1
    assert matmul(Qm, R) == R and Qm != identity(5)
    P2 = matmul(Qm, P5)
    assert matmul(P2, A) == R and P2 != P5 and rank(P2) == 5
    # swapping the last two rows of P5 also works
    P3 = [r[:] for r in P5]
    P3[3], P3[4] = P3[4], P3[3]
    assert matmul(P3, A) == R

    # ---- P is unique when the rows of A are independent -------------
    C = Q([[1, 0, 2, 1], [0, 1, 1, 3], [1, 1, 0, 5]])
    assert rank(C) == 3
    assert len(null_basis(transpose(C))) == 0   # unique P

    # ---- count of P over F_2: |F|^(r(m-r)) |GL_{m-r}| = 4 -----------
    A2 = [[1, 0, 1], [0, 1, 1], [1, 1, 0]]    # rank 2 over F_2
    R2 = [[1, 0, 1], [0, 1, 1], [0, 0, 0]]
    assert fp_matmul(A2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == \
        tuple(tuple(r) for r in A2)
    n2 = count_P_over_F2([[x for x in r] for r in A2], R2)
    assert n2 == 2**(2*1) * 1
    print('invertible P over F_2 with PA = R:', n2)

    # ---- boundary: non-invertible P gives strict inclusion ----------
    I2 = Q([[1, 0], [0, 1]])
    P1 = Q([[1, 0]])
    Bm = matmul(P1, I2)
    assert in_row_space(I2, Bm[0]) and rank(Bm) == 1 < rank(I2)
    # boundary: column operations change the row space
    E = Q([[1, 0], [0, 0]])
    Ecol = Q([[0, 1], [1, 0]])
    EQ = matmul(E, Ecol)
    assert rank(EQ) == rank(E) and not same_row_space(E, EQ)

    # ---- same row space <=> same RREF <=> row equivalent -------------
    A3 = Q([[1, 1, 0], [0, 1, 1]])
    B3 = Q([[1, 0, -1], [1, 2, 1]])
    assert same_row_space(A3, B3)
    assert rref(A3)[0] == rref(B3)[0]
    # RREF with the same row space is unique (pivot set = K(W))
    print('RREF of A3 =', [[str(x) for x in r] for r in rref(A3)[0]])

    # ---- implicit versus parametric description; intersection -------
    # W = row space of A:  b2 = 2 b1, b4 = 3 b1 + 4 b3
    eqs = Q([[2, -1, 0, 0, 0], [3, 0, 4, -1, 0]])
    for row in A:
        assert mat_vec(eqs, row) == [0, 0]
    assert len(null_basis(eqs)) == 3
    assert same_row_space(A, Q([list(v) for v in null_basis(eqs)]))
    # U = span{(1,0,0,0,0),(0,1,0,0,0),(0,0,0,0,1)}  (original)
    U = Q([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]])
    # W cap U: vectors c.U with eqs satisfied
    M = matmul(eqs, transpose(U))          # 2 x 3
    sol = null_basis(M)
    inter = [[sum((ci*u for ci, u in zip(cvec, col)), Fraction(0))
              for col in zip(*U)] for cvec in sol]
    assert len(inter) == 1
    assert inter[0] == Q([[0, 0, 0, 0, 1]])[0]
    dimsum = rank(A + U)
    assert dimsum == 3 + 3 - len(inter)
    print('dim(W) + dim(U) - dim(W cap U) =', 3, '+', 3, '-', len(inter),
          '=', dimsum)

    # ---- example of Section 3 (5 columns, 3 pivots) ------------
    R9 = Q([[1, 0, 3, 0, 2], [0, 1, 2, 0, -1], [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0]])
    assert is_rref(R9)
    c1, c2, c3 = Fraction(2), Fraction(-1), Fraction(5)
    vec = [c1*R9[0][j] + c2*R9[1][j] + c3*R9[2][j] for j in range(5)]
    assert coordinates(vec, R9, [0, 1, 3]) == [c1, c2, c3]
    print('R9 combination:', [int(x) for x in vec])
    print('all assertions passed')


if __name__ == '__main__':
    main()
