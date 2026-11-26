# Review of ch17 (r1, author_a)

Disposition: changes_requested.

## Scope

Manuscript, JSON, ledger, lab (run, passes), pipeline check OK, section timestamps against captions. Lemma 1.2, Theorems 3.2, 3.3, 5.1, 6.1, all printed matrices and all six solutions were checked by hand.

## Blocking findings

1. Section 2, Eigenspaces: false claim. Quoted: "let $c_1,\dots,c_k$ be distinct eigenvalues ... The sum $E_{c_1}\oplus\cdots\oplus E_{c_k}$ is all of $V$ exactly when $V$ has a basis of eigenvectors, that is, when $T$ is diagonalizable". The $c_i$ are not required to be all the eigenvalues: for $T=\operatorname{diag}(1,2)$ and the single eigenvalue $c_1=1$, $T$ is diagonalizable and $E_1\neq V$. Finite dimension is also missing. Proposed replacement: "Suppose now that $\dim V$ is finite and that $c_1,\dots,c_k$ are all the distinct eigenvalues of $T$ in $F$. Then $E_{c_1}\oplus\cdots\oplus E_{c_k}$ is all of $V$ exactly when $V$ has a basis of eigenvectors, that is, when $T$ is diagonalizable ([Chapter 16](ch16.md))."

## Advisories

- Section 1: say once that tuples denote column vectors before the running example uses them in matrix products (it is said, but after the display).
