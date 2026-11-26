# Review ch14, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory output, source ledger, six exercises and solutions; captions at positions 28 and 29; lab run; pipeline check.

Checked and correct: Theorem 1.1 (adjugate proof), the 2 by 2 check, Theorems 5.1, 5.2 and 5.3, the AB/BA example, the matrix G with its rank, minimal polynomial and inverse, and all six solutions (G^2 and G^{-1} recomputed). Section timestamps and caption topics agree.

## Findings

### F1 (blocking: under-hypothesized claim), Section 6, last paragraph

Quoted: "A matrix whose characteristic polynomial has distinct roots in $F$ is diagonalizable over $F$; the converse is false, since the running example is diagonalizable with a repeated root $0$."

The hypothesis that the polynomial splits is missing. The chapter's own matrix $G$ over $\mathbb{R}$ has a characteristic polynomial whose only root in $\mathbb{R}$ is $2$, with no repetition, and $G$ is not diagonalizable over $\mathbb{R}$.

Proposed replacement:

"If $c_M$ for $M\in F^{n\times n}$ is a product of $n$ distinct linear factors over $F$, then $M$ is diagonalizable over $F$, since it has $n$ distinct eigenvalues ([Chapter 12](ch12.md)). It is not enough that the roots which lie in $F$ are simple: over $\mathbb{R}$ the matrix $G$ has only the simple root $2$ and is not diagonalizable. The converse is also false, since the running example is diagonalizable with the repeated root $0$."

## Advisories

None.

## Recheck

Accepted. F1 repaired: the criterion now requires n distinct linear factors over F. Check passes.
