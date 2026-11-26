# Review ch12, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory output, source ledger, six exercises and solutions; captions at positions 24 and 25; lab run; pipeline check.

Checked and correct: Theorems 2.2, 2.3, 3.1, 6.1, 6.2, 7.2, 7.3, the two infinite-dimensional examples, all eigenvectors, P and P^{-1}, the growth ratios at k=10 and 20 (recomputed), and all six solutions.

## Findings

### F1 (blocking: unproved central claim), Section 8

Quoted: "The matrix $A$ has no such formula. It has only one eigendirection for the double eigenvalue $2$, and its powers contain terms $k\,2^k$, which explains the linear growth observed at the start. [Chapter 20](ch20.md) derives this shape from the decomposition of an operator into a diagonalizable and a nilpotent part."

This answers the opening question for $A$ by assertion and a forward reference. One line of eigenvectors does not by itself show a $k\,2^k$ term. The proof is short.

Proposed replacement:

"The matrix $A$ behaves differently. Let $u=(1,1,2)$, the eigenvector for $2$, and $w=(1,1,1)$. Direct multiplication gives $Aw=(3,3,4)=2w+u$. By induction, $A^kw=2^kw+k\,2^{k-1}u$: it holds for $k=1$, and $A^{k+1}w=2^kAw+k\,2^{k-1}Au=2^{k+1}w+(k+1)2^ku$. Hence $A^kw/2^k=w+\tfrac k2u$, which grows linearly in $k$ because $u\neq0$. So the entries of $A^k$ divided by $2^k$ cannot stay bounded, and no formula of the form $a+b\,2^k$ holds for all of them. [Chapter 20](ch20.md) explains the general shape through a diagonalizable and a nilpotent part."

## Advisories

None.

## Recheck

Accepted. F1 repaired: Section 8 now proves A^k w = 2^k w + k 2^(k-1) u by induction. Verified Aw = 2w+u and the induction step. Check passes.
