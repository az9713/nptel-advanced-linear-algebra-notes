# Review ch08, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; captions at positions 14, 15 and 16; lab run; pipeline check.

Verdict: accepted. Theorems 1.2, 2.1, 3.1, 4.1, 5.1, 6.2, 7.2, 7.3, 8.2 and Corollary 8.3 were checked. All printed matrices, the inverse of C, the one-sided inverse G and all six solutions were recomputed. Section timestamps and caption topics agree.

## Findings

None blocking.

## Advisories

1. Section 5: "bilinear" is not defined. Proposed replacement: "A vector space with an associative multiplication that is linear in each factor, as in parts 3 and 4, and has an identity is called a *linear algebra* over $F$ with identity."
2. Section 6, one-sided inverses: "Thus $T$ has a right inverse and no left inverse" does not follow from $GT\neq I_3$ alone. Proposed replacement: "Thus $T$ has a right inverse. It has no left inverse, because $ST=I$ would force $T$ to be injective."
3. `labs/ch08.py`, last block comment "Exercise 2: matrix of 2T - U" matches no exercise.
