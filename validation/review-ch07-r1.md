# Review ch07, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; captions at positions 12 and 13; lab run; pipeline check.

Verdict: accepted. Theorems 3.2, 4.1, 5.1, 6.1, Propositions 3.3 and 4.2, Corollary 5.2 and the infinite-dimensional boundary were checked. Running-example numbers (condition a=b+c, kernel (-2,4,3), forced value (1,0,-1)) and all six solutions were recomputed. Caption topics are covered.

## Findings

None blocking.

## Advisories

1. `labs/ch07.py` docstrings cite "Theorem 7.4" for the prescription criterion; the manuscript calls it Proposition 4.2.
2. Section 7 timestamp [29:44]: at that point the caption is still finishing the row-rank proof; the worked example starts later.
3. Solution 4, last sentence: "the kernel vector looks different" follows a statement that the two agree. Proposed replacement: "The rank is the same as over $\mathbb{Q}$ here."
