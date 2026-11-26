# Review ch04, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; captions at positions 7 and 8; lab run; pipeline check.

Verdict: accepted. Proofs of Proposition 2.2, Theorems 5.1, 6.2, 7.1, 8.1, Corollaries 5.2 and 6.3, Lemma 6.1 and Example 3.4 were checked. Running-example numbers and all six solutions were recomputed and agree. Caption topics are covered.

## Findings

None blocking.

## Advisories

1. Section 6, last paragraph: "Both hypotheses of the corollary matter in infinite dimension." Corollary 6.3 has one standing hypothesis. Proposed replacement: "The hypothesis of finite dimension matters in part (c)."
2. Section 9: "every minimal spanning set and every maximal independent set is a basis size" is garbled. Proposed replacement: "every minimal spanning set and every maximal independent set has three elements."
3. `labs/ch04.py`, section 10: the assertion `1 + 1 + 1 - ... == 3 != 2` is arithmetic on literals and tests nothing about the three lines; compute the pairwise intersection dimensions instead.
