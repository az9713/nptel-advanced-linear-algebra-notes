# Review ch02, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; caption at position 4; lab run; pipeline check.

Checked and correct: Theorems 3.1, 5.1, 5.2, 6.1, 6.2, 8.1, Corollaries 8.2 and 8.3, the nine-operation table, R, P, Z, both certificates, and all six solutions (recomputed by hand). Section timestamps agree with caption topics.

## Findings

### F1 (blocking: false claim), Section 7, first sentence

Quoted: "The smallest inconsistent system with nonzero coefficients is $x_1+2x_2=1$, $x_1+2x_2=4$."

This is false. The system $x_1=1$, $x_1=4$ is inconsistent, has nonzero coefficients, and has one unknown instead of two. No minimality is proved or needed.

Proposed replacement: "A simple inconsistent system in two unknowns is $x_1+2x_2=1$, $x_1+2x_2=4$."

## Advisories

1. The word "invertible" is used from the conventions paragraph on and is defined only in Exercise 1. Add the definition to the conventions paragraph.
2. `labs/ch02.py`, docstring of `rref_from_kernel`: "proof of Thm 6.1" should read "proof of Theorem 8.1".
3. References: the Axler chapter number is a weak pointer for echelon forms; consider dropping the chapter for that entry.

## Recheck

Accepted. F1 repaired: Section 7 now reads 'A simple inconsistent system in two unknowns'. The invertibility definition and lab docstring advisories were also applied. Check passes.
