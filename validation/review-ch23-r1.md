# Review of ch23 (r1, author_a)

Disposition: changes_requested.

## Scope

Manuscript, JSON, ledger, lab (run, passes), pipeline check OK, section timestamps against captions. Theorems 1.1, 5.2, 6.1, 9.1, Propositions 4.1 and 9.2, the general-basis formula, the infinite-dimensional arguments, all printed numbers and all six solutions were checked by hand.

## Blocking findings

1. Section 2: finite dimension is missing. Quoted: "Let $f\ne0$ and $N=\ker f$. Then $V=N\oplus N^{\perp}$ ... and $\dim N^{\perp}=\dim V-\dim N=1$." The conventions allow any dimension, and Section 4 itself gives a nonzero functional with no representing vector. Proposed replacement of the first words: "Let $V$ be finite-dimensional, let $f\ne0$ and put $N=\ker f$."
2. Proposition 9.2 is stated without hypothesis: "**Proposition 9.2.** $\ker T^*=(\operatorname{im}T)^{\perp}$." For a general operator $T^*$ need not exist (Section 8). Proposed replacement: "**Proposition 9.2.** If $T$ has an adjoint (for instance if $V$ is finite-dimensional), then $\ker T^*=(\operatorname{im}T)^{\perp}$."

## Advisories

- Section 3 uses the third recording before Section 4 returns to the second; the order is logical but not lecture order.
