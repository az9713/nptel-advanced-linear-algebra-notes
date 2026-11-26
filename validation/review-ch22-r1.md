# Review of ch22 (r1, author_a)

Disposition: changes_requested.

## Scope

Manuscript, JSON, ledger, lab (run, passes), pipeline check OK, section timestamps against captions. Theorems 2.1, 5.2, 6.1, Proposition 4.2, the infinite-dimensional example, all printed numbers and all six solutions were checked by hand.

## Blocking findings

1. Section 8, last sentence: under-hypothesized claim. Quoted: "Completeness of the space, a topic of functional analysis, is what restores existence." Completeness of $V$ alone does not give best approximations: the subspace must also be closed (in the complete space of square-summable sequences, the non-closed subspace of finitely supported sequences gives no best approximation to a vector outside it). Proposed replacement: "Existence is restored when the space is complete and the subspace $W$ is closed; both notions belong to functional analysis and are not developed here."

## Advisories

- Section 8 links back to the first recording after the sections on the second, and Section 5 carries inline timestamps 23:11 and then 12:08.
- Section 3: the existence of an orthogonal basis follows from Gram-Schmidt (Chapter 21), not from Theorem 2.1.
