# Review of ch20 (r1, author_a)

Disposition: changes_requested.

## Scope

Manuscript, JSON, ledger, lab (run, passes), pipeline check OK, section timestamps against captions. Lemma 2.2, Theorems 3.1, 4.1, 4.2, 5.2, 6.1, 7.1, 8.1, 9.1, Corollary 6.2, all printed matrices, the worked examples and all six solutions were checked by hand.

## Blocking findings

1. Section 6, undefined term and unproved step in the central proof. Quoted: "The residue classes $F[x]/(p)$ form a field $K$, and $[g]\beta=g(T)\beta$ makes $\ker p(T)$ a vector space over $K$." Neither "residue class" nor $F[x]/(p)$ nor $[g]$ is defined here or in an earlier chapter, and the field property is asserted. Theorem 6.1 depends on it. Proposed replacement: "For $g\in F[x]$ let $[g]$ be its *residue class* modulo $p$, the set of all $g+hp$ with $h\in F[x]$. Sums and products of classes, computed on any representatives, do not depend on the choice; the classes form a commutative ring $F[x]/(p)$. It is a field $K$: if $[g]\ne0$, then $p$ does not divide $g$, so $g$ and the irreducible $p$ are relatively prime, hence $ag+bp=1$ for some $a,b\in F[x]$ (the ideal argument of [Chapter 19](ch19.md)), and then $[a][g]=1$. The rule $[g]\beta=g(T)\beta$ is well defined on $\ker p(T)$ because $p(T)\beta=0$, and it makes $\ker p(T)$ a vector space over $K$."

## Advisories

- Section 6: "The recording states this result and leaves the proof to the listener. Here it is." narrates the source; delete it.
- Section 9: give the two lines that exclude a factorization of $y^4-8$ into two rational quadratics.
- Sections 5 to 8 do not follow the order of the recordings; the logical order is fine. The Section 9 timestamp 25:46 points a little before the rational form is stated.
