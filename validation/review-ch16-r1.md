# Review of ch16 (r1, author_a)

Disposition: changes_requested.

## Scope

Manuscript, JSON, ledger, lab (run, passes), pipeline check OK, all section timestamps against captions. Proofs of Lemmas 1.1, 1.2, 8.1 and Theorems 2.1, 5.1, 8.2, Proposition 6.1, all printed matrices and bases, and all six solutions were checked by hand; no mathematical error found.

## Blocking findings

1. Section 1: one term is used in two senses and the second sense is never defined. Quoted: "The **$T$-conductor of $\beta$ into $W$** is the set", followed by "A conductor of degree one, $g_\beta=x-c$", "The conductor is $g_\beta=x^2$", and in Section 4 "the conductor of $e_1$ is $x^3$". Chapter 15, Definition 8.2, gives the name to the monic polynomial. Proposed addition after the proof of Lemma 1.1: "As in [Chapter 15](ch15.md), the monic generator $g_\beta$ is also called the *conductor* of $\beta$ into $W$; below, a conductor with a degree always means this polynomial."
2. Sections 6 and 7 and Exercises 4 and 5 use "$\operatorname{char}F$", "characteristic $0$" and "characteristic $p$" without defining the characteristic of a field. Proposed addition before first use in Section 6: "The *characteristic* of $F$, written $\operatorname{char}F$, is the least positive integer $p$ for which the sum of $p$ copies of $1$ equals $0$; it is $0$ if there is no such integer. When positive, it is a prime."

## Advisories

- Introduction: "It is the worked matrix of the first recording used here (the captions ... read out its entries)" narrates the source; the house rule asks for no narration about the lecture. Delete the sentence and keep the timestamp on Section 4.
- Section 2: the parenthesis after Corollary 2.2 states Schur's theorem without proof; say that it is not proved here.
- Section 7, last sentence: triangulation is defined only in finite dimension; say that no flag of invariant subspaces can start.
