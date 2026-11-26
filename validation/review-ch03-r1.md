# Review ch03, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; captions at positions 5 and 6; lab run; pipeline check.

Checked and correct: Theorems 2.2, 3.1, 3.2, 4.1, 5.2, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, Propositions 9.2 and 9.3, the root-count induction, all running-example numbers, the finite-field counts against lab output, section timestamps, caption topic coverage, and solutions 2 to 6.

## Findings

### F1 (blocking: gap that makes the proof circular), solution to Exercise 1

Quoted: "Cancel $u$ on the left and $v$ on the right; cancellation uses only (A1), (A2) and (A3), as in the proof of Theorem 2.2(1) (the right-cancellation form is proved there by adding $-w$ on the right, without (A4))."

Cancelling $u$ on the LEFT is needed. The proof of Theorem 2.2(1) obtains left cancellation "by (A4)", which is the statement being proved. Axioms (A2) and (A3) give only $u+\mathbf 0=u$ and $u+(-u)=\mathbf 0$. Left cancellation needs $(-u)+u=\mathbf 0$ and $\mathbf 0+u=u$, and neither is proved.

Proposed replacement for the quoted sentence (and place it before "This leaves $v+u=u+v$"):

"Cancellation on both sides is needed, and (A4) may not be used. First we show that $-u$ is also a left inverse. Put $z=(-u)+u$. By (A1), (A3) and (A2), $z+z=(-u)+\big(u+(-u)\big)+u=\big((-u)+\mathbf 0\big)+u=(-u)+u=z$. Adding $-z$ on the right gives $z+\big(z+(-z)\big)=z+(-z)$, that is, $z=\mathbf 0$. Next, $\mathbf 0$ is also a left identity: $\mathbf 0+u=\big(u+(-u)\big)+u=u+\big((-u)+u\big)=u+\mathbf 0=u$. Now cancel $u$ on the left by adding $-u$ on the left, and cancel $v$ on the right by adding $-v$ on the right."

## Advisories

None.

## Recheck

Accepted. F1 repaired: the solution to Exercise 1 now proves that -u is a left inverse and 0 a left identity from (A1)-(A3) before cancelling. Each step was checked. Check passes.
