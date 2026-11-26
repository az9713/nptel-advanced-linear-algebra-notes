# Book-wide notation and conventions

Every chapter must still define each symbol locally at first substantive use. This file fixes the choices so that chapters agree.

## Fields, spaces, maps

- $F$ is a field. When a result needs more, say so: $\mathbb{R}$ or $\mathbb{C}$ for inner products; algebraically closed (or a split characteristic polynomial) for triangulation; characteristic zero only when it is actually used.
- $V, W, U$ are vector spaces over $F$. $\dim V$ is dimension. Finite dimension is stated as a hypothesis, not assumed silently.
- $F^n$ holds column vectors. $F^{m\times n}$ is the space of $m\times n$ matrices. $A^{\mathsf T}$ is transpose, $A^{*}=\overline{A}^{\mathsf T}$ is conjugate transpose.
- $\mathcal{L}(V,W)$ is the space of linear maps $V\to W$; $\mathcal{L}(V)=\mathcal{L}(V,V)$ holds linear operators. $I$ (or $I_V$) is the identity.
- $\ker T$ (null space) and $\operatorname{im} T$ (range). $\operatorname{rank} T=\dim\operatorname{im}T$, $\operatorname{nullity} T=\dim\ker T$.
- An ordered basis is written $\mathcal{B}=(v_1,\dots,v_n)$. Coordinates are $[v]_{\mathcal{B}}\in F^n$. The matrix of $T:V\to W$ relative to $\mathcal{B}$ and $\mathcal{C}$ is $[T]_{\mathcal{C}\leftarrow\mathcal{B}}$ (columns are $[Tv_j]_{\mathcal{C}}$); for an operator write $[T]_{\mathcal{B}}$. Change of basis: $[v]_{\mathcal{C}}=P_{\mathcal{C}\leftarrow\mathcal{B}}[v]_{\mathcal{B}}$.
- $\operatorname{span}(S)$; $W_1+W_2$; direct sum $W_1\oplus W_2$ (internal unless stated).

## Duality

- $V^{*}=\mathcal{L}(V,F)$ is the dual space; $f\in V^*$ is a linear functional. The dual basis of $\mathcal{B}$ is $\mathcal{B}^*=(f_1,\dots,f_n)$ with $f_i(v_j)=\delta_{ij}$.
- $S^{0}=\{f\in V^*: f(s)=0\ \forall s\in S\}$ is the annihilator. $V^{**}$ is the double dual, $L_v(f)=f(v)$.
- $T^{\mathsf t}:W^*\to V^*$, $T^{\mathsf t}g=g\circ T$, is the transpose of $T:V\to W$.

## Polynomials and spectra

- $F[x]$ is the polynomial ring. $c_T(x)=\det(xI-T)$ is the characteristic polynomial (monic convention). $m_T(x)$ is the minimal polynomial. $p(T)$ is polynomial evaluation.
- $E_\lambda=\ker(T-\lambda I)$ is the eigenspace; algebraic and geometric multiplicities are named in words.
- The $T$-conductor of $v$ into $W$ is $S_T(v;W)=\{g\in F[x]: g(T)v\in W\}$.
- Primary decomposition: $m_T=p_1^{r_1}\cdots p_k^{r_k}$ with distinct monic irreducible $p_i$, $W_i=\ker p_i(T)^{r_i}$.

## Inner products

- $\langle\cdot,\cdot\rangle$ is linear in the **first** argument and conjugate-linear in the second (Hoffman and Kunze convention; the course's reference style). State this once per chapter that uses it. $\|v\|=\sqrt{\langle v,v\rangle}$.
- Standard inner product on $\mathbb{C}^n$: $\langle x,y\rangle=y^{*}x=\sum_i x_i\overline{y_i}$.
- $W^{\perp}$ is the orthogonal complement. $T^{*}$ is the adjoint: $\langle Tv,w\rangle=\langle v,T^{*}w\rangle$.
- Unitary: $T^*T=TT^*=I$; normal: $T^*T=TT^*$; self-adjoint (Hermitian): $T^*=T$.
- The spectral theorem for normal operators requires $\mathbb{C}$ (or real self-adjoint/symmetric operators over $\mathbb{R}$). Never state it over $\mathbb{R}$ for a general normal operator.

## Typesetting rules (tool-enforced in practice)

- Use `$...$` inline and `$$...$$` display. Do **not** use `\boldsymbol`, `\bm`, `\cancel`, or mathtools-only macros (`\coloneqq`, `\mathclap`, `\DeclarePairedDelimiter`): the bundled MathJax has no extension folder and the whole page fails to render. Use `\mathbf` for bold.
- Use `\operatorname{...}` for rank, nullity, im, span, tr.
- Matrices use `\begin{pmatrix}...\end{pmatrix}` or `bmatrix`. Write `<` and `>` in math as `\lt` and `\gt` when adjacent to letters, to protect the HTML.
