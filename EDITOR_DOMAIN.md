# Domain: Advanced Linear Algebra

Audience: a master's student in mathematics or physics, or an engineering senior, who can compute with matrices and now needs proof-level command of vector spaces, linear maps, duality, canonical decompositions, and inner-product spaces. Manuscripts are `chapters/chNN.md` with metadata `chapters/chNN.json`; laboratories are `labs/chNN.py`. The outline and notation are in `plan/`. Recording identities and caption provenance are in `source-manifest.json`. Raw captions stay in the ignored `work/acquisition/` folder and are never published.

Apply the science-editor five-part standard (precise statement, every term defined, proof in the smallest setting, rescue or limit case, tie to something concrete), the rigorous-explainer prose rules, and the transcript-grounded coverage contract. Each chapter needs one current independent acceptance bound to the exact SHA-256 digests of its Markdown, JSON, and lab. A reviewer never accepts a chapter whose canonical author is that reviewer. Run the lab; do not infer success from plausible output. Review every full exercise solution.

## Distinctions a chapter must never blur

- A field-dependent statement versus a field-free one: say where a result needs $\mathbb{C}$, an algebraically closed field, a split characteristic polynomial, or characteristic not 2.
- Finite versus infinite dimension: dual-basis, double-dual, and rank-nullity statements need finite dimension; say so and give the infinite-dimensional failure when it is instructive.
- A linear map versus its matrix: a matrix depends on ordered bases. Similar matrices represent one operator.
- Row equivalence versus similarity versus equivalence ($PAQ$).
- Algebraic versus geometric multiplicity; characteristic versus minimal polynomial; diagonalizable versus triangulable.
- A direct sum versus a sum; a projection versus an orthogonal projection; an invariant subspace versus a reducing pair.
- Transpose of a linear map ($T^{\mathsf t}$ on duals) versus adjoint on an inner-product space ($T^*$); they agree with matrix transpose/conjugate transpose only in suitable (orthonormal / dual) bases.
- Unitary versus normal versus self-adjoint; real orthogonal diagonalization of symmetric matrices versus complex unitary diagonalization of normal matrices.
- Inner-product convention: linear in the first argument, conjugate-linear in the second. State it where used.
- Singular values versus eigenvalues; forms (sesquilinear/bilinear) versus operators, and the role of the inner product in identifying them.

## Evidence rules

- Captions do not carry the board matrices. Every worked numerical example is an original construction and must be labeled as such or presented as the chapter's own example. Do not claim the lecturer used a specific matrix unless the caption text states it.
- Exact arithmetic: labs use `fractions.Fraction` (stdlib) for row reduction, ranks, null spaces, minimal polynomials, and decompositions. `numpy` is allowed for floating-point illustrations (SVD, spectral checks) with explicit tolerances. Do not add sympy or other dependencies; `requirements.txt` is pinned.
- Primary references: the course lists Steven Roman, *Advanced Linear Algebra* (3rd ed., Springer GTM 135) and K. B. Datta, *Matrix and Linear Algebra* (PHI). The lecture sequence closely follows Hoffman and Kunze, *Linear Algebra* (2nd ed., Prentice-Hall, 1971). The books were not retrieved for this edition, so cite them by chapter only, never by section, theorem, or page number. Online sources (e.g. Axler's *Linear Algebra Done Right*, 4th ed., open access at https://linear.axler.net/) follow the same chapter-level rule.
- No claim of human expert certification, formal proof verification, or a comprehensive audiovisual audit.
