## How to use this guide

The chapters follow the lecture order. The table records which earlier chapters supply the definitions and theorems that a chapter uses most. Read a chapter by first attempting its opening problem, then the theorems with their proofs, then the six exercises before opening the solutions. Run the chapter laboratory after the exercises and change its inputs.

## Four routes

**Linear systems to dimension (Chapters 1–6).** Row operations, row-reduced echelon form, vector spaces, bases, coordinates, and subspace computations. This route needs only matrix arithmetic and fields.

**Maps and duality (Chapters 7–11).** Rank and nullity, the algebra of linear maps, matrices relative to ordered bases, dual spaces, annihilators, and the transpose.

**Canonical forms (Chapters 12–20).** Eigenvalues, annihilating and minimal polynomials, Cayley–Hamilton, invariant subspaces, direct sums, primary decomposition, and the diagonalizable-plus-nilpotent splitting. The field matters throughout: record at each theorem whether the characteristic or minimal polynomial must split.

**Inner-product spaces (Chapters 21–26).** Cauchy–Schwarz, orthogonal projection, adjoints, unitary and normal operators, spectral resolution, singular values, and forms. This route can be started after Chapter 11.

## Dependency map

| Chapter | Topic | Builds on |
|---|---|---|
| [01](ch01.md) | Linear systems as one matrix equation, and the operations that preserve its solutions | — |
| [02](ch02.md) | Row-reduced echelon form: one canonical answer to every linear system | 1 |
| [03](ch03.md) | Vector spaces and subspaces: the axioms that make linear combinations work | 1 |
| [04](ch04.md) | Linear independence, basis, and the invariance of dimension | 2, 3 |
| [05](ch05.md) | Coordinates and change of ordered basis | 4 |
| [06](ch06.md) | Row space and computations with subspaces | 2, 4, 5 |
| [07](ch07.md) | Linear transformations, null space, range, and the rank-nullity theorem | 4 |
| [08](ch08.md) | The algebra of linear transformations: sums, composition, invertibility, and isomorphism | 7 |
| [09](ch09.md) | Matrix representation of linear transformations and similarity | 5, 8 |
| [10](ch10.md) | Linear functionals, dual spaces, and dual bases | 7, 9 |
| [11](ch11.md) | Annihilators, the double dual, and the transpose of a linear map | 10 |
| [12](ch12.md) | Eigenvalues, eigenvectors, and the characteristic polynomial | 9 |
| [13](ch13.md) | Diagonalizable operators and annihilating polynomials | 12 |
| [14](ch14.md) | The Cayley-Hamilton theorem and the minimal polynomial in computation | 13 |
| [15](ch15.md) | Invariant subspaces, conductors, and triangulability | 14 |
| [16](ch16.md) | Applications of invariant subspaces: diagonalizability criteria and simultaneous triangulation | 15 |
| [17](ch17.md) | Direct-sum decompositions and projections | 4, 7 |
| [18](ch18.md) | Invariant direct sums and block-diagonal operators | 15, 17 |
| [19](ch19.md) | The primary decomposition theorem | 14, 18 |
| [20](ch20.md) | Applications of primary decomposition: diagonalizable plus nilpotent parts | 16, 19 |
| [21](ch21.md) | Inner products, norms, and the Cauchy-Schwarz inequality | 3, 4 |
| [22](ch22.md) | Orthogonal projection, Gram-Schmidt, and best approximation | 21 |
| [23](ch23.md) | Linear functionals and adjoint operators on inner product spaces | 11, 22 |
| [24](ch24.md) | Inner-product-space isomorphisms and unitary operators | 23 |
| [25](ch25.md) | Normal operators and spectral resolution | 12, 24 |
| [26](ch26.md) | Singular value decomposition and forms on inner product spaces | 23, 25 |

## Three capstone questions

1. Take one real $4\times4$ matrix with a repeated eigenvalue. Compute its characteristic and minimal polynomials, its primary decomposition, and its diagonalizable-plus-nilpotent splitting. Say which steps would fail over a field where the polynomial does not split.
2. For one linear map between finite-dimensional inner-product spaces, compare the transpose on dual spaces with the adjoint. State the basis conditions under which their matrices are the transpose and the conjugate transpose.
3. Solve an inconsistent system in the least-squares sense twice: once by orthogonal projection and once by the singular value decomposition. Explain which solution has minimum norm when the columns are dependent.
