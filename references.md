# References and evidence limits

## Source course

- Premananda Bera, *Advanced Linear Algebra*, NPTEL, IIT Roorkee. [Official course page](https://nptel.ac.in/courses/111107164) · [Supplied playlist](https://www.youtube.com/playlist?list=PLLy_2iUCG87ChONRU6EeO6A2rIzM3HUvG). The [recording map](source-map.html) lists every video, its playlist position, its lecture number, and the chapter that develops it. Lecture 44 (Inner Product Spaces I) is absent from the supplied playlist; this edition uses the video that the official course outline lists for it.

## Books cited in the chapters

- Steven Roman, *Advanced Linear Algebra*, 3rd ed., Graduate Texts in Mathematics 135, Springer, 2008. The course lists it as a reference text.
- K. B. Datta, *Matrix and Linear Algebra*, Prentice Hall of India. The course lists it as a reference text.
- Kenneth Hoffman and Ray Kunze, *Linear Algebra*, 2nd ed., Prentice-Hall, 1971. The order of the course topics (row reduction, vector spaces, linear transformations, duality, elementary canonical forms, inner-product spaces, operators on inner-product spaces) follows this book closely; chapters cite it by chapter.
- Sheldon Axler, *Linear Algebra Done Right*, 4th ed., Springer, 2024, open access at [linear.axler.net](https://linear.axler.net/). Cited for determinant-free spectral arguments.

Chapter-specific citations appear beside the claims they support. The books above were not retrieved while preparing this edition, so citations point to a book or a chapter as further reading, not to verified section, theorem, or page numbers. Results used in the text are proved or stated with their hypotheses in the text itself.

## What the evidence in this edition establishes

- **Caption reading.** Chapter authors read the English caption text of each recording completely and in order and recorded the topics in `validation/source-read-chNN.json`, bound to the caption digest in `source-manifest.json`. Captions do not contain board writing, so the matrices and vectors of the worked examples are original constructions. No comprehensive audit of the video, boards, or audio was made.
- **Proofs and examples.** Proofs were written for the edition. Where a lecture states a result without proof, the chapter supplies one. One AI reviewer that did not write the chapter checked each chapter, its six worked solutions, and its laboratory; the receipt `validation/review-chNN-r1.json` binds the exact file digests it accepted. Earlier editions in this series used two reviewers. This review is not human expert certification or formal proof verification.
- **Laboratories.** Each `labs/chNN.py` checks the chapter's finite examples, mostly in exact rational arithmetic with Python's `fractions.Fraction`, and fails loudly if an assertion breaks. A passing laboratory establishes those finite computations only.
- **Reproducibility.** `tools/release_check.py` rebuilds every HTML page from the public inputs alone and requires byte-identical output. This establishes reproducibility, not correctness.
- **Timestamps.** A timestamp beside a heading locates related discussion in the recording. It does not attribute every proof, example, or extension in that section to the lecturer.
