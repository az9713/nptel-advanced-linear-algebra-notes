# Review ch06, r1 (reviewer: author_b)

Scope: complete manuscript, metadata, laboratory, source ledger, six exercises and solutions; captions at positions 10 and 11; lab run; pipeline check.

Checked and correct: Theorems 2.1, 3.1, 4.1, 5.1, 8.1, Corollaries 2.2 to 3.3, the row reduction, Z, P, all running-example numbers, and all six solutions.

## Findings

### F1 (blocking: false provenance claim), `chapters/ch06.json`, `validation/source-read-ch06.json`

Quoted (JSON `source_scope`): "the worked 5 by 5 matrix is the chapter's own". Quoted (ledger, position 11 `clarifications`): "Board entries are not in the caption; the matrix used is the chapter's own and matches the spoken conditions". Quoted (ledger, position 10 `original_constructions`): "small 4 x 5 example".

All three are false. The caption at position 11 reads out the rows of the running matrix: "one two zero three zero" (row 1), "one two minus one minus one zero" (row 2), "zero zero one four zero" (row 3), "two four one ten one" (row 4), "zero zero zero zero one" (row 5), which is the chapter's $A$. The caption at position 10 reads out "one zero three zero two", "zero one two zero minus one", "zero zero zero one zero" and a zero row, which is the "small example" of Section 3.

Proposed replacements:

- JSON `source_scope`: "The captions of Lectures 9 and 10 were read in full; the 5 by 5 running matrix and the small example of Section 3 are the matrices read out in those captions, and the extra theorems and exercises are original."
- Ledger position 11 `clarifications`: ["The entries of the 5 by 5 matrix are read out in the caption and are used as the running example"].
- Ledger position 10 `original_constructions`: [] (move "small 4 x 5 example" to `topics` as "echelon example of Section 3").
- `chapters/ch06.md`, opening, replace "using $A$ as the running example." with "using $A$, the worked matrix of the second recording, as the running example."
- `chapters/ch06.md`, Section 3: replace "**A small example.**" with "**A small example (from the first recording).**"

## Advisories

1. The `labs/ch06.py` comment "Lecture-9 style example" can then read "example of Section 3, from the recording".

## Author response

F1 is correct for the metadata and the ledger; both are repaired. The two proposed sentences in `chapters/ch06.md` were not added: the house rule asks for no narration about the lecture in the chapter text, and the chapter text claims no originality for either matrix.

## Recheck

Accepted. F1 repaired: JSON source_scope and ledger now attribute both matrices to the captions. The two proposed attribution phrases for the chapter body were not added; the body makes no originality claim, so this is not blocking. Check passes.
