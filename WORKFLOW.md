# Production, review, and publication workflow

This independent edition follows the transcript-grounded NPTEL/Stanford textbook method used for the earlier editions in this series. The source recordings establish the topic and timing spine; the exposition independently develops definitions, complete proofs, worked constructions, primary-source references, executable laboratories, and fully reasoned exercise solutions.

1. Inventory the supplied playlist by stable video identity. Preserve playlist positions and official lecture numbers. The supplied playlist has a promo and 59 lectures; Lecture 44 is absent from it, so the official NPTEL outline (course 111107164) supplied its video, appended as position 61. Positions 46/47 and 54/55 are uploaded out of lecture order and are taught in lecture order inside their chapters.
2. Acquire English caption evidence (the NPTEL-verified manual YouTube track for 59 lectures; automatic captions for the promo and for Lecture 29, position 30, whose manual track carries the Lecture 22 text) with `yt-dlp`, without downloading video. Normalize with `tools/acquire_sources.py`, which records digests in `source-manifest.json`. Raw captions stay in the ignored `work/` folder.
3. Each chapter author reads its recordings' caption text in full and in order and records topics, clarifications, and original constructions in `validation/source-read-chNN.json`, bound to the caption digest.
4. Fix book-wide notation (`plan/notation.md`), the outline (`plan/outline.md`), the specification (`plan/book-spec.json`), and the review criteria (`EDITOR_DOMAIN.md`). Complete one pilot chapter through authorship, independent review, build, checks, and browser rendering before production.
5. Develop each chapter through coverage, prose, mathematics, research, implementation, exercise, parity, and validation passes. Put definitions before use; prove central results; interpret and test them on a running example; state boundaries with counterexamples.
6. Obtain one independent acceptance for each chapter's manuscript, metadata, and laboratory, bound to exact SHA-256 digests. Earlier editions in this series used two; this run was limited to two agents, so each chapter is reviewed only by the agent that did not author it. A blocking finding is repaired by the author and rechecked on the changed lines by the same reviewer.
7. Run all laboratories in fresh processes, check source coverage, timestamps, links, exercise structure, reading attestation, code style, and a source-free byte-identical rebuild (`tools/check.py`, `tools/audit_book.py`, `tools/release_check.py`).
8. Render every page offline in a browser at desktop and phone widths (`tools/browser_check.py`): math, navigation, solution disclosures, code copy, theme, print expansion, and the interactive laboratory.
9. Commit the public edition, push to a new repository, deploy with the checked GitHub Actions workflow, and verify the served files against the commit (`tools/verify_live.py`).

The chapter craft draws on the transcript-grounded-textbook and rigorous-explainer methods, with the science-editor claim standard. Book structure follows STEM-textbook dependency and notation practice.

## Evidence boundaries

Caption reading does not certify every board equation or spoken word. AI editorial review is not human expert certification or formal proof verification. Laboratory assertions establish their stated finite computations. A byte-identical rebuild establishes reproducibility, not mathematical correctness. Timestamps locate related discussion without attributing every added proof or example to the lecturer.
