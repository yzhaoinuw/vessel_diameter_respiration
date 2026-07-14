# Archive

This directory preserves analyses and documents that are not part of the retained workflow in `correlation_analysis_writeup.md`.

## Contents

- `docs/` — superseded pipeline guidance, planning notes, methodology drafts, and the earlier coupling report.
- `scripts/respiration/` — respiration preprocessing, coupling tests, sketches, and diagnostics.
- `scripts/legacy_coupling/` — early HR/respiration/vessel coupling utilities.
- `scripts/f153_prototypes/` — F153-only fine-resolution extraction and plotting prototypes.
- `scripts/superseded_hr_vessel/` — earlier five-minute, sliding, cross-correlation, and diagnostic workflows.
- `scripts/alternative_correlation_workflows/` — Pearson and non-Yue Spearman variants, with their local generated outputs preserved but ignored by Git.
- `scripts/hr_quality_superseded/` — unfiltered HR comparison scripts replaced by confidence-filtered versions.

## Archive Rule

Archived code is retained for provenance and may rely on old file layouts, schemas, or output locations. Do not use it as the implementation path for the current write-up and do not update it alongside active analysis changes. If an archived analysis is revived, move it back into an active location, document the scientific reason, verify it against current data, and update the project overview and work log.

Git history retains the former root and directory locations of tracked files. Generated CSV, XLSX, MAT, and PNG artifacts within archived workflows remain ignored.
