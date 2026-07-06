# Work Log

## 2026-07-06

### Multiple-Comparison Caveat Clarification

- Expanded the Results caveat to explain that unadjusted probability values do not use a stricter threshold for multiple related feature-period tests.
- Clarified why the caveat matters for interpreting nominally significant findings as exploratory.
- Model: Codex, GPT-5.
- Effort: standard interactive editing.
- Token budget: not specified.
- Verification: reviewed the Results paragraph diff and checked branch status before staging.

## 2026-06-23

### All-Recording Sensitivity Write-Up

- Added the sliding-window sensitivity result that removes the heart-rate-agreement exclusion while retaining the filtered-ECG coverage requirement.
- Regenerated the quality-curated sliding summary plot with descriptive vessel-feature labels and added the broader sensitivity summary plot for the write-up.
- Moved interpretation notes into the Results text next to the figures they qualify and removed Wilcoxon reporting from the manuscript prose.
- Model: Codex, GPT-5.
- Effort: standard interactive editing.
- Token budget: not specified.
- Verification: ran MATLAB batch generation for the curated and all-recording overlap Spearman summaries and visually inspected both updated write-up figures.

### Figure Hierarchy Revision

- Moved the across-mouse sliding-window Spearman summary into the Results section as the sole main figure.
- Moved the heart-rate validation and older non-sliding comparison figures into a supplementary section.
- Removed the representative one-mouse scatter figure because it did not add a necessary group-level conclusion.
- Model: Codex, GPT-5.
- Verification: confirmed the remaining Markdown image paths resolve and reviewed the revised figure numbering and placement.

## 2026-06-22

### Sliding-Window Correlation Write-Up

- Narrowed the manuscript write-up to Spearman correlations using quality-filtered one-minute ECG heart rate and overlapping ten-minute vessel windows advanced in one-minute steps.
- Added the primary sliding-window result interpretation and a brief comparison with the older non-sliding ten-minute analysis.
- Added representative heart-rate validation, within-mouse correlation, group-summary, and legacy comparison figures for direct rendering in the write-up.
- Clarified that the primary displayed analysis uses the five quality-curated recordings and that the findings remain exploratory.
- Model: Codex, GPT-5.
- Verification: checked the active overlap extraction and correlation scripts, confirmed the group statistics against the current CSV summaries, visually inspected the included figures, and verified all Markdown image paths.
