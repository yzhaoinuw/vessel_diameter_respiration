# Work Log

## 2026-07-14

### Treaty Badge and Repository Delivery (Codex GPT-5, standard delivery effort, token budget not specified)

- Verified that `README.md` displays the official centrally hosted tricolor Agent Collab Treaty adoption badge.
- Prepared the validated treaty refresh and archive cleanup for delivery on `main`.

- Verification:
  - Confirmed the badge image and adoption link use the official `agent_collab_treaty` repository URLs.
  - Re-ran treaty validation and Git staged-diff checks before committing.

## 2026-07-13

### Repository Cleanup and Treaty Alignment (Codex GPT-5, extended effort, token budget not specified)

- Reduced the supported MATLAB surface to the 15 scripts that reproduce the analyses retained by `correlation_analysis_writeup.md`.
- Preserved 47 obsolete or superseded MATLAB scripts under categorized archive paths, including the existing locally modified archive copies and the missing tracked pipeline sketch reconstructed from Git history.
- Archived four superseded planning/methodology documents and added an archive index that distinguishes historical code from supported entry points.
- Replaced placeholder treaty documentation with repository-specific runtime, workflow, data-contract, verification, and branch-handoff guidance.
- Canonicalized `work_log.md`, added the treaty badge and Copier metadata to the tracked documentation surface, and revised `.gitignore` so active governance/source files are tracked while raw data and generated outputs remain ignored.

- Verification:
  - Live-checked the central `agent_collab_treaty` remote: `main` and `dev` both resolved to `c8a4138909fdf5c53116cc968af06e25483c67a3`.
  - Ran MATLAB Code Analyzer over all 15 retained scripts; it completed without parser errors and reported only existing performance/style findings (`AGROW`, `NASGU`, `SAGROW`, `SVFIGC`, `MSNU`, and `LAXES`).
  - Confirmed that all four write-up figures exist and are byte-identical to their current source outputs.
  - Confirmed that every formerly tracked deleted MATLAB filename has a corresponding archived copy and that the archive contains 47 `.m` files.
  - Ran `git diff --check`, the local Markdown-link check, active-to-archive reference checks, and `treaty validate .`.

## 2026-07-06

### Multiple-Comparison Caveat Clarification (Codex GPT-5, standard interactive editing, token budget not specified)

- Expanded the Results caveat to explain that unadjusted probability values do not use a stricter threshold for multiple related feature-period tests.
- Clarified why the caveat matters for interpreting nominally significant findings as exploratory.

- Verification:
  - Reviewed the Results paragraph diff and checked branch status before staging.

## 2026-06-23

### All-Recording Sensitivity Write-Up (Codex GPT-5, standard interactive editing, token budget not specified)

- Added the sliding-window sensitivity result that removes the heart-rate-agreement exclusion while retaining the filtered-ECG coverage requirement.
- Regenerated the quality-curated sliding summary plot with descriptive vessel-feature labels and added the broader sensitivity summary plot for the write-up.
- Moved interpretation notes into the Results text next to the figures they qualify and removed Wilcoxon reporting from the manuscript prose.

- Verification:
  - Ran MATLAB batch generation for the curated and all-recording overlap Spearman summaries and visually inspected both updated write-up figures.

### Figure Hierarchy Revision (Codex GPT-5, effort not specified, token budget not specified)

- Moved the across-mouse sliding-window Spearman summary into the Results section as the sole main figure.
- Moved the heart-rate validation and older non-sliding comparison figures into a supplementary section.
- Removed the representative one-mouse scatter figure because it did not add a necessary group-level conclusion.

- Verification:
  - Confirmed the remaining Markdown image paths resolve and reviewed the revised figure numbering and placement.

## 2026-06-22

### Sliding-Window Correlation Write-Up (Codex GPT-5, effort not specified, token budget not specified)

- Narrowed the manuscript write-up to Spearman correlations using quality-filtered one-minute ECG heart rate and overlapping ten-minute vessel windows advanced in one-minute steps.
- Added the primary sliding-window result interpretation and a brief comparison with the older non-sliding ten-minute analysis.
- Added representative heart-rate validation, within-mouse correlation, group-summary, and legacy comparison figures for direct rendering in the write-up.
- Clarified that the primary displayed analysis uses the five quality-curated recordings and that the findings remain exploratory.

- Verification:
  - Checked the active overlap extraction and correlation scripts, confirmed the group statistics against the current CSV summaries, visually inspected the included figures, and verified all Markdown image paths.
