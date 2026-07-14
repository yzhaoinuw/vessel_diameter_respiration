# Guidelines and Tips for Agents

Read this file first when joining the repository. It identifies the active analysis path, the archive boundary, and the documentation that must move with substantive changes.

## Startup Rule

Do not read every Markdown file automatically. Start here, then use the documentation map below. Treat [`correlation_analysis_writeup.md`](correlation_analysis_writeup.md) as the authoritative statement of the analyses retained for reporting.

## Runtime Environment

Use the `ecg` conda environment:

```powershell
conda activate ecg
```

Run MATLAB scripts from the repository root so their `pwd`-relative data and output paths resolve correctly:

```powershell
matlab -batch "script_name"
```

The repository has no automated unit-test suite. Verification is script-specific and depends on the local, ignored `data/` tree.

## Active Analysis Workflow

The retained workflow supports the primary overlapping-window Spearman analysis, its all-recording sensitivity analysis, heart-rate quality checks, and the older Yue-HR Spearman supplement described in the write-up.

1. `PlotResults_continuous_v5.m` produces legacy `results.mat` variables used for Yue HR and vessel-derived cardiac frequency comparisons.
2. `extract_hr_amp_data_all_subjects.m` produces `data/<subject>/hr_amp_data.mat`, including quality-filtered one-minute ECG HR.
3. `extract_vessel_features_10min_overlap.m` produces overlapping ten-minute vessel features advanced in one-minute steps.
4. `hr_hr_correlation/` contains the retained filtered ECG/Yue/vessel HR quality checks.
5. `spearman_correlation_overlap/` computes, tests, diagnoses, and plots the primary curated and all-recording overlap analyses. Set `overlap_analysis_mode` to `"curated"` or `"all_subjects"` before running the scripts.
6. `spearman_correlation_yue/` reproduces the older non-sliding Yue-HR Spearman supplement.
7. `plot_three_hr_comparison.m` produces the representative three-HR quality-control overlay.

Do not revive archived Pearson, five-minute sliding, cross-correlation, respiration, or F153-only pilot workflows without an explicit request to reopen those analyses.

## Common Tasks

Primary overlap extraction and analysis:

```powershell
matlab -batch "extract_hr_amp_data_all_subjects"
matlab -batch "extract_vessel_features_10min_overlap"
matlab -batch "run('spearman_correlation_overlap/compute_hr_diam_corr_spearman_overlap.m')"
matlab -batch "run('spearman_correlation_overlap/compute_hr_diam_corr_spearman_overlap_ttest.m')"
matlab -batch "run('spearman_correlation_overlap/plot_hr_diam_corr_spearman_overlap_ttest_summary.m')"
```

Repository pre-flight:

```powershell
git diff --check
C:\Users\yzhao\python_projects\agent_collab_treaty\.venv\Scripts\treaty.exe validate .
```

For analysis changes, run the affected MATLAB scripts and inspect their CSV/PNG outputs. Record only commands actually run in `work_log.md`.

## When To Update Treaty Docs

At the end of substantive work, prepend an entry to `work_log.md`. Update `next_steps.md` when future work changes, and update `project_overview.md` when the active/archived boundary changes.

The live work log holds at most five unique calendar dates. Before adding a dated entry, verify the local date with `Get-Date -Format yyyy-MM-dd`. When a sixth date is introduced, rotate the oldest five dates together into `work_log_archive/work_log_<earliest>_to_<latest>.md`.

## Branch Handoff Discipline

`main` is the integration branch. Before switching branches, confirm that current work is committed, intentionally parked, or otherwise handed off:

```powershell
git status --short --branch
git log --oneline --left-right --cherry-pick main...HEAD
git merge-base --is-ancestor main HEAD
```

Do not switch away while substantive local analysis or archive changes are unresolved.

## Updating The Treaty

Upstream treaty changes are a maintainer action. `treaty update` requires a clean working tree and performs a three-way merge. When asked to update it:

1. Confirm the central `agent_collab_treaty` checkout and remote refs are current.
2. Commit or stash destination changes.
3. Run `treaty update`.
4. Resolve merge markers, review the diff, and run `treaty validate .`.

Do not edit `.copier-answers.yml` manually; Copier owns its source revision metadata.

## Documentation

- `README.md` — concise repository purpose, retained workflow, and setup.
- `correlation_analysis_writeup.md` — authoritative methods, results, interpretation, and figure set.
- `project_overview.md` — active pipeline, structure, data contracts, and archive boundary.
- `next_steps.md` — only current or genuinely paused work; remove completed analysis plans.
- `work_log.md` and `work_log_archive/` — recent and historical implementation/verification record.
- `archive/README.md` — archive categories and rules.

## Git And Generated Data

The repository intentionally ignores raw `data/`, generated result tables, diagnostic figures, MATLAB autosaves, and other local artifacts. The four curated figures under `writeup_figures/` are tracked because the write-up embeds them.

If Git reports dubious ownership, mark this repository safe:

```powershell
git config --global --add safe.directory C:/Users/yzhao/matlab_projects/vessel_diameter_respiration
```

## Commit Message Guidelines

Use a short title line. For a multi-part commit, add short flat bullets describing high-level behavior or organization changes. Do not mention tests or internal mechanics unless they are the purpose of the commit.

## Project-Specific Reminders

- The word "sliding" in the write-up refers to overlapping ten-minute vessel windows advanced by one minute, not the archived five-minute `hr_amp_data_sliding.mat` workflow.
- The primary overlap analysis pairs each vessel-window start with the quality-filtered one-minute ECG HR value at that same start time.
- `results.hr` is vessel-derived cardiac frequency in hertz; it is not ECG HR. Convert it to BPM before HR-quality comparisons.
- Keep the mouse as the unit of group-level inference. Overlapping time windows are repeated observations, not independent animals.
- Treat the reported association with slow vasomotion amplitude as exploratory because multiple feature-period tests were not adjusted.
- Archived scripts are historical references. Do not edit them as the implementation path for current analyses.
