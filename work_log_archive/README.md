# Work Log Archive

Older session notes rotated out of [`../work_log.md`](../work_log.md). Archive files use the same newest-first structure as the live log.

## Rotation Policy

The live `work_log.md` holds at most the five most recent unique calendar dates. When a new date would push it past five, move the oldest five dates together into a new archive file.

## File Naming

Use `work_log_<earliest>_to_<latest>.md`, where the dates describe the entries contained in the file. Each archive file holds exactly five unique calendar dates.

## Search All History

```powershell
rg -n '^## [0-9]{4}-[0-9]{2}-[0-9]{2}' work_log.md work_log_archive/
```
