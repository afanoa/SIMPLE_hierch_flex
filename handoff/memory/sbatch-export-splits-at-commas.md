---
name: sbatch-export-splits-at-commas
description: "TRAP: sbatch --export=ALL,VAR=a,b cuts VAR at the comma (b becomes a bogus entry) - pass comma values through the environment (VAR='a,b' sbatch --export=ALL ...) or inside --wrap"
metadata:
  type: feedback
---

Hit twice on 2026-10-20: GTX="--gt A,B" for the synthetic merge-up and CODES=swing:200:1.0,breath:300:1.0 for the fast factor
readout - both arrived truncated at the first comma, so the second truth code silently vanished (swing values up to 103 A
because breathing particles counted as swing).
**Why:** sbatch parses --export as a comma-separated list of NAME=value entries.
**How to apply:** never put a value containing commas in --export; set it in the submitting shell (`CODES='a,b' sbatch
--export=ALL,OTHER=x job.sbatch`) or pass it on the command line inside --wrap; check the job's first log lines for the value.
