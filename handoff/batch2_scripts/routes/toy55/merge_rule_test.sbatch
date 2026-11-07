#!/bin/bash
#SBATCH -J klem_merge_rule
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merge rule test on a synthetic run with truth (env RTAG=synth_55_200k, MODE=55_200k): the data-derived merge of the
# leaves with the fixed 12 A floor vs the auto floor (difference resolved at the pooled state's own resolution x MARGIN);
# each consolidated (all particles, reconstructed) and read against the truth: states, sizes, purity, coverage, strict recovery
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; RTAG=${RTAG:?}; MODE=${MODE:?}; GT=$W/sim$MODE/gt_particles.txt
for FL in 12 auto; do
  echo "######## floor $FL (margin ${MARGIN:-1.0})"; export RTAG LEAFGLOB="region*" FLOOR=$FL MARGIN=${MARGIN:-1.0} MASK=region WHITEN=0 UNCAP=1 TAG=rule_$FL
  bash $W/pooled_dd.sbatch 2>&1 | grep -a "leaf states\|pairs at the start\|merges ->\|^cut\|dropped\|identified\|recovered\|FSC median\|merge FAILED\|Traceback" | cut -c1-200
  C=$A/rnasep_klem/$RTAG/consolidated_${FL}_rule_$FL; [ "$FL" = auto ] && C=$A/rnasep_klem/$RTAG/consolidated_12_rule_auto
  echo "== truth purity of the merged states (floor $FL)"; $PY $W/consolidated_purity.py $C $GT 2>&1 | grep -v Warning
done
echo "MERGE RULE DONE"
