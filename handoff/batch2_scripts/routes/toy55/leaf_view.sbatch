#!/bin/bash
#SBATCH -J klem_leafview
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# an rflex run's leaves as viewing maps (lp 10, unsharpened, box 160) in view_rflex/<RT>_leaves/ + their closed-arm readout. env RT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/rflex_${RT:?}
V=$K_/view_rflex/${RT}_leaves; I=$K_/view_rflex/info; mkdir -p $I; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
awk '{printf "leaf%02d_%s %s %s\n", NR, $1, $3, $4}' $R/leaves.txt > $I/spec_${RT}_leaves.txt
$PY $T/make_view.py $V $I/spec_${RT}_leaves.txt --box 160
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $I/arm_${RT}_leaves.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1
awk 'NR == FNR {c[FNR] = $5; next} !/^#/ {n = split($1, p, "/"); f = p[n]; k = substr(f, 5, 2) + 0; printf "   %s (%d ptcls): closed-arm %.2f open-arm %.2f %s\n", f, c[k], $2, $3, $4}' $R/leaves.txt $I/arm_${RT}_leaves.txt | sort -k5,5gr
echo DONE
