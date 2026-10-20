#!/bin/bash
#SBATCH -J toy10_c18_tip
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
for n in true realigned consensus; do echo "== pure true-18 particles at the $n poses"; $PY $W/closed_arm_cells.py $W/lat/klem_c18_$n --top 1 2>&1 | grep -v Warning | tail -1; done
echo '== real closed subset boards (116 A mask run)'; $PY $W/closed_arm_cells.py /mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_zoom/board --top 3 2>&1 | grep -v Warning | tail -3
