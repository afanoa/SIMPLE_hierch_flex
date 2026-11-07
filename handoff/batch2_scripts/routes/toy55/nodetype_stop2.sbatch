#!/bin/bash
#SBATCH -J nodetype_stop2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Stop test on the HARD states of every node (nodetype_stop2.py, band 10-8 A, which the lp-10 refinement never used): s30a dev nodes and
# the real f692rf nodes.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for u in r11 r21 r22 r12 r1 r r2; do R=$K_/rflex_s30a; $PY -u $T/nodetype_stop2.py $(ls -td $R/node_$u/*_refine3D | head -1) $R/node_$u/refined_states.txt --bands 10:8,10:6 --tag "[s30a_$u hard]" 2>&1 | grep --line-buffered -v Warning; done
for u in r r1 r2 r21 r22; do R=$K_/rflex_f692rf; $PY -u $T/nodetype_stop2.py $(ls -td $R/node_$u/*_refine3D | head -1) $R/node_$u/refined_states.txt --bands 10:8,10:6 --mskdiam 232 --tag "[f692rf_$u hard]" 2>&1 | grep --line-buffered -v Warning; done
echo "[$(date '+%F %T') nodetype_stop2] DONE"
