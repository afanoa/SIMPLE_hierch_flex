#!/bin/bash
#SBATCH -J nodetype_stop3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Stop test v3 (nodetype_stop3.py: 2-shell margins, per-group 3-D shift fit) on every dev node (cheap children 13-9.5 A, hard states
# 9-7 A), the split session's s30split2 nodes and the real f692rf nodes.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
kid(){ $PY -u $T/nodetype_stop3.py $(ls -td $1/rec/*_reconstruct3D | head -1) $1/children.txt --bands 13:9.5 --mskdiam $3 --tag "[$2 children]" 2>&1 | grep --line-buffered -v Warning; }
hard(){ $PY -u $T/nodetype_stop3.py $(ls -td $1/*_refine3D | head -1) $1/refined_states.txt --bands 9:7 --mskdiam $3 --tag "[$2 hard]" 2>&1 | grep --line-buffered -v Warning; }
for m in null_50k_rc:root_null_50k_rc null_50k_r:root_null_50k_r two_50k_r:root_two_50k_r two_50k_rc:root_two_50k_rc; do kid $K_/cheapkids_${m#*:} ${m%%:*} 220; done
for u in r11 r21 r22 r12 r1 r r2; do kid $K_/cheapkids_s30a_$u s30a_$u 220; hard $K_/rflex_s30a/node_$u s30a_$u 220; done
for u in r1 r2 r3 r4 r5; do [ $u = r1 ] && D=$K_/rflex_s30split2/node_r1 C=$K_/cheapkids_s30split2_r1 || { D=$K_/rflex_s30split2_$u/node_$u; C=$K_/cheapkids_s30split2_${u}_$u; }
  [ -s $C/children.txt ] && kid $C s30split2_$u 220; [ -s $D/refined_states.txt ] && hard $D s30split2_$u 220; done
for u in r r1 r2 r21 r22; do kid $K_/cheapkids_f692rf_$u f692rf_$u 232; hard $K_/rflex_f692rf/node_$u f692rf_$u 232; done
echo "[$(date '+%F %T') nodetype_stop3] DONE"
