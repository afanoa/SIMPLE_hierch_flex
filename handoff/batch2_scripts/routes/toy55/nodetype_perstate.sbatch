#!/bin/bash
#SBATCH -J nodetype_perstate
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Per-state gate (nodetype_perstate.py: every hard state's cross-half deviation from the leave-one-out mean in the free band, max over
# states with Bonferroni) on the dev nodes (s30a, s30split2 r1-r5, the null_50k / two_50k panels; with truth) and the real f692rf nodes.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
G30=$T/simswing30/gt_particles.txt
ps(){ $PY -u $T/nodetype_perstate.py $1 $2 --bands 9:7,9:6 --mskdiam $3 --tag "[$4]" ${5:+--gt $5} 2>&1 | grep --line-buffered -v Warning; }
for u in r11 r12 r21 r22 r r1 r2; do D=$K_/rflex_s30a/node_$u; ps $(ls -td $D/*_refine3D | head -1) $D/refined_states.txt 220 s30a_$u $G30; done
for u in r1 r2 r3 r4 r5; do [ $u = r1 ] && D=$K_/rflex_s30split2/node_r1 || D=$K_/rflex_s30split2_$u/node_$u
  [ -s $D/refined_states.txt ] && ps $(ls -td $D/*_refine3D | head -1) $D/refined_states.txt 220 s30split2_$u $G30; done
for m in null_50k two_50k; do ps $K_/panel_$m/1_refine3D_states $K_/panel_$m/refined_states.txt 220 panel_$m $T/sim$m/gt_particles.txt; done
for u in r r1 r2 r21 r22; do D=$K_/rflex_f692rf/node_$u; ps $(ls -td $D/*_refine3D | head -1) $D/refined_states.txt 232 f692rf_$u; done
echo "[$(date '+%F %T') nodetype_perstate] DONE"
