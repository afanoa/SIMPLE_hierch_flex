#!/bin/bash
#SBATCH -J score_root_m79
#SBATCH -p norm
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is the m33 root motion useful? Reads the root replay stree_port/gate/m33_root_m79 node_r, writes nothing there:
# - fa_maps10.py --dock against the 55 deposited conformers: site 6's 5 bins, every other candidate's bins 1 and 5, the 3 stage-B
#   states (arm position s1, closed-arm occupancy, the paired arm / core test);
# - phaseE_maptests.py + phaseF_rel.py per candidate (effect size against half-map noise, in the node mask);
# - quality_bins.py per candidate and for stage B (corr / defocus / views per bin), and stage B against the A / middle / B seeds
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; FA=$K/full692_analysis
N=$K/stree_port/gate/m33_root_m79/1_hetero_tree/node_r; OUT=$K/scorecard/m33_root_m79; mkdir -p $OUT/series; L=$OUT/maplist.txt; : > $L; S=()
log(){ echo "[$(date "+%F %T") score_root] $*"; }
for D in $N/motion/c*/; do c=$(basename $D); W=$OUT/series/$c; mkdir -p $W/halves; [ -s $D/recvol_state05.mrc ] || continue
  for k in 1 2 3 4 5; do s=$(printf %02d $k); ln -sf $D/recvol_state$s.mrc $W/bin$k.mrc; ln -sf $D/recvol_state${s}_even.mrc $W/halves/bin${k}_even.mrc
    ln -sf $D/recvol_state${s}_odd.mrc $W/halves/bin${k}_odd.mrc; done
  awk '{n[$1]++} END {for (k = 1; k <= 5; k++) printf "bin %d: rows %d\n", k, n[k]}' $D/bins.txt > $W/bins.txt
  for k in 1 2 3 4 5; do [ $c = c6 ] || [ $k = 1 ] || [ $k = 5 ] || continue; s=$(printf %02d $k)
    echo "${c}_bin$k other $c $(grep -c "^$k\$" $D/bins.txt) $D/recvol_state$s.mrc $D/recvol_state${s}_even.mrc $D/recvol_state${s}_odd.mrc" >> $L; done
  S+=("$c")
done
for s in 1 2 3; do echo "stageB_s$s node stageB $(wc -l < $N/leaf_r_s${s}_rows.txt) $N/stageB/recvol_state0$s.mrc $N/stageB/recvol_state0${s}_even.mrc $N/stageB/recvol_state0${s}_odd.mrc" >> $L; done
log "${#S[@]} candidates, $(wc -l < $L) maps"
( OMP_NUM_THREADS=4 $PY -u $FA/fa_maps10.py $L $OUT/fa --dock --title "m33 root motion (m79 replay)" 2>&1 | grep -v "Warning\|warn(" > $OUT/fa.log ) &
for c in "${S[@]}"; do W=$OUT/series/$c
  ( $PY $T/phaseE_maptests.py root_$c:$W:$N/mask.mrc:225.456 > $W/maptests.txt 2>&1; $PY $T/phaseF_rel.py root_$c:$W:$N/mask.mrc > $W/relF.txt 2>&1
    cd $W && simple_exec prg=print_project_field projfile=$N/motion/$c/m.simple oritype=ptcl3D 2>/dev/null | $PY $T/quality_bins.py root_$c 1 - 3 > $W/quality.txt 2>&1 ) &
done
( cd $OUT && simple_exec prg=print_project_field projfile=$N/stageB/b.simple oritype=ptcl3D 2>/dev/null > $OUT/stageB_ptcl3D.txt
  $PY $T/quality_bins.py root_stageB 1 - max < $OUT/stageB_ptcl3D.txt > $OUT/stageB_quality.txt 2>&1
  grep -o "state=[0-9.]*" $OUT/stageB_ptcl3D.txt | cut -d= -f2 | cut -d. -f1 > $OUT/stageB_states.txt; rm -f $OUT/stageB_ptcl3D.txt
  $PY -c "import numpy as np, sys; st = np.loadtxt(sys.argv[1], dtype=int); r = np.loadtxt(sys.argv[2], dtype=int) - 1; ab = np.loadtxt(sys.argv[3], dtype=int)
print('stage-B state x seed (ab.txt: 0 middle, 1 A, 2 B), rows of the root:'); u, c = np.unique(np.c_[st[r], ab], axis=0, return_counts=True)
[print('  state %d seed %d: %d' % (a, b, n)) for (a, b), n in zip(u, c)]" $OUT/stageB_states.txt $N/rows.txt $N/motion/c6/ab.txt > $OUT/stageB_vs_ab.txt 2>&1 ) &
wait
grep -h "^SUMMARY" $OUT/series/*/maptests.txt > $OUT/effect_summary.txt; grep -h "^RELF" $OUT/series/*/relF.txt > $OUT/effect_relF.txt
grep -h "^QUAL\|^QEND" $OUT/series/*/quality.txt $OUT/stageB_quality.txt > $OUT/quality_summary.txt
log "DONE: fa $(grep -c . $OUT/fa/maps_table.txt 2>/dev/null) table lines, $(wc -l < $OUT/effect_summary.txt) effect, $(grep -c ^QUAL $OUT/quality_summary.txt) quality lines"
