#!/bin/bash
#SBATCH -J klem_lsamp_pqp
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the pipeline's own poles (stree_port poles_r1_p1 / p2, node r: three pole-seeded leaves of stageB): per run the score space against its
# three leaf maps at its stageB poses (lsamp_score2.sbatch), the poles named by overlap with the factorial poles, the three-pole coordinates
# (lsamp_pq_coords.py), the stage on the closed - open A axis (lsamp2p.py --coords), equal-N bins of every coordinate type; then the bins
# comparison of all five runs per type (with lsamp_pq_f.sbatch's factorial-pole runs), the stage comparison and the closed-child test.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=$K_/lsamp/repro_r1; O=$W/pq; mkdir -p $O
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
F=$K_/factorial/f692st_r1_polesets; REF=$F/script_z1lo/post_master.txt,$F/script_z1hi/post_master.txt,$F/script_z3hi/post_master.txt; TYPES="c_oA c_oB oA_oB c_rest"
log(){ echo "[$(date '+%F %T') pq] $*"; }
for p in p1 p2; do N=$K_/stree_port/poles_r1_$p/1_hetero_tree/node_r; D=$O/${p}_tree; mkdir -p $D
  if [ ! -s $D/leaves_list.txt ]; then
    ( for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done; B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
      cd $D && simple_exec prg=print_project_field projfile=$N/stageB/b.simple oritype=ptcl3D mkdir=no > dump.txt 2> /dev/null )
    $PY - <<PYEOF || { log "$p: tree setup FAILED"; continue; }
import numpy as np
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
L = open('$D/dump.txt').read().splitlines(); ch = np.loadtxt('$N/split/split_children.txt', dtype=int); assert len(L) == len(ch), 'dump %d lines vs %d labels' % (len(L), len(ch))
st = np.array([int(float(kv(l).get('state', 0))) for l in L]); lab = np.where(ch > 0, st, 0); np.savetxt('$D/final_labels.txt', lab, fmt='%d')
cnt = np.bincount(lab, minlength=4)[1:]; print('$p: %d node particles, stageB states %s (split children %s)' % ((lab > 0).sum(), ' '.join(map(str, cnt)), ' '.join(map(str, np.bincount(ch, minlength=4)[1:]))))
open('$D/leaves_list.txt', 'w').write(''.join('r_p%d $N/stageB/recvol_state%02d_even_unfil.mrc $N/stageB/recvol_state%02d_odd_unfil.mrc %d\n' % (k, k, k, cnt[k - 1]) for k in (1, 2, 3)))
PYEOF
  fi
  [ -s $O/${p}_ps/meta.npz ] || { TREE=$D OUT=$O/${p}_ps SIM=${H}_A,${H}_B MR=$MR CP=$CP CTFABS=1 EMLIM=12 MSKDIAM=232 POSES=$D/dump.txt bash $T/lsamp_score2.sbatch 2>&1 | grep "own leaf ranks\|DONE\|Error\|Traceback"; }
  $PY -u $T/lsamp_pq_coords.py --scores $O/${p}_ps --name $p --out $O --roles auto --ref_post $REF 2>&1 | grep -v "Warning\|warn("
  mkdir -p $O/stage_c_oA/$p; $PY -u $T/lsamp2p.py --scores $O/${p}_ps --sim ${H}_A,${H}_B --mr $MR --cache $CP --poses $D/dump.txt --ctfabs 1 --mskdiam 232 --out $O/stage_c_oA/$p --use r_ --axes 1 \
      --coords $O/c_oA/$p/coords.npz --ends centre --tag ${p}_c_oA 2>&1 | grep --line-buffered -v "Warning\|warn(" > $O/stage_c_oA/$p/lsamp.log
  grep "resolvable samples\|CHOSEN\|steps CV\|Error\|Traceback" $O/stage_c_oA/$p/lsamp.log | cut -c1-200
  for d in heldout_A heldout_B final_all; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/stage_c_oA/$p/arm_$d.txt $O/stage_c_oA/$p/$d 2>&1 | grep -v Warning > /dev/null
    echo "  $d: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $O/stage_c_oA/$p/arm_$d.txt)"; done
  for ty in $TYPES; do $PY -u $T/lsamp_repro_bins2.py --name $p --rundir $O/$ty/$p --scores $O/${p}_ps --poses $D/dump.txt --ref_rundir $O/$ty/script --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O/$ty/bins 2>&1 | grep -v "Warning\|warn("
    for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$ty/bins/$p/arm_bins_$f.txt $O/$ty/bins/$p/bins_$f 2>&1 | grep -v Warning > /dev/null; done; done
done
for ty in $TYPES; do log "bins comparison $ty"; $PY -u $T/lsamp_repro_bins2.py --cmp 1 --names script,port1,port2,p1,p2 --out $O/$ty/bins 2>&1 | grep -v "Warning\|warn(\|def lp" > $O/$ty/bins/bins_summary.txt; done
$PY -u $T/lsamp_repro_cmp2.py $O/stage_c_oA/repro_summary.txt script:$W/poleaxis_stage/script port1:$W/poleaxis_stage/port1 port2:$W/poleaxis_stage/port2 p1:$O/stage_c_oA/p1 p2:$O/stage_c_oA/p2 2>&1 | grep -v "Warning\|warn(\|def prep" > /dev/null
RT=""; for r in script port1 port2 port3; do RT="$RT,ax2_$r=$W/${r}_ax2/coords.npz"; done; for r in script port1 port2 p1 p2; do RT="$RT,oA_oB_$r=$O/oA_oB/$r/coords.npz,c_oA_$r=$O/c_oA/$r/coords.npz"; done
$PY -u $T/lsamp_pq_child.py --child p1=$K_/stree_port/childdisc/p1_r1/site,p2=$K_/stree_port/childdisc/p2_r1/site --root ${RT#,} 2>&1 | grep -v "Warning\|warn(" | tee $O/child_test.txt
log DONE
