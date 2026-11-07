#!/bin/bash
#SBATCH -J klem_maskref
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 05:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# data-driven site-mask refinement at the port tree's root: a site EM's two classes (tree junk excluded) -> 2-state reconstruct3D -> cross-half
# variance mask (arm_mask2.py) -> site EM (30 it) from the same seed with the new mask -> closed metrics; ROUNDS rounds. env O PT POST0 SEED ROUNDS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${O:-$K_/factorial/port_sites_cdisc}; PT=${PT:-$K_/stree_port/fulltoy_cdisc/1_hetero_tree}
D=$O/maskref; mkdir -p $D; H=$T/simsites; MEAN=$O/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc; SEED=${SEED:-$O/fits/r1_s4@r/hi/top_star.txt}; POST=${POST0:-$O/fits/r1_s4@r/hi/post_master.txt}
TJ=$K_/factorial/junk_test/simsites_truth_jt.txt; M=$T/maps160/simsites
$PY -c "
import numpy as np
s = np.loadtxt('$SEED', dtype=int) > 0
for h in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h).astype(int) - 1; np.savetxt('$D/start_%s.txt' % h, 1 + s[r].astype(int), fmt='%d')"
for k in $(seq 1 ${ROUNDS:-2}); do R=$D/round$k; mkdir -p $R; cd $R
  $PY -c "
import numpy as np
P = np.loadtxt('$POST'); TL = np.loadtxt('$PT/hetero_tree_labels.txt', dtype=int); L = np.where(TL > 0, 1 + (P >= 0.5).astype(int), 0); np.savetxt('$R/labels.txt', L, fmt='%d')
print('round $k: classes from %s: %d / %d (tree junk excluded: %d)' % ('$POST'.split('port_sites_cdisc/')[-1], (L == 1).sum(), (L == 2).sum(), (TL == 0).sum()))"
  if [ -z "$(ls $R/*_reconstruct3D/recvol_state02_odd_unfil.mrc 2>/dev/null)" ]; then cp -f $PT/n.simple r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$R/labels.txt mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=r.simple qsys_name=slurm nparts=4 nthr=8 mkdir=no > up.log 2>&1 && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=2 nparts=4 nthr=8 > rec.log 2>&1 \
    || { echo "round $k reconstruct FAILED: $(grep -a -m2 ERROR rec.log)"; exit 1; }; fi
  RD=$(ls -td $R/*_reconstruct3D | head -1); n1=$(awk '$1 == 1' $R/labels.txt | wc -l); n2=$(awk '$1 == 2' $R/labels.txt | wc -l)
  printf "s1 %s %s %d\ns2 %s %s %d\n" $RD/recvol_state01_even_unfil.mrc $RD/recvol_state01_odd_unfil.mrc $n1 $RD/recvol_state02_even_unfil.mrc $RD/recvol_state02_odd_unfil.mrc $n2 > $R/list.txt
  $PY $T/arm_mask2.py $R/mask.mrc --list $R/list.txt --mskdiam 232 --frac ${FRAC:-0.2} --dil 3 --falloff 6 2>&1 | grep -v Warning | tail -1
  $PY $T/mask_vs_truth.py $M/mask_A.mrc,$M/mask_B.mrc,$M/mask_C.mrc $R/mask.mrc 2>&1 | grep -v Warning
  $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/node_r/poses_A.txt,$O/node_r/poses_B.txt --labels $D/start_A.txt,$D/start_B.txt --cache $T/benchcache_sites --ctfabs 0 --cons_ref $MEAN \
    --mask $R/mask.mrc --soft 1 --iters ${ITERS:-30} --em_lim 15 --post 1 --nomaps 1 --out $R/em 2>&1 | grep -v Warning > $R/em.log; grep -a "^iteration" $R/em.log | tail -1 | cut -c1-110
  $PY - <<PYEOF
import numpy as np
T = np.loadtxt('$TJ'); hd = [l.split()[1:] for l in open('$TJ') if l.startswith('#')][0]; c = lambda n: T[:, hd.index(n)]; junk = c('junk') > 0; y = (c('C_closed') > 0) & ~junk; g = ~junk; NM = len(T)
def auc(s, yy): r_ = np.argsort(np.argsort(s)) + 1.; n1 = yy.sum(); n0 = len(yy) - n1; return (r_[yy].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
P = np.full(NM, np.nan)
for h in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h).astype(int) - 1; P[r] = np.loadtxt('$R/em/post_soft_%s.txt' % h)[:, 1]
np.savetxt('$R/post_master.txt', P, fmt='%.4f'); cl = P >= 0.5; o = np.argsort(-np.where(g, P, -np.inf))
print('round $k EM: class 2 = %d (junk %.3f); closed precision %.3f recall %.3f (base %.3f); closed AUC %.3f; top-2.5k purity %.3f, top-10k %.3f'
      % (cl.sum(), (cl & junk).sum() / max(cl.sum(), 1), (cl & y).sum() / max((cl & g).sum(), 1), (cl & y).sum() / y.sum(), y.sum() / g.sum(), auc(P[g], y[g]), y[o[:2500]].mean(), y[o[:10000]].mean()))
PYEOF
  POST=$R/post_master.txt; done
echo "[$(date '+%F %T') mask refine] DONE"
