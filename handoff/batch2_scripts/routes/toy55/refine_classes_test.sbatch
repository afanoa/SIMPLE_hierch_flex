#!/bin/bash
#SBATCH -J klem_refclasses
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# mask-refinement classes, option (a) = the site's pole label sets (top-N*) + rest vs option (b) = hard P >= 0.5 classes (maskref/round1):
# a3 = {lo pole, hi pole, rest}, a2 = {hi pole, rest}; tree junk excluded; 3- / 2-state reconstruct3D -> arm_mask2.py -> site EM (30 it,
# R1's seed) -> closed AUC / top-N / class-2 share and mask coverage of true C. env O PT SITE
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${O:-$K_/factorial/port_sites_cdisc}; PT=${PT:-$K_/stree_port/fulltoy_cdisc/1_hetero_tree}
D=$O/refclasses; mkdir -p $D; H=$T/simsites; MEAN=$O/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc; SITE=${SITE:-r1_s4}; TJ=$K_/factorial/junk_test/simsites_truth_jt.txt; M=$T/maps160/simsites
one(){ v=$1; R=$D/$v; mkdir -p $R; cd $R
  $PY -c "
import numpy as np
L = np.loadtxt('$O/labels/$SITE.txt', dtype=int); TL = np.loadtxt('$PT/hetero_tree_labels.txt', dtype=int); rows = np.loadtxt('$O/view/rows_r.txt') > 0
C = np.where(L == 1, 1, np.where(L == 2, 2, 3)) if '$v' == 'a3' else np.where(L == 2, 1, 2); C[~rows | (TL == 0)] = 0
np.savetxt('$R/classes.txt', C, fmt='%d'); print('$v classes (tree junk excluded):', np.bincount(C).tolist()[1:])"
  K=$(awk '$1 > 0' $R/classes.txt | sort -n | uniq | wc -l)
  if [ -z "$(ls $R/*_reconstruct3D/recvol_state0${K}_odd_unfil.mrc 2>/dev/null)" ]; then cp -f $PT/n.simple r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$R/classes.txt mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=r.simple qsys_name=slurm nparts=4 nthr=8 mkdir=no > up.log 2>&1 && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=$K nparts=4 nthr=8 > rec.log 2>&1 \
    || { echo "$v reconstruct FAILED: $(grep -a -m2 ERROR rec.log)"; return 1; }; fi
  RD=$(ls -td $R/*_reconstruct3D | head -1); : > $R/list.txt
  for s in $(seq 1 $K); do echo "s$s $RD/recvol_state0${s}_even_unfil.mrc $RD/recvol_state0${s}_odd_unfil.mrc $(awk -v s=$s '$1 == s' $R/classes.txt | wc -l)" >> $R/list.txt; done
  $PY $T/arm_mask2.py $R/mask.mrc --list $R/list.txt --mskdiam 232 --frac 0.2 --dil 3 --falloff 6 2>&1 | grep -v Warning | tail -1
  $PY $T/mask_vs_truth.py $M/mask_A.mrc,$M/mask_B.mrc,$M/mask_C.mrc $R/mask.mrc 2>&1 | grep -v Warning | sed "s/^/$v /"
  $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/node_r/poses_A.txt,$O/node_r/poses_B.txt --labels $O/maskref/start_A.txt,$O/maskref/start_B.txt --cache $T/benchcache_sites --ctfabs 0 \
    --cons_ref $MEAN --mask $R/mask.mrc --soft 1 --iters 30 --em_lim 15 --post 1 --nomaps 1 --out $R/em 2>&1 | grep -v Warning > $R/em.log; cd $D; }
one a3 & one a2 & wait
$PY - <<PYEOF
import numpy as np, os
T = np.loadtxt('$TJ'); hd = [l.split()[1:] for l in open('$TJ') if l.startswith('#')][0]; c = lambda n: T[:, hd.index(n)]; junk = c('junk') > 0; y = (c('C_closed') > 0) & ~junk; g = ~junk; NM = len(T)
def auc(s, yy): r_ = np.argsort(np.argsort(s)) + 1.; n1 = yy.sum(); n0 = len(yy) - n1; return (r_[yy].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
for v, d in (('b (R1, hard P >= 0.5 classes)', '$O/maskref/round1/em'), ('a3 {lo pole, hi pole, rest}', '$D/a3/em'), ('a2 {hi pole, rest}', '$D/a2/em')):
    if not os.path.exists(d + '/post_soft_B.txt'): print(v, 'EM missing'); continue
    P = np.full(NM, np.nan)
    for h in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h).astype(int) - 1; P[r] = np.loadtxt('%s/post_soft_%s.txt' % (d, h))[:, 1]
    o = np.argsort(-np.where(g, P, -np.inf)); print('%-32s closed AUC %.3f, top-2.5k %.3f, top-10k %.3f, class-2 share %.3f' % (v, auc(P[g], y[g]), y[o[:2500]].mean(), y[o[:10000]].mean(), (P[np.isfinite(P)] >= 0.5).mean()))
PYEOF
