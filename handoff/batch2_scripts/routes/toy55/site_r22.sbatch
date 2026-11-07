#!/bin/bash
#SBATCH -J klem_site_r22
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=110G
#SBATCH -t 02:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-SITE occupancy classification of real r22 (compositional route): 2 classes at fixed poses that differ only inside the closed-site
# mask (unit 7's cheap-children variation mask: closed-arm site 0.68, open site 0.01); the consensus outside the mask is shared (= signal
# subtraction), cross-half references, learned mixing weight, soft EM (minor_em3.py --post). Class 2 starts from the closed end of the
# 12 equal-width arm bins. Then the top-N particles by P(class 2) (N = 2500 5000 10000 20000) and the start bin are reconstructed with SIMPLE
# at r22's parent poses: closed-arm readout + the old closed states 15+16 traced into each subset. out rnasep_klem/site_r22/
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
MR=$K_/minor_em_r22; PP=$K_/rflex_f692rf/node_r2/2_refine3D/r.simple; MASK=${MASK:-$K_/rflex_f692u7/node_r/rsmask.mrc}; U=$K_/f692r/coordtree_fam1/units.txt
O=$K_/site_r22${TAG:-}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') site r22] $*"; }; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for h in A B; do [ -s $O/start_$h.txt ] || awk '{print ($1 == 0) ? 0 : ($1 == 12 ? 2 : 1)}' $MR/seed_$h.txt > $O/start_$h.txt; done
log "start: class 1 $(cat $O/start_A.txt $O/start_B.txt | awk '$1 == 1' | wc -l), class 2 (closed end bin) $(cat $O/start_A.txt $O/start_B.txt | awk '$1 == 2' | wc -l); mask $MASK"
[ -s $O/em/post_soft_B.txt ] || $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $MR/poses_A.txt,$MR/poses_B.txt --labels $O/start_A.txt,$O/start_B.txt --cache $T/benchcache_rnasep692 \
  --ctfabs 1 --cons_ref $T/maps160/toy_consensus.mrc --mask $MASK --soft 1 --iters ${ITERS:-10} --em_lim ${EMLIM:-15} --post 1 --nomaps 1 --out $O/em > $O/em.log 2>&1 || { log "EM FAILED: $(tail -2 $O/em.log)"; exit 1; }
grep -a '^start\|^iteration' $O/em.log | sed 's/deviation FSC.*//' | cut -c1-200
# master-row posterior of class 2, top-N subsets
$PY - <<PYEOF
import numpy as np
n = 932496; p = np.full(n, -1.); old = np.loadtxt('$U', dtype=int); cl = (old == 15) | (old == 16)
for h in 'AB':
    rows = np.loadtxt('$H' + '_%s/meta_rows.txt' % h).astype(int) - 1; P = np.loadtxt('$O/em/post_soft_%s.txt' % h)[:, 1]; s = np.loadtxt('$O/start_%s.txt' % h, dtype=int); m = s > 0; p[rows[m]] = P[m]
node = p >= 0; q = p[node]; print('posterior of class 2 over %d r22 particles: mean %.3f, > 0.5 %d, > 0.9 %d, > 0.99 %d; histogram (10 bins 0..1) %s' % (node.sum(), q.mean(), (q > 0.5).sum(), (q > 0.9).sum(), (q > 0.99).sum(), np.histogram(q, bins=10, range=(0, 1))[0].tolist()))
o = np.argsort(-np.where(node, p, -1.))
for N in (2500, 5000, 10000, 20000):
    lab = np.zeros(n, int); lab[o[:N]] = 1; np.savetxt('$O/top%d.txt' % N, lab, fmt='%d'); print('top %5d: P(class 2) >= %.3f, old closed 15+16 among them %d (%.1f %%; %.1f %% of the node\'s %d)' % (N, p[o[N - 1]], cl[o[:N]].sum(), 100 * cl[o[:N]].mean(), 100 * cl[o[:N]].sum() / max(cl[node].sum(), 1), cl[node].sum()))
start = np.zeros(n, int)
for h in 'AB':
    rows = np.loadtxt('$H' + '_%s/meta_rows.txt' % h).astype(int) - 1; s = np.loadtxt('$O/start_%s.txt' % h, dtype=int); start[rows[s == 2]] = 1
np.savetxt('$O/startbin.txt', start, fmt='%d'); print('start bin: %d particles, old closed %d (%.1f %%)' % (start.sum(), cl[start > 0].sum(), 100 * cl[start > 0].mean()))
PYEOF
: > $O/spec.txt
for sub in startbin top2500 top5000 top10000 top20000; do d=$O/rec_$sub; mkdir -p $d
  if [ -z "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $d && cp -f $PP r.simple
    simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/$sub.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1 \
      && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=1 nparts=4 nthr=16 > rec.log 2>&1 || log "$sub reconstruct FAILED"; cd $O; fi
  RR=$(ls -td $d/*_reconstruct3D 2>/dev/null | head -1); [ -n "$RR" ] && echo "$sub $RR/recvol_state01_even_unfil.mrc $RR/recvol_state01_odd_unfil.mrc" >> $O/spec.txt
done
V=$K_/view_rflex/site_r22${TAG:-}; $PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | sed 's/^/   /'
log DONE
