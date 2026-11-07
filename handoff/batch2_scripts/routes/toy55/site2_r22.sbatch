#!/bin/bash
#SBATCH -J klem_site2_r22
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=110G
#SBATCH -t 02:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# OPEN-site occupancy of real r22 + the four corners with the closed-site posterior (site_r22): open-site mask = r22's cheap-children
# variation mask (both arm sites) without unit 7's closed-site mask (mask_minus.py); 2 classes at fixed poses that differ only inside it,
# started as consensus vs consensus with the site removed (minor_em4.py --empty_init), cross-half, soft EM. Corners per particle
# (independent sites): closed-occupied x open-empty (the closed state proper), both occupied, open only, neither; populations; top-N of
# the closed corner and of open-empty alone reconstructed with SIMPLE at r22's parent poses -> closed-arm / open-arm readout.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
MR=$K_/minor_em_r22; S1=$K_/site_r22; PP=$K_/rflex_f692rf/node_r2/2_refine3D/r.simple; O=$K_/site2_r22; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') site2 r22] $*"; }
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
[ -s $O/open_mask.mrc ] || $PY $T/mask_minus.py $K_/focused_arm_f692rf_r22/arm_mask.mrc $K_/rflex_f692u7/node_r/rsmask.mrc $O/open_mask.mrc --cons $CONS --closed $CL 2>&1 | grep -v Warning
for h in A B; do [ -s $O/node_$h.txt ] || awk '{print ($1 > 0) ? 1 : 0}' $MR/seed_$h.txt > $O/node_$h.txt; done
[ -s $O/em_open/post_soft_B.txt ] || $PY -u $T/minor_em4.py --sim ${H}_A,${H}_B --poses $MR/poses_A.txt,$MR/poses_B.txt --labels $O/node_A.txt,$O/node_B.txt --cache $T/benchcache_rnasep692 \
  --ctfabs 1 --cons_ref $CONS --mask $O/open_mask.mrc --empty_init 1 --soft 1 --iters ${ITERS:-10} --em_lim ${EMLIM:-15} --post 1 --nomaps 1 --out $O/em_open > $O/em_open.log 2>&1 || { log "EM FAILED: $(tail -2 $O/em_open.log)"; exit 1; }
grep -a '^start\|^iteration' $O/em_open.log | sed 's/deviation FSC.*//' | cut -c1-160
$PY - <<PYEOF
import numpy as np
n = 932496; pc = np.full(n, -1.); po = np.full(n, -1.); old = np.loadtxt('$K_/f692r/coordtree_fam1/units.txt', dtype=int); cl = (old == 15) | (old == 16)
for h in 'AB':
    rows = np.loadtxt('$H' + '_%s/meta_rows.txt' % h).astype(int) - 1; m = np.loadtxt('$O/node_%s.txt' % h, dtype=int) > 0
    pc[rows[m]] = np.loadtxt('$S1/em/post_soft_%s.txt' % h)[m, 1]; po[rows[m]] = np.loadtxt('$O/em_open/post_soft_%s.txt' % h)[m, 1]
node = (pc >= 0) & (po >= 0); c, e = pc[node], po[node]
C = {'closed (closed site occupied, open site empty)': c * e, 'both sites occupied': c * (1 - e), 'open (open site occupied only)': (1 - c) * (1 - e), 'neither': (1 - c) * e}
print('r22 %d particles; P(closed site occupied) mean %.3f, P(open site empty) mean %.3f, correlation of the two %+.3f' % (node.sum(), c.mean(), e.mean(), np.corrcoef(c, e)[0, 1]))
for k, v in C.items(): print('   corner %-48s expected population %7.0f (%.1f %%), particles with P > 0.5: %6d' % (k, v.sum(), 100 * v.mean(), (v > 0.5).sum()))
print('   P(open site empty): > 0.5 %d, > 0.9 %d; histogram (10 bins) %s' % ((e > 0.5).sum(), (e > 0.9).sum(), np.histogram(e, bins=10, range=(0, 1))[0].tolist()))
score = np.full(n, -1.); score[node] = c * e; o = np.argsort(-score); so = np.full(n, -1.); so[node] = e; oo = np.argsort(-so)
for N in (2500, 5000, 10000):
    lab = np.zeros(n, int); lab[o[:N]] = 1; np.savetxt('$O/corner%d.txt' % N, lab, fmt='%d')
    print('closed corner top %5d: P >= %.3f, old closed 15+16 among them %d (%.1f %%)' % (N, score[o[N - 1]], cl[o[:N]].sum(), 100 * cl[o[:N]].mean()))
lab = np.zeros(n, int); lab[oo[:5000]] = 1; np.savetxt('$O/openempty5000.txt', lab, fmt='%d'); print('open-empty top 5000: P >= %.3f, old closed %d' % (so[oo[4999]], cl[oo[:5000]].sum()))
PYEOF
: > $O/spec.txt
for sub in corner2500 corner5000 corner10000 openempty5000; do d=$O/rec_$sub; mkdir -p $d
  if [ -z "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $d && cp -f $PP r.simple
    simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/$sub.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1 \
      && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=1 nparts=4 nthr=16 > rec.log 2>&1 || log "$sub reconstruct FAILED"; cd $O; fi
  RR=$(ls -td $d/*_reconstruct3D 2>/dev/null | head -1); [ -n "$RR" ] && echo "$sub $RR/recvol_state01_even_unfil.mrc $RR/recvol_state01_odd_unfil.mrc" >> $O/spec.txt
done
V=$K_/view_rflex/site2_r22; $PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $CONS $CL $O/arm.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | sed 's/^/   /'
log DONE
