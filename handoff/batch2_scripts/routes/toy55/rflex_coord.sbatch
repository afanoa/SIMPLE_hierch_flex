#!/bin/bash
#SBATCH -J klem_rflex_coord
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# continuous-coordinate leaves from labelled particles (a K-state classification): reconstruct the labels' maps -> score every
# labelled particle against them (refine3D, 1 iteration, prob_neigh geom; its pose/state output is discarded) -> per-particle best
# dist to each map -> PC1 of the centred score vectors (lrec_pc.py) -> KL quantile bins = next labels -> ... NIT times; the last
# bins are reconstructed as the leaves (final/). Poses stay those of BASE throughout.
# env SRC (lrec30 run label: BASE = its refined project, LAB0 = its final_states.txt), LABEL, BASE, LAB0, KL (8), NIT (2), LP (8),
# BOXC SMPDC (88 / 2.66181803: the classification's crop; full box needs >300 GB), XHALF (1: score each half against the OTHER half's maps = copies with _even/_odd swapped; no self-inclusion), EXTRA (refine3D keys)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; KD=$A/rnasep_klem
TAG=${TAG:-swing30}; GT=$T/sim$TAG/gt_particles.txt; O=${OUT:?}; LABEL=$(basename $O); KL=${KL:-8}; NIT=${NIT:-2}; LP=${LP:-8}
mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') lrec_loop $LABEL] $*"; }
if [ -n "${SRC:-}" ]; then S=$KD/lrec30/$SRC; RD=$(ls -td $S/*_refine3D_states | head -1); BASE=${BASE:-$RD/p.simple}; LAB0=${LAB0:-$S/final_states.txt}; fi
[ -s "${BASE:?}" ] && [ -s "${LAB0:?}" ] || { log "missing BASE $BASE or LAB0 $LAB0"; exit 1; }
cp -f $BASE base.simple; cp -f $LAB0 states_it0.txt; NROWS=$(wc -l < states_it0.txt)
recon(){ # $1 dir, $2 labels -> sets VOLS, K
  mkdir -p $1 && cd $1; cp -f $O/base.simple r.simple; K=$(awk '$1 > 0' $2 | sort -n | uniq | wc -l)
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$2 mkdir=no > rsel.log 2>&1 || { log "selection FAILED in $1"; exit 1; }
  simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > rup.log 2>&1
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=${MSKDIAM:-220} nstates=$K nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED in $1"; exit 1; }
  RR=$(ls -td $1/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*
  VOLS=""; for s in $(seq 1 $K); do VOLS="$VOLS vol$s=$RR/recvol_state$(printf %02d $s).mrc"; done; }
for it in $(seq 0 $NIT); do
  W=$O/it$it; recon $W $O/states_it$it.txt
  [ $it = 0 ] && [ -s $RR/r.simple ] && cp -f $RR/r.simple $O/base.simple   # keep the e/o partition the first reconstruction made (raw projects have none)
  cd $W; cp -f $O/base.simple s.simple
  simple_exec prg=selection projfile=s.simple oritype=ptcl3D infile=$O/states_it$it.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED it $it"; exit 1; }
  simple_exec prg=update_project projfile=s.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  if [ "${XHALF:-0}" = 1 ]; then XR=$W/xref; mkdir -p $XR; XV=""; for v in $VOLS; do f=${v#*=}; b=$(basename $f .mrc); cp -f $f $XR/$b.mrc; cp -f ${f%.mrc}_odd.mrc $XR/${b}_even.mrc; cp -f ${f%.mrc}_even.mrc $XR/${b}_odd.mrc; XV="$XV ${v%%=*}=$XR/$b.mrc"; done; VOLS=$XV; fi
  log "it $it: scoring all labelled particles against $K maps (lp $LP, cross-half ${XHALF:-0})"
  simple_exec prg=refine3D projfile=s.simple $VOLS nstates=$K refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=${MSKDIAM:-220} lp=$LP maxits=1 center=no box_crop=${BOXC:-88} smpd_crop=${SMPDC:-2.66181803} nparts=4 nthr=16 ${EXTRA:-} > score.log 2>&1 || { log "scoring refine3D FAILED: $(grep -a -m2 'ERROR' score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $W/*_refine3D | head -1)
  $PY $T/lrec_pc2.py $W/pc --tabs "$SD/dist_part_neigh_*.dat" --nstates $K --nrows $NROWS --k $KL $([ -n "${TAG:-}" ] && echo --gt $GT) --npc ${NPC:-2} $([ $it -gt 0 ] && echo --prev $O/it$((it - 1))/pc/coord.txt) 2>&1 | sed "s/^/  it $it /"
  [ -s $W/pc/bins.txt ] || { log "lrec_pc.py FAILED it $it"; exit 1; }
  rm -f $SD/*_part*_even.mrc $SD/*_part*_odd.mrc $SD/rho_* $SD/*.mrc; rm -rf $W/xref; cp -f $W/pc/bins.txt $O/states_it$((it + 1)).txt
done
recon $O/final $O/states_it$((NIT + 1)).txt; cp -f $O/states_it$((NIT + 1)).txt $O/final_states.txt; cp -f $O/it$NIT/pc/coord.txt $O/final_coord.txt
$PY - <<PYEOF
import numpy as np
W = '$O/it$NIT/pc'; Z = [np.loadtxt('%s/pc%d.txt' % (W, j)) for j in range(1, ${NPC:-2} + 1)]; m = np.isfinite(Z[0])
def gain(x):   # 2-component Gaussian mixture log-likelihood gain per particle over one Gaussian (EM, standardised)
    x = (x - x.mean()) / x.std(); mu = np.percentile(x, [25, 75]); sd = np.array([1., 1.]); pi = np.array([.5, .5])
    for it in range(200):
        p = pi * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / sd; r = p / p.sum(1, keepdims=True); nk = r.sum(0)
        pi = nk / len(x); mu = (r * x[:, None]).sum(0) / nk; sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(0) / nk) + 1e-6
    ll2 = np.log((pi * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))).sum(1)).mean(); ll1 = -0.5 * np.log(2 * np.pi) - 0.5
    return ll2 - ll1
g = [gain(z[m]) for z in Z]; j = int(np.argmin(g)); c = Z[j]
print('motion coordinate = PC%d (2-Gaussian gains %s: the smallest = the continuous axis)' % (j + 1, ' '.join('%.3f' % v for v in g)))
lo, hi = np.percentile(c[m], [${P:-1}, 100 - ${P:-1}]); e = np.linspace(lo, hi, ${NB:-10} + 1); lab = np.zeros(len(c), int)
lab[m] = np.clip(np.searchsorted(e, c[m]) - 1, 0, ${NB:-10} - 1) + 1; np.savetxt('$O/leaf_states.txt', lab, fmt='%d'); np.savetxt('$O/motion_coord.txt', c, fmt='%.6g')
print('leaf populations: %s' % ' '.join(str((lab == v).sum()) for v in range(1, ${NB:-10} + 1)))
PYEOF
recon $O/leaves $O/leaf_states.txt
if [ -n "${TAG:-}" ]; then
  $PY $T/scorecard.py $O/leaf_states.txt $GT --codes ${CODES:-swing:200:1} --tag "coordinate leaves" 2>&1 | grep -v Warning
  awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; ns[s]++}} END {for (s = 1; s <= 20; s++) if (n[s]) printf " %d:%.1f(%d,%.0f%%)", s, sw[s] / (ns[s] > 0 ? ns[s] : 1), n[s], 100 * ns[s] / n[s]; print ""}' $GT $O/leaf_states.txt | sed 's/^/per leaf mean swing (particles, swing share):/'
fi
log "DONE: leaves = $O/final (maps), labels final_states.txt, coordinate final_coord.txt"
