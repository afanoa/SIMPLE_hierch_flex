#!/bin/bash
#SBATCH -J lrec_leaves
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=96G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# leaves of a score-coordinate loop: NB equal-width bins of the coordinate over its P..100-P percentile range (the outer bins take the
# rest) -> reconstruct3D at the loop's poses -> bottom-up tree of the leaves by cross-half difference energy (state_merge.py, complete
# linkage) -> per-leaf truth readout. env LOOP (lrec30 loop label), NB (10), P (1), TAG (swing30)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; KD=$A/rnasep_klem
TAG=${TAG:-swing30}; GT=$T/sim$TAG/gt_particles.txt; L=$KD/lrec30/${LOOP:?}; NB=${NB:-10}; P=${P:-1}; O=$L/leaves${NB}; mkdir -p $O && cd $O
log(){ echo "[$(date '+%F %T') lrec_leaves $LOOP] $*"; }
$PY - <<PYEOF
import numpy as np
c = np.loadtxt('$L/final_coord.txt'); m = np.isfinite(c); lo, hi = np.percentile(c[m], [$P, 100 - $P]); e = np.linspace(lo, hi, $NB + 1)
lab = np.zeros(len(c), int); lab[m] = np.clip(np.searchsorted(e, c[m]) - 1, 0, $NB - 1) + 1; np.savetxt('$O/leaf_states.txt', lab, fmt='%d')
print('leaf populations: %s' % ' '.join(str((lab == v).sum()) for v in range(1, $NB + 1)))
PYEOF
cp -f $L/base.simple r.simple
simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/leaf_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=220 nstates=$NB nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td $O/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*
grep -a "POSTPROCESS: FSC" rec.log | sed 's/.*state\([0-9]*\).bin, 0.5\/0.143 at/ leaf \1 FSC 0.5\/0.143/' | tr '\n' ';'; echo
$PY $T/state_merge.py $O/smerge --rdir $RR --labels $O/leaf_states.txt --gt $GT --hp 40 --lp 10 --link complete --ncut 8,6,4,2 --unfil 2>&1 | grep -v Warning | tail -30
awk -v tag=leaves 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; sq[s] += x * x; ns[s]++; T1 += x; T2 += x * x; NT++}}
  END {for (s = 1; s <= '$NB'; s++) if (ns[s] > 0) {mu = sw[s] / ns[s]; W += sq[s] - sw[s] ^ 2 / ns[s]; printf " %d:%.1f+-%.1f(%d,%.0f%%)", s, mu, sqrt(sq[s] / ns[s] - mu ^ 2), n[s], 100 * ns[s] / n[s]}
    printf "\nLEAVES eta2 %.3f over %d swing particles\n", 1 - W / (T2 - T1 ^ 2 / NT), NT}' $GT $O/leaf_states.txt
log DONE
