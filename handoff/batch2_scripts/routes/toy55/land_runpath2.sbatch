#!/bin/bash
#SBATCH -J land_runpath2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_runpath2.%j.out
# data-derived path per run: the two 10 % tails of the run's own focused latent (PC1 of its posterior means) -> SIMPLE maps at the
# shared poses -> OT morph inside the node focus mask -> pose-marginalised r1 scores -> reweighting + basins. args: NAME COORDS (flex_pca_coordinates.txt)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 SIMPLE_EMAIL=nil
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; L=$K_/landscape; H=$T/simrnasep692
NAME=${1:?}; CO=${2:?}; O=$L/runpath2_$NAME; V=$K_/view_rflex/runpath2_r1_$NAME; PROJ=$K_/rflex_f692st/node_r1/fstageB/2_refine3D/b.simple
mkdir -p $O $V/halves && cd $O; log(){ echo "[$(date '+%F %T') runpath2 $NAME] $*"; }
$PY - <<PYEOF
import numpy as np
C = np.loadtxt('$CO', comments='#'); rows = C[:, 0].astype(int) - 1; z = C[:, 5:] - C[:, 5:].mean(0)
w, V = np.linalg.eigh(np.cov(z.T)); u = V[:, -1]; s = z @ u; q = np.quantile(s, [0.10, 0.90]); lab = np.zeros(932496, int)
lab[rows[s <= q[0]]] = 1; lab[rows[s >= q[1]]] = 2; np.savetxt('$O/tails.txt', lab, fmt='%d'); np.savetxt('$O/pc1.txt', np.c_[rows + 1, s], fmt='%d %.5f')
print('PC1 of the posterior means: %.1f %% of the latent variance; tails %d / %d particles' % (100 * w[-1] / w.sum(), (lab == 1).sum(), (lab == 2).sum()))
PYEOF
log "tails written"
if [ ! -s $V/state_02.mrc ]; then cp -f $PROJ $O/r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/tails.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=1 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=2 nparts=1 nthr=8 > rec.log 2>&1 || { log "reconstruct FAILED: $(grep -a -h -m2 ERROR sel.log up.log rec.log)"; exit 1; }
  RR=$(ls -td $O/*_reconstruct3D | head -1); : > spec.txt
  for s in 1 2; do b=$RR/recvol_state0$s; echo "state_0$s ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> spec.txt; cp -f ${b}_even_unfil.mrc $V/halves/state_0${s}_even.mrc; cp -f ${b}_odd_unfil.mrc $V/halves/state_0${s}_odd.mrc; done
  $PY $T/make_view.py $V spec.txt --lp 10 --box 160 > /dev/null 2>&1; fi
log "tail maps: $(ls $V/state_0?.mrc | tr '\n' ' ')"
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $V 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | sed 's#^#  #'
[ -s $O/morph/path_30.mrc ] || $PY -u $T/land_morph2.py $V/state_01.mrc $V/state_02.mrc $O/morph --mask $K_/rflex_f692st/node_r1/rsmask_variance.mrc --mask_thr 0.2 2>&1 | grep --line-buffered -v Warning
ls $O/morph/path_*.mrc > $O/maps.txt
[ -s $O/ps.npz ] || $PY -u $T/land_pscore.py --sim ${H}_A,${H}_B --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --cache $K_/lsamp/caches/rnasep692 \
   --poses $K_/lsamp/f692st_ps/dump_db61decd1108.txt --rowflags $K_/rflex_f692st/rows_r1.txt --ctfabs 1 --mskdiam 232 --maps $O/maps.txt --out $O/ps.npz 2>&1 | grep --line-buffered -v "Warning\|warn(\|block"
for lam in 0.001 0.01; do log "basins lam $lam"; $PY -u $T/land_basins.py --scores $O/ps.npz --lam $lam --thr 0.5,0.6,0.75 2>&1 | grep --line-buffered -v Warning; done
$PY -u $T/land_reweight.py --scores $O/ps.npz --out $O/rw_lmarg_split_0.001 --real --variant lmarg --tau split --lam 0.001 --nm 932496 --nboot 20 2>&1 | grep --line-buffered -v "Warning\|warn("
log "DONE: $O, tail maps $V"
