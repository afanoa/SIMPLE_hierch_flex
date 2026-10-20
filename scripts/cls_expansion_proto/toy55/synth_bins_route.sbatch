#!/bin/bash
#SBATCH -J synth_bins
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The direction-bin 2D route on a synthetic benchmark (env MODE): refined poses -> NBINS direction bins -> project with
# class = bin -> bin averages -> cls_expansion sub-classes per bin (pose mean) -> sub-class averages with KNOWN directions
# -> K 3D states by cross-half EM over the averages (avgs_states3d.py, no orientation search) -> particles inherit the
# state -> truth readouts (per-arm angle std ratio / per-site occupancy) for the raw sub-classes and for the K states.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_coupled/snap/9eaef62df; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=8   # cls_expansion pose-mean mode lives in the coupled snapshot
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; SP=$A/toyL; MODE=${MODE:?}; NBINS=${NBINS:-100}; NSUB=${NSUB:-6}; K=${K:-6}; MSK=230; SMPD=1.464
D=$W/binsroute_$MODE; mkdir -p $D; cd $D; log(){ echo "[$(date '+%F %T') binsroute $MODE] $*"; }
GT=$W/sim$MODE/gt_particles.txt; NPART=$(grep -vc '^#' $GT); CONS=$W/maps160/synth_$MODE/consensus_synth.mrc
[ -s $D/support.mrc ] || $PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import write_mrc
B = 160; c = B // 2; zz, yy, xx = np.mgrid[:B, :B, :B]; rr = np.sqrt((xx - c) ** 2. + (yy - c) ** 2. + (zz - c) ** 2.); R = $MSK / 2. / $SMPD
write_mrc('$D/support.mrc', np.clip((R + 3 - rr) / 6., 0, 1).astype(np.float32), $SMPD); print('support sphere written')
PYEOF
[ -s $D/all_binned.txt ] || $PY $SP/dir_bins.py $W/lat_$MODE/ptcl3D_r0.txt $NBINS $D/all_binned.txt $D/bins.npz || { log "binning FAILED"; exit 5; }
if [ ! -s $D/p/p.simple ]; then
  simple_exec prg=new_project projname=p > newproj.log 2>&1; cd $D/p
  simple_exec prg=import_particles projfile=p.simple stk=$W/sim$MODE/all.mrcs oritab=$D/all_binned.txt smpd=$SMPD kv=300 cs=2.7 fraca=0.1 ctf=yes mkdir=no > import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 5; }
fi
cd $D/p
if [ ! -s $D/p/cls_expansion_labels.txt ]; then
  log "make_cavgs: $NBINS bin averages at the refined in-plane poses"; simple_exec prg=make_cavgs projfile=p.simple ncls=$NBINS mskdiam=$MSK refs=cavgs_bins.mrc nparts=1 nthr=32 mkdir=no > MAKE_CAVGS_OUTPUT 2>&1 || { log "make_cavgs FAILED: $(grep -a -m2 'ERROR\|THROW' MAKE_CAVGS_OUTPUT | cut -c1-160)"; exit 6; }
  log "cls_expansion: pose mean, $NSUB sub-classes per bin, lp 20"; simple_exec prg=cls_expansion projfile=p.simple ncls=$NSUB lp=20 mskdiam=$MSK vol1=$CONS pcg_mskfile=$D/support.mrc nthr=32 mkdir=no > CLS_EXPANSION_OUTPUT 2>&1 || { log "cls_expansion FAILED: $(grep -a -m2 'ERROR\|THROW\|fopen' CLS_EXPANSION_OUTPUT | cut -c1-160)"; exit 7; }
fi
log "sub-classes: $(grep -vc '^#' cls_expansion_labels.txt) labelled particles, $(awk '!/^#/{print $3}' cls_expansion_labels.txt | sort -u | wc -l) sub-classes"
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$D/p/cls_expansion_labels.txt', comments='#').astype(int); pl = np.zeros($NPART, int); pl[L[:, 0] - 1] = L[:, 2]; np.savetxt('$D/subclass_labels.txt', pl, fmt='%d')
PYEOF
echo "== truth readout of the raw bin sub-classes"; $PY $W/label_truth_readout.py $D/subclass_labels.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -4
$PY $SP/avg_stack.py $D/p $W/lat_$MODE/ptcl3D_r0.txt $D/stitch || { log "avg_stack FAILED"; exit 8; }
for k in $K ${K2:-12}; do
  log "3D states from the averages: K=$k"; $PY $W/avgs_states3d.py $D/stitch $D/p $CONS $D/states_k$k --k $k --npart $NPART 2>&1 | grep -v Warning | grep -v "^it " ; grep "^it " $D/states_k$k/log.txt 2>/dev/null | tail -1
  echo "== truth readout of the K=$k states"; $PY $W/label_truth_readout.py $D/states_k$k/particle_states.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -$((k + 3))
done
log "DONE"
