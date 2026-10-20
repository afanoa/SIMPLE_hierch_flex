#!/bin/bash
#SBATCH -J rnasep_seed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The direction-bin route on REAL data, as a seed for a region board of f692e family 1 (env REG = region mask):
# bins80 sub-class averages (100 bins, 1,599 sub-classes, bin-frame stacks) -> per-bin principal contrast inside the
# projected region -> sign sync -> particle coordinate -> lattice init; then the region board seeded vs random on the
# family's particles at the refined poses; readouts without truth: cell reproducibility, supported states, closed-arm
# occupancy and deposited-conformer match of the reconstructed states.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; SP=$A/toyL; BX=$A/rnasep_coupled/bins80; R=$A/rnasep_klem/f692e; REG=${REG:?}; RT=$(basename $(dirname $REG))_$(basename $REG .mrc); O=$A/rnasep_klem/seed_rnasep/${RT}${TAG:-}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') rnasep seed $RT] $*"; }
NROWS=932496; export FB=$R/klem2_bin; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1)
[ -s $O/bins.npz ] || $PY - <<PYEOF
import numpy as np
e1, e2, cl = [], [], []
for l in open('$BX/all_binned80.txt'):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); e1.append(float(kv['e1'])); e2.append(float(kv['e2'])); cl.append(int(float(kv.get('class', 0))))
e1, e2, cl = np.radians(e1), np.radians(e2), np.array(cl); v = np.c_[np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)]; nb = cl.max()
C = np.array([v[cl == b].mean(0) for b in range(1, nb + 1)]); C /= np.linalg.norm(C, axis=1, keepdims=True); cnt = np.bincount(cl, minlength=nb + 1)[1:]
np.savez('$O/bins.npz', centres=C, label=cl, count=cnt); print('%d bins from all_binned80.txt, members %d-%d' % (nb, cnt.min(), cnt.max()))
PYEOF
[ -s $O/stitch/avgs.mrcs ] || $PY $W/avg_stack2.py $BX $BX/all_binned80.txt $O/stitch 2>&1 | tail -2
echo "== 2D-route seed of $RT on the real data"; $PY $W/region_contrast_seed.py $O/stitch $BX $O/bins.npz $BX/cavgs_bins.mrc $REG $W/maps160/toy_consensus.mrc $O --npart $NROWS --n1 7 --n2 3 2>&1 | grep -v Warning
PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple; SUB=$R/expand/family1_subset.txt; META=$R/fam1/meta_fb.npz
for I in ${INITS:-file random}; do
  log "region board on family 1 ($(grep -c . $SUB) particles), init $I, focus $RT"
  [ -s $O/board_$I/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=7 N2=3 NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$REG INIT=$I INITFILE=$O/init_positions.txt META=$META bash $W/fboard.sh $O/board_$I $PS $SUB "12 10 10 8 8 8" "1.2 1.0 1.0 0.8 0.7 0.7" > $O/board_$I.out 2>&1 || { log "board FAILED: $(tail -2 $O/board_$I.out | cut -c1-160)"; continue; }
  echo "-- init $I: $(awk '!/^#/ {n++; s+=$3; if($3>0.9) h++} END {printf "%d cells, mean cross-half reproducibility %.3f, %d > 0.9", n, s/n, h}' $O/board_$I/klem_states.txt)"
  L=$O/level1_$I; mkdir -p $L; $PY $W/node_states2.py $O/board_$I $L --nrows $NROWS --sim 0.7 --nmax 30000 --support --maxk 10 > $O/states_$I.out 2>&1; K=$(grep -vc "^#" $L/clusters.txt 2>/dev/null || echo 0); echo "-- init $I: $K supported states: $(grep -v '^#' $L/clusters.txt | awk '{printf "%s(%s) ", $4, $1}')"
  [ "$K" -ge 1 ] || continue; mkdir -p $L/rec $L/states_dir && cd $L/rec && cp -f $PS p.simple && simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$L/states_clusters.txt mkdir=no > sel.log 2>&1 && simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$K nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruction FAILED ($I)"; cd $O; continue; }
  RR=$(ls -td ./*_reconstruct3D | head -1); for s in $(seq -f "%02g" 1 $K); do cp -f $RR/recvol_state$s.mrc $L/states_dir/state_$s.mrc; done; rm -rf $RR/*_part*.mrc; cd $O
  echo "-- init $I: closed-arm occupancy"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm_$I.txt $L/states_dir 2>&1 | grep -v Warning | grep -v mask | sed "s#$L/states_dir/##" | sort -k3 -rn | head -4
  echo "-- init $I: deposited match"; $PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $L/matches.txt 2>&1 | grep -v Warning | grep identified
done
log DONE
