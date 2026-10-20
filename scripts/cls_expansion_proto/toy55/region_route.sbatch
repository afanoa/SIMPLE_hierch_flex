#!/bin/bash
#SBATCH -J region_route
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the route that recovered the closed state before (region board -> leaf -> zoom), inside one node of the real gold-standard tree,
# for ONE half (env H, TAG gs1, NODE r.2): (0) motion regions (derive_regions2) of the tree's own whole-particle board on the
# node, the region overlapping the arm-travel mask most; (1) board focused on that region, supported level-1 states,
# reconstructed, closed-arm readout + share of the earlier closed sub-state particles; (2) zoom board (7x7 + 1, same focus)
# on the most closed level-1 state, supported sub-states, reconstructed, same readouts. |CTF| and pscale throughout.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464 OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; H=${H:?}; TAG=${TAG:-gs1}; NODE=${NODE:-r.2}; MODE=rnasep692
R=$K/synth_${MODE}_$H; PS=$R/src.simple; META=$R/meta_toy.npz; NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); T=$K/dtree_${MODE}_$TAG
Z=$T/region_$(echo $NODE | tr '.' '_')_$H; mkdir -p $Z; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; CLOSED=$W/maps160/docked_all/conf18_emd74209_m160.mrc; ARM=$T/arm_mask.mrc
PROXY="$K/f692c/fam1/zoomregionmotion_1_s05/level1/states_clusters.txt 1:3 $W/simrnasep692_$H/meta_rows.txt"; export SERVER_LIST= NPARTS= NTHR=16 CTFABS=1 PSCALE=1
log(){ echo "[$(date '+%F %T') region $NODE $H] $*"; }
rec(){ d=$1; sf=$2; ns=$3; mkdir -p $d/rec $d/states_dir && cd $d/rec && cp -f $PS p.simple && simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=1 nthr=16 mkdir=no > up.log 2>&1 && \
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$sf mkdir=no > sel.log 2>&1 && simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=1 nthr=16 > rec.log 2>&1 || { log "rec FAILED in $d"; return 1; }
  RR=$(ls -td ./*_reconstruct3D | head -1); for s in $(seq -f "%02g" 1 $ns); do cp -f $RR/recvol_state$s.mrc $d/states_dir/state_$s.mrc; done; rm -f $RR/*part*.mrc; cd $Z; }
readout(){ d=$1; $PY $W/arm_readout.py $CONS $CLOSED $d/arm_readout.txt $d/states_dir 2>&1 | grep -v Warning | grep -v "^closed-arm mask" | sed "s#$d/states_dir/##"; $PY $W/proxy_in_labels.py $PROXY $d/states_clusters.txt "known-closed share" 2>&1 | grep -v Warning; }
# (0) motion regions from the tree's own board on the node
mkdir -p $Z/b0; for f in $R/$TAG/$NODE/board/*; do ln -sfn $(readlink -f $f) $Z/b0/$(basename $f); done
[ -s $Z/b0/klem_halves.npz ] || $PY $W/klem_fortran2npz.py $Z/b0 2>&1 | grep -v Warning | tail -2
[ -s $Z/regions_motion/region_1.mrc ] || $PY $W/derive_regions2.py $Z/b0 $Z/regions_motion --k 5 --ref $ARM 2>&1 | grep -v Warning | tail -8
MK=$($PY - <<PYEOF
import sys, glob; sys.path.insert(0, '$W'); from model_class import read_mrc
import numpy as np
arm = read_mrc('$ARM'); best = None
for f in sorted(glob.glob('$Z/regions_motion/region_*.mrc')):
    m = read_mrc(f); ov = float(np.minimum(m, arm).sum() / max(arm.sum(), 1e-9)); print('%s: covers %.2f of the arm-travel mask, %.1f %% of the box' % (f.split('/')[-1], ov, 100 * (m > 0.5).mean()), file=sys.stderr)
    if best is None or ov > best[1]: best = (f, ov)
print(best[0] if best else '')
PYEOF
); [ -n "$MK" ] || { log "no motion regions"; exit 1; }; log "focus region $(basename $MK)"
# (1) region board on the node, level-1 states
[ -s $Z/board/klem_latent.txt ] || NFREE=1 FOCUSMASK=$MK META=$META bash $W/fboard.sh $Z/board $PS $R/$TAG/$NODE/subset.txt "12 12 10 10 8 8 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" > $Z/board.out 2>&1 || { log "region board FAILED: $(tail -2 $Z/board.out | cut -c1-200)"; exit 1; }
$PY $W/node_states2.py $Z/board $Z/level1 --nrows $NROWS --sim 0.7 --nmax 30000 --maxk 30 --support > $Z/level1.out 2>&1; KS=$(grep -vc "^#" $Z/level1/clusters.txt 2>/dev/null || echo 0)
log "level 1: $KS supported states: $(grep -v '^#' $Z/level1/clusters.txt 2>/dev/null | awk '{printf "%s(%s) ", $4, $1}')"; [ "$KS" -ge 1 ] || exit 0
[ -s $Z/level1/states_dir/state_01.mrc ] || rec $Z/level1 $Z/level1/states_clusters.txt $KS || exit 1
echo "== level-1 states (region-focused board on $NODE, half $H)"; readout $Z/level1
# (2) zoom on the most closed level-1 state
S=$(grep -v "^#" $Z/level1/arm_readout.txt | sort -k2 -g -r | head -1 | awk '{print $1}' | grep -o "state_[0-9]*" | grep -o "[0-9]*$" | sed 's/^0*//'); log "zoom on level-1 state $S"
L=$Z/zoom_s$S; mkdir -p $L; awk -v s=$S '$1 == s {print NR}' $Z/level1/states_clusters.txt > $L/subset.txt
[ -s $L/board/klem_latent.txt ] || N1=7 N2=7 NFREE=1 FOCUSMASK=$MK META=$META bash $W/fboard.sh $L/board $PS $L/subset.txt "12 10 10 8 8 8" "1.2 1.0 1.0 0.8 0.7 0.7" > $L/board.out 2>&1 || { log "zoom board FAILED"; exit 1; }
$PY $W/node_states2.py $L/board $L/level1 --nrows $NROWS --sim 0.7 --nmax 10000 --maxk 10 --support > $L/level1.out 2>&1; KZ=$(grep -vc "^#" $L/level1/clusters.txt 2>/dev/null || echo 0)
log "zoom: $KZ supported sub-states: $(grep -v '^#' $L/level1/clusters.txt 2>/dev/null | awk '{printf "%s(%s) ", $4, $1}')"; [ "$KZ" -ge 1 ] || exit 0
[ -s $L/level1/states_dir/state_01.mrc ] || rec $L/level1 $L/level1/states_clusters.txt $KZ || exit 1
echo "== zoom sub-states of level-1 state $S (half $H)"; readout $L/level1
log DONE
