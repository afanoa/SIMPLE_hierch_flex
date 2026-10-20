#!/bin/bash
# f692 docked-node stage (run after the chain): the closed population is 7 conformers in the paper, so the docked
# node is enumerated as a multi-state node, not purified to one contact: (1) the region-board level-1 state of family
# FAM that matches conformer 18 best (identification readout; the state itself came from the data), (2) a within-state
# Fortran board with the inherited region as focus, (3) states = node_states clusters, reconstructed at the refined
# poses, identified against the deposited maps, scored for closed-arm occupancy, (4) the whole-reference local pose
# step on the largest closed state (Fortran board, 1x1, pose_every 2) and its re-imported reconstruction.
# env: R, FAM (1), META (meta_fb), SIM (0.7), NMAX (30000)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
FAM=${FAM:-1}; D=$R/fam$FAM; PS=$(ls -td $D/refine/*_refine3D_states | head -1)/p.simple; O=$D/docked; mkdir -p $O; log(){ echo "[$(date '+%F %T') dock] $*"; }
NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); LPSR="12 12 10 10 8 8 8 8 8"; SIGR="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"
# ---- (1) the docked state: best match to conformer 18 over the region boards' level-1 states
best=$($PY - <<PYEOF
import glob, numpy as np
rows = []
for f in glob.glob('$D/region*/level1/state_matches.txt'):
    for l in open(f):
        if l.strip() and not l.startswith('#'):
            w = l.split()
            if len(w) >= 3 and int(float(w[1])) == 18: rows.append((float(w[2]), f.replace('/level1/state_matches.txt', ''), int(float(w[0]))))
rows.sort(reverse=True); print('%s %d %.2f' % (rows[0][1], rows[0][2], rows[0][0]) if rows else 'none 0 0')
PYEOF
); set -- $best; RB=$1; ST=$2; log "docked state: region board $RB state $ST (corr with conformer 18 $3)"; [ "$RB" != none ] || exit 1
awk -v s=$ST '$1==s{print NR}' $RB/level1/states_clusters.txt > $O/subset.txt; n=$(grep -c . $O/subset.txt); log "$n particles"
MK=$(grep -a "^focus_mask" $RB/board/board.txt | awk '{print $2}')
# ---- (2) within-state board with the inherited region, (3) states
[ -s $O/board/klem_latent.txt ] || FOCUSMASK=$MK META=$META bash $W/fboard.sh $O/board $PS $O/subset.txt "$LPSR" "$SIGR" || exit 1
$PY $W/node_states.py $O/board $O/states --nrows $NROWS --sim ${SIM:-0.7} --nmax ${NMAX:-30000} 2>&1 | grep -v Warning | tail -3
KS=$(grep -vc "^#" $O/states/clusters.txt); [ "$KS" -ge 1 ] || { log "no states"; exit 1; }
mkdir -p $O/states/rec && cd $O/states/rec && cp -f $PS p.simple; export OMP_NUM_THREADS=8
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states/states_clusters.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$KS nparts=4 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $O/states/maps $O/states/states_dir; for s in $(seq -f "%02g" 1 $KS); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $O/states/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $O/states/states_dir/state_$s.mrc; done; cd $O
echo "== docked-node states vs the deposited conformers (readout)"; $PY $W/deposited_match.py $O/states/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $O/states/state_matches.txt 2>&1 | grep -v Warning | tail -$((KS + 2))
$PY - <<PYEOF
import numpy as np, glob
from scipy.ndimage import gaussian_filter
W = '$W'
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); cons = read_mrc(W + '/maps160/toy_consensus.mrc'); lp = lambda v: gaussian_filter(v, 2.5); L18, Lc = lp(c18), lp(cons)
dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
st = np.loadtxt('$O/states/states_clusters.txt', dtype=int); print('closed-arm occupancy per docked-node state (1 = the deposited closed map):')
for s in range(1, $KS + 1):
    V = read_mrc('$O/states/maps/recvol_state%02d.mrc' % s); Lv = lp(V); print('  state %d: %6d particles, arm occupancy %.2f' % (s, int((st == s).sum()), (Lv[arm].mean() / Lv[core].mean()) / ref))
PYEOF
log DONE
