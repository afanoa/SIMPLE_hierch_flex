#!/bin/bash
#SBATCH -J klem_view
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# viewing folder (lp 10, unsharpened, full maps only): rnasep_klem/view_rflex/<set>/ ; per-set notes in view_rflex/info/<set>.txt
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; V=$K_/view_rflex; I=$V/info; mkdir -p $I
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
# (1) f692gen: the 32 final states of the generic run
G=$K_/f692gen/coordtree_famall; RR=$(ls -td $G/rec/*_reconstruct3D | head -1); : > $I/spec_gen.txt
for s in $(seq 1 32); do b=$RR/recvol_state$(printf %02d $s); echo "state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $I/spec_gen.txt; done
$PY $T/make_view.py $V/f692gen_32states $I/spec_gen.txt
# (2) f692gen node 39 (the largest sampled branch, 203k): its 8 cells in axis order
N=$G/nodes2/node039; NR=$(ls -td $N/rec/*_reconstruct3D | head -1); : > $I/spec_n39.txt
for s in $(seq 1 8); do b=$NR/recvol_state$(printf %02d $s); echo "cell_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $I/spec_n39.txt; done
$PY $T/make_view.py $V/f692gen_node39_track $I/spec_n39.txt
# (3) rflex_f692rf (recursive, real): root states and node r1's leaf states (box 88)
for nd in r r1; do RD=$(ls -td $K_/rflex_f692rf/node_$nd/*_refine3D | head -1); : > $I/spec_rf_$nd.txt
  for s in $(seq 1 8); do b=$RD/recvol_state$(printf %02d $s); echo "${nd}_state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $I/spec_rf_$nd.txt; done
  $PY $T/make_view.py $V/f692rf_node_$nd $I/spec_rf_$nd.txt; done
# (4) synthetic realistic swing (rflex_s30a node r22): 10 coordinate leaves in coordinate order
CR=$(ls -td $K_/rflex_s30a/coord_r22/leaves/*_reconstruct3D | head -1); : > $I/spec_s30a.txt
for s in $(seq 1 10); do b=$CR/recvol_state$(printf %02d $s); echo "track_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $I/spec_s30a.txt; done
$PY $T/make_view.py $V/synthetic_swing_track $I/spec_s30a.txt
# references: consensus and the deposited closed conformer (lp 10)
$PY - <<PYEOF
import numpy as np, os
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
os.makedirs('$V/refs', exist_ok=True); import subprocess
open('$I/spec_refs.txt', 'w').write('consensus $T/maps160/toy_consensus.mrc $T/maps160/toy_consensus.mrc\nclosed_conf18_deposited $CL $CL\n')
PYEOF
$PY $T/make_view.py $V/refs $I/spec_refs.txt
# closed-arm readout of every set (written to info/, not into the map folders)
for d in f692gen_32states f692gen_node39_track f692rf_node_r f692rf_node_r1; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $I/arm_$d.txt $V/$d 2>&1 | grep -v Warning | tail -1 | sed "s/^/$d: /"; done
echo DONE
