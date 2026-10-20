#!/bin/bash
#SBATCH -J klem_unrefined_board
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Closed-state content of the never-refined quarter, with the chain's own machinery (the two-reference refinement drained
# 90 % into the sharper deposited reference in the unrefined set AND the control, so it measured nothing): one global
# single-reference refine3D (consensus) of the 240k rows for poses, then the arm-focused klem board (free cells, no
# reference), supported states reconstructed, closed-arm occupancy and population per state.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; O=$A/rnasep_klem/unrefined_closed/board; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') unrefboard] $*"; }
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); MK=$A/rnasep_klem/f692e/fam1/regions_motion/region_1.mrc; NROWS=932496
$PY - <<PYEOF
import numpy as np
fl = np.loadtxt('$A/rnasep_coupled/flags_refined.txt', dtype=int).ravel(); st = (fl == 0).astype(int); np.savetxt('$O/states_in.txt', st, fmt='%d'); np.savetxt('$O/subset.txt', np.where(st == 1)[0] + 1, fmt='%d'); print('%d never-refined rows' % st.sum())
PYEOF
if [ ! -s $O/refined.simple ]; then
  cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states_in.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  log "refine3D single reference (consensus), refine=prob global, lp 10, maxits 6, nparts 20 x 8"
  simple_exec prg=refine3D projfile=p.simple vol1=$W/maps160/toy_consensus.mrc refine=prob lp=10 maxits=6 mskdiam=232 pgrp=c1 nparts=20 nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
  cp -f $(ls -td $O/*_refine3D | head -1)/p.simple $O/refined.simple
fi
log "arm-focused board on the 240k at their new poses (7x7 + 1 free, focus $(basename $MK))"
[ -s $O/board/klem_latent.txt ] || FB=$A/rnasep_klem/f692e/klem2_bin SERVER_LIST= NPARTS= NTHR=16 N1=7 N2=7 NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$MK INIT=random bash $W/fboard.sh $O/board $O/refined.simple $O/subset.txt "12 10 10 8 8 8" "1.2 1.0 1.0 0.8 0.7 0.7" > $O/board.out 2>&1 || { log "board FAILED: $(tail -2 $O/board.out | cut -c1-160)"; exit 1; }
L=$O/level1; mkdir -p $L; $PY $W/node_states2.py $O/board $L --nrows $NROWS --sim 0.7 --nmax 100000 --support --maxk 10 > $O/states.out 2>&1; K=$(grep -vc "^#" $L/clusters.txt 2>/dev/null || echo 0); log "$K supported states: $(grep -v '^#' $L/clusters.txt | awk '{printf "%s(%s) ", $4, $1}')"
[ "$K" -ge 1 ] || { log "no supported states"; exit 0; }
mkdir -p $L/rec $L/maps $L/states_dir && cd $L/rec && cp -f $O/refined.simple p.simple && simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$L/states_clusters.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$K nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); for s in $(seq -f "%02g" 1 $K); do cp -f $RR/recvol_state$s.mrc $L/states_dir/state_$s.mrc; done
echo "== closed-arm occupancy of the never-refined quarter's states (population in clusters.txt)"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm_states.txt $L/states_dir 2>&1 | grep -v Warning | grep -v mask
$PY - <<PYEOF
import numpy as np, re
occ = {int(re.search(r'state_(\d+)', l.split()[0]).group(1)): float(l.split()[1]) for l in open('$O/arm_states.txt') if not l.startswith('#')}
pop = {int(l.split()[0]): int(l.split()[2]) for l in open('$L/clusters.txt') if not l.startswith('#')}
tot = sum(pop.values()); cl = sum(pop[s] for s in pop if occ.get(s, 0) >= 0.8); print('closed-arm occupancy >= 0.8: %d of %d particles in supported states (%.1f %%) of the never-refined quarter' % (cl, tot, 100 * cl / max(tot, 1)))
PYEOF
log DONE
