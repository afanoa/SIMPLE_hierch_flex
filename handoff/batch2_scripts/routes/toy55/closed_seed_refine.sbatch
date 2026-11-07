#!/bin/bash
#SBATCH -J klem_closed_seed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Diagnostic: is a closed population separable with full-SNR references? Two-state refine3D on the 150k dev subset seeded
# with vol1 = consensus, vol2 = the deposited closed conformer 18 (docked, same frame), then the arm readout of both maps.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; O=$A/rnasep_klem/closed_seed; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') closedseed] $*"; }
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; SUB=$A/rnasep_klem/f692dev/subset_all.txt; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1)
$PY - <<PYEOF
import numpy as np
n = $(grep -c . $SUB); sub = np.loadtxt('$SUB', dtype=int).ravel(); N = int(open('$A/rnasep_klem/f692c/consolidated_0.7/states.txt').read().count('\n')); st = np.zeros(N, int); st[sub - 1] = 1
np.savetxt('$O/states_in.txt', st, fmt='%d'); print('subset: %d particles, all in state 1 to start' % (st > 0).sum())
PYEOF
cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states_in.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "refine3D nstates=2 vol1=consensus vol2=closed(18), refine=prob, 10->8 A, maxits ${MAXITS:-6}, nparts ${NP:-20}"
simple_exec prg=refine3D projfile=p.simple nstates=2 vol1=$W/maps160/toy_consensus.mrc vol2=$CL refine=${REFINE:-prob} lp=${LP:-10} maxits=${MAXITS:-6} mskdiam=232 pgrp=c1 nparts=${NP:-20} nthr=16 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $O/*_refine3D | head -1); mkdir -p $O/maps; cp -f $RD/recvol_state0?.mrc $O/maps/ 2>/dev/null
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); print a[2]}}' > $O/states_out.txt
$PY -c "import numpy as np; b=np.loadtxt('$O/states_out.txt',dtype=int); print('final populations:', np.bincount(b)[1:])"
echo "== arm readout"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm_states.txt $O/maps 2>&1 | grep -v Warning
log DONE
