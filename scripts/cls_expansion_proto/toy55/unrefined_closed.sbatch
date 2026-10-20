#!/bin/bash
#SBATCH -J klem_unrefined_closed
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Does the never-refined quarter of the 932k set (e3 = x = y = 0 in 24_refine3D_auto) hold the missing closed particles?
# Two-state global refine3D (refine=prob) seeded with vol1 = consensus, vol2 = deposited closed conformer 18 (docked),
# on SET=unrefined (the 240k rows with flags 0) or SET=control (a random 240k of the refined rows); closed population + arm readout.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; SET=${SET:?}; O=$A/rnasep_klem/unrefined_closed/$SET; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') unrefined $SET] $*"; }
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1)
$PY - <<PYEOF
import numpy as np
fl = np.loadtxt('$A/rnasep_coupled/flags_refined.txt', dtype=int).ravel(); rng = np.random.default_rng(0)
if '$SET' == 'unrefined': st = (fl == 0).astype(int)
else: idx = rng.choice(np.where(fl == 1)[0], int((fl == 0).sum()), replace=False); st = np.zeros(len(fl), int); st[idx] = 1
np.savetxt('$O/states_in.txt', st, fmt='%d'); print('$SET: %d of %d rows selected, all in state 1 to start' % (st.sum(), len(fl)))
PYEOF
cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states_in.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "refine3D nstates=2 vol1=consensus vol2=closed(18), refine=prob (global), lp ${LP:-10}, maxits ${MAXITS:-6}, nparts ${NP:-20} x 8"
simple_exec prg=refine3D projfile=p.simple nstates=2 vol1=$W/maps160/toy_consensus.mrc vol2=$CL refine=prob lp=${LP:-10} maxits=${MAXITS:-6} mskdiam=232 pgrp=c1 nparts=${NP:-20} nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $O/*_refine3D | head -1); mkdir -p $O/maps; cp -f $RD/recvol_state0?.mrc $O/maps/ 2>/dev/null
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s=0; for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); s=a[2]}; print int(s)}' > $O/states_out.txt
$PY -c "import numpy as np; b=np.loadtxt('$O/states_out.txt',dtype=int); n=np.bincount(b, minlength=3); print('$SET final populations: consensus-seeded %d, closed-seeded %d (%.1f %% closed)' % (n[1], n[2], 100*n[2]/max(n[1]+n[2],1)))"
echo "== arm readout of the two maps"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm_states.txt $O/maps 2>&1 | grep -v Warning
log DONE
