#!/bin/bash
#SBATCH -J toy10_flex_node
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# flex_pca on ONE tree node of toy10 (family FAM of lat/expand) = the locally-linear chart; its latent evaluated against
# the truth on the node, next to the GLOBAL flex latent (lat/flex12) restricted to the same particles. env: FAM (1), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; B=$HOME/SIMPLE_kernel/snap/5e1e6aeee
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464 OMP_NUM_THREADS=8
FAM=${FAM:-1}; NP=${NP:-10}; R=$W/lat/flex_node_fam$FAM; SUB=$W/lat/expand/family${FAM}_subset.txt; GT=$W/sim/gt_particles.txt; CONS=$W/proj/1_refine3D/recvol_state01.mrc
STATES=1,10,18,27,28,37,46,48,49,51; AXIS=-0.246,-0.916,-0.442,1.000,0.972,-0.685,-0.489,0.446,0.499,-0.139; log(){ echo "[$(date '+%F %T') fnode$FAM] $*"; }
mkdir -p $R && cd $R && cp -f $W/proj/proj.simple p.simple
$PY - <<PY
import numpy as np
n = max(int(l.split()[0]) for l in open('$GT') if l.strip() and not l.startswith('#')); st = np.zeros(n, int); st[np.loadtxt('$SUB', dtype=int) - 1] = 1
np.savetxt('$R/states.txt', st, fmt='%d'); print('node: %d of %d particles' % (st.sum(), n))
PY
log "selection + flex_pca on the node (neigs 8, lp 20, nparts $NP)"
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$R/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=flex_pca projfile=p.simple vol1=$CONS mskdiam=200 lp=20 neigs=8 npreimages=2 rec_states=no nparts=$NP nthr=8 > flex.log 2>&1 || { log "flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; exit 1; }
FD=$(ls -td $R/*_flex_pca | head -1); log "flex done -> $FD"
echo "== node-local chart (flex on the node), evaluated on the node"; $PY $W/ptcl_latent_eval.py $FD $GT --states $STATES --axis $AXIS --subset $R/states.txt --tag node_local 2>&1 | grep -v Warning | tail -14
echo "== global chart (lat/flex12), evaluated on the same node"; $PY $W/ptcl_latent_eval.py $(ls -td $W/lat/flex12/*_flex_pca | head -1) $GT --states $STATES --axis $AXIS --subset $R/states.txt --tag global_on_node 2>&1 | grep -v Warning | tail -14
log DONE
