#!/bin/bash
#SBATCH -J subcls_maps
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Child maps of a per-row labelling at the parent's poses (no refinement): selection (state = child) -> reconstruct3D nstates = K with
# even/odd halves -> subcls_mapscore.py (truth: per child resolution / swing; end and neighbour difference maps vs the true difference).
# REFS scores existing reconstructions the same way (+-separated name:recdir:labels). env PP, LAB, OUT, TAG (swing30), MSKDIAM (220), REFS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55
PP=${PP:?}; LAB=${LAB:?}; O=${OUT:?}; TAG=${TAG:-swing30}; MD=${MSKDIAM:-220}; GT=$T/sim$TAG/gt_particles.txt; TM=$T/sim$TAG/maps
mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') subcls_maps $(basename $O)] $*"; }
$PY -c "import numpy as np; l = np.loadtxt('$LAB').astype(int); u = np.unique(l[l > 0]); m = np.zeros(l.max() + 1, int); m[u] = np.arange(1, len(u) + 1); np.savetxt('$O/labels.txt', m[l], fmt='%d')"
LAB=$O/labels.txt; K=$(awk '$1 > 0' $LAB | sort -n | uniq | wc -l); t0=$(date +%s)   # children relabelled 1..K in order
if [ -z "$(ls $O/*_reconstruct3D/recvol_state$(printf %02d $K)_odd_unfil.mrc 2>/dev/null)" ]; then
  cp -f $PP r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$LAB mkdir=no > rsel.log 2>&1 \
    && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=12 mkdir=no > rup.log 2>&1 || { log "project FAILED"; exit 1; }
  log "reconstruct3D: $K children"; simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$K nparts=4 nthr=12 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }
fi
RR=$(ls -td $O/*_reconstruct3D | head -1); log "reconstruction $(( $(date +%s) - t0 )) s"
$PY $T/subcls_mapscore.py $RR $LAB $GT $TM --mskdiam $MD --tag "$(basename $O)" 2>&1 | grep -v Warning
for r in $(echo ${REFS:-} | tr "+" " "); do nm=${r%%:*}; rest=${r#*:}; rd=${rest%%:*}; rl=${rest#*:}
  $PY $T/subcls_mapscore.py $rd $rl $GT $TM --mskdiam $MD --tag "$nm" 2>&1 | grep -v Warning; done
log DONE
