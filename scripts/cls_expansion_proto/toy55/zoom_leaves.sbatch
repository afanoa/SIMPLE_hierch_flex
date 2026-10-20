#!/bin/bash
#SBATCH -J klem_zoom
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# zoom on every level-1 leaf state of one family of a finished run (env RTAG, FAM): 4 zooms at a time on this node
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464 OMP_NUM_THREADS=${NTHR:-12}
export W=/mnt/beegfs/elmlund/afan/routes/toy55 PY=/mnt/nasapps/production/python/3.12/bin/python3 FB=${FB:-$HOME/SIMPLE_klem2/build} NTHR=${NTHR:-12} NPZ=${NPZ:-2} SRC=${SRC:-/mnt/beegfs/elmlund/afan/rnasep_coupled/full692k_2d/rnasep.simple} CTFABS=${CTFABS:-1} PSCALE=${PSCALE:-1}
A=/mnt/beegfs/elmlund/afan; R=$A/rnasep_klem/${RTAG:?}; FAM=${FAM:?}; D=$R/fam$FAM; cd $D; log(){ echo "[$(date '+%F %T') zoomfam$FAM] $*"; }
for RB in $D/region*/; do [ -s $RB/level1/clusters.txt ] || continue; for s in $(grep -v "^#" $RB/level1/clusters.txt | awk '{print $1}'); do echo "$R $FAM ${RB%/} $s"; done; done > $D/zoom_list.txt
log "$(wc -l < $D/zoom_list.txt) leaf states to zoom, ${PAR:-4} at a time x $NTHR threads"
xargs -P ${PAR:-4} -L 1 bash $W/zoom_one.sh < $D/zoom_list.txt
log "zooms finished: $(ls -d $D/zoomregion*/ 2>/dev/null | wc -l) dirs, $(ls $D/zoomregion*/level1/states_dir/state_*.mrc 2>/dev/null | wc -l) sub-state maps"; touch $D/ZOOM_DONE
