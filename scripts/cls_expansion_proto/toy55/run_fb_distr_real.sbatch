#!/bin/bash
#SBATCH -J rnasep_fb_distr
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Real-size validation + timing of the distributed board: the f692 family-1 feedback board (307k particles, focus mask,
# latent init, 14 -> 8 A, 17 min sequential at 16 threads) rerun with NP workers x NTHRP threads; identification vs the
# deposited conformers compared with the sequential board's. env: NP (6), NTHRP (8), RTAG (f692)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
FB=$HOME/SIMPLE_klem2/build; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692}; S=$R/fam1/klem_fb; T=$R/fam1/klem_fb_distr${NP:-6}; mkdir -p $T; cd $T; log(){ echo "[$(date '+%F %T') fbd] $*"; }
grep -v "^nthr_part" $S/board.txt > $T/board.txt; echo "nthr_part ${NTHRP:-8}" >> $T/board.txt; PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple
log "distributed feedback board: $(grep -c . $R/expand/family1_subset.txt) particles, ${NP:-6} workers x ${NTHRP:-8} threads"; t0=$(date +%s)
simple_exec prg=klem_board projfile=$PS infile=$T/board.txt dir=$T/board nparts=${NP:-6} nthr=${NTHRP:-8} mkdir=no > $T/board.log 2>&1 || { log "FAILED: $(grep -a -m3 'ERROR\|STOP\|rror' $T/board.log | cut -c1-200)"; ids=$(squeue -u $USER -h -o "%i %j" | awk '$2=="simple_proj"{print $1}'); [ -n "$ids" ] && scancel $ids; exit 1; }
grep -a ">>> it\|>>> states" $T/board.log | cut -c1-160; log "done in $(( $(date +%s) - t0 )) s (sequential: 1012 s at 16 threads)"
$PY $W/klem_fortran2npz.py $T/board 2>&1 | grep -v Warning
echo "== identification (distributed)"; $PY $W/deposited_match.py $T/board/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -12
echo "== identification (sequential, same board)"; $PY $W/deposited_match.py $S/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -2
$PY - <<PYEOF
import numpy as np
a = np.loadtxt('$S/klem_latent.txt', comments='#'); b = np.loadtxt('$T/board/klem_latent.txt', comments='#')
ia = {int(p): i for i, p in enumerate(a[:, 0])}; ib = {int(p): i for i, p in enumerate(b[:, 0])}; common = [p for p in ia if p in ib]
A = a[[ia[p] for p in common]]; B = b[[ib[p] for p in common]]
print('sequential vs distributed on %d particles: argmax agreement %.3f, posterior-mean distance median %.3f cells (p90 %.3f)' % (len(common), np.mean(A[:, 3] == B[:, 3]), np.median(np.hypot(A[:, 1] - B[:, 1], A[:, 2] - B[:, 2])), np.percentile(np.hypot(A[:, 1] - B[:, 1], A[:, 2] - B[:, 2]), 90)))
PYEOF
log DONE
