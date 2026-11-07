#!/bin/bash
#SBATCH -J toy10_klem_board_distr
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The DISTRIBUTED Fortran board on the 20k toy (prg=klem_board nparts=NP: persistent workers via SLURM), same config
# as the sequential toy test, compared with the sequential board's outputs (lat/klem_fboard20k): argmax agreement,
# posterior-mean distance, truth purity. env: NP (2), NTHRP (8), BDIR (build)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${BDIR:-$HOME/SIMPLE_klem2/build}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; S=$W/lat/klem_fboard20k; T=$W/lat/klem_fboard20k_distr${NP:-2}; mkdir -p $T; log(){ echo "[$(date '+%F %T') fbd] $*"; }
cp -f $S/subset.txt $T/subset.txt; sed -e "s#$S/subset.txt#$T/subset.txt#" $S/board.txt > $T/board.txt; echo "nthr_part ${NTHRP:-8}" >> $T/board.txt; cat $T/board.txt | tr '\n' ' '; echo
PJ=$(ls -td $W/proj/*_refine3D | head -1)/proj.simple
log "distributed klem_board, nparts ${NP:-2} x ${NTHRP:-8} threads, on $(grep -c . $T/subset.txt) particles"; t0=$(date +%s)
simple_exec prg=klem_board projfile=$PJ infile=$T/board.txt dir=$T/board nparts=${NP:-2} nthr=${NTHRP:-8} mkdir=no > $T/board.log 2>&1 || { log "klem_board FAILED: $(grep -a -m3 'ERROR\|STOP\|Error' $T/board.log | cut -c1-200)"; tail -8 $T/board.log; tail -5 $T/board/distr/klem_workers.log 2>/dev/null; exit 1; }
grep -a ">>>" $T/board.log | cut -c1-200; log "distributed board done in $(( $(date +%s) - t0 )) s (sequential reference: see $W/toy10_klem_board_loader.*.out)"
$PY $W/klem_fortran2npz.py $T/board 2>&1 | grep -v Warning
$PY - <<PYEOF
import numpy as np
a = np.loadtxt('$S/board/klem_latent.txt', comments='#'); b = np.loadtxt('$T/board/klem_latent.txt', comments='#')
ia = {int(p): i for i, p in enumerate(a[:, 0])}; ib = {int(p): i for i, p in enumerate(b[:, 0])}; common = [p for p in ia if p in ib]
A = a[[ia[p] for p in common]]; B = b[[ib[p] for p in common]]
print('sequential vs distributed on %d particles: argmax agreement %.3f, posterior-mean distance median %.3f cells (p90 %.3f), grid-mass corr %.3f' % (len(common), np.mean(A[:, 3] == B[:, 3]), np.median(np.hypot(A[:, 1] - B[:, 1], A[:, 2] - B[:, 2])), np.percentile(np.hypot(A[:, 1] - B[:, 1], A[:, 2] - B[:, 2]), 90), np.corrcoef(A[:, 5], B[:, 5])[0, 1]))
PYEOF
echo "== truth readout (distributed)"; $PY $W/state_coverage.py $T/board/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "distributed board 20k" 2>&1 | grep -v Warning | tail -1
echo "== truth readout (sequential)";  $PY $W/state_coverage.py $S/board/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "sequential board 20k" 2>&1 | grep -v Warning | tail -1
log DONE
