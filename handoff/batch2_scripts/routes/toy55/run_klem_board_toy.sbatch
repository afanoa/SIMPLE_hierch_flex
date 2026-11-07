#!/bin/bash
#SBATCH -J toy10_klem_board_test
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The Fortran klem board (prg=klem_board) on the toy: the same 20k particles, schedule and conventions as the Python
# reference dump board (lps 12 12 10, random init, no balance), outputs converted to npz, read out against the truth
# (state_coverage.py) and the deposited conformers (deposited_match.py); wall time per round printed by the program.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${BDIR:-$HOME/SIMPLE_klem2/build}; DONE=${DONEFILE:-make4.done}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; D=$W/lat/klem_dump20k; T=$W/lat/klem_fboard20k; mkdir -p $T; log(){ echo "[$(date '+%F %T') fb] $*"; }
n=0; while [ $n -lt 90 ] && [ ! -f $B/$DONE ]; do sleep 60; n=$((n+1)); done; cat $B/$DONE; grep -q "_EXIT=0" $B/$DONE || { log "build failed"; exit 1; }
$PY - <<PYEOF
import numpy as np
Z = np.load('$D/klem_estep_ref.npz'); np.savetxt('$T/subset.txt', Z['pinds'].astype(int), fmt='%d')
open('$T/board.txt', 'w').write('n1 7\nn2 7\nnfree 2\nlps 12 12 10\nsigmas 1.0 0.8 0.7\nshrink 0.02\neps0 0.05\nbalance 0\nfocus_mask none\nsubset $T/subset.txt\nnsub 0\nseed 1\nss %g\nfs %g\nconv %d\nsign %g\nctfabs 0\nwrite_halves 0\n' % (float(Z['ss']), float(Z['fs']), int(Z['conv']), float(Z['sign'])))
print(open('$T/board.txt').read().replace('\n', ' | '))
PYEOF
PJ=$(ls -td $W/proj/*_refine3D | head -1)/proj.simple
log "simple_exec prg=klem_board on $(grep -c . $T/subset.txt) particles"; t0=$(date +%s)
simple_exec prg=klem_board projfile=$PJ infile=$T/board.txt dir=$T/board nthr=16 nparts=1 mkdir=no > $T/board.log 2>&1 || { log "klem_board FAILED: $(grep -a -m3 'ERROR\|STOP\|Error' $T/board.log | cut -c1-200)"; tail -8 $T/board.log; exit 1; }
grep -a ">>>" $T/board.log | cut -c1-200; log "klem_board done in $(( $(date +%s) - t0 )) s"
$PY $W/klem_fortran2npz.py $T/board 2>&1 | grep -v Warning
echo "== truth readout (Fortran board)"; $PY $W/state_coverage.py $T/board/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "fortran board 20k" 2>&1 | grep -v Warning | tail -3
echo "== truth readout (Python reference board, same particles/schedule)"; $PY $W/state_coverage.py $D/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "python board 20k" 2>&1 | grep -v Warning | tail -3
echo "== identification against the deposited conformers (Fortran cells)"; $PY $W/deposited_match.py $T/board/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -13
log DONE
