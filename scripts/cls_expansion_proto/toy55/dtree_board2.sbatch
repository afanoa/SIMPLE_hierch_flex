#!/bin/bash
#SBATCH -J klem_dtboard
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# one node board of the data-derived tree (env R = half run dir, ND = node dir, DEPTH, FB build): 7x7 lattice + free components
# on the node's particles (ND/subset.txt), random lattice start, consensus poses, band finer with depth; local 16 threads
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; R=${R:?}; ND=${ND:?}; DEPTH=${DEPTH:-0}
FB=${FB:-$HOME/SIMPLE_klem2/build}; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
case $DEPTH in 0) LPS="20 20 20 14 14 14 10 10 10"; SIG="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7"; NF=2;; 1) LPS="14 14 12 12 10 10 8 8 8"; SIG="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"; NF=1;; *) LPS="12 12 10 10 8 8 8 8 8"; SIG="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"; NF=1;; esac
SEED=${SEED:-1}; O=$ND/board_s$SEED; mkdir -p $O   # SEED: random lattice start (EM restarts; the node keeps the most confident board)
$PY -c "import numpy as np; p = np.loadtxt('$ND/subset.txt', dtype=int, ndmin=1); g = np.random.default_rng($SEED).uniform(0, 6, (len(p), 2)); np.savetxt('$O/init_positions.txt', np.c_[p, g], fmt='%d %.4f %.4f')"
cat > $O/board.txt <<CFG
n1 7
n2 7
nfree $NF
lps $LPS
sigmas $SIG
shrink 0.02
eps0 0.05
balance 0.5
subset $ND/subset.txt
nsub 0
seed $SEED
conv 0
ss -1
fs 1
sign 1
ctfabs 0
write_halves 0
devw 0.05
pnoise 1
init_file $O/init_positions.txt
CFG
cd $O; t0=$(date +%s)
simple_exec prg=klem_board projfile=$R/src.simple infile=$O/board.txt dir=$O nthr=16 nparts=1 mkdir=no > $O/board.out 2>&1 || { echo "board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.out | cut -c1-200)"; exit 1; }
echo "[$(date '+%F %T') dtboard] $ND depth $DEPTH seed $SEED: $(grep -c . $ND/subset.txt) particles, board done in $(( $(date +%s) - t0 )) s"; grep -a ">>> it 9 (" $O/board.out | cut -c1-200
