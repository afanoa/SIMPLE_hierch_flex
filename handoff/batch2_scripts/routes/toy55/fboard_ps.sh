#!/bin/bash
# Fortran klem board wrapper (prg=klem_board from the klem2 worktree build), conventions of the validated real-data v2
# run: conv 0, ss -1, fs 1, sign 1, ctfabs 0, shrink 0.02, eps0 0.05, devw 0.05, pnoise 1, lattice init from the meta
# latent (klem_init_positions.py). Outputs in OUT: board.txt, board.log, klem_latent.txt, klem_states.txt, cells/,
# klem_resp.npz + klem_halves.npz (converted), klem_poses.txt when the pose step is on.
# usage: fboard.sh OUT PROJFILE SUBSET "LPS" "SIGMAS"
# env: FB (build dir), META (npz with the latent, for the init), N1 N2 (7), NFREE (0), BAL (0.5), FOCUSMASK (file or
#      empty), INIT (latent|random), POSE_EVERY POSE_STEP POSE_PRIOR (0 3 4), NTHR (16), SEED (1), WRITE_HALVES (0),
#      NPARTS (unset = one process; N = the distributed board: N persistent SLURM workers with NTHR threads each, this
#      process only merges; needs SIMPLE_QSYS=slurm and the project's compenv)
set -uo pipefail
OUT=$1; PROJ=$2; SUBSET=$3; LPS=$4; SIGMAS=$5
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
FB=${FB:-$HOME/SIMPLE_klem2/build}; NTHR=${NTHR:-16}; mkdir -p $OUT
INITLINE=""
if [ "${INIT:-latent}" = latent ]; then
  $PY $W/klem_init_positions.py ${META:?} $SUBSET $OUT/init_positions.txt --n1 ${N1:-7} --n2 ${N2:-7} 2>&1 | grep -v Warning
  INITLINE="init_file $OUT/init_positions.txt"
fi
FOCUSLINE=""; [ -n "${FOCUSMASK:-}" ] && FOCUSLINE="focus_mask $FOCUSMASK"
POSELINE=""; [ "${POSE_EVERY:-0}" != 0 ] && POSELINE="pose_every $POSE_EVERY
pose_step ${POSE_STEP:-3}
pose_prior ${POSE_PRIOR:-4}"
cat > $OUT/board.txt <<EOF
n1 ${N1:-7}
n2 ${N2:-7}
nfree ${NFREE:-0}
lps $LPS
sigmas $SIGMAS
shrink 0.02
eps0 0.05
balance ${BAL:-0.5}
subset $SUBSET
nsub 0
seed ${SEED:-1}
conv 0
ss -1
fs 1
sign 1
ctfabs ${CTFABS:-0}
write_halves ${WRITE_HALVES:-0}
devw 0.05
pnoise 1
$FOCUSLINE
$POSELINE
$INITLINE
EOF
[ -n "${NPARTS:-}" ] && echo "nthr_part $NTHR" >> $OUT/board.txt
[ -n "${PSCALE:-}" ] && echo "pscale $PSCALE" >> $OUT/board.txt
sed -i '/^$/d' $OUT/board.txt
echo "[$(date '+%F %T') fboard] $(grep -c . $SUBSET) particles -> $OUT ($(tr '\n' ' ' < $OUT/board.txt))"
t0=$(date +%s)
( export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH OMP_NUM_THREADS=$NTHR TOY_SMPD=1.464 SIMPLE_QSYS=${SIMPLE_QSYS:-slurm} SIMPLE_QSYS_PARTITION=${SIMPLE_QSYS_PARTITION:-csbdevel}
  simple_exec prg=klem_board projfile=$PROJ infile=$OUT/board.txt dir=$OUT nthr=$NTHR nparts=${NPARTS:-1} mkdir=no > $OUT/board.log 2>&1 )
rc=$?
if [ $rc -ne 0 ] || ! [ -s $OUT/klem_latent.txt ]; then echo "[fboard] klem_board FAILED (rc $rc): $(grep -a -m3 'ERROR\|STOP\|Error' $OUT/board.log | cut -c1-200)"; tail -5 $OUT/board.log; exit 1; fi
grep -a ">>>" $OUT/board.log | tail -4 | cut -c1-200
echo "[$(date '+%F %T') fboard] done in $(( $(date +%s) - t0 )) s"
$PY $W/klem_fortran2npz.py $OUT 2>&1 | grep -v Warning
