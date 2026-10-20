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
FB=${FB:-$HOME/SIMPLE_klem2/build}; NTHR=${NTHR:-16}; NTHRM=${NTHR_MASTER:-$NTHR}; mkdir -p $OUT
INITLINE=""
if [ "${INIT:-latent}" = file ]; then INITLINE="init_file ${INITFILE:?}"   # INIT=file: given "pind g1 g2" start positions
elif [ "${INIT:-latent}" = latent ]; then
  $PY $W/klem_init_positions.py ${META:?} $SUBSET $OUT/init_positions.txt --n1 ${N1:-7} --n2 ${N2:-7} 2>&1 | grep -v Warning
  INITLINE="init_file $OUT/init_positions.txt"
fi
FOCUSLINE=""; [ -n "${FOCUSMASK:-}" ] && FOCUSLINE="focus_mask $FOCUSMASK"
SOLVLINE=""; [ -n "${SOLVENTMASK:-}" ] && SOLVLINE="solvent_mask $SOLVENTMASK"
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
$SOLVLINE
$POSELINE
$INITLINE
EOF
SERVED=0; [ -n "${SERVER_LIST:-}" ] && [ "${SERVER_START:-0}" != 1 ] && SERVED=1   # served: queue on a running board server, no process of our own
[ -n "${NPARTS:-}" ] && [ $SERVED = 0 ] && echo "nthr_part $NTHR" >> $OUT/board.txt
[ "${SERVER_START:-0}" = 1 ] && { echo "server_list $SERVER_LIST" >> $OUT/board.txt; echo "kmax_load ${KMAX_LOAD:-29}" >> $OUT/board.txt; }
[ -n "${PSCALE:-}" ] && echo "pscale $PSCALE" >> $OUT/board.txt
[ -n "${GRIDDING:-}" ] && echo "gridding $GRIDDING" >> $OUT/board.txt
sed -i '/^$/d' $OUT/board.txt
echo "[$(date '+%F %T') fboard] $(grep -c . $SUBSET) particles -> $OUT ($(tr '\n' ' ' < $OUT/board.txt))"
t0=$(date +%s); rm -f $OUT/BOARD_DONE
waitdone(){ n=0; while [ ! -f $OUT/BOARD_DONE ]; do [ -f $SERVER_LIST.dead ] && { echo "[fboard] board server died"; return 1; }; sleep 10; n=$((n+1)); [ $n -gt ${SERVER_TIMEOUT:-2160} ] && { echo "[fboard] timed out waiting for the board server"; return 1; }; done; return 0; }
if [ $SERVED = 1 ]; then
  echo "$OUT" >> $SERVER_LIST; echo "[fboard] queued on the board server $SERVER_LIST"; waitdone; rc=$?
elif [ "${SERVER_START:-0}" = 1 ]; then
  : > $SERVER_LIST; rm -f $SERVER_LIST.dead $SERVER_LIST.pid
  ( cd $OUT; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH OMP_NUM_THREADS=$NTHR TOY_SMPD=1.464 SIMPLE_QSYS=${SIMPLE_QSYS:-slurm} SIMPLE_QSYS_PARTITION=${SIMPLE_QSYS_PARTITION:-csbdevel}
    simple_exec prg=klem_board projfile=$PROJ infile=$OUT/board.txt dir=$OUT nthr=$NTHRM nparts=${NPARTS:-1} mkdir=no > $OUT/board.log 2>&1 || touch $SERVER_LIST.dead ) &
  echo $! > $SERVER_LIST.pid; echo "[fboard] board server started (pid $!, list $SERVER_LIST, ${NPARTS:-1} workers)"; waitdone; rc=$?; [ $rc -ne 0 ] && echo STOP >> $SERVER_LIST
else
( cd $OUT; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH OMP_NUM_THREADS=$NTHR TOY_SMPD=1.464 SIMPLE_QSYS=${SIMPLE_QSYS:-slurm} SIMPLE_QSYS_PARTITION=${SIMPLE_QSYS_PARTITION:-csbdevel}
  simple_exec prg=klem_board projfile=$PROJ infile=$OUT/board.txt dir=$OUT nthr=$NTHRM nparts=${NPARTS:-1} mkdir=no > $OUT/board.log 2>&1 )
rc=$?
fi
if [ $rc -ne 0 ] || ! [ -s $OUT/klem_latent.txt ]; then echo "[fboard] klem_board FAILED (rc $rc): $(grep -a -m3 'ERROR\|STOP\|Error' $OUT/board.log | cut -c1-200)"; tail -5 $OUT/board.log; exit 1; fi
grep -a ">>>" $OUT/board.log | tail -4 | cut -c1-200
echo "[$(date '+%F %T') fboard] done in $(( $(date +%s) - t0 )) s"
$PY $W/klem_fortran2npz.py $OUT 2>&1 | grep -v Warning
