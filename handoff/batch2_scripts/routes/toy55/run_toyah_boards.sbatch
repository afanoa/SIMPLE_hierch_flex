#!/bin/bash
#SBATCH -J toyah_boards
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# apo/holo x open/closed on the "ah" toy, four questions on one simulated set (waits for run_toyah_sim):
#  A  plain 7x7 board (Fortran, refined poses): does any cell or axis separate holo from apo, closed from open, and
#     does the arm swing come out (state_coverage per conformer)?
#  B  oracle site: focal 1x2 board inside the truth rnpA mask (readout upper bound for occupancy classification)
#  C  truth-free site discovery: positive-excess candidates of the plain board (derive_excess.py), overlap with the
#     truth site, then the focal 1x2 board on the best candidate
#  D  product model: Python board with two tied lattice layers (base + ONE positive compact component), layer mass
#     = occupancy posterior; AUC vs the truth
# env: NT (16), STEPS (ABCD)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; S=$W/simah; L=$W/lat_ah; mkdir -p $L; export TOY_SMPD=1.464; log(){ echo "[$(date '+%F %T') ahb] $*"; }
n=0; while [ $n -lt 600 ] && [ ! -s $L/ptcl3D_r0.txt ]; do sleep 60; n=$((n+1)); done; [ -s $L/ptcl3D_r0.txt ] || { log "no toy"; exit 1; }
RD=$(ls -td $W/projah/*_refine3D | head -1); PROJ=$RD/projah.simple; CONS=$RD/recvol_state01.mrc; GT=$S/gt_particles.txt; CL=$S/closed_ids.txt; SITE=$W/maps160/holo/rnpa_mask.mrc
NALL=$(grep -vc "^#" $GT); seq 1 $NALL > $L/subset_all.txt; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $W/rnasep_conformers_ah.csv); AXIS=$(echo $STATES | sed 's/[0-9]*/0/g')
score(){ d=$1; tag=$2; $PY $W/score_ah.py $d/klem_latent.txt $GT $CL --tag "$tag" 2>&1 | grep -v Warning; $PY $W/state_coverage.py $d/klem_latent.txt $GT --argmax --tag "$tag" 2>&1 | grep -v Warning | tail -1; }
export FB=$HOME/SIMPLE_klem2/build NTHR=${NT:-16}
if [[ "${STEPS:-ABCD}" == *A* ]]; then
  log "A: plain 7x7 Fortran board at the refined poses"; [ -s $L/fA/klem_latent.txt ] || INIT=random NFREE=2 bash $W/fboard.sh $L/fA $PROJ $L/subset_all.txt "20 20 20 14 14 14 10 10 10" "1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" || exit 1
  score $L/fA "A plain board"
fi
if [[ "${STEPS:-ABCD}" == *B* ]]; then
  log "B: oracle focal 1x2 board inside the truth rnpA site"; [ -s $L/fB/klem_latent.txt ] || INIT=random NFREE=0 N1=1 N2=2 FOCUSMASK=$SITE bash $W/fboard.sh $L/fB $PROJ $L/subset_all.txt "14 14 12 12 10 10 8 8 8" "0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7" || exit 1
  score $L/fB "B oracle-site focal 1x2"
fi
if [[ "${STEPS:-ABCD}" == *C* ]]; then
  log "C: truth-free site candidates from the plain board"; $PY $W/derive_regions.py $L/fA $L/regA --maxreg 2 2>&1 | grep -v Warning | tail -3
  $PY $W/derive_excess.py $L/fA $L/exA --maxreg 4 --ref $SITE --exclude $L/regA/region_1.mrc 2>&1 | grep -v Warning
  $PY $W/derive_excess.py $L/fA $L/exA_noexcl --maxreg 4 --ref $SITE 2>&1 | grep -v Warning | tail -5
  if [ -s $L/exA/excess_1.mrc ]; then log "C: focal 1x2 board on candidate 1"; [ -s $L/fC/klem_latent.txt ] || INIT=random NFREE=0 N1=1 N2=2 FOCUSMASK=$L/exA/excess_1.mrc bash $W/fboard.sh $L/fC $PROJ $L/subset_all.txt "14 14 12 12 10 10 8 8 8" "0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7" && score $L/fC "C candidate-site focal 1x2"; fi
fi
if [[ "${STEPS:-ABCD}" == *D* ]]; then
  log "D: tied two-layer board (base + one positive compact component), Python"; export OMP_NUM_THREADS=1; mkdir -p $L/pD
  [ -s $L/pD/klem_latent.txt ] || $PY $W/klem.py $S/all.mrcs $L/ptcl3D_r0.txt $GT $L/pD --cons $CONS --states $STATES --axis=$AXIS --n1 7 --n2 7 --nfree 2 --nlayer 2 --tie_layers --tie_support 0.5 --tie_support_after 2 \
      --lps "20 20 20 14 14 14 10 10 10" --sigmas "1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" --nproc ${NT:-16} --init random --nsub 0 --devw 0.05 --pnoise --balance 0.5 --mskdiam 200 --tag "ah tied" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "tied\|final:\|states: cross-half\|layer\|Error\|Traceback" | cut -c1-220
  score $L/pD "D tied layers"
fi
log DONE
