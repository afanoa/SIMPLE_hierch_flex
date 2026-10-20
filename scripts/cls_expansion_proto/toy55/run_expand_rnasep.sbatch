#!/bin/bash
#SBATCH -J rnasep_expand
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Genetic expansion on real RNase P (the recipe that worked on the toy as v4): parent run -> families by spectral
# bisection of the state graph -> a child klem per family on its particles (no free comps, focus on, 14 -> 8 A, pose
# marginalisation + deviation weights + per-particle noise) -> per family: survival, identification vs the deposited
# conformers, 5 bins at full box. env: PARENT (klem_full_npz), NFAM (2), ETAG, META (pose npz for the children, default meta_all), BAL (balanced E-step lambda), PMARG (children pose marginalisation, 3; 0 = off)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; PARENT=${PARENT:-$A/rnasep_klem/klem_full_npz}; NFAM=${NFAM:-2}; ETAG=${ETAG:-}; E=$A/rnasep_klem/expand$ETAG; mkdir -p $E
META=${META:-$A/rnasep_coupled/basis_span_all/meta_all.npz}; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//'); export OMP_NUM_THREADS=1
$PY $W/expand_families.py split $PARENT $E --nfam $NFAM
for j in $(seq 1 $NFAM); do
  D=$E/fam$j; mkdir -p $D; echo "[$(date '+%F %T')] child run for family $j"
  $PY $W/klem.py $STK none none $D --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $E/family${j}_subset.txt --ctf auto \
      --lps "${LPS:-14 14 12 12 10 10 8 8 8}" --sigmas "${SIGMAS:-1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7}" --focus 1 --pose_marg ${PMARG:-3} --devw 0.05 --pnoise ${BAL:+--balance $BAL} --tag "rnasep expand$ETAG fam$j" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "subset\|focus mask\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j] /"
  echo "[fam $j] identification:"; $PY $W/deposited_match.py $D/cells $W/maps160/docked $W/maps160/toy_consensus.mrc | sed "s/^/[fam $j] /"
done
echo "[$(date '+%F %T')] expand DONE (bins per family: env KDIR=$E/famN NB=5 MSKDIAM=220 sbatch run_klem_cells_rec.sbatch)"
