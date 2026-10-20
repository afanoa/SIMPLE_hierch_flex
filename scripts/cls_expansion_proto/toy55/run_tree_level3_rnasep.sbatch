#!/bin/bash
#SBATCH -J rnasep_tree_level3
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Level 3 on a real family (the toy's best lever, 0.62 -> 0.66): the family's feedback child board (at its refined poses)
# is bisected by the spectral split of its state graph; a grandchild board per sub-family at the refined poses (focus,
# pose marginalisation, 12->8 A, balanced E-step), identification vs the deposited conformers per grandchild.
# env: E (expansion dir), FAM (1), CB (child board, default the family pipeline's feedback board), META (refined-pose
#      npz, default the pipeline's), BAL (0.5), L3 (_l3), PMARG (3; 0 = off)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=${E:-$A/rnasep_klem/expand1}; FAM=${FAM:-1}; P=$E/fam${FAM}_pipeline
CB=${CB:-$P/klem_fb}; META=${META:-$P/meta_fb.npz}; BAL=${BAL:-0.5}; L3=${L3:-_l3}; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//'); export OMP_NUM_THREADS=1
E3=$E/fam${FAM}${L3}; mkdir -p $E3; log(){ echo "[$(date '+%F %T') l3-$FAM] $*"; }
$PY $W/expand_families.py split $CB $E3 --nfam 2
for s in 1 2; do
  D=$E3/sub$s; mkdir -p $D; log "grandchild board: family $FAM sub-family $s ($(grep -c . $E3/family${s}_subset.txt) particles)"
  $PY $W/klem.py $STK none none $D --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $E3/family${s}_subset.txt --ctf auto \
      --lps "12 12 10 10 8 8 8 8 8" --sigmas "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" --focus 1 --pose_marg ${PMARG:-3} --devw 0.05 --pnoise ${BAL:+--balance $BAL} --tag "rnasep fam$FAM l3 sub$s" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "subset\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[sub $s] /"
  echo "[sub $s] identification:"; $PY $W/deposited_match.py $D/cells $W/maps160/docked $W/maps160/toy_consensus.mrc | sed "s/^/[sub $s] /"
done
log DONE
