#!/bin/bash
#SBATCH -J igm2_bagboard
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%A_%a.out
# one focused Python klem board with the 2D sub-class bag prior per array task: line t of CFG = OUT MASK META CONS BETA BAGS
# PARENTS (5x5 lattice + 1 free, random start, sub-class bags centred on the parent class, beta 0 = per-particle control).
# env CFG, MODE
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=1
read OUT MASK META CONS BETA BAGS PARS <<< "$(sed -n ${SLURM_ARRAY_TASK_ID:?}p ${CFG:?})"
[ -s $OUT/klem_latent.txt ] || $PY $W/klem.py $W/sim${MODE:?}/all.mrcs none none $OUT --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --init random --meta $META \
  --bags $BAGS --bag_parents $PARS --bag_beta $BETA --n1 5 --n2 5 --nfree 1 --focus_mask $MASK --lps "20 14 12 10 10 8 8 8" --sigmas "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" \
  --nproc 16 --mskdiam 220 --ctf auto --balance 0.5 --pnoise --pscale --devw 0.05 > $OUT.log 2>&1 || { echo "board FAILED $OUT: $(tail -2 $OUT.log)"; exit 1; }
echo "task done $OUT"
