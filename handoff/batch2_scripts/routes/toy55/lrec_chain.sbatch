#!/bin/bash
#SBATCH -J lrec_chain
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 05:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# the leaf recipe as one job: K-state family classification (lrec_cls2.sbatch, env as there: LABEL, K, SEED, NSUB, TAG, FAMFROM,
# FAMLAB, MAXITS) then the score-coordinate loop on its output (lrec_loop.sbatch: SRC=LABEL, LABEL=loop_LABEL, KL, NIT, XHALF 1)
set -uo pipefail
T=/mnt/beegfs/elmlund/afan/routes/toy55; L=${LABEL:?}
bash $T/lrec_cls2.sbatch || { echo "[lrec_chain $L] classification FAILED"; exit 1; }
grep -q . /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/$L/final_states.txt || { echo "[lrec_chain $L] no final states"; exit 1; }
SRC=$L LABEL=loop_$L XHALF=${XHALF:-1} bash $T/lrec_loop.sbatch || { echo "[lrec_chain $L] loop FAILED"; exit 1; }
