#!/bin/bash
#SBATCH -J klem_readout
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Sanity readout of a real-data klem run (nuisance / sub-class / state NMIs, figure) and the run-to-run match against a
# reference run on the same particles. env: KDIR (klem dir), REF (reference klem dir, default klem_gate_npz), META (default meta_all.npz)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export OMP_NUM_THREADS=4
KDIR=${KDIR:?}; REF=${REF:-$A/rnasep_klem/klem_gate_npz}; META=${META:-$A/rnasep_coupled/basis_span_all/meta_all.npz}; LAB=$A/rnasep_softem/states.txt
echo "== readout of $KDIR"; $PY $W/klem_real_readout.py $KDIR $META --labels $LAB --out $KDIR/fig_real_readout.png 2>&1 | grep -v Warning
echo "== match vs $REF"; $PY $W/klem_match.py $KDIR $REF
echo "== identification vs the deposited conformers"; $PY $W/deposited_match.py $KDIR/cells $W/maps160/docked $W/maps160/toy_consensus.mrc
