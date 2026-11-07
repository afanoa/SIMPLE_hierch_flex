#!/bin/bash
#SBATCH -J dk_debug
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; C=/mnt/beegfs/elmlund/afan/rnasep_klem/synth_null_200k/consolidated_12_vm
DOCK=1 FLOOR=12 MASK=region WHITEN=0 $PY $W/merge_data_derived.py $C $W/maps160/toy_consensus.mrc $C/merge_dk_12 2>&1 | grep -v Warning | tail -40
