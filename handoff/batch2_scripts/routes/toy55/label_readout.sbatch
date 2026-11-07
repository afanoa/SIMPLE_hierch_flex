#!/bin/bash
#SBATCH -J klem_label_readout
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
W=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3
for M in arms tethered; do echo "== 2D route $M: kernel sub-classes of the 5 expanded classes vs truth"; $PY $W/label_truth_readout.py $W/route2d_$M/subclass_labels.txt $W/sim$M/gt_particles.txt $W/truth_$M.txt $M 2>&1 | grep -v Warning | tail -14; done; echo LABEL DONE
