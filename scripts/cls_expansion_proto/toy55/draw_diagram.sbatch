#!/bin/bash
#SBATCH -J klem_diagram
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 00:05:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
mkdir -p /mnt/beegfs/elmlund/afan/routes/toy55/figs; /mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/draw_method_diagram.py
