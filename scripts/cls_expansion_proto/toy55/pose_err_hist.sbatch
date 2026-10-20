#!/bin/bash
#SBATCH -J synth_poseerr
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# refined-pose error shape of a synthetic set (env MODE)
W=/mnt/beegfs/elmlund/afan/routes/toy55; M=${MODE:?}; echo "== $M"
/mnt/nasapps/production/python/3.12/bin/python3 $W/pose_err_hist.py $W/sim$M/all.txt $W/lat_$M/ptcl3D_r0.txt $W/sim$M/all_perturbed.txt $W/sim$M/gt_particles.txt 2>&1 | grep -v Warning
