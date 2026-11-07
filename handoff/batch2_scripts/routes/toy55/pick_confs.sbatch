#!/bin/bash
#SBATCH -J pick_confs
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
$PY $W/pick_confs.py $W/truthtree_55_200k.npz $W/sim55_700k_A/gt_particles.txt $W/sim55_700k_B/gt_particles.txt ${EXTRA:-} 2>&1 | grep -v Warning
