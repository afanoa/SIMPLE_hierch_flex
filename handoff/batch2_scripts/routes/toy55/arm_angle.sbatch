#!/bin/bash
#SBATCH -J klem_arm_angle
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
A=/mnt/beegfs/elmlund/afan; /mnt/nasapps/production/python/3.12/bin/python3 $A/routes/toy55/arm_angle_readout.py $A/rnasep_klem/synth_arms $A/routes/toy55/simarms/gt_particles.txt $A/routes/toy55/truth_arms.txt 2>&1 | grep -v Warning; echo ARM ANGLE DONE
