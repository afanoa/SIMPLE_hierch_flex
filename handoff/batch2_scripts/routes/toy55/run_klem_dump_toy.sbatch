#!/bin/bash
#SBATCH -J toy10_klem_dump
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Reference run for the Fortran port: a short klem board on 20k toy particles at the refined poses, no pose
# marginalisation, no focus; KLEM_DUMP=1 writes the last E-step's band/noise/scale/log-likelihoods and the half cell
# volumes (klem_estep_ref.npz, cells_half/). The Fortran E-step must reproduce ll on these particles.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; export OMP_NUM_THREADS=1 KLEM_DUMP=1; export TOY_SMPD=1.464
RP=$(ls -td $W/proj/*_refine3D | head -1); CONS=$RP/recvol_state01.mrc; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $W/rnasep_conformers10.csv); AXIS=$(echo $STATES | sed 's/[0-9][0-9]*/0/g')
D=$W/lat/klem_dump20k; mkdir -p $D
$PY $W/klem.py $W/sim/all.mrcs $W/lat/ptcl3D_r0.txt $W/sim/gt_particles.txt $D --cons $CONS --states $STATES --axis=$AXIS --n1 7 --n2 7 --n3 1 --nfree 2 --lps "12 12 10" --sigmas "1.0 0.8 0.7" --nproc 16 --init random --nsub 20000 --mskdiam 232 --tag "dump20k" 2>&1 | grep -v "^   [0-9]* / " | tail -12
ls -la $D/klem_estep_ref.npz $D/cells_half | head -5
