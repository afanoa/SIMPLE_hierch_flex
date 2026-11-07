#!/bin/bash
#SBATCH -J toy10_grid2d
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# 2-D deformation latent on the supervised-seed windows: waits for the arm's deform_proto flows. env: TAG (_sup2), SIMDIR, CSV, PROJ
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=8
TAG=${TAG:-_sup2}; SIMDIR=${SIMDIR:-sim}; PJ=${PROJ:-proj}; OUT=$W/lat$TAG; T=$W/$SIMDIR; CSV=$W/${CSV:-rnasep_conformers10.csv}; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $CSV)
n=0; while [ $n -lt 300 ] && [ ! -s $OUT/deform/flow_scores.npy ]; do sleep 60; n=$((n+1)); done; [ -s $OUT/deform/flow_scores.npy ] || { echo "no flows"; exit 1; }
LOG=$(ls -t $W/toy10${TAG}_battery.*.out | head -1); AXIS=$(grep -o "truth axis (PC1 of the conformers, STATES order): .*" $LOG | head -1 | sed 's/.*: //')
R=$(grep -o "2 done: [0-9]* rounds" $LOG | tail -1 | awk '{print $3}'); TABLE=$(grep -o "poses /mnt[^ ]*" $LOG | tail -1 | awk '{print $2}'); [ -s "$TABLE" ] || TABLE=$OUT/ptcl3D_r0.txt
$PY $W/grid2d_estep.py $OUT/deform $OUT/maps_r$R $OUT/grid2d --nw ${NW:-8} --info $OUT/windows_r$R.txt --stk $T/all.mrcs --table $TABLE --gt $T/gt_particles.txt --states $STATES --axis=$AXIS --n1 7 --n2 5 --lp 12 --nsub 30000 --tag "8 grid2d$TAG" 2>&1 | grep -v "^   [0-9]* / "
echo DONE
