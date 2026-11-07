#!/bin/bash
#SBATCH -J toy10_oracle
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Oracle E-steps against the TRUTH conformer maps (rerun with the per-shell reference amplitude scale). env: SIMDIR, TAG, AXIS, STATES
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=8
SIMDIR=${SIMDIR:-sim}; TAG=${TAG:-}; T=$W/$SIMDIR; OUT=$W/lat$TAG; GT=$T/gt_particles.txt; STK=$T/all.mrcs
STATES=${STATES:-1,10,18,27,28,37,46,48,49,51}; AXIS=${AXIS:--0.246,-0.916,-0.442,1.000,0.972,-0.685,-0.489,0.446,0.499,-0.139}
TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc $W/maps160/variants/conf$(printf %02d $c)_*.mrc 2>/dev/null | head -1; done | tr "\n" , | sed 's/,$//')
for LP in 20 14 10 8; do
  $PY $W/latent_windows.py estep $STK $OUT/ptcl3D_r0.txt $GT --refs $TRUTH --lp $LP --out $OUT/oracle_lp$LP.txt --states $STATES --axis=$AXIS --nsub 40000 --tag "oracle lp$LP" 2>&1 | grep -v "^   [0-9]* / " | sed "s/^/[toy10$TAG] /"
done
echo DONE
