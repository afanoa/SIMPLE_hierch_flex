#!/bin/bash
#SBATCH -J synth_propagate
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# neighbouring-bin propagation of the direction-bin sub-classes into 3D states (env MODE), truth readout of the result
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; D=$W/binsroute_$MODE; O=$D/propagate; GT=$W/sim$MODE/gt_particles.txt; NP=$(grep -vc '^#' $GT)
echo "== $MODE: neighbour-bin propagation"; $PY $W/neighbor_propagate.py $D/stitch $D/p $D/bins.npz $D/p/cavgs_bins.mrc $O --nn ${NN:-6} --maxang ${MAXANG:-30} --npart $NP 2>&1 | grep -v Warning
echo "== local check: do matched sub-class pairs agree in truth?"; $PY $W/propagate_check.py $O $D/p $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning
echo "== truth readout of the propagated states"; $PY $W/label_truth_readout.py $O/particle_states.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -9; echo "PROPAGATE DONE $MODE"
