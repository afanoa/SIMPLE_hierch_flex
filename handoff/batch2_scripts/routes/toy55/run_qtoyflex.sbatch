#!/bin/bash
#SBATCH -J qtoyflex
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# positive control for flex <-> refinement: qtoy2 with its six arm stops replaced by a 30-level s1 continuum (qtoyflex_maps.py), the
# closed and holo states discrete as in qtoy2; then qtoy's own end-to-end simulation (the run_qtoy_sim.sbatch body, QT=qtoyflex) in this
# job -> simqtoyflex/src.simple and simqtoyflex/truth_qtoy.txt (conformer 201-230 = arm, 18 closed, 101 holo, 0 junk).
set -uo pipefail
W=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3; S=$W/simqtoyflex
[ -s $W/qtoyflex.csv ] || $PY $W/qtoyflex_maps.py $W/qtoy2.csv $W/maps160 armflex $W/qtoyflex.csv 30 || exit 1
mkdir -p $S
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py /mnt/beegfs/elmlund/afan/rnasep_coupled/basis_span_all/meta_all.npz $W/qtoyflex.csv $S ${N:-60000} 0.05 7 || exit 1
# simulate_particles refuses an odd number of particles in a class (build_refspiral): an odd class drops its last row
for f in $S/oris_in_*.txt; do [ $(( $(wc -l < $f) % 2 )) = 0 ] || sed -i '$ d' $f; done
QT=qtoyflex bash $W/run_qtoy_sim.sbatch
