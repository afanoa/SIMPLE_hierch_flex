#!/bin/bash
#SBATCH -J klem_landgraph
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# sample a run's state landscape evenly along its main path (landscape_graph.py), then reconstruct every sampled state from all
# its particles (mu_rec.sbatch) and read out (pooled_analysis.sbatch). env RTAG, SRCC (consolidated_12_fullw_zoom), TARGET (8),
# NMIN (data-derived if unset), MODE (fold|fps), NPC (3), OUT (landscape_<MODE>_n<NMIN>), PERI (prefix, default the run's mergeup/periphery_fam{F})
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; R=$K/${RTAG:?}; export OMP_NUM_THREADS=8
OUT=${OUT:-graph_s${NSEG:-10}_b${BUDGET:-5000}}; O=$R/$OUT; mkdir -p $O
$PY $W/landscape_graph.py $R/${SRCC:-consolidated_12_fullw_zoom} ${PARENTS:-$R/parents}/fam{F} ${PERI:-$R/mergeup/periphery_fam{F}} $O --nseg ${NSEG:-10} --budget ${BUDGET:-5000} 2>&1 | grep -v Warning | tee $O/sample.log
RTAG=$RTAG CDIR=$OUT bash $W/mu_rec.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^" | tail -3
RTAG=$RTAG CDIR=$OUT bash $W/pooled_analysis.sbatch 2>&1
bash $W/tidy_halves.sh $O; echo "LANDSCAPE GRAPH DONE $OUT"
