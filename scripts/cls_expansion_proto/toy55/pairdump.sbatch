#!/bin/bash
#SBATCH -J klem_pairdump
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# dump + render the extreme pairs of one node's shuffled (LG regionshuf1) or real (regioncand) groups (env MODE, TAG, ND, LG, STAT)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:-dtree1}; ND=${ND:-r}; LG=${LG:-regionshuf1}
NA=$K/synth_${MODE}_A/$TAG/$ND; NB=$K/synth_${MODE}_B/$TAG/$ND; O=$K/pairdump_${MODE}_${ND}_$LG
env FAMS=0 MASK=sphere NPROC=16 MINREC=100 FSCMAX=30 NUIS=2 STAT=${STAT:-z} LEAFGLOB=$LG DUMPPAIR=$O $PY $W/merge_cross_tree.py $NA $NB $W/maps160/toy_consensus.mrc $O.merge 2>&1 | grep -a "dumped\|pair statistics" | cut -c1-200
for t in max min; do $PY $W/pairdump_plot.py ${O}_$t.npz ${O}_$t.png 2>&1 | grep -v Warning; done
