#!/bin/bash
#SBATCH -J toyah_readouts
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# (a) truth purity of every leaf state of the finished no-refinement toy chain (18 conformers incl. holo, junk);
# (b) pose error vs truth of the chain's baseline refinement (local 20->8, maxits 6) that the matrix job did not print.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; GT=$W/simah/gt_particles.txt; R=$A/rnasep_klem/toyah_tree_nr
echo "######## (a) leaf states of toyah_tree_nr vs truth"
for f in $(ls $R/fam*/region*/level1/states_clusters.txt $R/fam*/region*/level2_s*/states_top.txt 2>/dev/null); do tag=$(echo $f | sed "s#$R/##; s#/states_[a-z]*.txt##"); $PY $W/score_states.py $f $GT x "$tag" 2>&1 | grep -v Warning | head -14; done
echo "######## (b) baseline refinement pose error (toyah_tree fam1, local 20->8 maxits 6)"
f=$A/rnasep_klem/toyah_tree/fam1/refined_ptcl3D.txt; [ -s $f ] && $PY $W/pose_error.py $f $W/simah/all.txt --states $A/rnasep_klem/toyah_tree/fam1/states.txt --tag baseline 2>&1 | tail -3
echo "ALL DONE"
