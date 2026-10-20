#!/bin/bash
#SBATCH -J klem_embed_now
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# embeddings available NOW: the corrected real run's boards so far (parent, children, feedback, finished regions) and the
# finished no-refinement toy chain with its truth
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
echo "######## f692c (partial tree, corrected conventions)"; $PY $W/global_embedding.py $A/rnasep_klem/f692c $A/rnasep_coupled/basis_span_all/meta_all.npz $A/rnasep_klem/f692c/fig_global_embedding_partial.png --nsub 80000 2>&1 | grep -v Warning | tail -14
echo "######## toyah_tree_nr (finished, truth)"; $PY $W/global_embedding.py $A/rnasep_klem/toyah_tree_nr $A/rnasep_klem/toyah_tree/meta_toy.npz $A/rnasep_klem/toyah_tree_nr/fig_global_embedding_final.png --nsub 80000 --gt $W/simah/gt_particles.txt --holo0 101 2>&1 | grep -v Warning | tail -24
echo "ALL DONE"
