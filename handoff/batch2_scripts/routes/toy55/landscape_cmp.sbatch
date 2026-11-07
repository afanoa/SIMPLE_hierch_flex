#!/bin/bash
#SBATCH -J klem_landscape_cmp
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape-level reproducibility of node f692st r1 across independent runs (landscape_cmp.py): script / port1 / port2 (/ port3) leaf maps,
# plus node r2's script leaves as a different-landscape control -> rnasep_klem/repro_r1/landscape_cmp.txt
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; K=/mnt/beegfs/elmlund/afan/rnasep_klem; P=$K/stree_port
A="script=$K/rflex_f692st/node_r1/fstageB/2_refine3D:$K/repro_r1/labels/S_B.txt port1=$P/tree_f692r1/1_hetero_tree/node_r/stageB:$K/repro_r1/labels/P_B.txt"
A="$A port2=$P/tree_f692r1_seed2/1_hetero_tree/node_r/stageB:$K/repro_r1/dscore_port2/labels.txt"
[ -s $K/repro_r1/dscore_port3/labels.txt ] && A="$A port3=$P/tree_f692r1_rep3/1_hetero_tree/node_r/stageB:$K/repro_r1/dscore_port3/labels.txt"
A="$A r2ctrl=$K/rflex_f692st/node_r2/fstageB/2_refine3D:$K/rflex_f692st/node_r2/refined_f.txt"
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/landscape_cmp.py $A --lp 15 > $K/repro_r1/landscape_cmp.txt 2>&1; cat $K/repro_r1/landscape_cmp.txt
