#!/bin/bash
#SBATCH -J klem_plot_tree_st
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# construction-tree figure of rflex_<RT> (plot_tree_st.py) with the delivered maps' arm readouts -> rnasep_klem/viz_<RT>/<RT>_tree.png. env RT
K_=/mnt/beegfs/elmlund/afan/rnasep_klem; mkdir -p $K_/viz_${RT:?}
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/plot_tree_st.py $K_/rflex_$RT $K_/tw_units/${RT}_plain/arm_emp2.txt $K_/viz_$RT/${RT}_tree.png 2>&1 | grep -v Warning
