#!/bin/bash
#SBATCH -J klem_score_embed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# score-space embedding of a finished rflex tree (score_embed.py): env RT (rflex_<RT>), POSES_DIR (dir with poses_A.txt / poses_B.txt =
# each particle's own leaf-node pose per half row; default rflex_<RT>/merge13)
export OMP_NUM_THREADS=16 NUMBA_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/rflex_${RT:?}
PD=${POSES_DIR:-$R/merge13}; H=$T/simrnasep692
$PY -u $T/score_embed.py --sim ${H}_A,${H}_B --cache $T/benchcache_rnasep692 --poses $PD/poses_A.txt,$PD/poses_B.txt --leaves $R/leaves_list.txt \
    --labels $R/final_labels.txt --out $K_/score_embed_$RT --em_lim ${EMLIM:-15} 2>&1 | grep -v "Warning\|warn("
