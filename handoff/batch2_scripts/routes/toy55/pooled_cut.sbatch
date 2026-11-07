#!/bin/bash
#SBATCH -J klem_pooled_cut
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# one job: consolidation of the region-local zoom tree at CUT (env RTAG, CUT, MINP) then the pooled-state readouts
set -uo pipefail
W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan/rnasep_klem; export RTAG=${RTAG:-f692c} CUT=${CUT:?} MINP=${MINP:-1000} CTAG=_zoomreg TREE=$A/${RTAG:-f692c}/merge_tree_zoom_reg.npz LEAFGLOB='zoomregion*'
bash $W/consolidate.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^\|ValueError\|allp = "
unset LEAFGLOB; export CDIR=consolidated_${CUT}${CTAG}; bash $W/pooled_analysis.sbatch
