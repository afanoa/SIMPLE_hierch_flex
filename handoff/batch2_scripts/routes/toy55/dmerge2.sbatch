#!/bin/bash
#SBATCH -J klem_dmerge2
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# confusion-driven merge of node f692st r1's leaves per run, then cross-run agreement at equal K (decisive_merge2.py) -> repro_r1/decisive_merge2.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/decisive_merge2.py ${RUNS:-script port1 port2} --thr 3 > /mnt/beegfs/elmlund/afan/rnasep_klem/repro_r1/decisive_merge2.txt 2>&1
cat /mnt/beegfs/elmlund/afan/rnasep_klem/repro_r1/decisive_merge2.txt
