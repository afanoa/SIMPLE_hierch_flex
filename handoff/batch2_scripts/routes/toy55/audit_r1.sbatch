#!/bin/bash
#SBATCH -J klem_audit_r1
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# audit of node f692st r1's partitions: confounders and decisiveness (audit_r1.py) -> rnasep_klem/repro_r1/audit.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/audit_r1.py > /mnt/beegfs/elmlund/afan/rnasep_klem/repro_r1/audit.txt 2>&1
cat /mnt/beegfs/elmlund/afan/rnasep_klem/repro_r1/audit.txt
