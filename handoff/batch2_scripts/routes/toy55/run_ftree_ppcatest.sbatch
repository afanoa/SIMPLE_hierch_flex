#!/bin/bash
#SBATCH -J toyah_ftree_ppca
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Fortran tree, stage=region on the toy's family 1 (region_1 mask) with level-2 ranking by the in-house cell PPCA
# (rank ppca, BIC rank <= 4) and reconstructions on; the log prints the PPCA/Gram-axis agreement and the identification.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
FB=$HOME/SIMPLE_klem2/build; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
A=/mnt/beegfs/elmlund/afan; R=$A/rnasep_klem/toyah_ftree; T=$A/rnasep_klem/toyah_tree; cd $R; log(){ echo "[$(date '+%F %T') ppcatest] $*"; }
sed -e 's/^region_tag .*/region_tag ppcatest/' -e 's/^norec .*/norec 0/' -e 's/^nl2 .*/nl2 1/' $R/fam1/tree_region_1.txt > $R/fam1/tree_ppcatest.txt; echo "rank ppca" >> $R/fam1/tree_ppcatest.txt; echo "qmax 4" >> $R/fam1/tree_ppcatest.txt
grep -q "^stage region" $R/fam1/tree_ppcatest.txt || { log "unexpected cfg"; cat $R/fam1/tree_ppcatest.txt; exit 1; }
log "klem_tree stage region (rank ppca) on family 1 / region_1"; t0=$(date +%s)
simple_exec prg=klem_tree projfile=$T/src.simple infile=$R/fam1/tree_ppcatest.txt nthr=8 mkdir=no > $R/fam1/tree_ppcatest.log 2>&1 || { log "klem_tree FAILED: $(grep -a -m3 'ERROR\|STOP\|rror' $R/fam1/tree_ppcatest.log | cut -c1-200)"; }
log "done in $(( ($(date +%s) - t0) / 60 )) min"; grep -a "PPCA\|cell score\|identified\|states\b" $R/fam1/tree_ppcatest.log | tail -12 | cut -c1-200; log DONE
