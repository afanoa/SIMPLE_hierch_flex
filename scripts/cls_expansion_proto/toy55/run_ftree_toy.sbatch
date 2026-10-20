#!/bin/bash
#SBATCH -J toyah_ftree
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The Fortran tree driver (prg=klem_tree) on the apo/holo toy, same settings as the no-refinement bash chain
# (toyah_tree_nr): compare families, states, identification and the landscape count. env: RTAG (toyah_ftree), SKIP (1), NOREC (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
FB=$HOME/SIMPLE_klem2/build; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-toyah_ftree}; T=$A/rnasep_klem/toyah_tree; mkdir -p $R; cd $R; log(){ echo "[$(date '+%F %T') ftree] $*"; }
n=0; while [ $n -lt 60 ] && [ ! -f $FB/make19.done ] && [ ! -f $FB/make20.done ]; do sleep 60; n=$((n+1)); done
ls $W/maps160/docked/conf??_*_m160.mrc | awk '{ match($0, /conf([0-9]+)_/, m); print m[1], $0 }' > $R/docked_list.txt; head -3 $R/docked_list.txt
cat > $R/tree.txt <<CFG
run $R
stage all
nfam 2
nparts 4
nthr_part 8
np_refine 10
maxits 6
skip_refine ${SKIP:-1}
norec ${NOREC:-1}
maxreg 3
nl2 1
parmax 2
docked_list $R/docked_list.txt
cons $W/maps160/toy_consensus.mrc
meta_init none
CFG
log "klem_tree on $T/src.simple"; t0=$(date +%s)
simple_exec prg=klem_tree projfile=$T/src.simple infile=$R/tree.txt nthr=8 mkdir=no > $R/tree.log 2>&1 || { log "klem_tree FAILED: $(grep -a -m3 'ERROR\|STOP\|rror' $R/tree.log | cut -c1-200)"; tail -8 $R/tree.log; exit 1; }
grep -a ">>>" $R/tree.log | cut -c1-200; log "done in $(( $(date +%s) - t0 )) s"
echo "== reference: the bash chain without refinement (toyah_tree_nr)"; grep -a "similarity 0\.\|identified [0-9]* /" $(ls -t $A/rnasep_klem/toyah_tree_nr.*.out | head -1) | tail -8 | cut -c1-200
