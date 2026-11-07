#!/bin/bash
#SBATCH -J lrec_trueproj
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# the benchmark's particles with their TRUE poses (simulation table, shifts 0) as a separate project, for the pose-error ceiling
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
T=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:-swing30}; S=$T/sim$TAG; PJ=proj${TAG}true
[ -e $T/$PJ ] && { echo "$T/$PJ exists"; exit 1; }
cd $T && simple_exec prg=new_project projname=$PJ > $T/np_$PJ.log 2>&1; cd $T/$PJ || { echo "new_project FAILED"; exit 1; }
sed -E 's/ state=[^ ]*//' $S/all.txt | sed -E 's/$/ state=1/' > $S/all_true_import.txt
simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_true_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { echo "import FAILED"; tail -3 import.log; exit 1; }
simple_exec prg=print_project_field projfile=$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null | head -2 | cut -c1-300
echo "DONE $T/$PJ/$PJ.simple"
