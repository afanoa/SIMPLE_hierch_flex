#!/bin/bash
#SBATCH -J klem_rflex13s_igm
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The state tree (rflex13s = rflex13 with synthetic half-stack inputs) on the IgM-like simarmsigm1 (5 independent arms x 2 angles), in ONE
# job: private copy of the refined master project (projarmsigm1, refine3D prob from 10-deg perturbed poses) -> meta_rows of the existing
# alternate-row halves simarmsigm1_A / _B (split_halves.py: A = master rows 1, 3, 5, ...) -> rflex13s inline (RT armsigm1s).
# env NP (8), PAR (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K_=$A/rnasep_klem; SDIR=$K_/armsigm1_src; mkdir -p $SDIR; cd $SDIR
NR=$(grep -c . $W/simarmsigm1/all.txt); NA=$(grep -c . $W/simarmsigm1_A/all.txt); NB=$(grep -c . $W/simarmsigm1_B/all.txt)
[ $NA = $(( (NR + 1) / 2 )) ] && [ $NB = $(( NR / 2 )) ] || { echo "halves $NA / $NB do not split $NR rows by parity"; exit 1; }
[ -s $SDIR/meta_rows_B.txt ] || { seq 1 2 $NR > $SDIR/meta_rows_A.txt; seq 2 2 $NR > $SDIR/meta_rows_B.txt; }
[ -s $SDIR/src.simple ] || { cp -f $W/projarmsigm1/1_refine3D/projarmsigm1.simple $SDIR/src.simple && simple_exec prg=update_project projfile=$SDIR/src.simple qsys_name=slurm qsys_partition=csbdevel \
  job_memory_per_task=24000 walltime=43200 nparts=8 nthr=8 mkdir=no > $SDIR/up.log 2>&1 || { echo "source project FAILED"; exit 1; }; }
echo "armsigm1 source $SDIR/src.simple: $NR rows, halves $NA / $NB"
RT=armsigm1s TAG=armsigm1 SRC=$SDIR/src.simple H=$W/simarmsigm1 MRA=$SDIR/meta_rows_A.txt MRB=$SDIR/meta_rows_B.txt CP=$W/benchcache_armsigm1 CTFABS=0 MSKDIAM=220 \
  KFLOOR=5000 MINP=15000 MAXD=4 BORROW=none PAR=${PAR:-1} NP=${NP:-8} bash $W/rflex13s.sbatch
