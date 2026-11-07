#!/bin/bash
#SBATCH -J compenv_repair
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# restore the queue settings of a project whose compenv was overwritten by a concurrent job-record write (arg: PROJ); keeps a backup
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm; P=$1; cp -p $P $P.bak_compenv; cd $(dirname $P)
simple_exec prg=update_project projfile=$P qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 mkdir=no > $P.repair.log 2>&1 || { echo "update_project FAILED"; tail -3 $P.repair.log; exit 1; }
cd $(mktemp -d); simple_exec prg=print_project_field projfile=$P oritype=compenv mkdir=no 2>&1 | grep -v "^$" | head -3
