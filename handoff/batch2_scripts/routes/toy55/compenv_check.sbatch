#!/bin/bash
#SBATCH -J compenv_check
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# print the compenv segment of the given projects (args)
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
cd $(mktemp -d); for p in "$@"; do echo "== $p ($(stat -c '%y %s' $p))"; simple_exec prg=print_project_field projfile=$p oritype=compenv mkdir=no 2>&1 | grep -v "^$" | head -5; done
