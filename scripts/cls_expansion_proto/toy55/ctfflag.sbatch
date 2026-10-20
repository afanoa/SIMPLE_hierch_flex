#!/bin/bash
#SBATCH -J rnasep_ctfflag
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$PATH; cd /tmp
echo "== stk segment ctf flags of the 692k project"; simple_exec prg=print_project_field projfile=/mnt/beegfs/elmlund/afan/rnasep_coupled/full692k_2d/rnasep.simple oritype=stk mkdir=no 2>/dev/null | grep -o "ctf=[a-z]*" | sort | uniq -c
echo "== toy project (simulated)"; simple_exec prg=print_project_field projfile=/mnt/beegfs/elmlund/afan/rnasep_klem/toyah_tree/src.simple oritype=stk mkdir=no 2>/dev/null | grep -o "ctf=[a-z]*" | sort | uniq -c
