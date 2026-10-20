#!/bin/bash
#SBATCH -J klem_tree
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH --mail-type=NONE
#SBATCH -o %x.%j.out
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${KLEM_BUILD:?set KLEM_BUILD to the SIMPLE build directory}; C=${CFG:?set CFG to the klem_tree config}
R=$(awk '$1 == "run" {print $2}' $C); [ -n "$R" ] || { echo "no run directory in $C"; exit 1; }
[ ! -d $R ] || [ -e $R/klem_bin ] || [ -z "$(ls -A $R)" ] || { echo "$R holds another run: use a new run directory"; exit 1; }; mkdir -p $R
[ -x $R/klem_bin/bin/simple_exec ] || { mkdir -p $R/klem_bin && cp -r $B/bin $B/scripts $R/klem_bin/; }
export SIMPLE_PATH=$R/klem_bin PATH=$R/klem_bin/bin:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${SIMPLE_QSYS_PARTITION:-csbdevel}
simple_private_exec prg=klem_tree projfile=${PROJ:?set PROJ to the project file} infile=$C nthr=8 mkdir=no
