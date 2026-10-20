#!/bin/bash
#SBATCH -J klem_parenthp
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rerun a tree's parent board with the low-frequency cut 'hp' (env R = tree dir, HP in A, FB = build), locally on 16 threads,
# same config, start and consensus poses; output R/parent_hp<HP>[_ps<PSC>] (PSC: per-particle contrast scale on/off)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
R=${R:?}; HP=${HP:?}; FB=${FB:-$HOME/SIMPLE_klem2/build}; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
O=$R/parent_hp$HP${PSC:+_ps$PSC}; mkdir -p $O; grep -v "^hp \|^pscale " $R/parent/board.txt > $O/board.txt; echo "hp $HP" >> $O/board.txt; [ -n "${PSC:-}" ] && echo "pscale $PSC" >> $O/board.txt; cd $O; t0=$(date +%s)
[ -s $O/klem_latent.txt ] || simple_exec prg=klem_board projfile=$R/src.simple infile=$O/board.txt dir=$O nthr=16 nparts=1 mkdir=no > $O/board.out 2>&1 || { echo "board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.out | cut -c1-200)"; exit 1; }
echo "[$(date '+%F %T') parenthp] $R hp $HP board done in $(( $(date +%s) - t0 )) s"; grep -a ">>> it 9 (" $O/board.out | cut -c1-200
