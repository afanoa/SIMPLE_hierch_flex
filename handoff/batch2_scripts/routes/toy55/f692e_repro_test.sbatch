#!/bin/bash
#SBATCH -J f692e_repro
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Reproduction test of the f692e build: the Oct 17 smoke parent board (f692smoke7/parent, made with the binaries the f692e
# run froze) rerun with the exact config (a) on the frozen f692e binaries and (b) on the current build; both compared
# with the Oct 17 output. Each run on a private copy of the smoke project (klem_board rewrites its projfile).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
S7=$K/f692smoke7; T=$K/f692e_repro_$(date +%m%d_%H%M); mkdir -p $T; echo "test dir $T"
for arm in frozen:$K/f692e/klem2_bin current:$HOME/SIMPLE_klem2/build; do name=${arm%%:*}; FB=${arm#*:}; O=$T/$name; mkdir -p $O
  cp $S7/src_local.simple $O/p.simple; cp $S7/parent/board.txt $O/board.txt
  ( cd $O; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH; simple_exec prg=update_project projfile=$O/p.simple qsys_name=local nparts=6 nthr=8 mkdir=no > $O/up.log 2>&1 ) || echo "$name: update_project FAILED"
  echo "== $name: $FB (simple_exec md5 $(md5sum $FB/bin/simple_exec | cut -c1-12))"; t0=$(date +%s)
  ( cd $O; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH OMP_NUM_THREADS=8 TOY_SMPD=1.464 SIMPLE_QSYS=local
    simple_exec prg=klem_board projfile=$O/p.simple infile=$O/board.txt dir=$O nthr=8 nparts=6 mkdir=no > $O/board.log 2>&1 ) || echo "$name: klem_board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.log | cut -c1-200)"
  echo "$name: board done in $(( $(date +%s) - t0 )) s"; grep -a ">>> it [0-9] (" $O/board.log | cut -c1-200
  $PY $W/repro_compare.py $S7/parent $O 2>&1 | grep -v Warning; rm -f $O/p.simple
done
echo "== frozen vs current"; $PY $W/repro_compare.py $T/frozen $T/current 2>&1 | grep -v Warning
echo "== Oct 17 reference iterations"; grep -a ">>> it [0-9] (" $S7/parent/board.log | cut -c1-200
