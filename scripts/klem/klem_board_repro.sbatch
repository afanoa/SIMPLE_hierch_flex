#!/bin/bash
#SBATCH -J klem_board_repro
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o %x.%j.out
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
K=/mnt/beegfs/elmlund/afan/rnasep_klem; S7=${REF:-$K/f692smoke7}; T=${OUT:-$K/klem_board_repro_$(date +%m%d_%H%M)}; mkdir -p $T; echo "test dir $T, reference $S7/parent"
ndiff_all=0
for arm in ${ARMS:-frozen:$K/f692e/klem2_bin current:${KLEM_BUILD:?set KLEM_BUILD}}; do name=${arm%%:*}; FB=${arm#*:}; O=$T/$name; mkdir -p $O
  cp $S7/src_local.simple $O/p.simple; cp $S7/parent/board.txt $O/board.txt
  ( cd $O; export SIMPLE_PATH=$FB PATH=$FB/bin:$PATH; simple_exec prg=update_project projfile=$O/p.simple qsys_name=local nparts=6 nthr=8 mkdir=no > $O/up.log 2>&1 ) || echo "$name: update_project FAILED"
  echo "== $name: $FB (simple_exec md5 $(md5sum $FB/bin/simple_exec | cut -c1-12))"; t0=$(date +%s)
  ( cd $O; export SIMPLE_PATH=$FB PATH=$FB/bin:$PATH OMP_NUM_THREADS=8 SIMPLE_QSYS=local
    simple_exec prg=klem_board projfile=$O/p.simple infile=$O/board.txt dir=$O nthr=8 nparts=6 mkdir=no > $O/board.log 2>&1 ) || echo "$name: klem_board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.log | cut -c1-200)"
  echo "$name: board done in $(( $(date +%s) - t0 )) s"
  n=0; d=0
  for f in klem_latent.txt klem_states.txt klem_resp.bin klem_halves.bin $(cd $S7/parent && ls cells/*.mrc); do
    n=$((n+1)); cmp -s $S7/parent/$f $O/$f || { d=$((d+1)); echo "$name: $f differs"; }
  done
  echo "$name: $n files compared, $d differ -> $([ $d = 0 ] && echo IDENTICAL || echo DIFFERENT)"; ndiff_all=$((ndiff_all + d)); rm -f $O/p.simple
done
echo "VERDICT: $([ $ndiff_all = 0 ] && echo IDENTICAL || echo DIFFERENT)"
