#!/bin/bash
#SBATCH -J klem_p3junk
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# is junk what spreads the realistic-SNR member masks over the molecule (problem 3)? rflex_s30a flex-state labels of node ND reconstructed at the
# node's poses from all their particles and from their non-junk particles only (truth: conformer 0 = junk), then p3_blur.py masks of both with truth.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$T/simswing30; O=$K_/p3junk
for ND in r22 r2; do N=$K_/rflex_s30a/node_$ND; W=$O/$ND; mkdir -p $W && cd $W
  awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; print ($1 > 0 && c[r] > 0) ? $1 : 0}' $S/gt_particles.txt $N/flex_states.txt > clean.txt
  awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; n[s]++; if (c[r] == 0) j[s]++} END {printf "'$ND' flex states (n, junk %%):"; for (s = 1; s <= 8; s++) printf " %d:%d/%.1f", s, n[s], 100 * j[s] / n[s]; print ""}' $S/gt_particles.txt $N/flex_states.txt
  for m in all clean; do lab=$N/flex_states.txt; [ $m = clean ] && lab=$W/clean.txt
    if [ -z "$(ls $W/rec_$m/*_reconstruct3D/recvol_state08_odd_unfil.mrc 2>/dev/null)" ]; then mkdir -p $W/rec_$m && cd $W/rec_$m; cp -f $N/n.simple r.simple
      simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$lab mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=8 mkdir=no > up.log 2>&1
      simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=220 nstates=8 nparts=4 nthr=8 > rec.log 2>&1 || { echo "$ND $m reconstruct3D FAILED"; continue; }; cd $W; fi
    RR=$(ls -td $W/rec_$m/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*
    awk -v d=$RR 'NR == FNR {if ($1 > 0) c[$1]++; next} END {for (s = 1; s <= 8; s++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n", s, d, s, d, s, c[s]}' $lab /dev/null > $W/list_$m.txt
    $PY $T/p3_blur.py $O --name ${ND}_$m --list $W/list_$m.txt --labels $lab --gt $S/gt_particles.txt --truth $S/maps --mskdiam 220 2>&1 | grep -v Warning | sed "s/^/$ND $m: /"
  done
done
echo DONE
