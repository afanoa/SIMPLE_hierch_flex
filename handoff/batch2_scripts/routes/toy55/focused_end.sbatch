#!/bin/bash
#SBATCH -J klem_focused_end
#SBATCH -p norm
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# (1) half-map reproducibility of the arm-focused latent's equal-width end bins (closed-arm readout of each bin's even and odd half-map);
# (2) finer end sampling: 16 equal-width bins over the 0.25..99.75 % range of the focused latent's view-centred PC1 -> reconstruct -> readout.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=$K_/focused_flex_f692rf_r22; R=$K_/rflex_f692rf; CK=$K_/cheapkids_f692rf_r22
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); set -- $(awk '$1 == "r22"' $R/queue.txt); PP=$3; ROWS=$4; FD=$(ls -td $O/*_flex_pca | head -1)
echo "== (1) half-map readout of the 12 equal-width bins (even / odd)"
RR=$(ls -td $O/rec_w1/*_reconstruct3D | head -1); mkdir -p $O/halfview
for h in even odd; do : > $O/spec_w1_$h.txt; for s in $(seq 1 12); do b=$RR/recvol_state$(printf %02d $s); echo "bin_$(printf %02d $s)_$h ${b}_${h}_unfil.mrc ${b}_${h}_unfil.mrc" >> $O/spec_w1_$h.txt; done
  $PY $T/make_view.py $O/halfview/$h $O/spec_w1_$h.txt > /dev/null 2>&1; $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_w1_$h.txt $(ls $O/halfview/$h/*.mrc) > /dev/null 2>&1; done
paste <(awk '!/^#/ {printf "%s %.2f %.2f\n", $1, $2, $3}' $O/arm_w1_even.txt) <(awk '!/^#/ {printf "%.2f %.2f\n", $2, $3}' $O/arm_w1_odd.txt) | awk '{n = split($1, p, "/"); printf "   %s: even closed %.2f open %.2f | odd closed %.2f open %.2f\n", p[n], $2, $3, $4, $5}'
echo "== (2) 16 equal-width bins over 0.25..99.75 %"
[ -s $O/children_w16.txt ] || { $PY $T/cheapkids_rf2.py --oris $CK/ptcl3D_parent.txt --rows $ROWS --flex $FD/flex_pca_coordinates.txt --out $O/w16 --k 16 --axis w1 --p 0.25 2>&1 | grep -v Warning; cp $O/w16/children_w1.txt $O/children_w16.txt; }
mkdir -p $O/rec_w16 && cd $O/rec_w16
if [ -z "$(ls $O/rec_w16/*_reconstruct3D/recvol_state16.mrc 2>/dev/null)" ]; then cp -f $PP r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/children_w16.txt mkdir=no > rsel.log 2>&1 \
  && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > rup.log 2>&1 && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=16 nparts=4 nthr=16 > rec.log 2>&1 || { echo "reconstruct FAILED"; exit 1; }; fi
RR=$(ls -td $O/rec_w16/*_reconstruct3D | head -1); V=$K_/view_rflex/focused_flex_f692rf_r22_w16; : > $O/spec_w16.txt
for s in $(seq 1 16); do b=$RR/recvol_state$(printf %02d $s); echo "bin_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec_w16.txt; done
$PY $T/make_view.py $V $O/spec_w16.txt > /dev/null 2>&1; $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_w16.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1
awk 'NR == FNR {n[$1]++; next} !/^#/ {k = $1; sub(/.*bin_/, "", k); sub(/\.mrc/, "", k); printf "   bin %s (%d ptcls): closed-arm %.2f open-arm %.2f %s\n", k, n[k + 0], $2, $3, $4}' $O/children_w16.txt $O/arm_w16.txt
echo DONE
