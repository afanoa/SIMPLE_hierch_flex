#!/bin/bash
#SBATCH -J klem_p2prep
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# prep of the E-step test (problem 2): one project with every particle at its leaf node's refined pose (rflex_s30a/post/
# poses_combined.txt) and its leaf label (final_labels.txt, 34 leaves), the leaves' half-maps at box 88 (reconstruct3D box_crop 88) =
# the it1 maps of every arm, and the truth readout of the starting point (p2_metrics.py).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_s30a; S=$T/simswing30
O=$K_/p2estep/prep; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') p2prep] $*"; }
if [ ! -s $O/p2/p2.simple ]; then rm -rf $O/p2; simple_exec prg=new_project projname=p2 > np.log 2>&1; cd $O/p2
  sed -E 's/ state=[^ ]*//' $R/post/poses_combined.txt | sed -E 's/$/ state=1/' > $O/poses_import.txt
  simple_exec prg=import_particles projfile=p2.simple stk=$S/all.mrcs oritab=$O/poses_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; tail -3 import.log; exit 1; }
  simple_exec prg=selection projfile=p2.simple oritype=ptcl3D infile=$R/final_labels.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p2.simple qsys_name=local nparts=4 nthr=12 mkdir=no > up.log 2>&1; fi
cd $O/p2
if [ -z "$(ls $O/p2/*_reconstruct3D/recvol_state34_odd_unfil.mrc 2>/dev/null)" ]; then log "reconstructing the 34 leaves at box 88"
  simple_exec prg=reconstruct3D projfile=p2.simple pgrp=c1 mskdiam=220 nstates=34 box_crop=88 nparts=4 nthr=12 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }; fi
RR=$(ls -td $O/p2/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*
RP=$RR/p2.simple; [ -s $RP ] || RP=$O/p2/p2.simple
simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null > $O/dump_it1.txt
awk '{for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); print a[2] + 0}}' $O/dump_it1.txt > $O/eo_new.txt
awk '{for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); print a[2] + 0}}' $R/post/poses_combined.txt > $O/eo_old.txt
log "rows $(wc -l < $O/dump_it1.txt); eo agreement with the node refinements' partition: $(paste $O/eo_new.txt $O/eo_old.txt | awk '$1 == $2 {a++} END {printf "%.3f", a / NR}')"
$PY $T/p2_metrics.py $O/dump_it1.txt --gt $S/gt_particles.txt --true $S/all.txt --ref0 $R/final_labels.txt --out $O/labels.txt --tag "start (combined node poses, 34 leaves)" 2>&1 | grep -v Warning
awk -v d=$RR 'NR == FNR {n[$1]++; next} END {for (k = 1; k <= 34; k++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n", k, d, k, d, k, n[k]}' $O/labels.txt /dev/null > $O/list.txt
cp -f $RP $O/proj_it1.simple && touch $O/prep.done; log DONE
