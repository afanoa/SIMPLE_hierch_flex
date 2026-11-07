#!/bin/bash
#SBATCH -J klem_mseries
#SBATCH -p csbdevel
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# motion series: per SERIES entry 'name:tree_dir:site:angle:variant[:box_crop]', the tree's root stage-B project (own poses, contrast
# field kept) with each row's state = its bin along the site posterior (motion_bins.py; 0 outside the site's level), the tree's own
# recB0 reconstruct3D line with nstates=5, postprocess=no (unsharpened) [and box_crop]; maps OUTROOT/name/bin1..5.mrc, halves in
# halves/, scored by motion_score.py. env SERIES OUTROOT NP NT [SNAP, default m8c52t] [TRUTH: IgM truth table for the angle stats]
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/${SNAP:-m8c52t}
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=${NT:?} SIMPLE_SEED=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
log(){ echo "[$(date "+%F %T") mseries] $*"; }
for e in ${SERIES:?}; do IFS=: read nm D site angle var crop <<< "$e"
  W=${OUTROOT:?}/$nm; [ -e $W ] && mv $W ${W}_old_$(date +%s); mkdir -p $W/halves $W/work; cd $W/work
  cp $D/node_r/stageB/b.simple p.simple
  simple_exec prg=update_project projfile=p.simple qsys_name=local mkdir=no > up.log 2>&1
  $PY $T/motion_bins.py $D/factorial/post/$site.txt states.txt $var ${TRUTH:+$TRUTH $angle} > ../bins.txt 2>&1 || { log "$nm bins FAILED: $(tail -1 ../bins.txt)"; continue; }
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=states.txt mkdir=no > sel.log 2>&1 || { log "$nm selection FAILED"; continue; }
  REC=$(grep -a -o "prg=reconstruct3D[^>]*" $D/node_r/recB0/hetero_step.sh | sed -e 's/projfile=[^ ]*/projfile=p.simple/' -e 's/nstates=[^ ]*/nstates=5/' -e "s/nparts=[^ ]*/nparts=$NP/" -e "s/nthr=[^ ]*/nthr=$NT/")
  REC="$REC postprocess=no${crop:+ box_crop=$crop}"; t0=$(date +%s)
  simple_exec $REC > rec.log 2>&1 || { log "$nm reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log)"; continue; }
  log "$nm reconstruct3D $(( $(date +%s) - t0 )) s: $REC"
  for k in 1 2 3 4 5; do s=$(printf %02d $k); cp recvol_state$s.mrc ../bin$k.mrc; cp recvol_state${s}_even.mrc ../halves/bin${k}_even.mrc; cp recvol_state${s}_odd.mrc ../halves/bin${k}_odd.mrc; done
  $PY $T/motion_score.py $W $D/node_r/mask.mrc $(echo "$REC" | grep -o "mskdiam=[0-9.]*" | cut -d= -f2) > ../score.txt 2>&1
  grep -a "POSTPROCESS: FSC\|RESOLUTION\|FSC AREA" rec.log | sed 's/^/   simple: /' >> ../score.txt
  { echo "######## $nm ($site, $angle, $var${crop:+, box_crop $crop})"; cat ../bins.txt ../score.txt; } | tee -a $OUTROOT/summary.txt
done
log DONE
