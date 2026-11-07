# sourced by the calE_*.sbatch jobs (set NP NT TAG first): SIMPLE on snapshot m8c52t with the local qsys, and
#   newdir NAME    a fresh OUT/NAME (an old one is moved aside) with halves/ and work/, entered at work/
#   recon NAME PROJECT RECB0_SH MASK [TRUTHKEY]   PROJECT with the states in work/states.txt, RECB0_SH's reconstruct3D with nstates=5
#                  and postprocess=no; maps to NAME/bin1..5.mrc and halves/, then phaseE_maptests.py in the background
#   finish         waits, then collects every OUT/*/maptests.txt into OUT/maptests.txt and the SUMMARY lines into OUT/summary.txt
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=${NT:?} SIMPLE_SEED=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
log(){ echo "[$(date "+%F %T") ${TAG:-calE}] $*"; }
newdir(){ W=${OUT:?}/$1; [ -e $W ] && mv $W ${W}_old_$(date +%s); mkdir -p $W/halves $W/work; cd $W/work; }
recon(){ local nm=$1 pj=$2 rs=$3 mk=$4 key=${5:-} REC md t0 k s; cd $W/work
  cp $pj p.simple; simple_exec prg=update_project projfile=p.simple qsys_name=local mkdir=no > up.log 2>&1
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=states.txt mkdir=no > sel.log 2>&1 || { log "$nm selection FAILED"; return 1; }
  REC=$(grep -a -o "prg=reconstruct3D[^>]*" $rs | sed -e 's/projfile=[^ ]*/projfile=p.simple/' -e 's/nstates=[^ ]*/nstates=5/' -e "s/nparts=[^ ]*/nparts=$NP/" -e "s/nthr=[^ ]*/nthr=$NT/")
  REC="$REC postprocess=no"; md=$(echo "$REC" | grep -o "mskdiam=[0-9.]*" | cut -d= -f2); t0=$(date +%s)
  simple_exec $REC > rec.log 2>&1 || { log "$nm reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log)"; return 1; }
  for k in 1 2 3 4 5; do s=$(printf %02d $k); cp recvol_state$s.mrc ../bin$k.mrc; cp recvol_state${s}_even.mrc ../halves/bin${k}_even.mrc; cp recvol_state${s}_odd.mrc ../halves/bin${k}_odd.mrc; done
  log "$nm reconstruct3D $(( $(date +%s) - t0 )) s: $REC"
  ( $PY $T/phaseE_maptests.py $nm:$W:$mk:$md${key:+:$key} > $W/maptests.txt 2>&1 ) & }
finish(){ wait; for d in $(ls -d $OUT/*/ | grep -v _old_ | sort); do [ -s $d/maptests.txt ] && grep -v Warning $d/maptests.txt; done > $OUT/maptests.txt
  grep "^SUMMARY" $OUT/maptests.txt > $OUT/summary.txt; log "DONE: $(wc -l < $OUT/summary.txt) series scored"; }
