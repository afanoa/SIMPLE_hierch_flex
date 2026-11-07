#!/bin/bash
#SBATCH -J klem_murec
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merged-up landscape of a run (env RTAG, CDIR with states.txt + groups.txt from merge_up.sbatch) -> reconstruct3D
# of every merged state from all its particles at the family-1 refined poses -> gallery. Output $R/consolidated_$CUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; CUT=${CUT:-mu}; C=$R/${CDIR:?}; log(){ echo "[$(date '+%F %T') consolidate] $*"; }
NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/fam1/meta_fb.npz')['flags']))")
[ -s $C/states.txt ] && [ -s $C/groups.txt ] || { log "needs $C/states.txt and $C/groups.txt (merge_up.sbatch)"; exit 1; }
K=$(grep -vc "^#" $C/groups.txt); L=$C/fam0/regionall/level1; mkdir -p $L/states_dir $L/maps $C/rec && cd $C/rec
# poses: each family's particles are reconstructed with that family's refined project (same frame, local refinements from
# the consensus); a merged state with particles in both families keeps the map of the family holding more of them
CH=${CHUNK:-60}; NCH=$(( (K + CH - 1) / CH ))   # reconstruct3D takes at most ~100 states at once: chunks of CH
for F in 1 2; do
  PS=$(ls -td $R/fam$F/refine/*_refine3D_states 2>/dev/null | head -1)/p.simple; [ -f "$PS" ] || PS=$R/fam$F/refine/p.simple; [ -f "$PS" ] || PS=${SRC:-$A/rnasep_coupled/full692k_2d/rnasep.simple}; log "family $F poses from $PS"
  for q in $(seq 0 $((NCH - 1))); do
    $PY - <<PYEOF
import numpy as np
st = np.loadtxt('$C/states.txt', dtype=int); sub = np.loadtxt('$R/expand/family${F}_subset.txt', dtype=int).ravel() - 1; m = np.zeros(len(st), bool); m[sub] = True
lo, hi = $q * $CH, ($q + 1) * $CH; loc = np.where(m & (st > lo) & (st <= hi), st - lo, 0); np.savetxt('$C/states_f${F}_q$q.txt', loc, fmt='%d')
if $q == 0: print('family $F: %d particles in merged states' % int((st[m] > 0).sum()))
PYEOF
    NS=$(( K - q * CH )); [ $NS -gt $CH ] && NS=$CH; mkdir -p $C/rec$F/q$q && cd $C/rec$F/q$q && cp -f $PS p.simple
    EO=$R/fam$F/eo_rows.txt; [ -s $EO ] || { simple_exec prg=print_project_field projfile=$PS oritype=ptcl3D mkdir=no 2>/dev/null | awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}' > $EO.$$ && mv -f $EO.$$ $EO; }
    h=${MINHALF:-10}; E=$EO; [ "$(grep -c . $EO 2>/dev/null)" = "$(grep -c . $C/states_f${F}_q$q.txt)" ] || { h=0; E=/dev/null; }   # a state of a few particles or an (almost) empty even/odd half crashes reconstruct3D: no map for it
    awk -v m=${MINREC:-50} -v h=$h 'FILENAME == ARGV[1] {e[FNR] = $1; next} {l[FNR] = $1; n = FNR; c[$1]++; if (e[FNR] == 0) ce[$1]++; else if (e[FNR] == 1) co[$1]++} END {for (i = 1; i <= n; i++) {s = l[i]; print ((s > 0 && (c[s] < m || ce[s] < h || co[s] < h)) ? 0 : s)}}' $E $C/states_f${F}_q$q.txt > $C/rec$F/q$q/states_rec.txt
    awk '$1 > 0 {f = 1; exit} END {exit !f}' $C/rec$F/q$q/states_rec.txt || { log "family $F chunk $q: no particles in its states, skipped"; continue; }
    simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$C/rec$F/q$q/states_rec.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED (fam $F chunk $q)"; exit 1; }
    log "reconstruct3D family $F: states $((q * CH + 1))-$((q * CH + NS)) of $K (cut $CUT)"; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$NS nparts=${NPREC:-6} nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED (fam $F chunk $q)"; exit 1; }
  done
done
$PY - <<PYEOF
import numpy as np, glob, shutil, os
st = np.loadtxt('$C/states.txt', dtype=int); K = $K; CH = $CH
for F in (1, 2):
    sub = np.loadtxt('$R/expand/family%d_subset.txt' % F, dtype=int).ravel() - 1; m = np.zeros(len(st), bool); m[sub] = True; globals()['n%d' % F] = np.bincount(st[m], minlength=K + 1)
for s in range(1, K + 1):
    F = 1 if n1[s] >= n2[s] else 2; q = (s - 1) // CH; sl = s - q * CH; RR = (sorted(glob.glob('$C/rec%d/q%d/*_reconstruct3D' % (F, q))) or [''])[-1]
    for h in ('', '_even', '_odd'):
        src = '%s/recvol_state%02d%s.mrc' % (RR, sl, h)
        if os.path.isfile(src): shutil.copy(src, '$L/maps/recvol_state%02d%s.mrc' % (s, h))
    if not os.path.isfile('%s/recvol_state%02d.mrc' % (RR, sl)): print('merged state %d: no map (switched off: too few particles or an empty half)' % s); continue
    shutil.copy('%s/recvol_state%02d.mrc' % (RR, sl), '$L/states_dir/state_%02d.mrc' % s); print('merged state %d: family %d map (%d + %d particles)' % (s, F, n1[s], n2[s]))
PYEOF
awk '!/^#/{print $1, 0, $2, $2}' $C/groups.txt > $L/clusters.txt
$PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $L/state_matches.txt 2>&1 | grep -v Warning | tail -2
$PY $W/volume_gallery.py $C $W/maps160/toy_consensus.mrc $W/figs/gallery_${RTAG}_consolidated_$CUT${CTAG:-} --nrow 12 --lp ${LP:-7} 2>&1 | grep -v Warning
log "DONE: $K merged states, maps in $L/states_dir, groups $C/groups.txt"
