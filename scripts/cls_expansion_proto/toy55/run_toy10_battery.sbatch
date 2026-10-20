#!/bin/bash
#SBATCH -J toy10_battery
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Latent-variable battery on the 10-conformer RNase P toy (truth-scored), run after the simulation chain:
#   0 oracle E-steps against the TRUTH conformer maps per band (the per-particle inference ceiling at this SNR)
#   1 flex_pca linear latent per particle: ordering, separability, closed recall, delivered states
#   2 coarse-to-fine latent windows: PC1 -> 8 windows -> window reconstructions -> per-window local pose refinement
#     (rounds 2+) -> soft cross-half E-step with a smoothness prior -> posterior-mean coordinate, 3 rounds (20, 14, 10 A)
#   3 deformation-latent prototype on the final window maps (flow PCA, flow-corrected consensus, grid E-step)
#   4 refine3D_states from 4 coarse windows (local poses) + mixture of local linear models (flex_pca per coarse group)
#   5 2D arm (solve2D + pose-mean expansion) for the sub-class-pooled readout of the particle latent
# env: SIMDIR (sim) PROJ (proj) CSV TAG NW (8) NP (10) K (4) LPS ("20 14 10") REFINE_WINDOWS (1) DO2D (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; BC=$HOME/SIMPLE_coupled/snap/9eaef62df; BF=$HOME/SIMPLE_master1013/snap/c6e344c6b_reseed27
usebuild(){ export SIMPLE_PATH=$1 PATH=$1/bin:$1/scripts:/mnt/nasapps/production/gcc/15.2.0/bin:/usr/bin:/bin SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8; }
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/toyexp/pv
export TOY_SMPD=1.464 OMP_NUM_THREADS=8
SIMDIR=${SIMDIR:-sim}; PJ=${PROJ:-proj}; CSV=$W/${CSV:-rnasep_conformers10.csv}; TAG=${TAG:-}; NW=${NW:-8}; NP=${NP:-10}; K=${K:-4}; LPS=${LPS:-"20 14 10"}; SEEDCOL=${SEEDCOL:-1}; FLEXSRC=${FLEXSRC:-}; SKIPORACLE=${SKIPORACLE:-0}
T=$W/$SIMDIR; OUT=$W/lat$TAG; mkdir -p $OUT; GT=$T/gt_particles.txt; STK=$T/all.mrcs; TRUTHTAB=$T/all.txt
STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $CSV); CLOSED=18,19,30,50,51,52,53
log(){ echo "[$(date '+%F %T') toy10$TAG] $*"; }
score(){ $PY - $1 $GT $2 <<'PYEOF'
import sys, numpy as np
from itertools import permutations
lab = {}
if sys.argv[3] == 'proj':
    for i, l in enumerate(open(sys.argv[1])):
        kv = dict(t.split('=') for t in l.split() if '=' in t); s = int(float(kv.get('state', 0)))
        if s > 0: lab[i + 1] = s
else:
    for l in open(sys.argv[1]):
        if l.startswith('#'): continue
        w = l.split(); lab[int(w[0])] = int(w[2])
gt = {}; closed = {}
for l in open(sys.argv[2]):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1]); closed[int(w[0])] = int(w[2])
p = np.array(sorted(lab)); L = np.array([lab[i] for i in p]); t = np.array([gt[i] for i in p]); c = np.array([closed[i] for i in p]); ks = np.unique(L); keep = t > 0
ts = np.unique(t[keep]); Lk, tk = L[keep], t[keep]
best = max(((sum(np.sum((Lk == k) & (tk == s)) for k, s in zip(pm, ts[:len(pm)])), pm) for pm in permutations(ks, min(len(ks), len(ts)))), key=lambda x: x[0]) if len(ks) <= len(ts) else (0, ())
cl = [(np.sum((L == k) & (c == 1)), k) for k in ks]; nc, kc = max(cl)
print('labelled %d (junk among them %d = %.3f) | accuracy vs %d conformers (K=%d labels) %.3f | closed: best label %d holds %.2f of closed particles at purity %.2f | sizes %s' %
      (len(p), int((t == 0).sum()), (t == 0).mean(), len(ts), len(ks), best[0] / max(keep.sum(), 1), kc, nc / max((c == 1).sum(), 1), nc / max(np.sum(L == kc), 1), ' '.join('%d:%d' % (k, np.sum(L == k)) for k in ks)))
PYEOF
}
RP=$(ls -td $W/$PJ/*_refine3D | head -1); PROJFILE=$(ls $RP/*.simple | head -1); CONS=$RP/recvol_state01.mrc; [ -s $CONS ] || CONS=$(ls $RP/recvol_state01*.mrc | head -1)
TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc $W/maps160/variants/conf$(printf %02d $c)_*.mrc 2>/dev/null | head -1; done | tr "\n" , | sed 's/,$//')
NPTCL=$(grep -vc "^#" $GT)
log "START: project $PROJFILE, consensus $CONS, states $STATES, $NPTCL particles, truth maps $(echo $TRUTH | tr , '\n' | wc -l)"
AXIS=$($PY - <<PYEOF
import numpy as np
def rd(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(fn, np.float32, offset=1024 + int(h[23]), count=nx*ny*nz).reshape(nz, ny, nx).astype(np.float64)
fns = '$TRUTH'.split(','); V = [rd(f) for f in fns]; B = V[0].shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx-B/2, yy-B/2), zz-B/2)
msk = r < 100/1.464; F = [np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))) for v in V]; band = (np.rint(r) >= B*1.464/60) & (np.rint(r) <= B*1.464/10)
X = np.array([np.concatenate([f[band].real, f[band].imag]) for f in F]); X -= X.mean(0); U, s, Vt = np.linalg.svd(X, full_matrices=False); pc = U[:, 0] * s[0]; pc = pc / np.abs(pc).max()
import sys; sys.stderr.write('conformer PCA variance shares %s; PC1 scores %s\n' % (np.array2string(s**2/np.sum(s**2), precision=2), np.array2string(pc, precision=2)))
print(','.join('%.3f' % v for v in pc))
PYEOF
); log "truth axis (PC1 of the conformers, STATES order): $AXIS"
LW="$PY $W/latent_windows.py"; SA="--states $STATES --axis=$AXIS"
usebuild $BK
simple_exec prg=print_project_field projfile=$PROJFILE oritype=ptcl3D mkdir=no 2>/dev/null > $OUT/ptcl3D_r0.txt
$PY $W/pose_error.py $OUT/ptcl3D_r0.txt $TRUTHTAB --tag "refined consensus poses" | sed "s/^/[toy10$TAG] /"
# ---- 0 oracle: truth references, per band (the ceiling of per-particle state inference at this SNR)
for LP in 20 14 10 8; do [ "$SKIPORACLE" = 1 ] && break
  [ -s $OUT/oracle_lp$LP.txt ] || $LW estep $STK $OUT/ptcl3D_r0.txt $GT --refs $TRUTH --lp $LP --out $OUT/oracle_lp$LP.txt $SA --nsub 40000 --tag "oracle lp$LP" 2>&1 | grep -v "^   [0-9]* / " | sed "s/^/[toy10$TAG] /"
done
# ---- 1 flex_pca linear latent on the particles
D1=$OUT/flex; [ -n "$FLEXSRC" ] && D1=$FLEXSRC
if [ ! -e $D1/flex.done ]; then mkdir -p $D1 && cd $D1 && cp -f $PROJFILE proj.simple; usebuild $BF
  simple_exec prg=flex_pca projfile=proj.simple vol1=$CONS mskdiam=200 lp=20 neigs=16 npreimages=$K state_placement=gmm rec_states=yes nparts=$NP nthr=8 > flex.log 2>&1 || {
    log "flex_pca (gmm placement) FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200); retrying with defaults"
    simple_exec prg=flex_pca projfile=proj.simple vol1=$CONS mskdiam=200 lp=20 neigs=16 npreimages=$K rec_states=yes nparts=$NP nthr=8 > flex2.log 2>&1 || { log "flex_pca FAILED again: $(grep -a -m2 'ERROR\|THROW' flex2.log | cut -c1-200)"; }; }
  touch $D1/flex.done; fi
FD=$(ls -td $D1/*_flex_pca 2>/dev/null | head -1)
if [ -n "$FD" ] && [ -s $FD/flex_pca_embedding.bin ]; then
  log "1 flex_pca lp20 neigs16: $FD"
  $PY $W/ptcl_latent_eval.py $FD $GT $SA --out $OUT/coord_r1.txt --tag "flex lp20" 2>&1 | sed "s/^/[toy10$TAG] /"
else log "no flex_pca embedding; the windows are seeded from the oracle lp20 posterior mean instead (flagged)"
  awk '!/^#/ { print $1, $2, 0 }' $OUT/oracle_lp20.txt > $OUT/coord_r1.txt; fi
# ---- 2 coarse-to-fine latent windows
usebuild $BK; r=0; COORD=$OUT/coord_r1.txt; TABLE=$OUT/ptcl3D_r0.txt; PROJCUR=$PROJFILE; COL=$SEEDCOL
for LP in $LPS; do r=$((r+1))
  log "2 round $r at $LP A: windows from $COORD"
  $LW assign $COORD $GT --nptcl $NPTCL --nw $NW --out $OUT/states_r$r.txt --info $OUT/windows_r$r.txt $SA --col $COL --tag "assign r$r (col $COL)" 2>&1 | sed "s/^/[toy10$TAG] /"; COL=1
  DR=$OUT/rec_r$r
  if [ ! -e $DR/rec.done ]; then mkdir -p $DR && cd $DR && cp -f $PROJCUR proj.simple
    simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$OUT/states_r$r.txt mkdir=no > sel.log 2>&1 || { log "selection r$r FAILED"; tail -2 sel.log; }
    simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $DR/ptcl3D.txt
    simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=200 nstates=$NW nparts=$NP nthr=8 > rec.log 2>&1 || { log "reconstruct3D r$r FAILED: $(grep -a -m2 'ERROR\|THROW' rec.log | cut -c1-200)"; break; }
    RD=$(ls -td $DR/*_reconstruct3D | head -1); mkdir -p $OUT/maps_r$r; for w in $(seq -f %02g 1 $NW); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $OUT/maps_r$r/ 2>/dev/null; done; done
    rm -f "${RD:?}"/recvol_state*_part* "${RD:?}"/*_unfil.mrc 2>/dev/null; touch $DR/rec.done; fi
  TABLE=$DR/ptcl3D.txt
  log "   window maps r$r: $(ls $OUT/maps_r$r/recvol_state??.mrc 2>/dev/null | wc -l) maps; cFAR $(grep -a 'cFAR' $DR/rec.log | awk '{print $NF}' | tr '\n' ' ' | cut -c1-100)"
  $LW maps $OUT/maps_r$r --nw $NW --info $OUT/windows_r$r.txt --truth $TRUTH $SA --cons $W/maps160/toy_consensus.mrc --tag "maps r$r" 2>&1 | sed "s/^/[toy10$TAG] /"
  # per-window local pose refinement against the window maps (rounds 2+), time-guarded
  if [ $r -ge 2 ] && [ "${REFINE_WINDOWS:-1}" = 1 ]; then TABS=""; t_first=0
    for w in $(seq 1 $NW); do DW=$OUT/ref_r${r}_w$w; W2=$(printf %02d $w)
      if [ ! -s $DW/ptcl3D.txt ]; then mkdir -p $DW && cd $DW && cp -f $DR/proj.simple p.simple; awk -v w=$w '{ print ($1 == w) ? 1 : 0 }' $OUT/states_r$r.txt > members.txt
        t0=$(date +%s); simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=members.txt mkdir=no > sel.log 2>&1
        simple_exec prg=refine3D projfile=p.simple vol1=$OUT/maps_r$r/recvol_state$W2.mrc nstates=1 refine=shc objfun=cc pgrp=c1 mskdiam=200 lp=$LP maxits=2 nparts=$NP nthr=8 > refine.log 2>&1 || { log "   refine3D window $w FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-160)"; continue; }
        RW=$(ls -td $DW/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RW/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $DW/ptcl3D.txt; rm -f "${RW:?}"/recvol_state01_part* "${RW:?}"/*_unfil.mrc 2>/dev/null
        dt=$(( $(date +%s) - t0 )); [ $w = 1 ] && t_first=$dt; log "   refine3D window $w: $dt s"
        if [ $t_first -gt 2400 ]; then log "   window refinement too slow ($t_first s per window): the remaining windows keep their poses"; break; fi
      fi
      [ -s $DW/ptcl3D.txt ] && TABS="$TABS,$w:$DW/ptcl3D.txt"; done
    if [ -n "$TABS" ]; then $LW merge $DR/ptcl3D.txt $OUT/states_r$r.txt --tables ${TABS#,} --out $OUT/ptcl3D_r$r.txt | sed "s/^/[toy10$TAG] /"
      $PY $W/pose_error.py $DR/ptcl3D.txt $TRUTHTAB --tag "   poses before window refinement r$r" | sed "s/^/[toy10$TAG] /"; $PY $W/pose_error.py $OUT/ptcl3D_r$r.txt $TRUTHTAB --tag "   poses after window refinement r$r" | sed "s/^/[toy10$TAG] /"
      TABLE=$OUT/ptcl3D_r$r.txt
      DP=$OUT/proj_r$r; mkdir -p $DP && cd $DP && simple_exec prg=new_project projname=p > np.log 2>&1 && cd p && simple_exec prg=import_particles projfile=p.simple stk=$STK oritab=$OUT/ptcl3D_r$r.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 \
        && simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1 && PROJCUR=$DP/p/p.simple || log "   re-import with refined poses FAILED; next round keeps the previous project"
    fi
  fi
  $LW estep $STK $TABLE $GT --mapdir $OUT/maps_r$r --nw $NW --lp $LP --prev $OUT/states_r$r.txt --prior_sigma 1.5 --out $OUT/estep_r$r.txt $SA --tag "estep r$r lp$LP" 2>&1 | grep -v "^   [0-9]* / " | sed "s/^/[toy10$TAG] /"
  [ -s $OUT/estep_r$r.txt ] || { log "E-step r$r produced no coordinate; stopping the rounds"; break; }
  COORD=$OUT/estep_r$r.txt
done
RLAST=$r; log "2 done: $RLAST rounds; final coordinate $COORD, poses $TABLE"
# ---- 3 deformation-latent prototype on the final window maps
[ -s $OUT/deform/coord_deform.txt ] || $PY $W/deform_proto.py $OUT/maps_r$RLAST $OUT/deform --nw $NW --info $OUT/windows_r$RLAST.txt --cons $W/maps160/toy_consensus.mrc --truth $TRUTH $SA \
    --stk $STK --table $TABLE --gt $GT --nsub 20000 --lp 12 --ngrid 9 --nproc 8 --sigma 3 --tag "3 deform" 2>&1 | grep -v "^   [0-9]* / " | sed "s/^/[toy10$TAG] /"
# ---- 4 refine3D_states from K coarse windows of the final coordinate (local poses) + mixture of local linear models
$LW assign $COORD $GT --nptcl $NPTCL --nw $K --out $OUT/states_K.txt --info $OUT/windows_K.txt $SA --tag "4 coarse groups" 2>&1 | sed "s/^/[toy10$TAG] /"
D4=$OUT/r3s
if [ ! -e $D4/done ]; then mkdir -p $D4 && cd $D4 && cp -f $PROJCUR p.simple; usebuild $BK
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$OUT/states_K.txt mkdir=no > sel.log 2>&1
  simple_exec prg=reconstruct3D projfile=p.simple mskdiam=200 pgrp=c1 nstates=$K nparts=$NP nthr=8 > rec0.log 2>&1 || log "   reconstruct3D (K groups) FAILED"
  VOLS=""; for s in $(seq 1 $K); do VOLS="$VOLS vol$s=$D4/recvol_state0$s.mrc"; done
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=200 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1   # this build derives the state volumes from the labels (vol1..N rejected) || log "   refine3D_states FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-200)"
  RD=$(ls -td $D4/*_refine3D_states 2>/dev/null | head -1); [ -n "$RD" ] && simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D4/ptcl3D_final.txt; touch $D4/done; fi
awk '{print "state="$1}' $OUT/states_K.txt > $OUT/states_K_proj.txt; log "4 seed = the $K coarse windows: $(score $OUT/states_K_proj.txt proj)"
[ -s $D4/ptcl3D_final.txt ] && { log "4 refine3D_states from the $K coarse windows: $(score $D4/ptcl3D_final.txt proj)"
  RD=$(ls -td $D4/*_refine3D_states | head -1); MS=$(ls $RD/recvol_state0?.mrc 2>/dev/null | tr "\n" , | sed 's/,$//'); [ -n "$MS" ] && $PY $P/pv_readout.py maps 1.464 15 200 $TRUTH $MS 2>&1 | grep -v " - " | head -6 | cut -c1-300 | sed "s/^/[toy10$TAG] /"; }
if [ -n "$FD" ]; then usebuild $BF
  for g in $(seq 1 $K); do DG=$OUT/flex_g$g; awk -v g=$g '{ print ($1 == g) ? 1 : 0 }' $OUT/states_K.txt > $OUT/members_g$g.txt
    if [ ! -e $DG/done ]; then mkdir -p $DG && cd $DG && cp -f $PROJCUR proj.simple
      simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$OUT/members_g$g.txt mkdir=no > sel.log 2>&1
      simple_exec prg=flex_pca projfile=proj.simple vol1=$D4/recvol_state0$g.mrc mskdiam=200 lp=20 neigs=8 npreimages=2 rec_states=no nparts=$NP nthr=8 > flex.log 2>&1 || log "   flex_pca group $g FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-160)"
      touch $DG/done; fi
    FG=$(ls -td $DG/*_flex_pca 2>/dev/null | head -1)
    [ -n "$FG" ] && [ -s $FG/flex_pca_embedding.bin ] && { $PY $W/ptcl_latent_eval.py $FG $GT $SA --subset $OUT/members_g$g.txt --tag "4 local flex g$g" 2>&1 | grep -E "axis:|knn:|closed" | sed "s/^/[toy10$TAG] /"
      $PY $W/ptcl_latent_eval.py $FD $GT $SA --subset $OUT/members_g$g.txt --tag "4 global flex on g$g" 2>&1 | grep -E "axis:|knn:|closed" | sed "s/^/[toy10$TAG] /"; }
  done; fi
# ---- 5 2D arm: solve2D + pose-mean expansion -> sub-class-pooled readout of the particle latent
if [ "${DO2D:-1}" = 1 ] && [ -n "$FD" ]; then usebuild $BK; D2=$OUT/2D
  if [ ! -e $D2/solve2D.done ]; then mkdir -p $D2 && cd $D2 && cp -f $PROJFILE proj.simple
    simple_exec prg=solve2D projfile=proj.simple ncls=${NCLS2D:-100} mskdiam=200 nparts=$NP nthr=8 mkdir=no > solve2D.log 2>&1 || { log "5 solve2D FAILED"; grep -a -m3 "ERROR\|THROW" solve2D.log; }
    simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl2D mkdir=no 2>/dev/null > ptcl2D.txt; simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > ptcl3D.txt; touch solve2D.done; fi
  D5=$OUT/exppm
  if [ ! -e $D5/expansion.done ] && [ -s $D2/ptcl2D.txt ]; then mkdir -p $D5/in && cd $D5/in; usebuild $BC
    simple_exec prg=new_project projname=proj > np1.log 2>&1 && cd proj && simple_exec prg=import_particles projfile=proj.simple stk=$STK oritab=$D2/ptcl2D.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1
    cd $D5/in && simple_exec prg=new_project projname=poses > np2.log 2>&1 && cd poses && simple_exec prg=import_particles projfile=poses.simple stk=$STK oritab=$D2/ptcl3D.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1
    cd $D5/in/proj && simple_exec prg=replace_project_field projfile=proj.simple projfile_target=$D5/in/poses/poses.simple oritype=ptcl3D > replace.log 2>&1
    simple_exec prg=update_project projfile=proj.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
    cd $D5 && cp -f $D5/in/proj/proj.simple proj.simple
    env SIMPLE_FLEXCLS_HALF_CONSISTENT=yes SIMPLE_FLEXCLS_DIRECT_AVGS=yes simple_exec prg=cls_expansion projfile=proj.simple ncls=16 lp=20 mskdiam=200 vol1=$CONS nparts=$NP nthr=8 > clsexp.log 2>&1 || { log "5 pose-mean cls_expansion FAILED"; grep -a -m3 "ERROR\|THROW" clsexp.log; }
    EX=$(ls -td $D5/*_cls_expansion 2>/dev/null | head -1); [ -n "$EX" ] && simple_exec prg=print_project_field projfile=$EX/proj.simple oritype=ptcl2D mkdir=no 2>/dev/null > $D5/ptcl2D.txt; touch $D5/expansion.done; fi
  [ -s $D5/ptcl2D.txt ] && { log "5 pose-mean expansion: $($PY $P/pv_readout.py exp $T $D5/ptcl2D.txt 2>&1 | tail -1 | cut -c1-300)"
    $PY $W/ptcl_latent_eval.py $FD $GT $SA --sublab $D5/ptcl2D.txt --tag "5 flex lp20 pooled" 2>&1 | grep -E "pooled|knn:" | sed "s/^/[toy10$TAG] /"; }
fi
# ---- 6 finer linear latent (lp 12): does a finer basis order the particles better than lp 20?
D6=$OUT/flex12
if [ ! -e $D6/flex.done ]; then mkdir -p $D6 && cd $D6 && cp -f $PROJFILE proj.simple; usebuild $BF
  simple_exec prg=flex_pca projfile=proj.simple vol1=$CONS mskdiam=200 lp=12 neigs=16 npreimages=$K state_placement=gmm rec_states=no nparts=$NP nthr=8 > flex.log 2>&1 || log "6 flex_pca lp12 FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"
  touch $D6/flex.done; fi
F6=$(ls -td $D6/*_flex_pca 2>/dev/null | head -1); [ -n "$F6" ] && [ -s $F6/flex_pca_embedding.bin ] && $PY $W/ptcl_latent_eval.py $F6 $GT $SA --out $OUT/coord_lp12.txt --tag "6 flex lp12" 2>&1 | grep -E "axis:|ladder|knn:|closed|states:" | sed "s/^/[toy10$TAG] /"
log "DONE toy10$TAG battery -> $OUT"
