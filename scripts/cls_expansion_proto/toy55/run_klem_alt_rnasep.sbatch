#!/bin/bash
#SBATCH -J rnasep_klem_alt
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Pose alternation on real RNase P: klem latent -> NB bins -> refine3D_states (local 20->8 A, poses start from the project) ->
# poses into a copy of the metadata -> klem on the refined poses -> sanity readout. env: KDIR (start latent dir), NB (5),
# SUBSET (klem_latent.txt restricting the EM to those particles, e.g. a 50k gate; "" = all), KTAG, INIT (npz), MAXITS (6), NP (10), NROUND (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; NB=${NB:-5}; NP=${NP:-10}; NROUND=${NROUND:-1}; KTAG=${KTAG:-_alt}
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; META=$A/rnasep_coupled/basis_span_all/meta_all.npz; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
KDIR=${KDIR:-$A/rnasep_klem/klem_full_npz}; SUBSET=${SUBSET:-}; ALT=$A/rnasep_klem/alt$KTAG; mkdir -p $ALT; log(){ echo "[$(date '+%F %T') alt$KTAG] $*"; }
for R in $(seq 1 $NROUND); do D=$ALT/round$R; mkdir -p $D && cd $D; export OMP_NUM_THREADS=8
  $PY $W/klem_bins.py $KDIR $META $D/states.txt --nb $NB
  cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; tail -2 sel.log; exit 1; }
  simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/seed_ptcl3D.txt
  NR=$(wc -l < $D/seed_ptcl3D.txt); NS=$(wc -l < $D/states.txt); [ "$NR" = "$NS" ] || { log "project rows $NR != metadata rows $NS: row mapping unsafe, abort"; exit 1; }
  log "round $R: $(awk '{for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); if(a[2]+0>0) n++}} END{print n}' $D/seed_ptcl3D.txt) particles in $NB bins enter refine3D_states (local 20->8 A, maxits ${MAXITS:-6})"
  simple_exec prg=refine3D_states projfile=p.simple nstates=$NB pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=200 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
  RD=$(ls -td $D/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt
  $PY - <<PYEOF
import numpy as np
M = dict(np.load('$META')); keys = ('e1', 'e2', 'e3', 'x', 'y'); new = {k: M[k].copy() for k in keys}; n = 0; dang = []; dsh = []
def mat(e1, e2, e3):
    a, b, g = np.radians([e1, e2, e3]); ca, sa, cb, sb, cg, sg = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    return np.array([[cg*cb*ca - sg*sa, cg*cb*sa + sg*ca, -cg*sb], [-sg*cb*ca - cg*sa, -sg*cb*sa + cg*ca, sg*sb], [sb*ca, sb*sa, cb]])
for i, l in enumerate(open('$D/refined_ptcl3D.txt')):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t)
    if float(kv.get('state', 0)) <= 0: continue
    v = [float(kv[k]) for k in keys]
    R0 = mat(M['e1'][i], M['e2'][i], M['e3'][i]); R1 = mat(*v[:3]); dang.append(np.degrees(np.arccos(np.clip((np.trace(R0.T @ R1) - 1) / 2, -1, 1)))); dsh.append(np.hypot(v[3] - M['x'][i], v[4] - M['y'][i]))
    for k, val in zip(keys, v): new[k][i] = val
    n += 1
for k in keys: M[k] = new[k]
np.savez('$D/meta_alt.npz', **M)
print('poses updated for %d particles: angular change median %.2f deg (90th pct %.2f), shift change median %.2f px (90th %.2f)' % (n, np.median(dang), np.percentile(dang, 90), np.median(dsh), np.percentile(dsh, 90)))
PYEOF
  export OMP_NUM_THREADS=1; mkdir -p $D/klem
  $PY $W/klem.py $STK none none $D/klem --meta $D/meta_alt.npz --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 2 --nproc 16 --init ${INIT:-npz} --nsub 0 ${SUBSET:+--subset $SUBSET} --ctf auto --tag "rnasep alt$KTAG r$R" 2>&1 | grep --line-buffered -v "^   [0-9]* / "
  $PY $W/klem_real_readout.py $D/klem $D/meta_alt.npz --labels $LAB --out $D/klem/fig_real_readout.png 2>&1 | grep -v Warning
  KDIR=$D/klem
done
log DONE
