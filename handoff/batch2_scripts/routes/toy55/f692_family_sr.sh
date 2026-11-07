#!/bin/bash
# f692 family stage with the Fortran board: child board on the family's particles at the consensus poses (focus = the
# parent node's own moving region), k-means K on its posterior means -> refine3D_states (K states, local 20->8 A) ->
# refined poses -> meta_fb.npz -> feedback board at the refined poses (focus = the child's own moving region) ->
# identification of the feedback cells against the deposited conformers (readout only).
# env: R (run root), FAM, SRC (project with the consensus poses), META (meta_all.npz), K (10), NP (10), MAXITS (6), LPSC, SIGC
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
E=$R/expand; D=$R/fam$FAM; mkdir -p $D; SUB=$E/family${FAM}_subset.txt; K=${K:-10}; NP=${NP:-10}; log(){ echo "[$(date '+%F %T') fam$FAM] $*"; }
LPSC=${LPSC:-"14 14 12 12 10 10 8 8 8"}; SIGC=${SIGC:-"1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"}; NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))")
# ---- 1 child board at the consensus poses
if [ ! -s $D/child/klem_latent.txt ]; then
  FM=""; [ -s $R/regions_parent/region_1.mrc ] && FM=$R/regions_parent/region_1.mrc
  FOCUSMASK=$FM META=$META bash $W/fboard.sh $D/child $SRC $SUB "$LPSC" "$SIGC" || exit 1
fi
echo "== child board of family $FAM vs the deposited conformers (readout)"; $PY $W/deposited_match.py $D/child/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -12
# ---- 2 k-means K on the child board's posterior means -> states.txt (one row per project particle)
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$D/child/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z = L[:, 1:3]; gm = L[:, 5]; ok = gm > 0.5; rng = np.random.default_rng(0); X = z[ok]; K = $K; best = None
for seed in range(10):
    C = X[rng.choice(len(X), K, replace=False)]
    for it in range(100):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2); a = d.argmin(1); Cn = np.array([X[a == k].mean(0) if (a == k).any() else C[k] for k in range(K)])
        if np.allclose(Cn, C): break
        C = Cn
    cost = d[np.arange(len(X)), a].sum()
    if best is None or cost < best[0]: best = (cost, a.copy())
st = np.zeros($NROWS, int); st[pind[ok] - 1] = best[1] + 1; np.savetxt('$D/states.txt', st, fmt='%d'); print('family $FAM: k-means K=%d on the child board posterior means: sizes %s' % (K, np.bincount(best[1])))
PYEOF
# ---- 3 refine3D_states: K states, local poses, 20 -> 8 A (SKIPREFINE=1: the consensus poses stand in as the "refined" project)
if [ "${SKIPREFINE:-0}" = 1 ] && [ ! -s $D/refined_ptcl3D.txt ]; then
  mkdir -p $D/refine/1_refine3D_states && cp -f $SRC $D/refine/1_refine3D_states/p.simple && cd $D/refine/1_refine3D_states
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt; mkdir -p $D/maps; cd $D; log "SKIPREFINE: consensus poses used as the refined project"
fi
if [ ! -s $D/refined_ptcl3D.txt ]; then
  mkdir -p $D/refine && cd $D/refine && cp -f $SRC p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  export OMP_NUM_THREADS=8; log "refine3D_states nstates=$K local 20->8 A maxits ${MAXITS:-6} nparts $NP"
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
  RD=$(ls -td $D/refine/*_refine3D_states | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 $K); do cp -f $RD/recvol_state$w.mrc $D/maps/ 2>/dev/null; done
  simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt; cd $D
fi
RD=$(ls -td $D/refine/*_refine3D_states | head -1); PS=$RD/p.simple; log "refined project $PS ($(grep -c . $D/refined_ptcl3D.txt) rows)"
if [ "${SKIPREFINE:-0}" != 1 ]; then echo "== refined state maps of family $FAM vs the deposited conformers (readout)"; mkdir -p $D/states_dir; for w in $(seq -f %02g 1 $K); do cp -f $D/maps/recvol_state$w.mrc $D/states_dir/state_$w.mrc 2>/dev/null; done
$PY $W/deposited_match.py $D/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -12; fi
# ---- 4 refined poses into the meta (meta_fb.npz: the init latent and CTF stay)
$PY - <<PYEOF
import numpy as np
M = dict(np.load('$META')); keys = ('e1', 'e2', 'e3', 'x', 'y'); new = {k: M[k].astype(float).copy() for k in keys}; n = 0
for i, l in enumerate(open('$D/refined_ptcl3D.txt')):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t)
    if float(kv.get('state', 0)) <= 0: continue
    for k in keys: new[k][i] = float(kv[k])
    n += 1
for k in keys: M[k] = new[k]
np.savez('$D/meta_fb.npz', **M); print('refined poses for %d particles -> meta_fb.npz' % n)
PYEOF
# ---- 5 feedback board at the refined poses, focus = the child board's own moving region
[ -s $D/regions_child/region_1.mrc ] || $PY $W/derive_regions.py $D/child $D/regions_child --maxreg 2 2>&1 | grep -v Warning | tail -4
if [ ! -s $D/klem_fb/klem_latent.txt ]; then
  FM=""; [ -s $D/regions_child/region_1.mrc ] && FM=$D/regions_child/region_1.mrc
  FOCUSMASK=$FM META=$D/meta_fb.npz bash $W/fboard.sh $D/klem_fb $PS $SUB "$LPSC" "$SIGC" || exit 1
fi
echo "== feedback board of family $FAM vs the deposited conformers (readout)"; $PY $W/deposited_match.py $D/klem_fb/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $D/klem_fb/cell_matches.txt 2>&1 | grep -v Warning | tail -12
log DONE
