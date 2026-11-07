#!/bin/bash
#SBATCH -p norm
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# A node's own pure-quality reference for the blur screen. env TAG, PROJ (the project the judged maps came from), RECSH (a hetero_step.sh
# holding its reconstruct3D line), MASK, CROP (optional box_crop, e.g. 88 for stage-B maps), NP, NT. In PROJ's rows with state > 0: the
# top and bottom corr quintiles inside dfx deciles (matched defocus) become states 1 (H) and 2 (L), reconstructed with RECSH's line
# (nstates=2, postprocess=no[, box_crop]); dB_ref +- SE = quality_screen3.py's blur of L vs H from the unfiltered halves in the stable region
# Corr and poses come from PROJ alone (mixing alignments flips dB_ref); where selection refuses PROJ, set_states_binoris.py sets the states
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
NP=${NP:-4}; NT=${NT:-8}; B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=$NT
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; OUT=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/${TAG:?}
[ -e $OUT ] && mv $OUT ${OUT}_old_$(date +%s); mkdir -p $OUT; cd $OUT; log(){ echo "[$(date "+%F %T") quality_ref $TAG] $*"; }
cp ${PROJ:?} p.simple; simple_exec prg=update_project projfile=p.simple qsys_name=local mkdir=no > up.log 2>&1
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py dump.npz > dump.log 2>&1
# RANK=ptcl2D: rank by the 2D-classification corr (no 3D poses) instead of ptcl3D's; states and reconstruction stay on ptcl3D
[ "${RANK:-}" = ptcl2D ] && simple_exec prg=print_project_field projfile=p.simple oritype=ptcl2D 2>/dev/null | $PY $T/dump_ptcl.py dump2d.npz > dump2d.log 2>&1
RANK=${RANK:-} $PY - <<'EOF' > groups.log 2>&1
import os, numpy as np
D = np.load('dump.npz'); st, c, f = np.nan_to_num(D['state']), D['corr'], D['dfx']
if os.environ['RANK'] == 'ptcl2D': c = np.load('dump2d.npz')['corr']; assert len(c) == len(st), (len(c), len(st))
node = (st > 0) & np.isfinite(c) & np.isfinite(f); out = np.zeros(len(st), int)
dd = np.searchsorted(np.quantile(f[node], np.linspace(0.1, 0.9, 9)), f)
for j in range(10):
    m = np.nonzero(node & (dd == j))[0]; r = np.argsort(np.argsort(c[m])) / max(len(m) - 1, 1); out[m[r >= 0.8]] = 1; out[m[r <= 0.2]] = 2
print('node rows %d: top quintile %d, bottom %d; median corr %.4f / %.4f' % (node.sum(), (out == 1).sum(), (out == 2).sum(), np.median(c[out == 1]), np.median(c[out == 2])))
np.savetxt('states.txt', out, fmt='%d')
EOF
cat groups.log
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=states.txt mkdir=no > sel.log 2>&1 || { log "selection refused p.simple: $(grep -a -m1 Inconst sel.log | cut -c1-80); setting the states in a copy"
  $PY $T/set_states_binoris.py p.simple p_states.simple states.txt dump.npz > setst.log 2>&1 && mv p_states.simple p.simple && log "$(tail -1 setst.log | cut -c1-200)" || { log "set_states FAILED: $(tail -1 setst.log)"; exit 1; }; }
REC=$(grep -a -o "prg=reconstruct3D[^>]*" ${RECSH:?} | sed -e 's/projfile=[^ ]*/projfile=p.simple/' -e 's/nstates=[^ ]*/nstates=2/' -e "s/nparts=[^ ]*/nparts=$NP/" -e "s/nthr=[^ ]*/nthr=$NT/" -e 's/ postprocess=[^ ]*//')
REC="$REC postprocess=no${CROP:+ box_crop=$CROP}"; t0=$(date +%s); simple_exec $REC > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log)"; exit 1; }
log "reconstruct3D $(( $(date +%s) - t0 )) s: $REC"
nH=$(grep -c "^1$" states.txt); nL=$(grep -c "^2$" states.txt)
echo "$TAG ref $OUT/recvol_state01.mrc $OUT/recvol_state01_even.mrc $OUT/recvol_state01_odd.mrc $OUT/recvol_state02.mrc $OUT/recvol_state02_even.mrc $OUT/recvol_state02_odd.mrc ${MASK:?} dump.npz states.txt" > case.txt
echo "QSCREEN $TAG ref | L = B ($nL) vs H ($nH)" | $PY $T/quality_screen3.py case.txt > dBref.txt 2>&1; cat dBref.txt; log DONE
