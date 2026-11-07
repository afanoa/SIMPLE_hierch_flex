#!/bin/bash
#SBATCH -J klem_fastcomp
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# composition stage test: the families of a fast run (spectral clusters of its whole-particle board's cells, expand/family*_subset.txt)
# seed a K-state refine3D_states (local poses: it takes local|global only, low resolution, local parts on this node) that re-classifies every particle against the
# state maps; purity of the seeds and of the re-classified states against the truth. env TAG, RT (ff2_<TAG>), LPSTART (20), LPSTOP (12), MAXITS (6)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; TAG=${TAG:?}; R=$K/${RT:-ff2_$TAG}
O=$R/comp${KSEED:+_k$KSEED}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') fastcomp $TAG] $*"; }; NF=$(cat $R/expand/nfam.txt); GT=$T/sim$TAG/gt_particles.txt
if [ -n "${KSEED:-}" ]; then NF=$KSEED; $PY - <<PYEOF   # KSEED: over-cluster the whole board's cells into KSEED spectral clusters as the seeds
import glob, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23])).reshape(n, n, n).astype(float)
fs = sorted(glob.glob('$R/whole/board/cells/cell_*.mrc')); V = np.array([rd(f) for f in fs]); D = V - V.mean(0); X = np.array([np.fft.fftn(d).ravel() for d in D]); X /= np.linalg.norm(X, axis=1)[:, None]
S = np.maximum(np.real(X @ X.conj().T), 0); np.fill_diagonal(S, 0); dg = S.sum(1) + 1e-9; ev, U = np.linalg.eigh(np.eye(len(S)) - S / np.sqrt(np.outer(dg, dg))); k = $KSEED
E = U[:, :k] / np.maximum(np.linalg.norm(U[:, :k], axis=1)[:, None], 1e-12); rng = np.random.default_rng(0); best = None
for rep in range(20):
    C = E[rng.choice(len(E), k, replace=False)]
    for it in range(100): lab = np.argmin(((E[:, None] - C[None]) ** 2).sum(-1), 1); C = np.array([E[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
    sse = ((E - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
lab = best[1]; L_ = np.loadtxt('$R/whole/board/klem_latent.txt', comments='#'); pind, cell = L_[:, 0].astype(int), L_[:, 3].astype(int); ncell = len(fs)
n = sum(1 for l in open('$T/lat_$TAG/ptcl3D_r0.txt') if l.strip()); out = np.zeros(n, int); on = (cell >= 1) & (cell <= ncell); out[pind[on] - 1] = lab[cell[on] - 1] + 1
np.savetxt('$O/seed_states.txt', out, fmt='%d'); print('KSEED %d seed labels: %s' % (k, ' '.join(str(int((out == j).sum())) for j in range(1, k + 1))))
PYEOF
else
$PY - <<PYEOF
import numpy as np
n = sum(1 for l in open('$T/lat_$TAG/ptcl3D_r0.txt') if l.strip()); lab = np.zeros(n, int)
for F in range(1, $NF + 1): lab[np.loadtxt('$R/expand/family%d_subset.txt' % F, dtype=int).ravel() - 1] = F
np.savetxt('$O/seed_states.txt', lab, fmt='%d'); print('seed labels: %s' % ' '.join(str(int((lab == F).sum())) for F in range(1, $NF + 1)))
PYEOF
fi
purity(){ awk -v tag="$1" 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; cls = (k == 0) ? "junk" : (k >= 400) ? "grid" : (k >= 300) ? "breath" : (k >= 200) ? "swing" : (k >= 101) ? "holo" : "apo"; n[s]++; m[s, cls]++; N++}
  END {printf "%s:", tag; tot = 0; for (s in n) {best = 0; line = ""; for (q in m) {split(q, p, SUBSEP); if (p[1] == s) {line = line sprintf(" %s %.0f%%", p[2], 100 * m[q] / n[s]); if (m[q] > best) best = m[q]}}; tot += best; printf " | state %s (%d):%s", s, n[s], line}; printf " || purity %.2f\n", tot / N}' $GT $2; }
purity seeds $O/seed_states.txt
cp -f $T/proj$TAG/proj$TAG.simple p.simple
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
log "refine3D_states: $NF states seeded by the board families, local poses, ${LPSTART:-20} -> ${LPSTOP:-12} A, maxits ${MAXITS:-6}"
simple_exec prg=refine3D_states projfile=p.simple nstates=$NF pose_policy=local pgrp=c1 mskdiam=220 lpstart=${LPSTART:-20} lpstop=${LPSTOP:-12} maxits=${MAXITS:-6} ml_reg=no nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR\|error' ref.log | cut -c1-150)"; exit 1; }
RD=$(ls -td $O/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/refined_states.txt
purity refined $O/refined_states.txt
log DONE
