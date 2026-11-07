#!/bin/bash
#SBATCH -J rnasep_closed_focal_fam2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Focal classification of the closed state (Degenhardt: closed = TL-TLR docking of P10.1 onto P12, ~13 % of the apo
# ensemble, resolved with focal masks on the periphery). Focal mask = where the deposited closed conformer (18) differs
# from the open ones (27/28/48/49), dilated and softened; a klem board whose E-step references differ only inside it
# (consensus outside), at the refined poses, pose marginalisation, balanced; cells scored by the focal correlation with
# closed vs open deposited maps and by arm occupancy; the docked cells reconstructed at full box with NO radial mask.
# env: SUBSET (family1 subset or the closed-family subset), ZTAG, NTOP (3), BAL (0.5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=${PIPE:-$E/fam1_pipeline}; KD=$P/klem_fb; META=$P/meta_fb.npz; METAALL=$A/rnasep_coupled/basis_span_all/meta_all.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
SUBSET=${SUBSET:-$E/family1_subset.txt}; D=$KD/closed_focal${ZTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cf] $*"; }
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
W = '$W'; SMPD = 1.464; g = lambda v: gaussian_filter(v, 2.0)
C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in (18, 27, 28, 48, 49)}; cons = read_mrc(W + '/maps160/toy_consensus.mrc')
Lc = g(C[18]); Lo = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0)); diff = Lc - Lo; dens = Lc > Lc.mean() + 1.5 * Lc.std()
core = (diff > diff.mean() + 2.5 * diff.std()) & dens; m = binary_dilation(core, iterations=4); soft = gaussian_filter(m.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
n = 160; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD
write_mrc('$D/focal_mask_p12.mrc', soft.astype(np.float32), SMPD); print('focal mask: %d core voxels (closed - open > 2.5 sd inside the closed density), %d after dilation; radius of the core median %.0f A, max %.0f A' % (core.sum(), m.sum(), np.median(r[core]), r[core].max()))
PYEOF
export OMP_NUM_THREADS=1; log "focal board at the refined poses on $(grep -c . $SUBSET) particles"
$PY $W/klem.py $STK none none $D/board --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $SUBSET --ctf auto \
    --lps "12 12 10 10 8 8 8 8 8" --sigmas "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" --focus_mask $D/focal_mask_p12.mrc --pose_marg 3 --devw 0.05 --pnoise --balance ${BAL:-0.5} --tag "rnasep closed focal${ZTAG:-}" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "subset\|focus mask\|final:\|states: cross-half\|Error\|Traceback" | head -6
echo "== cells: focal correlation closed vs open, arm occupancy"; $PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464
fs = sorted(glob.glob('$D/board/cells/cell_*.mrc')); V = [read_mrc(f) for f in fs]; nc = V[0].shape[0]; n = 160
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in (18, 27, 28, 48, 49)}; mk = read_mrc('$D/focal_mask_p12.mrc')
def crop(v):
    Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); c = n // 2; h = nc // 2; sub = Fv[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(sub)))) * (nc / n) ** 3
mkc = np.clip(crop(mk), 0, 1); closed = crop(C[18] - cons); openm = crop(np.mean([C[i] for i in (27, 28, 48, 49)], 0) - cons); Vm = np.mean(V, 0)
def fc(x, y): return float(np.sum(mkc * x * y) / max(np.sqrt(np.sum(mkc * x * x) * np.sum(mkc * y * y)), 1e-30))
sc = []
for k, v in enumerate(V):
    d = v - Vm; sc.append((fc(d, closed) - fc(d, openm), fc(d, closed), fc(d, openm)))
sc = np.array(sc); order = np.argsort(-sc[:, 0]); np.savetxt('$D/board/focal_scores.txt', np.c_[np.arange(1, len(V) + 1), sc], fmt='%d %.4f %.4f %.4f', header='cell closed_minus_open corr_closed corr_open (inside the focal mask, deviation from the mean cell)')
print('focal score (corr with closed - corr with open inside the P12/P10.1 mask): median %+.2f, max %+.2f, min %+.2f' % (np.median(sc[:, 0]), sc[:, 0].max(), sc[:, 0].min()))
for k in order[:8]: print('  %s: closed-open %+.2f (closed %+.2f, open %+.2f)' % (fs[k].split('/')[-1], *sc[k]))
L = np.loadtxt('$D/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]; N = len(np.load('$METAALL')['flags']); st = np.zeros(N, int)
top = [int(k) + 1 for k in order[:${NTOP:-3}]]; bot = [int(k) + 1 for k in order[-${NTOP:-3}:]]
for c in top: st[pind[(gm > 0.5) & (am == c)] - 1] = 1
for c in bot: st[pind[(gm > 0.5) & (am == c)] - 1] = 2
np.savetxt('$D/states.txt', st, fmt='%d'); print('state 1 = top %d docked cells %s: %d particles; state 2 = bottom %d cells %s: %d particles' % (len(top), top, int((st == 1).sum()), len(bot), bot, int((st == 2).sum())))
PYEOF
export OMP_NUM_THREADS=8; PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=2 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in 01 02; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY $W/sharpen_matched.py $D/maps $D/for_download --rmask 0 --bfac -150 2>/dev/null | tail -3 || true
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; n = 160; g = lambda v: gaussian_filter(v, 2.0); mk = read_mrc('$D/focal_mask_p12.mrc'); cons = read_mrc(W + '/maps160/toy_consensus.mrc')
C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in (18, 27, 28, 48, 49)}; closed = g(C[18] - cons); openm = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0) - cons)
def fc(x, y): return float(np.sum(mk * x * y) / max(np.sqrt(np.sum(mk * x * x) * np.sum(mk * y * y)), 1e-30))
zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
for w, name in ((1, 'docked cells'), (2, 'undocked cells')):
    v = read_mrc('$D/maps/recvol_state%02d.mrc' % w); d = g(v - cons); print('FOCAL MAP %d (%s): corr with closed-18 inside the P12/P10.1 mask %+.2f, with open %+.2f | mask-region density relative to the closed map %.2f' % (w, name, fc(d, closed), fc(d, openm), float(np.sum(mk * lp(v, 8)) / max(np.sum(mk * lp(C[18], 8)), 1e-30) * (np.sum(mk * lp(C[18], 8)) / max(np.sum(mk * lp(cons, 8)), 1e-30)) / (np.sum(mk * lp(v, 8)) / max(np.sum(mk * lp(cons, 8)), 1e-30)) if False else np.sum(mk * lp(v, 8)) / max(np.sum(mk * lp(C[18], 8)), 1e-30))))
    bsharp = v.copy(); write_mrc('$D/for_download/focal_state%d_%s_raw.mrc' % (w, name.split()[0]), bsharp.astype(np.float32), SMPD)
PYEOF
log DONE
