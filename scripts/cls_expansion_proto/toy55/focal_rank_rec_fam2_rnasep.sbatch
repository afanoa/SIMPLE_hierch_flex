#!/bin/bash
#SBATCH -J rnasep_focal_rank_rec_fam2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# From a focal P12/P10.1 board: (1) arm / tip / base occupancy and FSC of its docked and undocked cell maps; (2) every
# particle ranked by its responsibility-weighted focal score (sum_m r_im x (corr closed - corr open)_m), the top N1 and
# the next N2 reconstructed at the refined poses (no sharpening), scored the same way. env: FD (focal dir), N1 (15000), N2 (10000)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=${PIPE:-$A/rnasep_klem/expand1/fam1_pipeline}; FD=${FD:-$P/klem_fb/closed_focal_closedfam}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
D=$FD/rank; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') fr] $*"; }
cat > $D/score.py <<'PYEOF'
import numpy as np, glob, sys
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; g = lambda v: gaussian_filter(v, 2.5)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD; tip = arm & (r > 95); base = arm & (r < 80)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
def score(path, name):
    v = read_mrc(path); Lv = g(v); c = Lv[core].mean(); occ = (Lv[arm].mean() / c) / ref; ot = (Lv[tip].mean() / c) / (L18[tip].mean() / L18[core].mean()); ob = (Lv[base].mean() / c) / (L18[base].mean() / L18[core].mean())
    try: res = fsc_res(read_mrc(path.replace('.mrc', '_even.mrc')), read_mrc(path.replace('.mrc', '_odd.mrc')))
    except Exception: res = float('nan')
    print('  %-34s: arm %.2f | tip %.2f | base %.2f | FSC %.1f A' % (name, occ, ot, ob, res))
for p_, nm in zip(sys.argv[1::2], sys.argv[2::2]): score(p_, nm)
PYEOF
echo "== (1) the focal board's own cell maps"; $PY $D/score.py $FD/maps/recvol_state01.mrc "docked cells (top 3)" $FD/maps/recvol_state02.mrc "undocked cells (bottom 3)" $P/klem_fb/matched_rec/maps/recvol_state03.mrc "all 18-matched (40k, reference)"
echo "== (2) particles ranked by the responsibility-weighted focal score"; $PY - <<PYEOF
import numpy as np
F = np.loadtxt('$FD/board/focal_scores.txt', comments='#'); fs = F[:, 1]; R = np.load('$FD/board/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1'])
s = r[:, :ngrid] @ fs[:ngrid]; order = np.argsort(-s); N = len(np.load('$META')['flags']); st = np.zeros(N, int); N1 = ${N1:-15000}; N2 = ${N2:-10000}
st[pinds[order[:N1]] - 1] = 1; st[pinds[order[N1:N1 + N2]] - 1] = 2; st[pinds[order[-N1:]] - 1] = 3
np.savetxt('$D/states.txt', st, fmt='%d'); print('per-particle focal score: median %+.3f, top %d >= %+.3f, next %d >= %+.3f, bottom %d <= %+.3f' % (np.median(s), N1, s[order[N1 - 1]], N2, s[order[N1 + N2 - 1]], N1, s[order[-N1]]))
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=3 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in 01 02 03; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY $D/score.py $D/maps/recvol_state01.mrc "top ${N1:-15000} by focal score" $D/maps/recvol_state02.mrc "next ${N2:-10000}" $D/maps/recvol_state03.mrc "bottom ${N1:-15000} (undocked)"
$PY $W/sharpen_matched.py $D/maps $D/for_download --rmask 0 --bfac -150 > /dev/null 2>&1 || true
VIS=$P/klem_fb/closed_focal/lever/vis.py; [ -s $VIS ] && $PY $VIS $D/maps/recvol_state01.mrc "top ${N1:-15000} by focal score" $D/maps/recvol_state02.mrc "next ${N2:-10000}" $D/maps/recvol_state03.mrc "bottom ${N1:-15000}" $P/klem_fb/closed_focal/lever/maps/recvol_state01.mrc "family-1-wide docked top 30k (reference)"
log DONE
