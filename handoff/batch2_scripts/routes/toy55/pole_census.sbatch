#!/bin/bash
#SBATCH -J klem_pole_census
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=220G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# fixed-point census of node r1's two-class site EM (site_02 = 1f1, shared input poses, as the pole sets): where does it settle from
# starts that are NOT the z1/z3 tails? Starts: both tails of port1's fflex z4..z6, the middle 15k of z1..z3 (interior), three random
# 15k subsets, the 2.5k tails of z1 (start size). Endpoints grouped (Spearman >= 0.99 / 0.97, single linkage) together with the 18
# pole-set endpoints (script/port1/port2 z1..z3); a group = one fixed point. env OUT PARS
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; P=$K/stree_port
BASE=$K/factorial/f692st_r1_inposes; PS=$K/factorial/f692st_r1_polesets; O=${OUT:-$K/factorial/f692st_r1_census}; mkdir -p $O/starts
log(){ echo "[$(date '+%F %T') pole census] $*"; }
$PY - <<PYEOF || { log "starts FAILED"; exit 1; }
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0
C = np.loadtxt('$P/tree_f692r1/1_hetero_tree/node_r/fflex/flex_pca_coordinates.txt', comments='#'); r = C[:, 0].astype(int) - 1; Z = C[:, 5:]; k = use[r]; r, Z = r[k], Z[k]
def save(nm, ix): l = np.zeros(len(use), int); l[r[ix]] = 1; np.savetxt('$O/starts/%s.txt' % nm, l, fmt='%d')
for d in range(3, min(6, Z.shape[1])):
    o = np.argsort(Z[:, d]); save('p1z%dhi' % (d + 1), o[::-1][:15000]); save('p1z%dlo' % (d + 1), o[:15000])
for d in range(3):
    o = np.argsort(Z[:, d]); m = len(o) // 2; save('p1z%dmid' % (d + 1), o[m - 7500:m + 7500])
for s in (101, 102, 103): save('rand%d' % s, np.random.default_rng(s).choice(len(r), 15000, replace=False))
o = np.argsort(Z[:, 0]); save('p1z1hi2k', o[::-1][:2500]); save('p1z1lo2k', o[:2500])
print('port1 fflex latent: %d dims over %d rows' % (Z.shape[1], len(r)))
PYEOF
SEEDS=""; for f in $O/starts/*.txt; do n=$(basename $f .txt); SEEDS="$SEEDS c$n=$f@1"; done; log "starts:$SEEDS"
RT=f692st NODE=r1 SRC=$K/rflex_f692st/node_r1/n.simple BASE=$BASE J=2 MSKDIAM=232 OUT=$O SEEDS="$SEEDS" PARS=${PARS:-4} MAPN=2500 \
  TOPN=1000,2500,5000,10000,15000 bash $T/site_reseed.sbatch > $O/reseed.log 2>&1
$PY - <<PYEOF | tee $O/census.txt
import numpy as np, os, glob
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; old = np.loadtxt('$K/f692r/coordtree_fam1/units.txt', dtype=int); old = ((old == 15) | (old == 16))[use]
ents = sorted((os.path.basename(os.path.dirname(f)), f) for f in glob.glob('$PS/*_z*/post_master.txt')) + \
       sorted((os.path.basename(os.path.dirname(f)), f) for f in glob.glob('$O/c*/post_master.txt'))
Pm = {n: np.loadtxt(f)[use] for n, f in ents}; nm = list(Pm); R = {n: np.argsort(np.argsort(-Pm[n])) for n in nm}
S = np.array([[np.corrcoef(R[a], R[b])[0, 1] for b in nm] for a in nm])
arm = {}
for f in ('$PS/arm_top.txt', '$O/arm_top.txt'):
    if os.path.exists(f): arm.update({os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open(f)) if w and not w[0].startswith('#')})
def nstar(m):
    for d in ('$PS', '$O'):
        f = '%s/%s/nstar.txt' % (d, m)
        if os.path.exists(f): return open(f).read().strip()
    return '-'
for thr in (0.99, 0.97):
    g = list(range(len(nm)))
    for i in range(len(nm)):
        for j in range(i + 1, len(nm)):
            if S[i, j] >= thr: gi, gj = g[i], g[j]; g = [gi if x == gj else x for x in g]
    print('==== fixed points at Spearman >= %.2f (single linkage): %d groups over %d endpoints' % (thr, len(set(g)), len(nm)))
    for G in sorted(set(g), key=lambda G: -sum(x == G for x in g)):
        mem = [nm[i] for i in range(len(nm)) if g[i] == G]; ix = [nm.index(m) for m in mem]; o = np.argsort(-Pm[mem[0]])
        rd = ' '.join('%s %.2f/%.2f' % (m, *arm[m + '_top2500']) for m in mem if m + '_top2500' in arm)
        print('  group of %d (min Spearman %.3f): %s | N* %s | old-closed top-2.5k %.2fx | top-2.5k closed/open: %s' % (len(mem),
              S[np.ix_(ix, ix)].min(), ' '.join(mem), ' '.join(nstar(m) for m in mem), old[o[:2500]].mean() / old.mean(), rd or '-'))
print('==== each census start: its best-matching earlier pole-set endpoint (Spearman)')
pe = [n for n in nm if not n.startswith('c')]
for n in nm:
    if n.startswith('c'):
        j = max(pe, key=lambda e: S[nm.index(n), nm.index(e)]); print('  %-12s -> %-12s %+.3f   (N* %s)' % (n, j, S[nm.index(n), nm.index(j)], nstar(n)))
PYEOF
log DONE
