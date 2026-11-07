#!/bin/bash
#SBATCH -J klem_pole_sets
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# pole sets across runs: both tails of each run's top NDIR latent coordinates (COORDS 'run=flex_pca_coordinates.txt ..') start BASE's site J EM at BASE's poses; endpoints grouped across runs; env RT NODE SRC BASE J COORDS NDIR NPOLE EXTRA OUT MSKDIAM MAPN PARS
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; BASE=${BASE:?}; J=${J:-2}; O=${OUT:?}; mkdir -p $O/poles
log(){ echo "[$(date '+%F %T') pole sets] $*"; }
$PY - <<PYEOF || { log "pole starts FAILED"; exit 1; }
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0
for e in '${COORDS:?}'.split():
    run, f = e.split('=', 1); C = np.loadtxt(f, comments='#'); r = C[:, 0].astype(int) - 1; Z = C[:, 5:]; k = use[r]; r, Z = r[k], Z[k]
    for d in range(${NDIR:-3}):
        o = np.argsort(Z[:, d])
        for tag, ix in (('hi', o[::-1][:${NPOLE:-15000}]), ('lo', o[:${NPOLE:-15000}])):
            l = np.zeros(len(use), int); l[r[ix]] = 1; np.savetxt('$O/poles/%s_z%d%s.txt' % (run, d + 1, tag), l, fmt='%d')
    print('%s: %d rows of %s, coordinates 1..${NDIR:-3}, both tails of %d' % (run, len(r), f, ${NPOLE:-15000}))
PYEOF
SEEDS=""; for e in $COORDS; do for d in $(seq 1 ${NDIR:-3}); do for t in hi lo; do SEEDS="$SEEDS ${e%%=*}_z$d$t=$O/poles/${e%%=*}_z$d$t.txt@1"; done; done; done
RT=${RT:?} NODE=${NODE:-r} SRC=${SRC:?} BASE=$BASE J=$J MSKDIAM=${MSKDIAM:-220} OUT=$O SEEDS="$SEEDS" PARS=${PARS:-4} MAPN="${MAPN:-2500}" TOPN=1000,2500,5000,10000,15000 bash $T/site_reseed.sbatch > $O/reseed.log 2>&1
$PY - <<PYEOF | tee $O/pole_sets.txt
import numpy as np, os
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; old = np.loadtxt('$K_/f692r/coordtree_fam1/units.txt', dtype=int); old = ((old == 15) | (old == 16))[use]
ents = [(s.split('=')[0], '$O/%s/post_master.txt' % s.split('=')[0]) for s in '$SEEDS'.split()] + [tuple(e.split('=', 1)) for e in '${EXTRA:-}'.split()]
P = {n: np.loadtxt(f)[use] for n, f in ents if os.path.exists(f)}; nm = list(P); R = {n: np.argsort(np.argsort(-P[n])) for n in nm}
S = np.array([[np.corrcoef(R[a], R[b])[0, 1] for b in nm] for a in nm]); g = list(range(len(nm)))
for i in range(len(nm)):
    for j in range(i + 1, len(nm)):
        if S[i, j] >= 0.99: gi, gj = g[i], g[j]; g = [gi if x == gj else x for x in g]
arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open('$O/arm_top.txt')) if w and not w[0].startswith('#')} if os.path.exists('$O/arm_top.txt') else {}
ns = lambda m: open('$O/%s/nstar.txt' % m).read().strip() if os.path.exists('$O/%s/nstar.txt' % m) else '-'
runs = sorted(set(n.split('_z')[0] for n in nm if '_z' in n))
print('endpoints of %d starts over runs %s, grouped at Spearman >= 0.99 (single linkage); readout = SIMPLE top-2.5k closed / open; N* = hold-out sizing' % (sum('_z' in n for n in nm), ' '.join(runs)))
for G in sorted(set(g), key=lambda G: -sum(x == G for x in g)):
    mem = [nm[i] for i in range(len(nm)) if g[i] == G]; o = np.argsort(-P[mem[0]]); ix = [nm.index(m) for m in mem]; mn = S[np.ix_(ix, ix)].min()
    enr = ' '.join('%d: %.2fx' % (N, old[o[:N]].mean() / old.mean()) for N in (1000, 2500, 10000)); per = ' '.join('%s %d' % (r, sum(m.startswith(r + '_z') for m in mem)) for r in runs)
    rd = ' '.join('%.2f/%.2f' % arm[m + '_top2500'] for m in mem if m + '_top2500' in arm)
    print('group of %d (min Spearman %.3f; per run %s): %s | old-closed %s | N* %s | top-2.5k %s' % (len(mem), mn, per, ' '.join(mem), enr, ' '.join(ns(m) for m in mem), rd or '-'))
print('Spearman matrix (%s):' % ' '.join(nm))
for i, x in enumerate(nm): print('%-16s %s' % (x, ' '.join('%+.2f' % v for v in S[i])))
PYEOF
log DONE
