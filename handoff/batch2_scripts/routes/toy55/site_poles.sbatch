#!/bin/bash
#SBATCH -J klem_site_poles
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# site EM from both tails (NPOLE rows) of the top NDIR latent coordinates of FLEXDIR, endpoints grouped at Spearman >= 0.99; env RT NODE SRC BASE J FLEXDIR NDIR NPOLE EXTRA OUT MSKDIAM MAPN PARS
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; BASE=${BASE:?}; J=${J:-2}; O=${OUT:?}; mkdir -p $O/poles
log(){ echo "[$(date '+%F %T') site poles] $*"; }
$PY - <<PYEOF || { log "pole starts FAILED"; exit 1; }
import numpy as np, glob, os
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; fd = sorted(glob.glob('${FLEXDIR:?}/*_flex_pca'), key=os.path.getmtime)[-1]
C = np.loadtxt(fd + '/flex_pca_coordinates.txt', comments='#'); r = C[:, 0].astype(int) - 1; Z = C[:, 5:]; k = use[r]; r, Z = r[k], Z[k]
for d in range(${NDIR:-3}):
    o = np.argsort(Z[:, d])
    for tag, ix in (('hi', o[::-1][:${NPOLE:-15000}]), ('lo', o[:${NPOLE:-15000}])):
        l = np.zeros(len(use), int); l[r[ix]] = 1; np.savetxt('$O/poles/z%d%s.txt' % (d + 1, tag), l, fmt='%d')
print('pole starts from %s: %d rows, latent coordinates 1..${NDIR:-3}, both tails, %d rows each' % (fd, len(r), ${NPOLE:-15000}))
PYEOF
SEEDS=""; for d in $(seq 1 ${NDIR:-3}); do for t in hi lo; do SEEDS="$SEEDS z$d$t=$O/poles/z$d$t.txt@1"; done; done
RT=${RT:?} NODE=${NODE:-r} SRC=${SRC:?} BASE=$BASE J=$J MSKDIAM=${MSKDIAM:-220} OUT=$O SEEDS="$SEEDS" PARS=${PARS:-2} MAPN="${MAPN:-2500 10000}" TOPN=1000,2500,5000,10000,15000 bash $T/site_reseed.sbatch
$PY - <<PYEOF | tee $O/groups.txt
import numpy as np, os
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; old = np.loadtxt('$K_/f692r/coordtree_fam1/units.txt', dtype=int); old = ((old == 15) | (old == 16))[use]
ents = [('base_r1s5', '$BASE/site_$J/post_master.txt')] + [(s.split('=')[0], '$O/%s/post_master.txt' % s.split('=')[0]) for s in '$SEEDS'.split()] + [tuple(e.split('=', 1)) for e in '${EXTRA:-}'.split()]
P = {n: np.loadtxt(f)[use] for n, f in ents if os.path.exists(f)}; nm = list(P); R = {n: np.argsort(np.argsort(-P[n])) for n in nm}
S = np.array([[np.corrcoef(R[a], R[b])[0, 1] for b in nm] for a in nm]); g = list(range(len(nm)))
for i in range(len(nm)):
    for j in range(i + 1, len(nm)):
        if S[i, j] >= 0.99: gi, gj = g[i], g[j]; g = [gi if x == gj else x for x in g]
arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open('$O/arm_top.txt')) if w and not w[0].startswith('#')} if os.path.exists('$O/arm_top.txt') else {}
print('endpoints grouped at Spearman >= 0.99 (single linkage); old f692r closed enrichment of the top-N; SIMPLE top-N readouts closed / open')
for G in sorted(set(g)):
    mem = [nm[i] for i in range(len(nm)) if g[i] == G]; p = P[mem[0]]; o = np.argsort(-p)
    enr = ' '.join('%d: %.2fx' % (N, old[o[:N]].mean() / old.mean()) for N in (1000, 2500, 10000))
    rd = '; '.join('%s %s' % (m, ' '.join('top%s %.2f / %.2f' % (k.split('_top')[1], *arm[k]) for k in sorted(arm, key=lambda k: int(k.split('_top')[1])) if k.split('_top')[0] == m)) for m in mem if any(k.split('_top')[0] == m for k in arm))
    print('group of %d: %s | enrichment %s | %s' % (len(mem), ' '.join(mem), enr, rd or '-'))
print('Spearman matrix (%s):' % ' '.join(nm))
for i, a in enumerate(nm): print('%-12s %s' % (a, ' '.join('%+.3f' % x for x in S[i])))
PYEOF
log "poles DONE"
