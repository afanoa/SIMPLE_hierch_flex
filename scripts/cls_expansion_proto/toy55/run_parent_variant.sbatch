#!/bin/bash
#SBATCH -J rnasep_f692var_parent
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the real 692k parent board WITH the per-particle contrast scale (pscale 1), same settings as f692 otherwise, into
# f692ps/parent; readouts: Spearman(z1,z2 ; defocus), junk mass vs defocus, family split, identification vs the 10 and
# the 55 deposited apo maps. Frozen copy of the rebuilt klem2 build.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692${TAG:-ps}; mkdir -p $R; cd $R; log(){ echo "[$(date '+%F %T') f692${TAG:-ps}] $*"; }
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; METAALL=$A/rnasep_coupled/basis_span_all/meta_all.npz
if [ ! -x $R/klem2_bin/bin/simple_exec ]; then mkdir -p $R/klem2_bin && cp -r $HOME/SIMPLE_klem2/build/bin $R/klem2_bin/ && { [ -d $HOME/SIMPLE_klem2/build/scripts ] && cp -r $HOME/SIMPLE_klem2/build/scripts $R/klem2_bin/ || true; }; fi
export FB=$R/klem2_bin NTHR=8 NPARTS=10 PSCALE=${PSCALE:-1} CTFABS=${CTFABS:-0}
LPSP="20 20 20 14 14 14 10 10 10"; SIGP="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7"; cp -f $A/rnasep_klem/f692/subset_all.txt $R/subset_all.txt
log "parent board (pscale $PSCALE ctfabs $CTFABS) on $(grep -c . $R/subset_all.txt) particles"
[ -s $R/parent/klem_latent.txt ] || NFREE=2 META=$METAALL bash $W/fboard_ps.sh $R/parent $SRC $R/subset_all.txt "$LPSP" "$SIGP" || { log "parent board FAILED"; exit 1; }
grep -a "contrast scale\|noise factor" $R/parent/board.log $R/parent/klem_worker_01.log 2>/dev/null | tail -4
$PY - <<PY
import numpy as np
from scipy.stats import spearmanr
M = np.load('$METAALL'); dfa = 0.5 * (M['dfx'] + M['dfy'])
for tag, d in (('f692 (no pscale)', '$A/rnasep_klem/f692/parent'), ('f692$TAG', '$R/parent')):
    L = np.loadtxt(d + '/klem_latent.txt', comments='#'); p = L[:, 0].astype(int); df = dfa[p - 1]; R = np.load(d + '/klem_resp.npz'); r = R['r']; jm = r[:, 49:].sum(1) / np.maximum(r.sum(1), 1e-30)
    print('%-18s Spearman(z1, defocus) %+.3f  (z2 %+.3f) | junk mass: Spearman with defocus %+.3f, share of particles with junk > 0.5: %.3f | mean defocus of junk>0.5 vs rest: %.2f vs %.2f' % (tag, spearmanr(L[:, 1], df)[0], spearmanr(L[:, 2], df)[0], spearmanr(jm, df)[0], (jm > 0.5).mean(), df[jm > 0.5].mean() if (jm > 0.5).any() else 0, df[jm <= 0.5].mean()))
PY
$PY $W/expand_families.py split $R/parent $R/expand --nfam 2 2>&1 | grep -v Warning | tail -4
$PY - <<PY
import numpy as np
M = np.load('$METAALL'); dfa = 0.5 * (M['dfx'] + M['dfy'])
for f in (1, 2):
    s = np.loadtxt('$R/expand/family%d_subset.txt' % f, dtype=int); print('family %d: %d particles, mean defocus %.2f' % (f, len(s), dfa[s - 1].mean()))
PY
echo "== identification of the parent cells vs the 10 deposited (readout)"; $PY $W/deposited_match.py $R/parent/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -3
echo "== identification of the parent cells vs ALL 55 deposited apo"; $PY $W/deposited_match.py $R/parent/cells $W/maps160/docked_all $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -3
log DONE
