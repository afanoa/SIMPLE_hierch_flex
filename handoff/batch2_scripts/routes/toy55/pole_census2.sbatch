#!/bin/bash
#SBATCH -J klem_pole_census2
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=220G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# cross-run check of the r1 census's middle fixed point ("both arms weak", reached by port1's middle 15k of z1-z3): the same middle-band
# starts from the script's and port2's own fflex latents, the same site EM (site_02, shared input poses); grouped with census + pole sets
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; P=$K/stree_port
BASE=$K/factorial/f692st_r1_inposes; C1=$K/factorial/f692st_r1_census; O=$K/factorial/f692st_r1_census2; mkdir -p $O/starts
log(){ echo "[$(date '+%F %T') pole census2] $*"; }
$PY - <<PYEOF || { log "starts FAILED"; exit 1; }
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0
for run, f in (('s', '$K/rflex_f692st/node_r1/focus/1_flex_pca/flex_pca_coordinates.txt'), ('p2', '$P/tree_f692r1_seed2/1_hetero_tree/node_r/fflex/flex_pca_coordinates.txt')):
    C = np.loadtxt(f, comments='#'); r = C[:, 0].astype(int) - 1; Z = C[:, 5:]; k = use[r]; r, Z = r[k], Z[k]
    for d in range(3):
        o = np.argsort(Z[:, d]); m = len(o) // 2; l = np.zeros(len(use), int); l[r[o[m - 7500:m + 7500]]] = 1
        np.savetxt('$O/starts/%sz%dmid.txt' % (run, d + 1), l, fmt='%d')
PYEOF
SEEDS=""; for f in $O/starts/*.txt; do n=$(basename $f .txt); SEEDS="$SEEDS c$n=$f@1"; done; log "starts:$SEEDS"
RT=f692st NODE=r1 SRC=$K/rflex_f692st/node_r1/n.simple BASE=$BASE J=2 MSKDIAM=232 OUT=$O SEEDS="$SEEDS" PARS=${PARS:-4} MAPN=2500 \
  TOPN=1000,2500,5000,10000,15000 bash $T/site_reseed.sbatch > $O/reseed.log 2>&1
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for d in $O/rec_*_top2500; do m=$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null | head -1); n=$(basename $d); n=${n#rec_}
  [ -n "$m" ] && echo "$n $($PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $m 2>&1 | grep -v Warning | tail -1 | awk '{print $2, $3, $4, $5}')"; done | tee $O/arm_readouts.txt
$PY - <<PYEOF | tee $O/census2.txt
import numpy as np, os, glob
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0
ents = sorted((os.path.basename(os.path.dirname(f)), f) for f in glob.glob('$O/c*/post_master.txt')) + \
       sorted((os.path.basename(os.path.dirname(f)), f) for f in glob.glob('$C1/c*/post_master.txt')) + \
       sorted((os.path.basename(os.path.dirname(f)), f) for f in glob.glob('$K/factorial/f692st_r1_polesets/*_z*/post_master.txt'))
Pm = {n: np.loadtxt(f)[use] for n, f in ents}; nm = list(Pm); R = {n: np.argsort(np.argsort(-Pm[n])) for n in nm}
mid = [n for n in nm if n.endswith('mid')]
print('Spearman among the middle-band endpoints (port1 cp1*, script cs*, port2 cp2*):')
for a in mid: print('  %-10s %s' % (a, ' '.join('%s %+.3f' % (b, np.corrcoef(R[a], R[b])[0, 1]) for b in mid if b != a)))
print('each middle endpoint vs its best non-middle endpoint:')
for a in mid:
    o = [n for n in nm if not n.endswith('mid')]; b = max(o, key=lambda x: np.corrcoef(R[a], R[x])[0, 1]); print('  %-10s -> %-12s %+.3f' % (a, b, np.corrcoef(R[a], R[b])[0, 1]))
PYEOF
log DONE
