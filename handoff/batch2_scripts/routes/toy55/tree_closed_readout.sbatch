#!/bin/bash
#SBATCH -J klem_tree_readout
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Closed-state readout of EVERY finished node state of a real-data rflex tree (read-only on the tree; outputs in OUT): stage-B state maps
# (lp 10, box 160, view copies in VIEW) -> closed-arm / open-arm occupancy (arm_readout.py) -> per state: particles, the r22 closed-SITE
# posterior (site_r22 soft EM: mean over the state's r22 members, members among the site route's top 2.5k / 5k / 10k) and the old f692r
# closed states 15+16 (count, share of the state, share of all of them) -> radial extent (radial_extent.py) of the 6 most closed states.
# Also every state's level / node and whether it became a leaf. env RT (f692st), OUT (rnasep_klem/factorial/<RT>_readout)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:-f692st}; R=$K_/rflex_$RT
O=${OUT:-$K_/factorial/${RT}_readout}; V=$K_/view_rflex/factorial_${RT}_states; mkdir -p $O $V; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
: > $O/spec.txt; : > $O/states.txt
for D in $R/node_*; do u=${D##*/node_}; [ -e $D/node.done ] || continue; RD=$(ls -td $D/fstageB/*_refine3D 2>/dev/null | head -1); [ -n "$RD" ] || continue
  LAB=$D/refined_f.txt; [ -s $D/refined_s.txt ] && { LAB=$D/refined_s.txt; RD=$(ls -td $D/softrec/*_reconstruct3D | head -1); }
  for s in $(awk '$1 > 0' $LAB | sort -n | uniq); do b=$RD/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] || continue
    nm=${u}_s$(printf %02d $s); echo "$nm ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; echo "$nm $u $s $LAB" >> $O/states.txt; done; done
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $CONS $CL $O/arm.txt $(awk '{print "'$V'/" $1 ".mrc"}' $O/spec.txt) 2>&1 | grep -v Warning | tail -1
$PY - <<PYEOF
import numpy as np, os
H = '$T/simrnasep692'; n = 932496; p = np.full(n, np.nan)
for h in 'AB':
    rows = np.loadtxt(H + '_%s/meta_rows.txt' % h).astype(int) - 1; s = np.loadtxt('$K_/site_r22/start_%s.txt' % h, dtype=int); m = s > 0
    p[rows[m]] = np.loadtxt('$K_/site_r22/em/post_soft_%s.txt' % h)[m, 1]
top = {N: np.loadtxt('$K_/site_r22/top%d.txt' % N, dtype=int) > 0 for N in (2500, 5000, 10000)}
old = np.loadtxt('$K_/f692r/coordtree_fam1/units.txt', dtype=int); cl = (old == 15) | (old == 16)
arm = {os.path.basename(l.split()[0]).replace('.mrc', ''): (float(l.split()[1]), float(l.split()[2])) for l in open('$O/arm.txt') if not l.startswith('#')}
leaf = {}   # node state -> leaf / junk (from every finished node's decision)
for l in open('$O/states.txt'):
    u = l.split()[1]; g = '$R/node_%s/groups.txt' % u
    if os.path.exists(g):
        for w in (x.split() for x in open(g) if x.strip()):
            if w[0] in ('leaves', 'junk'):
                for s in w[2].split('+'): leaf['%s_s%02d' % (u, int(s))] = w[0] if '+' not in w[2] else w[0] + '(group)'
rows = []
for l in open('$O/states.txt'):
    nm, u, s, lab = l.split(); L = np.loadtxt(lab, dtype=int) == int(s); r22 = L & np.isfinite(p)
    rows.append((nm, int(L.sum()), arm.get(nm, (np.nan, np.nan)), np.nanmean(p[r22]) if r22.any() else np.nan, int(r22.sum()), [int((L & top[N]).sum()) for N in (2500, 5000, 10000)], int((L & cl).sum()), leaf.get(nm, '')))
rows.sort(key=lambda t: -t[2][0]); tc = cl.sum()
print('%-12s %7s %6s %6s | %8s %7s | %5s %5s %5s | %6s %6s %6s | leaf' % ('state', 'ptcls', 'closed', 'open', 'P(site)', 'in r22', 'top2.5', 'top5k', 'top10k', 'old15+16', '%state', '%all'))
for nm, n_, (c_, o_), pm, nr, tp, k, lf in rows:
    print('%-12s %7d %6.2f %6.2f | %8.3f %7d | %5d %5d %5d | %6d %5.1f%% %5.1f%% | %s' % (nm, n_, c_, o_, pm, nr, tp[0], tp[1], tp[2], k, 100 * k / max(n_, 1), 100 * k / tc, lf))
open('$O/top6.txt', 'w').write('\n'.join(r[0] for r in rows[:6]) + '\n')
PYEOF
$PY $T/radial_extent.py $CONS $CL --maps $(awk -v v=$V '{printf "%s:%s/%s.mrc ", $1, v, $1}' $O/top6.txt) 2>&1 | grep -v Warning
