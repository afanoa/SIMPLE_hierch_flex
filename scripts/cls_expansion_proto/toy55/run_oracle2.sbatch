#!/bin/bash
#SBATCH -J toy10_oracle2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Headroom measurements: (1) how much purity the HARD family splits throw away (per conformer, the share of its
# particles that sit in its majority family / sub-family); (2) the per-particle oracle (true conformer maps as
# references) at the REALIGNED poses at 10, 8 and 6 A, against the consensus-pose oracle (0.58-0.59 at 10-8 A).
# env: TABLE (pose table), LPS ("10 8 6")
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=8; cd $W
T=$W/sim; GT=$T/gt_particles.txt; STK=$T/all.mrcs; STATES=1,10,18,27,28,37,46,48,49,51; AXIS=-0.246,-0.916,-0.442,1.000,0.972,-0.685,-0.489,0.446,0.499,-0.139
TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc | head -1; done | tr "\n" , | sed 's/,$//'); TABLE=${TABLE:-$W/lat/realign_realign1/ptcl3D_realigned.txt}
echo "== purity of the hard family splits (share of each conformer's particles in its majority family)"
$PY - <<PYEOF
import numpy as np, glob
gt = {}
for l in open('$GT'):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
for E in ('expand2b', 'expand_realign1_bal05', 'expand_l3b_fam1', 'expand_l3b_fam2'):
    fs = sorted(glob.glob('$W/lat/%s/family*_subset.txt' % E)); fam = {}
    for j, f in enumerate(fs):
        for p in np.loadtxt(f, dtype=int): fam[int(p)] = j
    confs = sorted(set(gt[p] for p in fam if gt.get(p, 0) > 0)); tot = 0; inmaj = 0; per = []
    for c in confs:
        cnt = np.bincount([fam[p] for p in fam if gt.get(p, 0) == c], minlength=len(fs)); per.append('%d:%.2f(f%d)' % (c, cnt.max() / cnt.sum(), cnt.argmax() + 1)); tot += cnt.sum(); inmaj += cnt.max()
    print('  %-24s %d families: particles in their conformer majority family %.3f | %s' % (E, len(fs), inmaj / tot, ' '.join(per)))
PYEOF
for LP in ${LPS:-10 8 6}; do
  $PY $W/latent_windows.py estep $STK $TABLE $GT --refs $TRUTH --lp $LP --out $W/lat/oracle_realigned_lp$LP.txt --states $STATES --axis=$AXIS --nsub 40000 --tag "oracle realigned lp$LP" 2>&1 | grep -v "^   [0-9]* / " | sed "s/^/[realigned] /"
done
echo "[$(date '+%F %T')] oracle2 DONE"
