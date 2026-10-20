#!/bin/bash
#SBATCH -J klem_readout55
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# all-55 identification readout of a finished run (env RTAG): every leaf state map (level 1 + level 2) of every region
# against the 55 docked apo conformers (maps160/docked_all); prints the distinct conformers identified and the per-map table
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; S=$R/readout55${RTAGS:-}; rm -rf $S; mkdir -p $S/states
for d in $R/fam?/${LEAFGLOB:-region*}/level1/states_dir $R/fam?/${LEAFGLOB:-region*}/level2_s*/states_dir; do [ -d $d ] || continue; t=$(echo $d | sed "s#$R/##; s#/states_dir##; s#/#_#g"); for m in $d/state_*.mrc; do [ -f $m ] && ln -sf $m $S/states/state_${t}_$(basename $m .mrc | sed "s/^state_//").mrc; done; done
ls $S/states > $S/map_order.txt   # deposited_match numbers the states in this (sorted) order
echo "$RTAG: $(ls $S/states | wc -l) leaf state maps"
$PY $W/deposited_match.py $S/states $W/maps160/docked_all $W/maps160/toy_consensus.mrc --save $S/matches.txt 2>&1 | grep -v Warning | tail -4
$PY - <<PYEOF
names = [l.strip() for l in open('$S/map_order.txt')]; rows = [l.split() for l in open('$S/matches.txt') if l.strip() and not l.startswith('#')]
best = {}
for r in rows:
    i, c, cc = int(r[0]) - 1, r[1], float(r[2])
    if c not in best or cc > best[c][0]: best[c] = (cc, names[i])
print('$RTAG: %d of 55 deposited apo conformers are the best match of at least one leaf state map' % len(best))
with open('$S/identification55.txt', 'w') as o:
    o.write('# map  best_matching_deposited_conformer  corr (55 docked apo conformers, deviation correlation 60-20 A)\n')
    for r in rows: o.write('%-44s %4s  %s\n' % (names[int(r[0]) - 1], r[1], r[2]))
PYEOF
