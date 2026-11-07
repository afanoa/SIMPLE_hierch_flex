#!/bin/bash
#SBATCH -p norm
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# bshare.py on the option (a) validation sets (quality_a/<set>/cases.txt + qscreen.txt), joined with their ratios.txt labels and ratio_a
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; Q=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_a
for s in ${SETS:-m79 d80 cd77 igm77 ole10 ole11}; do ( cd $Q/$s && OMP_NUM_THREADS=4 $PY $T/bshare.py cases.txt < qscreen.txt > bshare.txt 2>&1 ) & done; wait
for s in ${SETS:-m79 d80 cd77 igm77 ole10 ole11}; do echo "=== $s"
  $PY - $Q/$s <<'PYEOF'
import re, sys
D = sys.argv[1]; R = {}
for l in open(D + '/ratios.txt'):
    if l.startswith('#'): continue
    w = l.split(' | ')[0].split(); R[w[0]] = (w[14], w[16], w[19], l.split(' | ')[1].strip() if ' | ' in l else '')
for l in open(D + '/bshare.txt'):
    m = re.match(r'BS (\S+)\s+\S+\s+\| (\d+) shells \| share (\S+) \| c0 (\S+) dBm (\S+) A\^2 \| rel (\S+)', l)
    if not m: print(l.rstrip()); continue
    r = R.get(m.group(1), ('-', '-', '-', ''))
    print('%-22s share %s c0 %s dBm %7s rel %s | ratio_a %6s dB/ref %6s | %-14s %s' % (m.group(1), m.group(3), m.group(4), m.group(5), m.group(6), r[0], r[1], r[2], r[3][:70]))
PYEOF
done
