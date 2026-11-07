#!/bin/bash
#SBATCH -J lrec_recursive
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# recursive hard classification of the oracle swing family: K1 states (lrec_cls2, random seeds), then K2 states inside each of
# them (random seeds), leaves = (state, substate); readout = per-leaf mean +- sd of the true swing and eta2 over the family.
# env LABEL (rec), K1 (4), K2 (4), TAG, FAMFROM, FAMLAB (as lrec_cls2)
set -uo pipefail
T=/mnt/beegfs/elmlund/afan/routes/toy55; KD=/mnt/beegfs/elmlund/afan/rnasep_klem; L=${LABEL:-rec}; K1=${K1:-4}; K2=${K2:-4}
LABEL=$L K=$K1 SEED=random MAXITS=8 bash $T/lrec_cls2.sbatch || { echo "[lrec_recursive] level 1 FAILED"; exit 1; }
for s in $(seq 1 $K1); do
  LABEL=${L}_s$s K=$K2 SEED=random MAXITS=8 FAMFROM=$KD/lrec30/$L/final_states.txt FAMLAB=$s bash $T/lrec_cls2.sbatch || echo "[lrec_recursive] level 2 state $s FAILED"
done
paste $(for s in $(seq 1 $K1); do echo $KD/lrec30/${L}_s$s/final_states.txt; done) | awk -v k2=$K2 '{l = 0; for (i = 1; i <= NF; i++) if ($i > 0) l = (i - 1) * k2 + $i; print l}' > $KD/lrec30/$L/leaf_states.txt
awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; if (k >= 200 && k <= 230) {x = k - 200; n[s]++; sw[s] += x; sq[s] += x * x; T1 += x; T2 += x * x; NT++}}
  END {m = 0; for (s in n) {mu[s] = sw[s] / n[s]; ord[++m] = s; W += sq[s] - sw[s] ^ 2 / n[s]}
    for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (mu[ord[j]] < mu[ord[i]]) {t = ord[i]; ord[i] = ord[j]; ord[j] = t}
    printf "RECURSIVE leaves (%d): eta2 %.3f; by mean:", m, 1 - W / (T2 - T1 ^ 2 / NT); for (i = 1; i <= m; i++) {s = ord[i]; printf " %d:%.1f+-%.1f(%d)", s, mu[s], sqrt(sq[s] / n[s] - mu[s] ^ 2), n[s]}; print ""}' $T/simswing30/gt_particles.txt $KD/lrec30/$L/leaf_states.txt
echo "[lrec_recursive $L] DONE"
