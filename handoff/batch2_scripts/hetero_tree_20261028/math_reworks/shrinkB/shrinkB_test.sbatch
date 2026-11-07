#!/bin/bash
#SBATCH -J shrinkBt
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/logs/%x_%j.out
# shrinkB.py on production's final stage-B maps of three gate nodes: copy must be byte-identical, dev1 / dev2 must keep the
# count-weighted node mean of each half; prints the per-band shrinkage of every state's deviation.
set -uo pipefail
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
counts() {
  awk -v K=$2 '/PCG DISTRIBUTED \| STATE=/ && /KIND=base/ { split($0, a, "|"); s = a[2]; gsub(/[^0-9]/, "", s); n = a[5]; gsub(/[^0-9]/, "", n)
    if (s + 0 > 0) { if (a[3] ~ /even/) e[s + 0] = n; else o[s + 0] = n } }
    END { for (k = 1; k <= K; k++) printf "%s%d", (k > 1 ? "," : ""), e[k] + o[k] }' $1
}
for n in mm2_qtoyflex/1_hetero_tree/node_r mm2_r131/1_hetero_tree/node_r2 mm2_qtoynull/1_hetero_tree/node_r; do
  N=$G/$n; K=$(grep -o 'nstates=[0-9]*' $N/stageB/hetero_step.sh | cut -d= -f2); MSK=$(grep -o 'mskdiam=[0-9.]*' $N/stageB/hetero_step.sh | cut -d= -f2)
  cnt=$(counts $N/stageB/hetero_step.log $K); T=$M/test/$(echo $n | cut -d/ -f1)_$(basename $N)
  echo "######## $n: K $K counts $cnt mskdiam $MSK"
  for mode in copy dev1 dev2; do
    $PY $M/shrinkB.py $N/stageB $T/$mode $mode --counts $cnt --mskdiam $MSK || echo "FAILED $mode"
  done
  for k in $(seq -f %02g 1 $K); do for s in "" _even _odd; do cmp -s $N/stageB/recvol_state$k$s.mrc $T/copy/recvol_state$k$s.mrc || echo "copy differs: state $k$s"; done; done
  $PY - $N/stageB $T $K $cnt <<'EOF'
import sys, numpy as np
src, T, K, cnt = sys.argv[1], sys.argv[2], int(sys.argv[3]), np.array([float(c) for c in sys.argv[4].split(',')])
w = cnt / cnt.sum()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0])
    return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).astype(np.float64)
for mode in ('dev1', 'dev2'):
    for h in ('_even', '_odd'):
        a = [rd('%s/recvol_state%02d%s.mrc' % (src, k + 1, h)) for k in range(K)]; b = [rd('%s/%s/recvol_state%02d%s.mrc' % (T, mode, k + 1, h)) for k in range(K)]
        ma = sum(w[k] * a[k] for k in range(K)); mb = sum(w[k] * b[k] for k in range(K))
        dev_a = np.sqrt(np.mean([((a[k] - ma) ** 2).mean() for k in range(K)])); dev_b = np.sqrt(np.mean([((b[k] - mb) ** 2).mean() for k in range(K)]))
        print('%s%s: node mean kept to %.2e (rms %.3e); rms deviation %.3e -> %.3e (x %.3f)' % (mode, h, np.abs(ma - mb).max(), np.sqrt((ma ** 2).mean()), dev_a, dev_b, dev_b / dev_a))
EOF
done
echo "done $(date)"
