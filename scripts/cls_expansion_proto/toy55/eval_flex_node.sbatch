#!/bin/bash
#SBATCH -J toy10_flex_node_eval
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# evaluates the node-local flex latents (lat/flex_node_fam1,2) against the truth next to the global latent on the same
# particles; waits for the flex runs (the in-job evaluation of run_flex_node died on the --axis argument parsing)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; GT=$W/sim/gt_particles.txt
STATES=1,10,18,27,28,37,46,48,49,51; AXIS=-0.246,-0.916,-0.442,1.000,0.972,-0.685,-0.489,0.446,0.499,-0.139
for i in $(seq 1 200); do n=0; for f in 1 2; do ls $W/lat/flex_node_fam$f/*_flex_pca/flex_pca_embedding.bin >/dev/null 2>&1 && n=$((n+1)); done; [ $n = 2 ] && break; grep -q FAILED $W/toy10_flex_node.*.out 2>/dev/null && { echo "a node flex run FAILED"; grep -h FAILED $W/toy10_flex_node.*.out; }; sleep 60; done
for f in 1 2; do FD=$(ls -td $W/lat/flex_node_fam$f/*_flex_pca 2>/dev/null | head -1); [ -n "$FD" ] || { echo "fam$f: no flex output"; continue; }
  echo "=================== node fam$f: local chart"; $PY $W/ptcl_latent_eval.py $FD $GT --states $STATES --axis=$AXIS --subset $W/lat/flex_node_fam$f/states.txt --tag node_local 2>&1 | grep -v Warning | tail -14
  echo "=================== node fam$f: global chart on the same particles"; $PY $W/ptcl_latent_eval.py $(ls -td $W/lat/flex12/*_flex_pca | head -1) $GT --states $STATES --axis=$AXIS --subset $W/lat/flex_node_fam$f/states.txt --tag global_on_node 2>&1 | grep -v Warning | tail -14
done; echo "ALL DONE"
