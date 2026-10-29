#!/bin/bash
# Run from a configured SIMPLE environment. First argument is the absolute Python launcher path.
#SBATCH --job-name=ht_state_maps
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
launcher=${1:?Pass the absolute path of hetero_tree_state_maps.py}
shift
[[ $launcher == /* && -f $launcher ]] || exit 2
exec python3 "$launcher" --cpus "${SLURM_CPUS_PER_TASK:?}" "$@"
