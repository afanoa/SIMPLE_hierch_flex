#!/bin/bash
#SBATCH -J rnasep_render_subclass
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# env: SUB (sub-class id)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs.mrc; [ -s $CAVGS ] || CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs_even.mrc
D=$KF/closed_focal/sub${SUB}; mkdir -p $D; export OMP_NUM_THREADS=4
$PY $W/render_subclass.py $SUB $META $CAVGS $D/sub${SUB}_vs_3d.png --poses $P/meta_fb.npz --scores $KF/closed_focal/lever/focal_score_per_particle.txt --top $KF/closed_focal/lever/states_top.txt,$KF/closed_focal_docked30k/rank/states.txt \
  --maps "docked within-10k=$KF/closed_focal_docked30k/rank/maps/recvol_state01.mrc" "docked top-30k=$KF/closed_focal/lever/maps/recvol_state01.mrc" "deposited 18 (closed)=$(ls $W/maps160/docked/conf18_*.mrc | head -1)" "deposited 27 (open)=$(ls $W/maps160/docked/conf27_*.mrc | head -1)" "consensus=$W/maps160/toy_consensus.mrc" 2>&1 | grep -v Warning
