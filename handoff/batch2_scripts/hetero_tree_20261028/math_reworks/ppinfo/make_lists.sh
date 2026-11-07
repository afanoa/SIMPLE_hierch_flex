#!/bin/bash
# the ppinfo run lists: make_lists.sh > list_<phase>.txt. Sites = the seven saved sites of the lever study (math_reworks/lever/sites,
# replays c0_* carry site_contrast.txt and the saved posteriors); levers = python args of ppinfo_site.py. Run names <config>_<site>.
H=/mnt/beegfs/elmlund/afan/hetero_tree_20261028; M=$H/math_reworks; L=$M/lever/sites; D2=$H/detlim2/runs; FX=$H/detlim/fix
CONF=/mnt/beegfs/elmlund/afan/routes/toy55/maps160/docked_all; TRUE=$H/detlim/pool2/all.txt
# site: key node sitedir truth rare
site() { case $1 in
  L065s1) echo "L065s1 $D2/new_s1/P_L065/1_hetero_tree/node_r $L/c0_L065s1_site5 $FX/P_L065/truth.txt 18";;
  L065s2) echo "L065s2 $D2/new_s2/P_L065/1_hetero_tree/node_r $L/c0_L065s2_site5 $FX/P_L065/truth.txt 18";;
  L08s1)  echo "L08s1 $D2/new_s1/P_L08/1_hetero_tree/node_r $L/c0_L08s1_site1 $FX/P_L08/truth.txt 18";;
  L08s2)  echo "L08s2 $D2/new_s2/P_L08/1_hetero_tree/node_r $L/c0_L08s2_site1 $FX/P_L08/truth.txt 18";;
  L10s1)  echo "L10s1 $D2/new_s1/P_L10/1_hetero_tree/node_r $L/c0_L10s1_site1 $FX/P_L10/truth.txt 18";;
  qtoynull) echo "qtoynull $M/gates/df_qtoynull/1_hetero_tree/node_r $L/c0_qtoynull_site1 - 0";;
  olenull)  echo "olenull $M/gates/df_olenull/1_hetero_tree/node_r $L/c0_olenull_site1 - 0 CB=96";;   # box 180 at 1.908 A: box_crop 70 (lp 10) / 92 (lp 7.5)
esac; }
# config: name -> python args
cfg() { case $1 in
  base)     echo "";;
  b64)      echo "--lp 7.5";;                              # box_crop 64 on pool2 (92 on OLE): hold-out to 7.6 A, EM band unchanged
  b64e12)   echo "--lp 7.5 --emlim 12";;                   # EM shells to 12 A, hold-out 22..31
  b64e10)   echo "--lp 7.5 --emlim 10";;                   # EM shells to 10 A, hold-out 25..31
  b88e8)    echo "--lp 5.5 --emlim 8";;                    # box 88 on pool2: EM shells to 8 A, hold-out 31..43
  core)     echo "--maskmode core";;                       # the site mask's core (>= 0.5) instead of its falloff
  sphere)   echo "--maskmode sphere";;                     # no site mask: the whole molecule (control)
  oracle)   echo "--maskmode oracle:$(ls $CONF/conf18_emd*_m160.mrc | head -1),$(ls $CONF/conf01_emd*_m160.mrc | head -1),0.3";;
  stageB)   echo "--poses stageB";;
  truepose) echo "--poses true:$TRUE";;
  stageB_core) echo "--poses stageB --maskmode core";;
  b64_core) echo "--lp 7.5 --maskmode core";;
  b64_stageB) echo "--lp 7.5 --poses stageB";;
  b64_stageB_core) echo "--lp 7.5 --poses stageB --maskmode core";;
  truth)    echo "--labels truth";;                        # oracle labels: the ceiling of the acceptance under the config
  truth_b64) echo "--labels truth --lp 7.5";;
  truth_core) echo "--labels truth --maskmode core";;
  truth_truepose) echo "--labels truth --poses true:$TRUE";;
  truth_b88) echo "--labels truth --lp 5.5";;
  seed1)    echo "--labels seedtruth --maxits 1";;             # one M + E step from the truth: the posterior of perfect references
  seed1_b64) echo "--labels seedtruth --maxits 1 --lp 7.5";;
  seed1_b64e10) echo "--labels seedtruth --maxits 1 --lp 7.5 --emlim 10";;
  seed1_core) echo "--labels seedtruth --maxits 1 --maskmode core";;
  seed1_truepose) echo "--labels seedtruth --maxits 1 --poses true:$TRUE";;
  seed1_sphere) echo "--labels seedtruth --maxits 1 --maskmode sphere";;
  seed1_stageB) echo "--labels seedtruth --maxits 1 --poses stageB";;
  seed1_oracle) echo "--labels seedtruth --maxits 1 --maskmode oracle:$(ls $CONF/conf18_emd*_m160.mrc | head -1),$(ls $CONF/conf01_emd*_m160.mrc | head -1),0.3";;
  seed1_truepose_core) echo "--labels seedtruth --maxits 1 --poses true:$TRUE --maskmode core";;
  rsph)     echo "--labels file:$M/ppinfo/runs/seed1_sphere_L065s1/site_post_truth1.txt";;   # the whole-molecule E-step posterior accepted under the site mask
  rorc)     echo "--labels file:$M/ppinfo/runs/seed1_oracle_L065s1/site_post_truth1.txt";;
  rtp)      echo "--labels file:$M/ppinfo/runs/truepose_L065s1/site_post_z1hi.txt,$M/ppinfo/runs/truepose_L065s1/site_post_z3lo.txt";;
  truepose_n15) echo "--poses true:$TRUE --posenoise 1.5";;    # the pose-accuracy requirement: true poses jittered
  truepose_n3)  echo "--poses true:$TRUE --posenoise 3";;
  seed1_truepose_n15) echo "--labels seedtruth --maxits 1 --poses true:$TRUE --posenoise 1.5";;
  seed1_truepose_n3)  echo "--labels seedtruth --maxits 1 --poses true:$TRUE --posenoise 3";;
  seedem)   echo "--labels seedtruth";;                      # the full EM started from the truth: does the truth survive the EM
esac; }
for spec in "$@"; do
  c=${spec%%:*}; sites=${spec#*:}
  for s in ${sites//,/ }; do echo "${c}_$s $(site $s) $(cfg $c)"; done
done
