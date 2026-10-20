#!/bin/bash
# Morning digest of the latent-variable night: key lines of every stream, in order. usage: collect_results.sh [TAG]
W=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${1:-}
echo "################ MODEL CLASS (truth maps): unexplained fraction per band 117-30 / 30-20 / 20-12 / 12-8 A"
for d in model_class_s3 model_class_s6 model_class_c3; do echo "=== $d"; grep -E "^SUMMARY|^   linear K=|^   shift_only|^   deform |^   deform_occ|^   \(K=all|^flow-field|^volume-difference" $W/$d/log.txt 2>/dev/null | cut -c1-300; grep -E "^   deform conf(18|51|61|62)" $W/$d/log.txt 2>/dev/null | cut -c1-260; done
echo; echo "################ SIMULATION"
for f in $W/toy55_sim.*.out $W/toy10c_sim.*.out; do [ -s $f ] && { echo "=== $(basename $f)"; grep -E "calibration:|merged|refined vs|DONE|FAILED" $f | cut -c1-220; }; done
echo; echo "################ BATTERY lat$TAG"
F=$(ls -t $W/toy10${TAG}_battery.*.out 2>/dev/null | head -1); [ -s "$F" ] || { echo "no battery log"; exit 0; }; echo "log: $F"
grep -E "\] (START|truth axis|refined consensus|0 |1 |2 |3 |4 |5 |DONE|.*FAILED)" $F | cut -c1-300
echo "--- oracle"; grep -E "oracle lp[0-9]+ oracle" $F | cut -c1-300
echo "--- flex latent"; grep -E "flex lp20 (axis|ladder|knn|closed|states|junk|pooled)" $F | cut -c1-300
echo "--- windows per round"; grep -E "assign r[0-9] windows|maps r[0-9] maps:|estep r[0-9] lp[0-9]+ E-step|estep r[0-9] lp[0-9]+ (free|smooth) windows \(|poses (before|after) window|window maps r" $F | cut -c1-300
echo "--- window maps vs truth (last round)"; R=$(grep -o "maps r[0-9]" $F | tail -1); grep -E "$R   window" $F | cut -c1-260
echo "--- deformation prototype"; grep -E "3 deform" $F | grep -E "canonical|displacement|flow PCA|flow-corrected|grid references|grid (free|E-step)|windows (free|E-step)" | cut -c1-320
echo "--- refine3D_states / local flex"; grep -E "4 (seed|refine3D_states|local flex|global flex|coarse groups windows)" $F | cut -c1-300
echo "--- pooled 2D"; grep -E "5 " $F | cut -c1-300
echo "--- lp12 latent"; grep -E "6 flex lp12" $F | cut -c1-300
