---
name: flex-tests-run-on-cluster-only
description: "verification/test jobs go on batch2 (slurm), never on the box fwl-c142339, even when the box has the binary and data"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d3cd5123-c794-474a-98bc-c9d5a47fe521
  modified: 2026-09-16T04:34:33.962Z
---

2026-09-16: the user stopped me the moment a verification chain started on the box: "dont run
things on the box, everything should run on the cluster."

**Why:** the box is the user's interactive machine (their own flex_pca runs, e.g. trpm4 under
`~/SIMPLE_fast`, start there at any time) and it is shared with elmlundho; batch2/csbdevel is where
batch work belongs and where the per-arm records (sacct, out files) live.

**How to apply:** build the branch under test on batch2 and submit through `~/flex_verify/driver.sh`
(sbatch masters that spawn their own parts). Use the box only to read/score results from beegfs or
for a warm plane cache the user explicitly asks for. Always list every submitted job (id, arm,
build, dataset) in the report. Related: [[box-shared-with-elmlundho]], [[flex-pca-fork-map-2026-09-15]].
