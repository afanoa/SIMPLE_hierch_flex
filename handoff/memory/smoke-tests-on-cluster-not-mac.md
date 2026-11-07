---
name: smoke-tests-on-cluster-not-mac
description: "User 2026-10-30 ('cant you just submit on cluster'): run smoke/pipeline tests as small sbatch jobs, not on the Mac (its disk was full)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 36dc059e-25bd-4e47-9f7c-87d5f98d1a8b
  modified: 2026-10-30T20:17:08.917Z
---

2026-10-30, after I smoke-tested the GRFP Stage 1 SIMPLE chain on the Mac build (~/SIMPLE_master/build-mac) with
synthetic data and hit a full Mac disk (372 MB free of 926 GB): "cant you just submit on cluter".

**Why:** the Mac has no room and is not the place for test runs; the cluster has the data, builds and space.

**How to apply:** test new pipelines with a small sbatch probe on batch2/batch (norm or csbdevel), e.g. a
600-frame subset, writing to the experiment folder; keep the Mac to editing scripts and reading small result
files. Related: [[flex-tests-run-on-cluster-only]], [[head-node-runs-nothing-but-sbatch]].
