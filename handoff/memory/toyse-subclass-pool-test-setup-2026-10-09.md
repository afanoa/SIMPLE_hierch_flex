---
name: toyse-subclass-pool-test-setup-2026-10-09
description: "Subclass-pooled placement test on toySE: scripts/cls_expansion_proto/subclass_pool/{pool_place.py,run_toyse_pool.sbatch} written 2026-10-09 off-site (batch2/box NXDOMAIN), scorer verified on a synthetic cache; paths TOYSE/FLEXDIR/B are assumptions to check before the first sbatch"
metadata:
  node_type: memory
  type: project
  originSessionId: a6c2687a-b3cd-4070-9eaa-6163aa554f8d
  modified: 2026-10-10T02:39:53.584Z
---

Set up 2026-10-09 while off the institute network (batch2.ncifcrf.gov and fwl-c142339 both NXDOMAIN), so nothing was
submitted. Files (in the repo, not the scratchpad): `scripts/cls_expansion_proto/subclass_pool/pool_place.py` (reads
flex_pca_embedding.bin SIMPLFXC v3 + cls_expansion_labels.txt + gt_particles.txt; arms ptcl / class / sub / sub_prec =
lowest-inertia k-means over particles, parent-class means, subclass means, precision-weighted subclass means; Hungarian
accuracy, per-state recall/purity, NMI, label-only ceiling) and `run_toyse_pool.sbatch` (one job: copy toySE project ->
update_project qsys_name=slurm -> abinitio2D ncls=50 nparts=10 nthr=8 -> cls_expansion ncls=4 lp=20 on $B build ->
pool_place.py on the EXISTING toySE flex embedding). Scorer verified on a synthetic cache written in the Fortran stream
layout (sub 0.975 vs ptcl 0.28 at K=4, trivially state-pure classes: a mechanics check only).
**Assumptions to verify on batch2 before submitting:** TOYSE=/mnt/beegfs/elmlund/afan/softem/toySE holds one .simple
project and gt_particles.txt (toy2p_concat format); FLEXDIR holds the toySE flex_pca delivery with
flex_pca_embedding.bin (default $TOYSE/flex is a guess); B=$HOME/SIMPLE_clsrender is a cls-expansion-solvent build;
scripts copied to /mnt/beegfs/elmlund/afan/softem/scripts. Expected readout: does `sub` recover state 13 (10 %) where
`ptcl` and the flex GMM fold lose it ([[flex-placement-not-latent-loses-small-states]]); `class` reproduces the old
cluster-tracking arm ([[pose-state-coupling-toy-verdict]]). Related: [[cls-expansion-vs-flex-pca-latent-diagnosis-2026-10-09]].
