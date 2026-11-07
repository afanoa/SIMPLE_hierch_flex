---
name: flex-pca-hierarchical-gmm-auto
description: "SIMPLE_COV_GMM_AUTO (default ON): island/continuum factorization inside the GMM; the two respawn bugs, the 5000 floor, and the validated Sinkhorn-EM continuum option"
metadata: 
  node_type: memory
  type: project
  originSessionId: fb1d9e33-0137-47d8-8d28-e3a32330f274
  modified: 2026-09-03T13:53:49.405Z
---

Implemented + validated 2026-09-02 on 10028/cryoDRGN (h0901-int2, commits
c0817fe33..a86bdd2d2; final run auto4_gmm_br360). Solves the factor-mixing
failure where one joint GMM cuts across composition x conformation (the
missing-40S minority polluted "ratcheting" states).

**The mechanism (knob-free, all inside the mixture):** over-fit tied-cov GMM at
kfit=2*npreimages -> merge pairs whose 1-D density along the line between means
is unimodal (EXACT under tied covariance: d<=2sigma unimodal for any weights,
else grid dip test with the fitted pis) -> transitive closure = macro-clusters
(continuum tiles chain, islands keep a dip) -> >=1 state per macro, seats by
mass with capacity minocc*(seats+1) -> respawn-free within-macro GMM refits
with pi floor minocc/N. minocc=5000 default (SIMPLE_COV_MIN_STATE). Falls back
to plain GMM when nmac<=1, nmac>=nstates, orient/deconv active.

**Result on 10028 (no pruning, full 105k):** macro-clusters=2; island 7,881
ptcls -> 1 state @7.2A (owns the whole UMAP satellite); continuum -> 5 states
18.1k/24.8k/21.7k/22.7k/10.0k @3.9-4.6A, centroids tracing an ordered arc.
Matches the manual selection(states=2,3,4,5)+prune_project reference arm.

**Three engine-vs-prototype bugs, all found in PYTHON ON THE CACHE first:**
1. Discovery needs ~60 CLEAN EM iters to separate the island; respawn resets
   ate the 60 cap -> maxit=300 for discovery+refits.
2. RESPAWN IS A BRIDGE-MAKER: it re-seeds a component on the worst-explained
   particle = an outlier in the gap BETWEEN clusters; a pi=0.0028 bridge at
   1.76sigma chained 2 macros into 1. Discovery + within-macro refits run
   respawn=.false.; trivial comps (pi<0.25/kfit) get no merge-graph vote.
   Respawn ALSO caused the 220/381-ptcl junk states in within-macro refits
   (respawn-free same fit: min occupancy 5,761).
3. Component-1 quantile init fails at some kfit (amplitude anomaly reorders
   per-axis variance); init along the max OBSERVED-variance axis passes at
   every kfit in the sweep.

**Workflow lesson (user-directed):** iterate placement approaches in python on
flex_pca_embedding.bin (seconds) before any build+run cycle (~35 min); the
engine run is only the final confirmation + maps. infile= replay skips the EM
fit for placement-only engine checks.

**Continuum readout findings:** per-particle rotational state is measurement-
noise-limited, NOT an embedding defect (mean 3.4 active states/ptcl, sub-sigma
axis separations; population motion solid at 0.99 basis reproducibility). The
stored per-particle precision in the cache is BROKEN on current runs (50%
asymmetric, indefinite for every particle) - blocks measurement-error
deconvolution and is now the top reason to finish the amplitude-anomaly
diagnosis. Axis-free equal-mass continuum placement VALIDATED in python:
Sinkhorn-EM (= GMM with mixing weights replaced by OT duals so every state
gets N/B mass; farthest-first init; no principal-component choice) gives occ
14.9-23.8k vs free GMM 5.8-33.2k with the same separations, states forming a
4.65sigma arc. WIRED into Fortran later that day (4ea25d7c8,
sinkhorn_state_weights): continuum macros deliver HARD argmax equal slices
(soft responsibilities of overlapping components made every state map the
consensus — neff ~= N for all states is the red flag in the occupancy table).

**AFM/RNA 932k campaign verdict (2026-09-02 evening):** poses must come from
24_refine3D_auto (23's give blob states; user confirmed). AUTO on 24: 3
macro-clusters (27.6k junk / 819.5k main / 85.3k junk-ish) — user judged all
three non-main states junk. Sinkhorn hard-balanced slices of the main 819k
(delivered: 164-232k disjoint quarters at 4.4-4.5 A, consensus-corr 0.90->0.85)
STILL look identical: the embedding is noise-limited on this data, placement
cannot fix it. Escalation (user-directed): multi-state abinitio3D (nstates=4,
39_abinitio3D_multi) for particle-level junk removal BEFORE flex; then the
basis-lever ladder (convergence-gated iters, frequency marching as band
schedule, rec_backend=pcg mask-inside-solve = solvent prior) is the next
battery — this dataset can rank those levers, 10028 could not. NOTE: the AFM
and 10028 auto configs ran the SINGLE full-N fit (no SIMPLE_COV_PAIRED in that
env lineage) — basis self-agreement 0.93->0.98 at it6, still climbing =
NOT converged; no cross-half table exists for these runs.

**Defaults flip + CNGA1 (2026-09-04):** the config-of-record IS now the engine
default (9eee2b6): bare `flex_pca projfile vol1 mskdiam [box_rec] nparts nthr`
runs polar E-step + ONLINE + EM_MIX16 + EM_DEFLATE4 + BG-deflation + relprior,
self-stop off, neigs=15 npre=6 it=6, band DERIVED (derive_flex_pca_band; it
chose 9.5 A on CNGA1 vs the hand-set 16). Bare-command 10028 regression PASSED
(island 8,136 vs auto4's 7,881 = derived-band boundary shift; sinkhorn slices
all >=5k). Single-continuum delivery re-routed through Sinkhorn (b866026)
after CNGA1's plain-GMM fallback hoarded 87% into one state and left 3 states
sub-floor. CNGA1 (245k, C4, box 288@0.844, mskdiam 180, poses 10_refine3D_auto)
= pure continuum (1 macro), balanced slices 33.6-50.7k @4.8-5.0 A. IO lesson
requantified: BeegFS-hosted stacks leave workers in D-state (31/80 cores,
2017 s); staging the stack on LOCAL disk = 924 s. reimport_particles REFUSES a
stack whose image count exceeds the project's (245,424 of 274,333 referenced);
the working repoint is the length-preserving BYTE-PATCH of the stk path
(pad with a double slash to match length). NOTE: paired halves (SIMPLE_COV_
PAIRED) is still NOT in the defaults - the RECOVAR-style cross-half covariance
regularization argues it should be; candidate next flip.

**RESTORE 2026-09-04 (user-directed):** the box build (~/SIMPLE-h0901-int2) is
checked out at branch `auto4-restore` = e60f95070 (local a86bdd2d2, "floor +
robust discovery init") — the binary that produced auto4_gmm_br360, with SOFT
within-macro GMM placement. Everything after it is preserved at tag
`tip-20260904` (f3d0d9c92): Sinkhorn hard-balanced continuum slices
(within-macro AND single-continuum), the defaults flip (bare command =
config-of-record + derived band), and the SIMPLE_COV_FOCUS_MSK protomer mask.
To return: `git checkout h0901-int2` (or the tag) + make install. Reason:
user wants the auto4-exact behavior as the working binary for now; Sinkhorn
placement is hardwired after e60f950 (no runtime switch back to soft GMM).
Open at restore time: star-import geometry defect (user is taking it over) —
imported projects' D5 dies at shell>=3 with a non-monotonic sh4 bounce
(aliasing-like) while all ptcl3D fields round-trip bit-identical; untested
discriminator: cache=no on the imported project.
