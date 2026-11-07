---
name: cavg3D-states-findings-2026-10-08
description: NOT/synthetic solve3D_cavgs states from cls_expansion rendered averages - plain docked split is chance-level, view-stratified hack makes streaky maps, the good "state 3" holds the most variable well-rendered parents with broad views; route = band-rule expansion -> solve3D single -> refine3D_states flex
metadata:
  type: project
---
2026-10-08, NOT 2000 view-balanced rendered averages (run 23, 24_cavg3D_states): plain docked 3 states split views
(hemisphere coverage 38/41/66%); the user judged state 3 (66%) the best map. View-stratified hack (b91898965,
SIMPLE_CAVGS_VIEW_STATES=yes: per-view-group standardised state scores) gave 60/58/60% coverage but EXTREMELY STREAKY maps
(forces every view into every state -> view-inconsistent mixtures). Keep it off.
Synthetic cavg harness (scripts/cls_expansion_proto/cavg3D_test, batch2 cavg3D_test/, SNAP=<commit> BOXD=96 LPSTOP=12 ->
~6 min/arm): EMD-8440/8445 60/40, uniform views, SNR 1: plain docked assignment accuracy 0.53 (chance 0.59), both maps
mixtures closer to the majority state. So solve3D_cavgs states do not separate structure even in the easy case.
Balancing in solve3D (balance=cavg) = equal quota per avglink view group (20), then equal per class - select_balanced.py
copies it. Snapshots: build_clsrender.sh copies bin/lib/scripts to ~/SIMPLE_clsrender/snap/<commit>; run jobs from there.
**How to apply:** don't trust state coverage as success; score states on the synthetic harness first.
Why state 3 was good (state_why.py, 24_cavg3D_states/state_why): states sort PARENT CLASSES by rendering regime, not
conformation (85% majority-state share per parent). s1 = collapsed-band parents (capture_mid 0.2-0.5, agree ~0.15) + the
worst-fitting averages; s2 = moderate-variation parents over-rendered (capture ~1.3), 12 view groups; s3 = the most
variable parents (real_mid 75-130) with broadest views. Only s3 fits its own map better than the single-state map (+0.13
corr). Flex_pca on class averages: synthetic 100% purity (true and estimated poses); NOT state-3 alone (661, cFAR ~0) gave
view-grouped flex states (NMI 0.50, 0/8 reproducible axes). Agreed route (user 2026-10-08): band-rule cls_expansion ->
12,600 rendered as particles (class=parent, parents.mrc = mean of rendered) -> solve3D single balance=cavg
force_lp_range lpstart=20 lpstop=8 -> refine3D_states flex=yes nstates=3 (29_render_states/scripts/run_chain.sh).
Import trap: import_particles oritab=<project>.simple oritype=cls3D leaves poses at zero; print cls3D to text first.
Longer-term idea: fit flex 3D covariance from the per-class 2D covariances of cls_expansion (Python toy first).
Recipe trace (2026-10-15, user asked): states3 = abinitio2D ncls=300 -> map_cavgs_selection (selected.spi, 126 classes,
3_selection) -> cls_expansion ncls=100 at fc252c2cc (run 23, MODEL-RENDERED averages cls_expansion_cavgs_rendered.mrc)
-> view groups from solve3D balance=cavg nclust=20 (killed after view_partition.txt) -> select_balanced.py 100/group =
2000 -> even/odd = byte copies of the subset -> solve3D_cavgs nstates=3 docked lpstart_ini3D=20 lpstop_ini3D=8 (639/700/661).
Master e232ed660 (2026-10-14) REMOVED the model/rendered route (kernel only), and no fc252c2cc snapshot exists
(~/SIMPLE_clsrender is now 8588389c7 dirty); run 23's rendered stack and 24_cavg3D_states/subset are still on disk.
The failed kernel retry (kernel100_expansion -> s3d420_lp30_12_k3): all 12,600, balance none, lp 30-12 -> 3144/767/8689.
2026-10-15 kernel rerun (user: "only change" kernel averages + 12.6k total, then real even/odd, more than 2000):
34_kernel_cavg3D_states, job 53886125: kernel100_expansion averages + real halves, select_nested.py = master solve3D
balance=cavg nested quota drawn once (target 4000 -> 4008; reproduces 100/group = 2000 on run 23), 3D at
snap/b91898965 (docked; master solve3D_cavgs rejects multivol_mode). Master cls_expansion ncls is a TOTAL budget
(groups get pop**0.5 shares, minpop caps), not per class: 33-182 per class on NOT.
Control (user: is balancing doing anything or is docked doing it all): 35_kernel_cavg3D_all, job 53886494 = same 3D
step on all 12,600 kernel averages (symlinked stacks + real halves), no balancing. Compare with 34 when both finish.
User cancelled the 10-part versions ("one part with 32 cores, it's faster"): reruns 53887035 (34, balanced 4008) and
53887036 (35, all 12,600), -c 32 --mem=96G, SIMPLE_QSYS=local, nthr=32, no nparts (10-14 one-node 12.6k run peaked 49 GB).
Master control (user: "a control solve3D_cavgs in the latest code"): 36_kernel_cavg3D_master, job 53887885, origin/master
b57538b7b built in ~/SIMPLE_m1015 (login node, snap/b57538b7b; clone of ~/SIMPLE_kernel, origin = GitHub) on the 34
subset: independent states, reseed_states default no. Docked was removed from master solve3D_cavgs by fbffebf73 (10-07).
Results: 34 (docked, 4008) 1133/1867/1008 in 1h56; 36 (master independent, 4008) 2321/1202/485 in 4h19 (pops swung
25/50/25 -> 3/91/6 at stage 3 -> 51/39/10). FSC 0.143 = 6.5 A everywhere = Nyquist of the 168-px crop at lp 8, not a
measurement. Timing (bench files): one-node cost is the prob pre-step (rows x nspace x nstates, every row every iteration)
+ 3-state assembly (98 s/iter one-node vs 9 s distributed); 10-part distributed was 16 s/iter vs 28 s one-node at stage 1
on 2026-10-15, so "one node is faster" did NOT hold that day (queue fast).
2026-10-15 CLEANUP (user): /mnt/beegfs/elmlund/afan/not holds ONLY 1_abinitio2D, 3_selection, kernel100_expansion
(1_cls_expansion -> 2_selection [selection_4008.txt = run 36's sub-classes] -> 3_solve3D_cavgs, run as 36: job 53894402)
and 36_kernel_cavg3D_master (reference; its project's input stack is in the archive's 34 dir). Everything else, incl.
24_cavg3D_states, 34, 35, the cls_expansion test dirs and root logs/scripts, moved to not/_archive_2026-10-15 (not
deleted). Keep that dir lean: new steps go inside kernel100_expansion as N_prg.
