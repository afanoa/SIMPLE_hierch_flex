# hetero_tree: the mathematical reworks (November 2026)

Branch `claude/tree-math-reworks` off `claude/tree-dist` (9005a8c77). Spec: sections 4, 6 and 7 of
`hetero_tree_algorithmic_audit_2026_11_01.md`, amended on 2026-11-01 (the junk screen stays, with two soft uses). Every
rework is a same-seed replay on the gate fixtures before any default changes: qtoy2 (maxdepth 1), the r131 closed
subtree (maxdepth 1, the closed corner must be retained), the nulls qtoynull / olenull / micnull (maxdepth 0, the largest
leaf must hold >= 80 % of the particles) and, for reworks 1 and 2, the OLE span. Replay runs live under
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/` on batch2; builds there follow the dist-build recipe (an
x6342 node, Release, ARCHOPT, fast gate), so the numerics match the baseline builds.

Gate keys are the recipe18d gate runs' (the latest accepted arms): qtoy2 `nparts=4 mskdiam=232 minpop=900 site_layer=yes
site_states=yes kway_maxdepth=2` at maxdepth 1; r131 `nparts=12 mskdiam=232 minpop=3000 site_layer=yes site_states=yes
kway_maxdepth=-1 split_ntests=460` at maxdepth 1, seed 2, from `ancestral/fixtures/r131`; the nulls at maxdepth 0 with their
recorded keys; all with `oneflex=yes sta0_npeaks=2 junk_screen=no` and `SIMPLE_SPLIT_CONTRAST=maps`, maps kept.

## 1. Endpoint reproducibility by one subspace test (audit 6.3)

(in progress)

### 1.1 What changed

`site_rank` keeps the four endpoint EM fits (the lo / hi tails of fflex directions 1 and 3) and drops every halves and
quarters refit. Each fit already yields the endpoint's cross-half deviation maps at N* (even and odd, inside the site
mask). The hold-out-band Gram matrices of these maps (even.even, odd.odd, even.odd; `dev_gram`) carry every decision:

- an endpoint is accepted when the cross-half correlation of its deviation is `HALVES_KSIG` (3) Fisher SEs above zero,
  with `n_eff` from `mask_neff` (the same n_eff the hold-out SE in `site_starts.txt` used);
- two accepted endpoints are one pole when they correlate positively and neither's residual after the other reproduces
  across halves (the principal-angle test, `gram_resid`: project the other's direction out in each half, correlate the
  residuals, Fisher z below 3); groups merge closest first (highest disattenuated correlation) when every cross pair is
  such a pair; the two ends of one motion are anti-correlated and stay distinct poles by construction;
- a pole whose deviation lies in the span of two or more other positively correlated poles is redundant and dropped;
- the rest (best pole posterior below 0.5) is ranked by one minus that posterior and scored once (no EM); it is a state
  when its deviation reproduces and is distinct from every pole it correlates with positively;
- the pole rule keeps a pole when its share is below 0.5 and its seed set, scored once by the pole's posterior, gives a
  deviation that reproduces and is distinct from the other poles' directions. `site_pole_rule.txt` keeps its columns:
  the two margins are now the Fisher z of reproducibility and of distinctness, so `strongest_pole` reads them unchanged.

Per site: 4 EM fits + 1 scoring pass for the rest + 1 per minority pole, against 4 x 7 + 6 + 6 x (kept poles) fits
before. `site_halves.txt`, `site_poles.txt`, `site_seeds.txt`, `site_starts.txt`, `site_neff.txt` keep their formats
(the controller's parsers are unchanged); `site_gram.txt` is new. The `rowmap` key is still accepted and now unused.

Prototype (`subspace_proto.py` on the Gram dumps of 32 saved sites, qtoy2 / qtoynull / olenull roots): on qtoy2 the
grouping reproduces the halves verdict where that verdict was decisive (sites 1, 2, 5) and keeps apart endpoints the
halves test merged at a disattenuated correlation of 0.85 with a residual that reproduces at z 8 (sites 3, 4); on the
nulls the Fisher-z test rejects every endpoint at n_eff below about 100 and accepts the same quality-driven deviations
the halves test accepted at high n_eff (those are left to the pole rule and the split test, as before). The eigenvalue
rank of the symmetrised cross-half correlation matrix has no power against near-collinear ends and is not used.
