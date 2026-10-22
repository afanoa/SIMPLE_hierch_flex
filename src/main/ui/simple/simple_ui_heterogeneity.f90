!@descr: module defining the user interfaces for heterogeneity-analysis programs in the simple_exec suite
module simple_ui_heterogeneity
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('heterogeneity', 'Heterogeneity Analysis', 65)
type(ui_program), target :: flex_pca
type(ui_program), target :: refine3D_states
type(ui_program), target :: refine3D_states_cavgs
type(ui_program), target :: classify3D_refs
type(ui_program), target :: ptcl3D_state_consensus
type(ui_program), target :: state_tree
type(ui_program), target :: state_mergeup
type(ui_program), target :: state_split
type(ui_program), target :: site_rank
type(ui_program), target :: hetero_tree
type(ui_program), target :: landscape_compare
type(ui_program), target :: site_discover
type(ui_program), target :: landscape_sample

contains

    subroutine construct_heterogeneity_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_flex_pca(prgtab)
        call new_refine3D_states(prgtab)
        call new_refine3D_states_cavgs(prgtab)
        call new_classify3D_refs(prgtab)
        call new_ptcl3D_state_consensus(prgtab)
        call new_state_tree(prgtab)
        call new_state_mergeup(prgtab)
        call new_state_split(prgtab)
        call new_site_rank(prgtab)
        call new_hetero_tree(prgtab)
        call new_landscape_compare(prgtab)
        call new_site_discover(prgtab)
        call new_landscape_sample(prgtab)
    end subroutine construct_heterogeneity_programs

    subroutine new_flex_pca( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call flex_pca%new(&
        &'flex_pca',&
        &'Projection-aware covariance heterogeneity states',&
        &'fits a sigma-whitened low-rank 3D covariance factor directly from fixed-pose particles, performs MAP latent inference, and reconstructs combined/even/odd kernel states without a diffusion graph',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_STANDARD, display_name='Covariance Heterogeneity States')
        call flex_pca%add_input(UI_IMG, 'vol1', 'file', &
            'Consensus mean volume', 'Fixed mean subtracted from every particle; when omitted, the &
            &project consensus map (out segment, state 1) at native sampling is used', &
            'e.g. vol1.mrc (consensus mean)', .false., '', &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'neigs', 'num', &
            'Covariance components (default 10)', 'Number of fitted low-rank covariance factors; capped at 48. &
            &Over-ranking dilutes cross-half reproducibility (cnga1: 9 of 10 reproducible at 10, 3 of 32 at 32; &
            &TRPM4 kept 10 of 10 at 10 and 2 of 20 at 20) and the fit cost grows as the square of the rank', &
            '# components', .false., 10.0, &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'npreimages', 'num', &
            'Max state volumes (default 4; 8 with preimage_auto)', &
            'Upper bound on the kernel-regression targets in latent space; with the default state_axis=0 &
            &and state_placement=kcenter these are diffusion k-centers over all retained components. The two-gate merge collapses &
            &indistinct states, so the recovered count is <= this. Left blank, it is 4, or 8 with preimage_auto=yes; &
            &a value given with preimage_auto=yes is honoured as the ceiling', &
            'max # states 3-32', .false., 4.0, &
        &visibility=UI_VIS_STANDARD, preserve_default=.true.)
        call flex_pca%add_input(UI_FILT, 'preimage_auto', 'binary', &
            'Determine the state count automatically (default no)', &
            'Raises the state ceiling to 8 (unless npreimages is given) and enables the two-gate merge, &
            &so the delivered state count is recovered from the data rather than requested(yes|no){no}', &
            '', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'rec_states', 'binary', &
            'Reconstruct the state volumes (default yes)', &
            'no delivers the latent embedding, the state labels and the tables in seconds, without &
            &reconstructing any map -- the maps are the dominant cost of a run(yes|no){yes}', &
            '', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'state_axis', 'num', &
            'Latent axis for state targets (default 0 = all components)', &
            'With 0 the state targets are placed over ALL retained covariance components by state_placement. &
            &A negative value places them along a density-spread path instead. A positive value cuts &
            &equal-occupancy slices along that single component, which discards the other components', &
            'component index', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'state_placement', 'multi', &
            'State target placement (default kcenter)', &
            'kcenter places the state targets at diffusion k-centers over all retained components, which &
            &covers a continuous reaction coordinate and branched compositional states with the same &
            &constants. equal_occ cuts a reliability-ordered path through the latent space into slices &
            &holding equal particle counts, so every state gets the same occupancy; use it when the latent &
            &clusters are not well separated. Applies with state_axis=0(kcenter|equal_occ){kcenter}', &
            '', .false., 'kcenter', &
        &choices=ui_choices([character(len=9) :: 'kcenter', 'equal_occ']), &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'nkern', 'num', &
            'Latent components entering the state kernel (default 0 = all)', &
            'Decouples the state stage from neigs. neigs sets how many eigenvolumes are estimated; &
            &nkern sets how many of them define "nearby" for target placement and kernel weighting. &
            &Components past the first few are usually dominated by fitting noise, and each one still &
            &contributes to the Mahalanobis distance, so leaving them in lets noise directions decide &
            &which particles pool at a target and drives every state map toward the consensus. Rank &
            &components by observed spread over posterior variance and keep those above ~1.5', &
            'leading components, 0=all', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'rec_states_backend', 'multi', &
            'Backend for the final state maps (default gridding)', &
            'Independent of rec_backend, which governs the coupled M-step. The two stages want &
            &different estimators: the PCG M-step with a support constraint gives a markedly more &
            &reproducible basis, while PCG state maps have been measured to spread their integrated &
            &density several-fold across states at identical placement, which makes cross-state &
            &difference maps mix real change with amplitude change. Gridding is therefore the default &
            &here even when rec_backend=pcg(gridding|pcg){gridding}', &
            '', .false., 'gridding', &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'infile', 'file', &
            'Cached embedding to resume from', &
            'Path to a flex_pca_embedding.bin written by an earlier run. Skips the covariance &
            &basis fit and the per-particle embedding (~77% of the runtime) and re-enters at the &
            &state-weighting stage, so a different npreimages/min_neff/state_axis can be tried &
            &without refitting. The cache is tied to the particle selection it was built from', &
            'e.g. flex_pca_embedding.bin', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'nbins', 'num', &
            'Cross-validated bandwidth bins (default 1 = off)', &
            'Number of kernel bandwidth bins swept per state. With 1 the &
            &chi-squared floor is used directly, which is the most local and noisiest choice. &
            &With >1 each bin is reconstructed on both halfsets and the bandwidth maximising &
            &even/odd agreement is kept per state. Costs 2*nbins extra reconstruction passes', &
            '# bins, 1=off', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! min_neff is not a flex_pca UI input; params%min_neff (default 2000) still sets the state-pruning
        ! and GMM-island floors of run_flex_pca.
        call flex_pca%add_input(UI_FILT, 'column_separation', 'num', &
            'Minimum grid separation between columns (default 2)', &
            'Selected frequencies closer than this are suppressed; also decorrelates the column noise', &
            'grid units', .false., 2.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_SRCH, 'n_probe_iters', 'num', &
            'Half-set EM iterations (default 4)', &
            'EM iterations refining the covariance basis. Two fits on disjoint particle halves each run &
            &this many iterations; after their merge, one joint EM iteration over all particles follows, &
            &which this count does not include. More iterations raise the likelihood but lowered cross-half &
            &reproducibility on low-SNR data. At least 2, as the half-set merge needs', &
            '# iterations >= 2', .false., 4.0, &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        ! lp/smpd_target/box_crop are derived from the project (derive_flex_pca_sampling); preserve_default
        ! keeps the declared values (lp 16 A = COV_LP_DEFAULT)
        call flex_pca%add_input(UI_FILT, 'lp', 'num', &
            'Variance resolution limit (default 16 A)', &
            'Resolution to which the variability is resolved. It sets the working sampling to lp/2.5 and the &
            &working box from it (magic-box autoscale, never finer than the data, floor 64). With an explicit &
            &box_crop, lp instead defaults to 2.5*smpd_crop', &
            'low-pass limit in Angstroms', .false., 16.0, &
        &group="regularization", visibility=UI_VIS_STANDARD, preserve_default=.true.)
        call flex_pca%add_input(UI_PARM, smpd_target, required_override=.false., &
            label_override='Working sampling distance (expert alias of lp; default lp/2.5)', &
            help_override='Sampling distance of the covariance lattice in Angstroms. Ignored when lp is given; &
            &the default is lp/2.5 = 6.4 A', &
            group="regularization", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call flex_pca%add_input(UI_PARM, 'box_crop', 'num', &
            'Working box size (override; default from lp)', 'Even box used for covariance fitting and &
            &the latent embedding. Normally derived from lp (sampling lp/2.5) through the magic-box autoscale (never finer &
            &than the data, floor 64); set it only to pin a box for tests', &
            'pixels', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'oritype', 'str', &
            'Particle orientation segment', 'Fixed to ptcl3D', 'ptcl3D', .false., 'ptcl3D', &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'umap', 'binary', &
            'UMAP plot of the final embedding', &
            'Write flex_pca_umap.txt: 2-D UMAP coordinates of a bounded subsample of the final all-N latent &
            &embedding, for plotting. Serial and not part of state placement', &
            '(yes|no){yes}', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_MASK, mskdiam, required_override=.false., &
            group="mask", visibility=UI_VIS_STANDARD)
        ! preserve_default: these defaults live in apply_flex_pca_defaults, which default_audit.py does not trace
        call flex_pca%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Backend of the coupled M-step basis solve (the state maps follow rec_states_backend); PCG solves &
        &the same weighted least-squares problem with the support constraint inside the solve(gridding|pcg){pcg}', &
        &'', .false., 'pcg', &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call flex_pca%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Iterations of each warm-started PCG basis solve (every M-step and the merge''s joint solve); &
        &used only when rec_backend=pcg', 'iterations{4}', .false., 4., &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true., &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call flex_pca%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Stop at this true L2 relative residual (a diminishing-returns stop on the update applies too); &
        &<=0, the default, runs exactly maxits_pcg iterations', 'tolerance{0}', &
        &.false., 0.0, visibility=UI_VIS_ADVANCED, preserve_default=.true., &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call flex_pca%add_input(UI_FILT, 'pcg_mskfile', 'file', 'PCG support-constraint mask volume', &
        &'Real-space [0,1] mask volume at the project box installed as the hard support of every PCG solve &
        &(state maps and basis); spherical mskdiam support when absent', &
        &'e.g. automask3D_state01.mrc', .false., '', group="mask", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call flex_pca%add_input(UI_COMP, nparts, required_override=.false., &
            group="compute", visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
call flex_pca%add_input(UI_FILT, 'min_state_frac', 'num', &
            'Minimum state population fraction (default 0 = off)', &
            'Population floor: every delivered state must hold at least this fraction of the embedded &
            &particles. Under-populated clusters are dropped, the targets re-placed on the retained &
            &particles with the count raised by the deficit until npreimages states qualify, and the &
            &dropped or unassigned particles receive a random label among the delivered states, which &
            &are then reconstructed from their hard labels. Incompatible with preimage_auto and the merge', &
            'fraction of particles 0-1', .false., 0.0, &
        &visibility=UI_VIS_STANDARD)
call flex_pca%add_input(UI_FILT, 'heldout', 'binary', &
            'Cross-halfset (held-out) embedding', &
            'Fit the covariance basis on one halfset and embed the other, then swap, so no particle is &
            &projected onto a basis estimated from it; removes in-sample bias and reports the halfset &
            &subspace principal angles. Costs two covariance estimations', &
            '(yes|no){no}', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call add_ui_program('flex_pca', flex_pca, prgtab, UI_CATEGORY)
    end subroutine new_flex_pca

    subroutine new_refine3D_states( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call refine3D_states%new(&
        &'refine3D_states',&
        &'Refine conformational states from a same-lineage particle orientation scaffold',&
        &'refines conformational states from existing particle poses. Input without state labels is split into states by flex PCA; the pose policy searches a local neighborhood or permits global matching',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='Conformational State Refinement')
        call refine3D_states%add_input(UI_SRCH, maxits, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, nstates, required_override=.false., group='search', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_SRCH, 'min_state_frac', 'num', 'Minimum flex state population fraction', &
        &'For state=0/1 input, initialized by flex PCA (nstates>=3), every initial state must hold at least this fraction of the particles; under-populated &
        &flex clusters are dropped and their particles randomized over the delivered states{0.1}', &
        &'fraction of particles', .false., 0.1, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'nsample', 'num', 'Particle sample target', &
        &'Particles sampled per iteration; set 0 to derive the automatic target from the number of states', &
        &'particles (0=automatic)', .false., 0., group='search', visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call refine3D_states%add_input(UI_SRCH, 'pose_policy', 'multi', 'Pose-search policy', &
        &'Pose-search policy(local|global){global}. Local uses the current geometric neighborhood; global permits full probabilistic matching', &
        &'local or global', .false., 'global', group='search', &
        &choices=ui_choices([character(len=6) :: 'local', 'global']), visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_SRCH, 'local_ang_bound', 'num', 'Local projection-angle bound', &
        &'Advanced override for the automatically derived local projection-direction neighborhood in degrees', &
        &'degrees (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, 'local_inpl_bound', 'num', 'Local in-plane angle bound', &
        &'Advanced upper bound for the automatically derived probabilistic in-plane neighborhood in degrees', &
        &'degrees (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, 'local_shift_bound', 'num', 'Local shift bound', &
        &'Advanced override for the automatically derived local translational half-width in pixels', &
        &'pixels (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', &
        &.false., 'no', group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', &
        &'Continue previous refinement(yes|no){no}', '', .false., 'no', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'overlap', 'num', 'State-overlap convergence target', &
        &'Required overlap of state assignments for convergence during frequency marching', 'overlap fraction', &
        &.false., .99, group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(fsc|nonuniform_lpset|none){nonuniform_lpset}', '', .false., 'nonuniform_lpset', group='filter', &
        &choices=ui_choices([character(len=16) :: 'fsc', 'nonuniform_lpset', 'none']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Starting low-pass limit for the common frequency schedule across states', 'low-pass limit in Angstroms', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Final low-pass limit for the common frequency schedule across states', 'low-pass limit in Angstroms', &
        &.false., 6., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, ml_reg, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, envfsc, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, envmsklp, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_MASK, 'automsk', 'multi', 'Refinement envelope mode', &
        &'Use the density envelope, or prefer the lag-one NU-evidence envelope with density fallback(yes|nu|no){no}', &
        &'', .false., 'no', group='mask', choices=ui_choices([character(len=3) :: 'yes', 'nu', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_MASK, nu_msk_sig, group='mask', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('refine3D_states', refine3D_states, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_states

    subroutine new_refine3D_states_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call refine3D_states_cavgs%new(&
        &'refine3D_states_cavgs',&
        &'Refine conformational states of class-expansion sub-class averages',&
        &'sorts the sub-class averages of cls_expansion into conformational states: flex PCA initializes the states on the even and odd class averages, the refine3D_states frequency march refines them with PCG reconstruction, and the states are handed back to the particles',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_DEVELOPER, display_name='Class-Average State Refinement')
        call refine3D_states_cavgs%add_input(UI_SRCH, maxits, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_SRCH, nstates, required_override=.false., group='search', visibility=UI_VIS_STANDARD)
        call refine3D_states_cavgs%add_input(UI_SRCH, 'min_state_frac', 'num', 'Minimum flex state population fraction', &
        &'For class averages without multi-state labels, initialized by flex PCA (nstates>=3), every initial state must hold &
        &at least this fraction of the class-average rows{0.1}', &
        &'fraction of rows', .false., 0.1, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_SRCH, 'pose_policy', 'multi', 'Pose-search policy', &
        &'Pose-search policy(local|global){global}. Local uses the current geometric neighborhood; global permits full probabilistic matching', &
        &'local or global', .false., 'global', group='search', &
        &choices=ui_choices([character(len=6) :: 'local', 'global']), visibility=UI_VIS_STANDARD)
        call refine3D_states_cavgs%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call refine3D_states_cavgs%add_input(UI_SRCH, 'overlap', 'num', 'State-overlap convergence target', &
        &'Required overlap of state assignments for convergence during frequency marching', 'overlap fraction', &
        &.false., .99, group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states_cavgs%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Starting low-pass limit of the frequency march, also the covariance band of the flex PCA initialization', &
        &'low-pass limit in Angstroms', .false., 10., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Final low-pass limit of the frequency march', 'low-pass limit in Angstroms', &
        &.false., 6., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_FILT, ml_reg, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_FILT, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend of the flex PCA initialization, the frequency march and the final maps(pcg|gridding){pcg}', &
        &'', .false., 'pcg', group='filter', choices=ui_choices([character(len=8) :: 'pcg', 'gridding']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_states_cavgs%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call refine3D_states_cavgs%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call refine3D_states_cavgs%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('refine3D_states_cavgs', refine3D_states_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_states_cavgs

    subroutine new_classify3D_refs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call classify3D_refs%new(&
        &'classify3D_refs',&
        &'Classify particles against supplied references and reconstruct the resulting state maps',&
        &'competitively classifies particles against a complete set of references that may have independent provenance, then reconstructs maps from the resulting hard assignments',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='Reference-guided 3D Classification')
        call classify3D_refs%add_input(UI_IMG, 'vol1', 'file', 'Reference volume', &
        &'First member of the complete vol1..volN fixed-reference set', 'input volume e.g. vol1.mrc', .false., 'vol1.mrc', &
        &visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, 'maxits', 'num', 'Maximum iterations', &
        &'Total number of frequency-marched classification iterations; must be >= 1', 'iterations{50}', &
        &.false., 50., group='search', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call classify3D_refs%add_input(UI_SRCH, nstates, required_override=.true., group='search', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, 'nsample', 'num', 'Particle sample target per state', &
        &'Particles sampled per iteration per state', 'particles{10000}', .false., 10000., group='search', &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call classify3D_refs%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) and map shifts back to particles(yes|no){no}', '', .false., 'no', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(nonuniform_lpset|none){nonuniform_lpset}', '', .false., 'nonuniform_lpset', group='filter', &
        &choices=ui_choices([character(len=16) :: 'nonuniform_lpset', 'none']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Starting low-pass limit for frequency-marched classification', 'low-pass limit in Angstroms', .false., 10., &
        &group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Final low-pass limit for frequency-marched classification', 'low-pass limit in Angstroms', .false., 6., &
        &group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, ml_reg, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, envfsc, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, envmsklp, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_MASK, 'automsk', 'multi', 'Refinement envelope mode', &
        &'Use the density envelope, or prefer the lag-one NU-evidence envelope with density fallback(yes|nu|no){no}', &
        &'', .false., 'no', group='mask', choices=ui_choices([character(len=3) :: 'yes', 'nu', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_MASK, nu_msk_sig, group='mask', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('classify3D_refs', classify3D_refs, prgtab, UI_CATEGORY)
    end subroutine new_classify3D_refs

    subroutine new_ptcl3D_state_consensus( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ptcl3D_state_consensus%new(&
        &'ptcl3D_state_consensus', &                                    ! name
        &'Build consensus particle-state assignments across projects',& ! summary
        &'is a program that builds a consensus particle state assignment from a file table of SIMPLE projects &
        &and writes it to the target project ptcl3D field', &           ! help
        &'simple_exec',&                                                ! executable
        &.true., &
        &display_name='Particle State Consensus')                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ptcl3D_state_consensus%add_input(UI_FILE, projtab,&
        &help_override        = 'Text file listing SIMPLE project files (*.simple) containing ptcl3D state assignments',&
        &placeholder_override = 'e.g. projtab.txt',&
        &required_override          = .true.,&
        &group="data", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call ptcl3D_state_consensus%add_input(UI_FILE, projfile,&
        &help_override        = 'Target SIMPLE project file that receives the consensus ptcl3D state assignment',&
        &required_override          = .true.,&
        &group="data", visibility=UI_VIS_STANDARD)
        call ptcl3D_state_consensus%add_input(UI_SRCH, nstates,&
        &help_override        = 'Number of state labels to match; inferred from projtab when omitted',&
        &required_override          = .false.,&
        &group="state", visibility=UI_VIS_ADVANCED)
        call ptcl3D_state_consensus%add_input(UI_PARM, prune,&
        &group="data", visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('ptcl3D_state_consensus', ptcl3D_state_consensus, prgtab, UI_CATEGORY)
    end subroutine new_ptcl3D_state_consensus

    subroutine new_state_tree( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call state_tree%new(&
        &'state_tree',&
        &'Bottom-up tree of states from their even/odd maps',&
        &'is a program that organises the states of a multi-state project bottom-up. For every pair of states the &
        &cross-half difference energy of their even/odd maps inside a common envelope, over a band, is an unbiased &
        &estimate of how different the states are; relative to their signal it orders agglomerative merging. Writes &
        &merge_tree.txt, tree.nwk, the leaf tables and the nested cuts cutN/groups.txt and cutN/states.txt',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='State Tree')
        call state_tree%add_input(UI_FILE, 'dir', 'dir', 'Directory of the state half-maps', &
        &'Directory holding recvol_stateXX_even/odd maps (unfiltered ones when present); by default the directory &
        &of the state maps registered in the project', 'e.g. 1_refine3D_states', .false., '', visibility=UI_VIS_STANDARD)
        call state_tree%add_input(UI_FILT, 'hp', 'num', 'High-pass limit of the band', &
        &'Lowest resolution of the band over which the states are compared{40}', 'in Angstroms{40}', .false., 40., &
        &group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_tree%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit of the band', &
        &'Highest resolution of the band over which the states are compared{10}', 'in Angstroms{10}', .false., 10., &
        &group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_tree%add_input(UI_PARM, 'linkage', 'multi', 'Linkage', &
        &'Cluster distance: the largest (complete) or the particle-weighted mean (average) of the leaf distances, or the &
        &distance of the merged maps (centroid)(complete|average|centroid){complete}', '', .false., 'complete', &
        &choices=ui_choices([character(len=8) :: 'complete', 'average', 'centroid']), visibility=UI_VIS_ADVANCED, &
        &preserve_default=.true.)
        call state_tree%add_input(UI_PARM, 'minpop', 'num', 'Particle budget', &
        &'With a budget, states with fewer particles first merge into their nearest states until every state holds &
        &the budget, then the tree is built over those states (0 = off){0}', 'particles{0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_tree%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('state_tree', state_tree, prgtab, UI_CATEGORY)
    end subroutine new_state_tree

    subroutine new_state_mergeup( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call state_mergeup%new(&
        &'state_mergeup',&
        &'Merge states up their tree by cross-half split tests',&
        &'is a program that decides which states of a multi-state project the particles support. The state tree is &
        &cut into nested groups (N/2, N/4, ... 1); in every group the states are tested together on the even and odd &
        &particles (state_split): states without high-resolution signal are set aside as junk, every pair is tested at &
        &the equal N of its smaller member with group maps docked to the other half, against a label-shuffle null &
        &stratified by view and defocus, then a cross-half EM and a re-test in the band beyond the EM band. States that &
        &do not survive merge by majority. Writes the merged states to the project and to merged_units.txt',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='State Merge-up')
        call state_mergeup%add_input(UI_FILE, 'infile', 'file', 'State tree', &
        &'merge_tree.txt written by state_tree for the states of this project', 'e.g. 2_state_tree/merge_tree.txt', &
        &.true., '', visibility=UI_VIS_STANDARD)
        call state_mergeup%add_input(UI_FILE, 'infile2', 'file', 'Per-state focus masks', &
        &'Text file of state and focus-mask volume pairs; every group is tested inside the union of its states'' masks', &
        &'e.g. leaf_masks.txt', .false., '', visibility=UI_VIS_ADVANCED)
        call state_mergeup%add_input(UI_FILT, 'hp', 'num', 'High-pass limit of the test band', &
        &'Lowest resolution of the band of the split statistic{40}', 'in Angstroms{40}', .false., 40., &
        &group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_mergeup%add_input(UI_FILT, 'lp', 'num', 'Resolution of the test', &
        &'The particles are Fourier-cropped to this Nyquist resolution; ignored with box_crop{10}', 'in Angstroms{10}', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_mergeup%add_input(UI_PARM, 'box_crop', 'num', 'Working box (override)', &
        &'Fourier-cropped box of the test; by default the smallest even box whose Nyquist reaches lp', 'pixels', &
        &.false., 0., visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call state_mergeup%add_input(UI_SRCH, pgrp, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call state_mergeup%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call state_mergeup%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('state_mergeup', state_mergeup, prgtab, UI_CATEGORY)
    end subroutine new_state_mergeup

    subroutine new_state_split( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call state_split%new(&
        &'state_split',&
        &'Test whether the states of a project are supported by its particles',&
        &'is a program that tests the states of a project as one node on the even and odd particles: states without &
        &high-resolution signal are set aside as junk, every pair of states is tested at the equal N of its smaller &
        &member with group maps docked to the other half, against a label-shuffle null stratified by view and defocus; &
        &the merged groups are refined by a cross-half EM and re-tested in the band beyond the EM band. With a focus mask &
        &every map is compared inside it. Writes the supported groups per particle to split_children.txt and the junk &
        &states to split_junk.txt',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_ADVANCED, display_name='State Split Test')
        call state_split%add_input(UI_FILE, 'pcg_mskfile', 'file', 'Focus mask', &
        &'Real-space mask with values in [0,1]: every map is compared only inside it', 'e.g. node_mask.mrc', &
        &.false., '', visibility=UI_VIS_ADVANCED)
        call state_split%add_input(UI_FILT, 'hp', 'num', 'High-pass limit of the test band', &
        &'Lowest resolution of the band of the split statistic{40}', 'in Angstroms{40}', .false., 40., &
        &group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_split%add_input(UI_FILT, 'lp', 'num', 'Resolution of the test', &
        &'The particles are Fourier-cropped to this Nyquist resolution; ignored with box_crop{10}', 'in Angstroms{10}', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call state_split%add_input(UI_PARM, 'box_crop', 'num', 'Working box (override)', &
        &'Fourier-cropped box of the test; by default the smallest even box whose Nyquist reaches lp', 'pixels', &
        &.false., 0., visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call state_split%add_input(UI_SRCH, pgrp, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call state_split%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call state_split%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('state_split', state_split, prgtab, UI_CATEGORY)
    end subroutine new_state_split

    subroutine new_site_rank( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call site_rank%new(&
        &'site_rank',&
        &'Rank the particles toward both poles of a direction',&
        &'is a program that finds the two extremes of a direction of variability at the fixed poses of a project. From &
        &each tail of the per-particle coordinate (min(15000, N/6) particles) a two-class soft EM is run on the even and &
        &odd particles, every particle scored against the other half''s references, the class deviation confined to the &
        &site mask; each pole is sized by the cross-half correlation of its top-N deviation in the band beyond the EM &
        &band. Writes per particle the pole posteriors (site_post_*.txt), the top-N* selections (site_top_*.txt) &
        &and the sizing curves. Every endpoint is also ranked in two disjoint particle halves: a pole counts when its &
        &half-versions beat their best other match by 3 standard errors of the cross-half quarter pairs (an endpoint whose &
        &best other match is one endpoint in both directions joins it) and endpoints merge unless 3 standard errors apart; &
        &a particle is seeded to the accepted pole of highest mean posterior when that posterior reaches 0.5, the others &
        &form a rest state kept only when it passes the same test (site_seeds.txt, site_poles.txt, site_halves.txt)',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_ADVANCED, display_name='Site Ranking')
        call site_rank%add_input(UI_FILE, 'infile', 'file', 'Direction', &
        &'Text file with one coordinate per particle of the project (e.g. a flex latent direction), or with site_dirs &
        &the flex_pca_coordinates.txt of a flex_pca run', 'e.g. z1.txt', .true., '', visibility=UI_VIS_STANDARD)
        call site_rank%add_input(UI_PARM, 'site_dirs', 'str', 'Latent directions', &
        &'Comma-separated latent components of a flex_pca_coordinates.txt infile whose two tails are ranked', &
        &'e.g. 1,3', .false., '', visibility=UI_VIS_ADVANCED)
        call site_rank%add_input(UI_FILE, 'infile2', 'file', 'Extra starts', &
        &'Text file with one more start per line: a name and a file with one value per particle of the project (positive = &
        &in the start); every start''s pole enters the grouping (site_starts.txt)', 'e.g. starts.txt', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call site_rank%add_input(UI_PARM, 'pole_merge', 'num', 'Pole merge line', &
        &'0: every endpoint is also ranked in two disjoint particle halves; a pole counts when its half-versions beat any &
        &other match by 3 standard errors and endpoints merge unless 3 standard errors apart (site_halves.txt). A &
        &Spearman value instead merges endpoints whose rankings correlate at least that much{0}', 'Spearman{0}', .false., &
        &0., visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call site_rank%add_input(UI_PARM, 'junk_screen', 'binary', 'Junk screen', &
        &'Every particle is scored at its pose against every seed state''s reference from the other half (its amplitude in &
        &units of a genuine particle, shells to 12 A); particles of a mixture component of the amplitudes at least 2 sd below &
        &the genuine one are seeded 0 (site_junk.txt, site_junk_report.txt)(yes|no){yes}', '', .false., 'yes', &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call site_rank%add_input(UI_FILE, 'pcg_mskfile', 'file', 'Site mask', &
        &'Real-space mask with values in [0,1] confining the class deviation', 'e.g. site.mrc', &
        &.true., '', visibility=UI_VIS_STANDARD)
        call site_rank%add_input(UI_FILT, 'lp', 'num', 'Resolution of the working box', &
        &'The particles are Fourier-cropped to this Nyquist resolution{10}', 'in Angstroms{10}', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call site_rank%add_input(UI_SRCH, pgrp, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call site_rank%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call site_rank%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('site_rank', site_rank, prgtab, UI_CATEGORY)
    end subroutine new_site_rank

    subroutine new_hetero_tree( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call hetero_tree%new(&
        &'hetero_tree',&
        &'Recursive tree of the heterogeneity of a refined project',&
        &'is a distributed workflow that organizes the heterogeneity of a consensus-refined project as a tree of nodes. &
        &Every node runs its mean, flex_pca, view-balanced cheap children (at the root the flex states) giving a &
        &falloff-prior focus mask, a focused flex_pca, the two-tailed site ranking of the chosen directions of its latent &
        &(site_rank), whose merged poles seed a two-stage focused refine3D with one state per pole plus one for the &
        &particles no pole claims, and the cross-half &
        &split test of the refined states (state_split). Supported groups of at least minpop particles become child nodes below maxdepth, &
        &the others are leaves. The node steps run as queue jobs, nodepar nodes at a time, with one part per 20000 &
        &particles of the node (at least 4, at most nparts) and the site ranking on 32 threads; rerun with mkdir=no in the &
        &run directory to resume. Writes hetero_tree_labels.txt, hetero_tree_leaves.txt and the tree of the leaves',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_ADVANCED, display_name='Heterogeneity Tree')
        call hetero_tree%add_input(UI_FILE, 'infile', 'file', 'Root particles', &
        &'One value per particle of the project: the rows with a positive value form the root node, a node of an earlier tree at its &
        &poses, masked like an inner node{all particles with a state}', &
        &'e.g. root_rows.txt', .false., '', visibility=UI_VIS_ADVANCED)
        call hetero_tree%add_input(UI_PARM, 'kfloor', 'num', 'Particles per flex state', &
        &'flex_pca states per node, which build its focus mask: the node particles divided by kfloor, between 2 and 8{5000}', &
        &'particles{5000}', .false., 5000., visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call hetero_tree%add_input(UI_PARM, 'site_dirs', 'str', 'Pole directions', &
        &'Comma-separated components of every node''s focused flex latent whose two tails start its site ranking{1,3}', &
        &'e.g. 1,3', .false., '1,3', visibility=UI_VIS_ADVANCED)
        call hetero_tree%add_input(UI_PARM, 'site_disc', 'binary', 'Discover every node''s site', &
        &'Every node first refines its focused flex states (stage A0), discovers the sites that vary between them and ranks &
        &the poles of the site carrying the largest share of the variance, instead of using pcg_mskfile(yes|no){yes}', '', &
        &.false., 'yes', visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call hetero_tree%add_input(UI_PARM, 'pole_merge', 'num', 'Pole merge line', &
        &'0: a node''s poles are accepted and merged by their versions in two disjoint particle halves; a Spearman value &
        &instead merges endpoints whose rankings correlate at least that much{0}', 'Spearman{0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call hetero_tree%add_input(UI_PARM, 'junk_screen', 'binary', 'Junk screen at every node', &
        &'At every node''s seeding, particles that no seed state explains get label 0 at that node and below &
        &(site/site_junk_report.txt)(yes|no){yes}', '', .false., 'yes', visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call hetero_tree%add_input(UI_FILE, 'pcg_mskfile', 'file', 'Site mask', &
        &'Real-space mask with values in [0,1] confining the pole classes of every node''s site ranking{the node''s &
        &focus mask}', 'e.g. site.mrc', .false., '', visibility=UI_VIS_ADVANCED)
        call hetero_tree%add_input(UI_PARM, 'minpop', 'num', 'Smallest child node', &
        &'A supported group with at least this many particles becomes a child node{15000}', 'particles{15000}', .false., &
        &15000., visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call hetero_tree%add_input(UI_PARM, 'maxdepth', 'num', 'Maximum depth', &
        &'Nodes at this depth are not split further{4}', 'depth{4}', .false., 4., visibility=UI_VIS_ADVANCED, &
        &preserve_default=.true.)
        call hetero_tree%add_input(UI_PARM, 'nodepar', 'num', 'Concurrent nodes', &
        &'Nodes whose steps run at the same time{2}', 'nodes{2}', .false., 2., visibility=UI_VIS_ADVANCED, &
        &preserve_default=.true.)
        call hetero_tree%add_input(UI_SRCH, pgrp, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call hetero_tree%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call hetero_tree%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call hetero_tree%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('hetero_tree', hetero_tree, prgtab, UI_CATEGORY)
    end subroutine new_hetero_tree

    subroutine new_landscape_compare( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call landscape_compare%new(&
        &'landscape_compare',&
        &'Reproducibility of a conformational landscape across runs',&
        &'is a program that compares the state maps of independent runs on one node. Inside the mask every map is low-&
        &passed and taken as its deviation from its run''s fraction-weighted mean, with its half-map reliability as the &
        &noise ceiling. Reports how well each map is explained by every other run''s maps, the cross-run deviation FSC &
        &of matched maps, the principal angles of the runs'' deviation subspaces, the occupancy along the pooled &
        &principal components, the extreme maps at both ends of those components matched across runs, and optionally &
        &how the pole rankings of the runs group. Writes landscape_compare.txt',&
        &'simple_exec',&
        &.false.,&
        &visibility=UI_VIS_ADVANCED, display_name='Landscape Compare')
        call landscape_compare%add_input(UI_FILE, 'infile', 'file', 'Runs', &
        &'Text file with one run per line: a name, the directory of its recvol_stateNN_even/odd_unfil.mrc maps and its &
        &per-particle state labels file (absolute paths; every run on the same particles and box)', 'e.g. runs.txt', &
        &.true., '', visibility=UI_VIS_STANDARD)
        call landscape_compare%add_input(UI_FILE, 'pcg_mskfile', 'file', 'Node mask', &
        &'Real-space mask with values in [0,1] inside which the maps are compared', 'e.g. node_mask.mrc', &
        &.true., '', visibility=UI_VIS_STANDARD)
        call landscape_compare%add_input(UI_FILE, 'infile2', 'file', 'Pole rankings', &
        &'Text file with one pole per line: a name and its per-particle posterior file (site_post_lo/hi.txt of &
        &site_rank, negative = not ranked)', 'e.g. poles.txt', .false., '', visibility=UI_VIS_ADVANCED)
        call landscape_compare%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit', &
        &'Gaussian low-pass applied to every map before the comparison{15}', 'in Angstroms{15}', .false., 15., &
        &group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call landscape_compare%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('landscape_compare', landscape_compare, prgtab, UI_CATEGORY)
    end subroutine new_landscape_compare

    subroutine new_site_discover( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call site_discover%new(&
        &'site_discover',&
        &'Discover the sites that vary between the states of a node',&
        &'is a program that finds the sites of a set of states from their even/odd maps: the cross-half variance of the &
        &low-passed, scaled state deviations, its core split into connected components and, where the voxels of a &
        &component vary independently across the states, into two factors; every candidate gets a falloff-prior mask and &
        &the cross-half reproducibility of the state deviations inside it. Writes the sites ranked by reproducibility &
        &(site_NN.mrc, sites.txt)',&
        &'simple_exec',&
        &.false.,&
        &visibility=UI_VIS_ADVANCED, display_name='Site Discovery')
        call site_discover%add_input(UI_FILE, 'dir', 'dir', 'Directory of the state half-maps', &
        &'Directory holding the recvol_stateNN_even/odd_unfil.mrc maps of the states', 'e.g. recB', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call site_discover%add_input(UI_FILE, 'infile', 'file', 'State labels', &
        &'Text file with the state of every particle (0 = none); it weighs the states', 'e.g. states.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call site_discover%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call site_discover%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('site_discover', site_discover, prgtab, UI_CATEGORY)
    end subroutine new_site_discover

    subroutine new_landscape_sample( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call landscape_sample%new(&
        &'landscape_sample',&
        &'Sample the landscape between the pole pairs of a node',&
        &'is a program that samples the path between every accepted pair of poles of a node at its stage-B poses: each &
        &particle is scored against the even/odd state maps (the other half''s map, 12 A band, two sets of Fourier pixels &
        &split by the parity of h+k); the pair coordinate is the score difference of the two poles. The number of samples is &
        &the largest candidate whose neighbour steps are resolved on the pixels the coordinate did not use, in both folds; &
        &the rest state never ends a pair but its particles enter every pair. Writes per pair the samples &
        &(landscape_A_B_samples.txt), the coordinates and the held-out and final sample maps; landscape_pairs.txt lists the pairs',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_ADVANCED, display_name='Landscape Sampling')
        call landscape_sample%add_input(UI_FILE, 'infile', 'file', 'Pole list', &
        &'site_poles.txt of the node: one line per state in state order, the line holding rest marks the rest state', &
        &'e.g. site_poles.txt', .true., '', visibility=UI_VIS_STANDARD)
        call landscape_sample%add_input(UI_FILE, 'infile2', 'file', 'Split labels', &
        &'split_children.txt of the node: poles in one child count as one, rows with child 0 are left out', &
        &'e.g. split_children.txt', .false., '', visibility=UI_VIS_ADVANCED)
        call landscape_sample%add_input(UI_FILE, 'dir', 'dir', 'Directory of the state half-maps', &
        &'Directory holding the recvol_stateNN_even/odd_unfil.mrc maps of the states', 'e.g. stageB', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call landscape_sample%add_input(UI_FILT, 'lp', 'num', 'Resolution of the working box', &
        &'The particles are Fourier-cropped to this Nyquist resolution{10}', 'in Angstroms{10}', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call landscape_sample%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call landscape_sample%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('landscape_sample', landscape_sample, prgtab, UI_CATEGORY)
    end subroutine new_landscape_sample

end module simple_ui_heterogeneity
