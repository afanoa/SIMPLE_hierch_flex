!@descr: supporting 3D orientation search
module simple_commanders_refine3D
use simple_commanders_api
use simple_pftc_srch_api
use simple_refine3D_fnames,                        only: refine3D_state_vol_fname, refine3D_fsc_fname
use simple_refine3D_stage_plan,                    only: refine3D_stage_plan_entry, plan_refine3D_frequency_stages
use simple_external_reference_pose_initialization, only: initialize_poses_against_external_references
use simple_gui_communicator,                       only: gui_communicator
implicit none
private :: exec_refine3D_workflow, prepare_refine3D_state_work_project
private :: configure_refine3D_sampling, configure_refine3D_grid, seed_refine3D_auto_nonuniform_lpset
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_refine3D_auto
  contains
    procedure :: execute      => exec_refine3D_auto
end type commander_refine3D_auto

type, extends(commander_base) :: commander_refine3D_local
  contains
    procedure :: execute      => exec_refine3D_local
end type commander_refine3D_local

type, extends(commander_base) :: commander_refine3D_states
  contains
    procedure :: execute      => exec_refine3D_states
end type commander_refine3D_states

type, extends(commander_base) :: commander_refine3D_states_cavgs
  contains
    procedure :: execute      => exec_refine3D_states_cavgs
end type commander_refine3D_states_cavgs

type, extends(commander_base) :: commander_classify3D_refs
  contains
    procedure :: execute      => exec_classify3D_refs
end type commander_classify3D_refs

type, extends(commander_base) :: commander_refine3D
  contains
    procedure :: execute      => exec_refine3D
end type commander_refine3D

type, extends(commander_base) :: commander_refine3D_distr_worker
  contains
    procedure :: execute      => exec_refine3D_distr_worker
end type commander_refine3D_distr_worker

type, extends(commander_base) :: commander_bootstrap_rec3D
  contains
    procedure :: execute      => exec_bootstrap_rec3D
end type commander_bootstrap_rec3D

type, extends(commander_base) :: commander_state_tree
  contains
    procedure :: execute      => exec_state_tree
end type commander_state_tree

type, extends(commander_base) :: commander_state_mergeup
  contains
    procedure :: execute      => exec_state_mergeup
end type commander_state_mergeup

type, extends(commander_base) :: commander_state_split
  contains
    procedure :: execute      => exec_state_split
end type commander_state_split

type, extends(commander_base) :: commander_site_rank
  contains
    procedure :: execute      => exec_site_rank
end type commander_site_rank

type, extends(commander_base) :: commander_site_split
  contains
    procedure :: execute      => exec_site_split
end type commander_site_split

type, extends(commander_base) :: commander_hetero_tree
  contains
    procedure :: execute      => exec_hetero_tree
end type commander_hetero_tree

type, extends(commander_base) :: commander_landscape_compare
  contains
    procedure :: execute      => exec_landscape_compare
end type commander_landscape_compare

type, extends(commander_base) :: commander_site_discover
  contains
    procedure :: execute      => exec_site_discover
end type commander_site_discover

type, extends(commander_base) :: commander_landscape_sample
  contains
    procedure :: execute      => exec_landscape_sample
end type commander_landscape_sample

contains

    subroutine exec_refine3D_auto( self, cline )
        class(commander_refine3D_auto), intent(inout) :: self
        class(cmdline), intent(inout)                 :: cline
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        call exec_refine3D_workflow(cline, .false.)
    end subroutine exec_refine3D_auto

    subroutine exec_refine3D_local( self, cline )
        class(commander_refine3D_local), intent(inout) :: self
        class(cmdline), intent(inout)                  :: cline
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'no')
        call exec_refine3D_workflow(cline, .true.)
    end subroutine exec_refine3D_local

    ! The commander owns workflow sequencing, sampling and final reconstruction.
    ! The controller prepares startup state on the caller's project.
    subroutine exec_refine3D_workflow( cline, l_local )
        use simple_refine3D_workflow_controller, only: refine3D_workflow_controller
        use simple_final_rec,                    only: calc_final_rec
        use simple_commanders_rec,               only: commander_rec3D
        use simple_commanders_euclid,            only: commander_calc_pspec
        use simple_gui_communicator,             only: gui_communicator
        use simple_refine3D_strategy,            only: strip_refine3D_search_only_args
        use simple_strategy3D_cont,              only: check_cont_seeds
        use simple_refine3D_fnames,              only: refine3D_state_vol_fname
        use simple_defs,                         only: logfhandle
        use simple_parameters,                   only: parameters
        use simple_sp_project,                   only: sp_project
        use simple_error,                        only: simple_exception
        class(cmdline), intent(inout) :: cline
        logical, intent(in)           :: l_local
        type(refine3D_workflow_controller) :: controller
        type(commander_refine3D)           :: xrefine3D
        type(commander_bootstrap_rec3D)    :: xbootstrap_rec3D
        type(commander_rec3D)              :: xrec3D
        type(commander_calc_pspec)         :: xcalc_pspec
        type(gui_communicator)             :: gui_comm
        class(parameters), allocatable     :: params
        class(sp_project), allocatable     :: spproj
        class(cmdline), allocatable        :: cline_boot
        type(string)                       :: pose_cont_arg, init_vol, pose_init_refs(1), pose_init_checkpoint(1)
        integer                            :: maxits_user, nactive
        logical                            :: l_polish, l_have_init_vol, l_external_input, l_maxits_defined
        integer                            :: box_crop
        real                               :: smpd_crop
        character(len=14)                  :: workflow_label

        workflow_label = 'REFINE3D_AUTO'
        if( l_local ) workflow_label = 'REFINE3D_LOCAL'
        l_external_input = cline%defined('vol1')
        l_maxits_defined = cline%defined('maxits')
        l_polish = .false.
        if( cline%defined('refine') ) THROW_HARD('The workflow derives refine; remove the refine argument')
        if( l_local )then
            call cline%set('prg', 'refine3D_local')
            if( cline%defined('pose_cont') ) THROW_HARD('refine3D_local selects continuous refinement; remove pose_cont')
            if( cline%defined('nspace') .or. cline%defined('nspace_sub') ) &
                THROW_HARD('refine3D_local uses no projection grid; remove nspace and nspace_sub')
            call require_refine3D_local_value(cline, 'ref_pose_init', 'none')
            call require_refine3D_local_value(cline, 'regpass', 'no')
            call require_refine3D_local_value(cline, 'objfun', 'euclid')
            call require_refine3D_local_value(cline, 'ml_reg', 'yes')
            call require_refine3D_local_value(cline, 'center', 'no')
            call require_refine3D_local_value(cline, 'combine_eo', 'no')
            if( cline%defined('lp') ) THROW_HARD('refine3D_local requires half-specific references; use lpstop instead of lp')
            if( cline%defined('filt_mode') )then
                if( cline%get_carg('filt_mode') == 'nonuniform_lpset' ) &
                    THROW_HARD('refine3D_local does not support filt_mode=nonuniform_lpset')
            endif
            if( cline%defined('nstates') )then
                if( cline%get_iarg('nstates') /= 1 ) THROW_HARD('refine3D_local refines one state; use state=X')
            endif
            if( cline%defined('trail_seed') ) THROW_HARD('trail_seed is reserved for refine3D_local')
            call cline%set('regpass', 'no')
        else
            if( .not. cline%defined('pose_cont') ) call cline%set('pose_cont', 'no')
            pose_cont_arg = cline%get_carg('pose_cont')
            select case(pose_cont_arg%to_char())
                case('no')
                case('yes')
                    l_polish = .true.
                case('only')
                    THROW_HARD('pose_cont=only is retired; use prg=refine3D_local without pose_cont')
                case DEFAULT
                    THROW_HARD('refine3D_auto supports pose_cont=no|yes')
            end select
        endif

        ! Set workflow defaults before the single typed parameter construction.
        call cline%set('balance',       'none') ! global lowest-updatecnt sampling, no sampling units
        call cline%set('greedy_sampling', 'no') ! only active with sampling units (balance=class|cavg)
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        if( l_local )then
            call cline%set('refine',        'cont') ! continuous Cartesian refinement of the stored poses
        else
            call cline%set('refine',  'prob_neigh') ! probabilistic neighborhood 3D refinement
        endif
        call cline%set('ml_reg',         'yes') ! ML regularization is on
        call cline%set('overlap',         0.99) ! discrete convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('incrreslim',      'no')

        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',            'yes')
        ! Density-envelope/phase-randomization-corrected FSC steers the
        ! matching LP, NU gating, and convergence.
        if( .not. cline%defined('envfsc')      ) call cline%set('envfsc', 'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',            'no')
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',        'no')
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',        'yes')
        if( .not. cline%defined('nsample')     ) call cline%set('nsample',          25000)
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode', 'nonuniform')
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',          'yes') ! automatic density-envelope masking

        ! Explicit maxits also sets minits; otherwise allow at least three
        ! iterations for the startup model to be superseded.
        if( l_maxits_defined )then
            maxits_user = cline%get_iarg('maxits')
            if( maxits_user < 1 ) THROW_HARD('maxits must be >= 1 for '//workflow_label)
            call cline%set('minits', maxits_user)
        else if( cline%defined('minits') )then
            call cline%set('minits', max(3, cline%get_iarg('minits')))
        else
            call cline%set('minits', 3)
        endif
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol', 'no') ! we do not keep volumes for each iteration by default

        if( .not. cline%defined('state') .and. cline%defined('m_estimator') )then
            if( cline%get_carg('m_estimator') == 'flex' ) &
                THROW_HARD(workflow_label//' m_estimator=flex refines one state of the project''s state weight set; give state=')
        endif
        call controller%new(l_local, l_external_input)
        ! Validate first: params%new can create missing halves from the merged map.
        call controller%preflight(cline)
        allocate(params, spproj, cline_boot)
        call params%new(cline)
        ! pose_cont is this workflow's switch; children receive the derived refine mode.
        call cline%delete('pose_cont')

        ! Run state selection before this workflow's GUI communicator starts.
        ! The selection child must own the process communicator while it runs.
        if( cline%defined('state') )then
            call prepare_refine3D_state_work_project(params, cline, spproj, workflow_label)
            call controller%publish_project(params, spproj)
        else
            call controller%load_project(params, spproj)
        endif
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        if( trim(params%ref_pose_init) == 'cc' .and. .not. l_external_input ) &
            THROW_HARD(workflow_label//' ref_pose_init=cc requires an external vol1 input')
        ! Check stored poses before any startup particle work.
        if( l_local ) call check_cont_seeds(spproj%os_ptcl3D)
        call configure_refine3D_sampling(params, cline, spproj, workflow_label, l_maxits_defined)
        call configure_refine3D_grid(params, cline, l_local, box_crop, smpd_crop)
        call controller%prepare_startup(params, cline, spproj, init_vol, l_have_init_vol, box_crop, smpd_crop)
        call gui_comm%new(params)

        if( trim(params%ref_pose_init) == 'cc' )then
            nactive = spproj%os_ptcl3D%count_state_gt_zero()
            if( nactive < 1 ) THROW_HARD(workflow_label//' pose initialization found no active particles')
            pose_init_refs(1) = init_vol
            call initialize_poses_against_external_references(params, cline, xrefine3D, xrec3D, &
                &nactive, pose_init_refs, pose_init_checkpoint, 1)
            init_vol = pose_init_checkpoint(1)
            call cline%set('startit',    2)
            call cline%set('which_iter', 2)
            call cline%set('extr_iter',  2)
            write(logfhandle,'(A)') '>>> REFINE3D_AUTO POSE INITIALIZATION CHECKPOINT RECONSTRUCTED; ENTERING EUCLIDEAN REFINEMENT'
        endif

        ! Auto builds startup references from particles; local uses its staged saved pair.
        if( .not. l_local )then
            cline_boot = cline
            call strip_refine3D_search_only_args(cline_boot)
            call cline_boot%set('prg', 'calc_pspec')
            call cline_boot%set('which_iter', 1)
            call xcalc_pspec%execute(cline_boot)
            cline_boot = cline
            call cline_boot%set('prg', 'reconstruct3D')
            call cline_boot%delete('trail_rec')
            call cline_boot%delete('update_frac')
            call cline_boot%delete('endit')
            ! Alignment-loop controls have no meaning for a reconstruction.
            call cline_boot%delete('refine')
            call cline_boot%delete('maxits')
            call cline_boot%delete('minits')
            call cline_boot%delete('continue')
            call cline_boot%set('postprocess', 'no')
            call cline_boot%set('which_iter', 1)
            ! The bootstrap pair is the lag-one reference of the first auto iteration.
            if( l_have_init_vol .and. .not. cline_boot%defined('vol1') ) call cline_boot%set('vol1', init_vol)
            call xrec3D%execute(cline_boot)
            call cline%set('vol1', refine3D_state_vol_fname(1))
            call cline_boot%kill
        endif

        call seed_refine3D_auto_nonuniform_lpset(params, cline, workflow_label)
        ! Global registration against the masked startup references.
        call run_refine3D_registration_pass(params, cline, xrefine3D, workflow_label)
        call cline%set('prg',                   'refine3D')
        if( l_local )then
            ! Continuous Cartesian refinement: no projection grid.
            call cline%set('refine',              'cont')
        else
            ! Neighbourhood refinement uses the main direction set.
            call cline%set('refine',        'prob_neigh')
            if( .not. cline%defined('nspace')     ) call cline%set('nspace',     20000)
            if( .not. cline%defined('nspace_sub') ) call cline%set('nspace_sub', 500)
            if( l_polish ) call cline%set('pose_cont', 'yes')
        endif
        call cline%set('maxits',             params%maxits)
        call xrefine3D%execute(cline)

        ! Reconstruct all eligible particles at native sampling. The shared helper
        ! owns sigma fallback, project registration and final products.
        call calc_final_rec(params, spproj, params%projfile, cline, xrec3D, xbootstrap_rec3D, &
            &l_postprocess=.true., lp_snapshot=params%res_target)
        call spproj%read_segment('cls3D',  params%projfile)
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj%read_segment('out', params%projfile)
        call gui_comm%add_metadata(spproj, oritype='cls3D', stage=0, selection=.true.)
        call controller%kill
        call spproj%kill
        call init_vol%kill
        call gui_comm%kill()
    end subroutine exec_refine3D_workflow


    subroutine configure_refine3D_grid( params, cline, l_local, box_crop, smpd_crop )
        use simple_magic_boxes, only: autoscale
        use simple_parameters,  only: parameters
        use simple_cmdline,     only: cmdline
        use simple_defs_conv,   only: AHELIX_WIDTH
        class(parameters), intent(in) :: params
        class(cmdline), intent(inout) :: cline
        logical, intent(in)           :: l_local
        integer, intent(out)          :: box_crop
        real, intent(out)             :: smpd_crop
        real               :: smpd_target, scale, trslim
        logical            :: l_autoscale
        logical, parameter :: DEBUG = .true.
        integer, parameter :: MINBOX = 256
        smpd_target = 1.3
        if( cline%defined('box_crop') )then
            box_crop    = params%box_crop
            smpd_crop   = params%smpd_crop
            smpd_target = smpd_crop
            scale       = real(box_crop) / real(params%box)
            l_autoscale = box_crop < params%box
        else if( .not. params%l_autoscale .or. params%box <= MINBOX )then
            smpd_target = params%smpd
            smpd_crop   = params%smpd
            box_crop    = params%box
            scale = 1.
            l_autoscale = .false.
        else
            call autoscale(params%box, params%smpd, smpd_target, box_crop, smpd_crop, scale, minbox=MINBOX)
            l_autoscale = box_crop < params%box
        endif
        trslim = min(8., max(2., AHELIX_WIDTH / smpd_crop))
        if( l_local .and. cline%defined('trs') ) trslim = params%trs
        if( DEBUG )then
            print *, 'smpd_target: ', smpd_target
            print *, 'box:         ', params%box
            print *, 'box_crop:    ', box_crop
            print *, 'smpd:        ', params%smpd
            print *, 'smpd_crop:   ', smpd_crop
            print *, 'scale:       ', scale
            print *, 'trslim:      ', trslim
            print *, 'l_autoscale: ', l_autoscale
        endif
        call cline%set('trs', trslim)
        if( l_autoscale )then
            call cline%set('box_crop',  box_crop)
        else
            call cline%delete('box_crop')
        endif
        call cline%delete('smpd_crop')
    end subroutine configure_refine3D_grid


    ! Refine state X on a separate single-state copy; never write the parent.
    ! Flex selects positive applied-weight rows and remaps the weight set.
    ! The selected state's map and FSC become state 1's references.
    subroutine prepare_refine3D_state_work_project( params, cline, work, workflow_label )
        use simple_refine3D_strategy,       only: select_refine3D_state_rows
        use simple_syslib,                  only: file_exists, del_file, simple_abspath
        use simple_fileio,                  only: simple_copy_file, get_fbody, fopen, fclose
        use simple_string_utils,            only: int2str_pad
        use simple_commanders_project_core, only: commander_selection
        use simple_state_weight_set,        only: state_weight_set
        use simple_parameters,              only: parameters
        use simple_sp_project,              only: sp_project
        use simple_cmdline,                 only: cmdline
        use simple_defs,                    only: logfhandle
        use simple_error,                   only: simple_exception
        character(len=*), intent(in)     :: workflow_label
        class(parameters), intent(inout) :: params
        class(cmdline), intent(inout)    :: cline
        class(sp_project), intent(inout) :: work
        type(commander_selection) :: xselection
        type(cmdline)             :: cline_sel
        type(sp_project)          :: parent
        type(state_weight_set)    :: parent_set, work_set
        type(string)              :: src_projfile, work_projfile, work_projname, selfile
        type(string)              :: vol_fname, fsc_fname, pstk, wstk
        integer,      allocatable :: parent_rows(:)
        logical,      allocatable :: l_sel(:)
        real                      :: vol_smpd
        integer                   :: state, nparent, nwork, i, funit, vol_box, fsc_box, nstates_out, pstkind, wstkind, pind, wind
        logical                   :: l_vol, l_fsc
        state = params%state
        if( state < 1 ) THROW_HARD(workflow_label//' state= must be at least 1')
        ! the copy and its single-state set live in this run's own directory, beside a copy of the parent
        ! that is removed below: never in the parent's directory, nor beside another state's work project
        if( trim(params%mkdir) /= 'yes' )then
            THROW_HARD(workflow_label//' state= runs in its own directory: mkdir=yes')
        endif
        src_projfile  = params%projfile
        work_projfile = params%prg%to_char()//'_state'//int2str_pad(state,2)//'.simple'
        work_projname = get_fbody(work_projfile, 'simple')
        ! the parent opens under its own name: its state weight set is tied to it
        call parent%read(src_projfile)
        nparent = parent%os_ptcl3D%get_noris()
        if( nparent < 1 .or. parent%is_virgin_field('ptcl3D') ) &
            &THROW_HARD(workflow_label//' state= requires existing ptcl3D orientations')
        if( file_exists(work_projfile) ) call del_file(work_projfile)
        call simple_copy_file(src_projfile, work_projfile)
        call cline_sel%set('prg',      'selection')
        call cline_sel%set('projfile', work_projfile)
        call cline_sel%set('projname', work_projname)
        call cline_sel%set('oritype',  'ptcl3D')
        call cline_sel%set('prune',    'yes')
        call cline_sel%set('append',   'no')
        call cline_sel%set('mkdir',    'no')
        call select_refine3D_state_rows(parent, state, params%l_m_estimator_flex, l_sel, parent_set)
        ! the rows, as the selection program's per-row selection file
        selfile = params%prg%to_char()//'_state'//int2str_pad(state,2)//'_rows.txt'
        call fopen(funit, file=selfile, status='REPLACE', action='WRITE')
        do i = 1, nparent
            write(funit,'(I1)') merge(1, 0, l_sel(i))
        enddo
        call fclose(funit)
        call cline_sel%set('infile', selfile)
        call xselection%execute(cline_sel)
        call work%read(work_projfile)
        nwork = work%os_ptcl3D%get_noris()
        if( nwork < 1 ) THROW_HARD(workflow_label//' the work project selected no particle')
        ! work row i is the i-th selected parent row: check it by stack and index in stack
        parent_rows = pack([(i, i = 1, nparent)], l_sel)
        if( size(parent_rows) /= nwork )then
            THROW_HARD(workflow_label//' work project rows disagree with the selection')
        endif
        do i = 1, nwork
            call parent%map_ptcl_ind2stk_ind('ptcl3D', parent_rows(i), pstkind, pind)
            call work%map_ptcl_ind2stk_ind('ptcl3D', i, wstkind, wind)
            pstk = parent%os_stk%get_str(pstkind, 'stk')
            wstk = work%os_stk%get_str(wstkind, 'stk')
            if( pind /= wind .or. pstk /= wstk ) &
                &THROW_HARD(workflow_label//' a work project row is not the parent row it was selected from')
        enddo
        ! state X's map and FSC (the parent's: the selection drops the entries of states none of its rows
        ! is labelled with) become state 1's; every other state's map and FSC go
        l_vol = parent%isthere_in_osout('vol', state)
        if( l_vol ) call parent%get_vol('vol', state, vol_fname, vol_smpd, vol_box)
        if( l_vol ) l_vol = file_exists(vol_fname)
        l_fsc = parent%isthere_in_osout('fsc', state)
        if( l_fsc ) call parent%get_fsc(state, fsc_fname, fsc_box)
        if( l_fsc ) l_fsc = file_exists(fsc_fname)
        nstates_out = 0
        if( work%os_out%get_noris() > 0 ) nstates_out = work%os_out%get_n('state')
        do i = 1, max(state, nstates_out)
            call work%remove_state_artifacts_from_osout(i)
        enddo
        if( l_vol ) call work%add_vol2os_out(vol_fname, vol_smpd, 1, 'vol', box=vol_box)
        if( l_fsc ) call work%add_fsc2os_out(fsc_fname, 1, fsc_box)
        ! the parent's set does not describe the copy's rows; the copy gets its own below
        call work%remove_state_weights_from_osout
        ! the copy's sigmas are seeded in this run, as in solve3D's state continuation
        if( work%projinfo%get_noris() == 1 )then
            if( work%projinfo%isthere(1, 'sigma2_state') ) call work%projinfo%delete_entry('sigma2_state')
        endif
        call work%update_projinfo(work_projfile)
        call work%write(work_projfile)
        if( params%l_m_estimator_flex )then
            call work_set%publish_work_state(parent_set, state, parent_rows, work, work%os_ptcl3D, work_projfile)
            write(logfhandle,'(A,I0,A,I0,A,F12.3,A,F12.3,A,I0)') '>>> '//workflow_label//' STATE ', state, &
                &' WEIGHT SET: particles ', nwork, ', applied mass ', work_set%get_mass(1), &
                &', effective sample size ', work_set%get_ess(1), ', hard population ', work_set%get_pop(1)
            call work_set%kill
            call parent_set%kill
        endif
        ! the run continues on the copy, whose only state is 1; this directory's copy of the parent goes,
        ! so that the work project is the run directory's only project
        call del_file(src_projfile)
        call cline%set('projfile', work_projfile)
        call cline%set('projname', work_projname)
        params%projfile = simple_abspath(work_projfile)
        params%projname = work_projname
        params%nptcls   = nwork
        call cline%delete('state')
        params%l_state_defined = .false.
        params%state           = 1
        write(logfhandle,'(A,I0,A,I0,A,A)') '>>> '//workflow_label//' STATE/PARTICLES: ', state, ' / ', nwork, &
            &' WORK PROJECT: ', params%projfile%to_char()
        call parent%kill
        call cline_sel%kill
    end subroutine prepare_refine3D_state_work_project

    subroutine configure_refine3D_sampling( params, cline, project, workflow_label, l_maxits_defined )
        use simple_parameters, only: parameters
        use simple_sp_project, only: sp_project
        use simple_cmdline,    only: cmdline
        use simple_defs,       only: logfhandle
        use simple_error,      only: simple_exception
        character(len=*), intent(in)     :: workflow_label
        logical, intent(in)              :: l_maxits_defined
        class(parameters), intent(inout) :: params
        class(cmdline), intent(inout)    :: cline
        class(sp_project), intent(inout) :: project
        integer            :: maxits_auto, nptcls_per_iter, nsample_target, nptcls_eff
        real               :: update_frac_auto
        integer, parameter :: MAXITS_REFINE3D_CAP = 50
        real, parameter    :: TARGET_UPDATES_PER_PARTICLE_REFINE3D = 4.0
        nsample_target = params%nsample
        if( nsample_target < 1 ) THROW_HARD('nsample must be >= 1 for '//workflow_label)
        nptcls_eff = project%count_state_gt_zero()
        if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//workflow_label)
        nptcls_per_iter = min(nptcls_eff, nsample_target)
        if( nptcls_eff <= nsample_target )then
            params%update_frac   = 1.0
            params%l_update_frac = .false.
            params%l_trail_rec   = .false.
            call cline%delete('update_frac')
            write(logfhandle,'(A,I0,A,I0,A)') &
                &'>>> '//workflow_label//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
        else
            update_frac_auto = real(nsample_target) / real(nptcls_eff)
            if( update_frac_auto <= 0.99 )then
                params%update_frac   = update_frac_auto
                params%l_update_frac = .true.
                params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                call cline%set('update_frac', update_frac_auto)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//workflow_label//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                    &nptcls_eff, '/', nsample_target, '/', update_frac_auto
            else
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)') &
                    &'>>> '//workflow_label//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
            endif
        endif
        if( .not. l_maxits_defined )then
            ! target updates x the iterations one sweep of the global lowest-updatecnt draw takes
            maxits_auto = ceiling(TARGET_UPDATES_PER_PARTICLE_REFINE3D * &
                &real(ceiling(real(nptcls_eff) / real(nptcls_per_iter))))
            maxits_auto = max(params%minits, min(MAXITS_REFINE3D_CAP, max(2, maxits_auto)))
            params%maxits = maxits_auto
            call cline%set('maxits', params%maxits)
            write(logfhandle,'(A,I0,A,F5.1,A,I0,A)') '>>> '//workflow_label//' MAXITS: ', &
                &params%maxits, ' FOR ~', TARGET_UPDATES_PER_PARTICLE_REFINE3D, &
                &' UPDATES/PARTICLE (MINIMUM: ', params%minits, ')'
        else
            write(logfhandle,'(A,I0)') &
                &'>>> '//workflow_label//' MAXITS COMMAND-LINE OVERRIDE: ', params%maxits
        endif
    end subroutine configure_refine3D_sampling

    subroutine seed_refine3D_auto_nonuniform_lpset( params, cline, workflow_label )
        use simple_refine3D_fnames, only: refine3D_fsc_fname
        use simple_syslib,          only: file_exists
        use simple_fileio,          only: file2rarr
        use simple_math_ft,         only: get_resarr, calc_fourier_index
        use simple_estimate_ssnr,   only: get_resolution
        use simple_parameters,      only: parameters
        use simple_sp_project,      only: sp_project
        use simple_cmdline,         only: cmdline
        use simple_defs,            only: TINY, logfhandle
        use simple_error,           only: simple_exception
        character(len=*), intent(in)     :: workflow_label
        class(parameters), intent(inout) :: params
        class(cmdline), intent(inout)    :: cline
        type(sp_project)  :: seed_proj
        type(string)      :: fsc_fname
        real, allocatable :: fsc(:), res(:)
        real              :: fsc05, fsc0143
        integer           :: fsc_box
        if( .not. params%l_nonuniform_lpset ) return
        if( cline%defined('lp') ) return
        fsc_box = 0
        call seed_proj%read_segment('out', params%projfile)
        if( seed_proj%isthere_in_osout('fsc', 1) )then
            call seed_proj%get_fsc(1, fsc_fname, fsc_box)
        else
            fsc_fname = refine3D_fsc_fname(1)
            fsc_box = params%box
        endif
        call seed_proj%kill
        if( .not. file_exists(fsc_fname) )then
            THROW_HARD(workflow_label//' filt_mode=nonuniform_lpset requires starting-volume FSC metadata or explicit lp')
        endif
        fsc = file2rarr(fsc_fname)
        if( size(fsc) < 2 )then
            THROW_HARD('starting-volume FSC has too few shells; '//workflow_label//' filt_mode=nonuniform_lpset')
        endif
        if( fsc_box < 1 ) fsc_box = params%box
        res = get_resarr(fsc_box, params%smpd)
        if( size(res) < size(fsc) ) THROW_HARD('starting-volume FSC/box size mismatch; '//workflow_label)
        call get_resolution(fsc, res, fsc05, fsc0143)
        if( fsc0143 <= TINY )then
            THROW_HARD(workflow_label//' nonuniform_lpset needs a positive starting resolution; set lp')
        endif
        params%lp         = fsc0143
        params%kfromto(2) = calc_fourier_index(params%lp, params%box, params%smpd)
        params%l_lpset    = .true.
        call cline%set('lp', params%lp)
        write(logfhandle,'(A,F8.3,A)') &
            &'>>> '//workflow_label//' nonuniform_lpset seeded matching low-pass from starting FSC 0.143: ', &
            &params%lp, ' A'
        if( allocated(fsc) ) deallocate(fsc)
        if( allocated(res) ) deallocate(res)
        call fsc_fname%kill
    end subroutine seed_refine3D_auto_nonuniform_lpset

    subroutine require_refine3D_local_value( cline, key, expected )
        use simple_error, only: simple_exception
        class(cmdline), intent(in)   :: cline
        character(len=*), intent(in) :: key, expected
        if( .not. cline%defined(key) ) return
        if( cline%get_carg(key) /= expected ) THROW_HARD('refine3D_local requires '//key//'='//expected)
    end subroutine require_refine3D_local_value

    ! One global greedy pass over every particle, banded by the startup FSC via lpstop.
    ! Keep the incumbent pose; log direction/shift reassignment fractions.
    ! The main run continues from this pass's reference and poses as iteration 2.
    subroutine run_refine3D_registration_pass( params, cline, xrefine3D, workflow_label )
        use simple_ori,             only: ori
        use simple_oris,            only: oris
        use simple_refine3D_fnames, only: refine3D_fsc_fname, refine3D_state_vol_fname
        use simple_syslib,          only: file_exists
        use simple_fileio,          only: file2rarr
        use simple_math_ft,         only: get_resarr
        use simple_linalg,          only: rad2deg
        use simple_estimate_ssnr,   only: get_resolution_at_fsc
        use simple_parameters,      only: parameters
        use simple_sp_project,      only: sp_project
        use simple_cmdline,         only: cmdline
        use simple_commander_base,  only: commander_base
        use simple_defs,            only: TINY, logfhandle
        use simple_error,           only: simple_exception
        class(parameters), intent(in)        :: params
        class(cmdline), intent(inout)        :: cline
        class(commander_base), intent(inout) :: xrefine3D
        character(len=*), intent(in)         :: workflow_label
        integer, parameter :: REGPASS_NSPACE = 5000
        type(cmdline)      :: cline_pass
        type(sp_project)   :: pass_proj
        type(oris)         :: os_before, os_after
        type(ori)          :: o_before
        type(string)       :: fsc_fname
        real, allocatable  :: fsc(:), res(:)
        real               :: res_pass, delta, dist, sh_before(2), sh_after(2), user_lpstop
        integer            :: fsc_box, i, n, n_active, n_moved1, n_moved2, n_shift
        logical            :: l_user_lpstop
        if( .not. params%l_regpass ) return
        ! A user lpstop can only make this pass coarser.
        l_user_lpstop = cline%defined('lpstop')
        user_lpstop   = 0.
        if( l_user_lpstop ) user_lpstop = params%lpstop
        ! Resolution at FSC=regpass_fsc of the startup pair.
        fsc_box = 0
        call pass_proj%read_segment('out', params%projfile)
        if( pass_proj%isthere_in_osout('fsc', 1) )then
            call pass_proj%get_fsc(1, fsc_fname, fsc_box)
        else
            fsc_fname = refine3D_fsc_fname(1)
            fsc_box   = params%box
        endif
        call pass_proj%kill
        if( .not. file_exists(fsc_fname) )then
            write(logfhandle,'(A)') '>>> '//workflow_label//' REGISTRATION PASS SKIPPED: no startup FSC on disk'
            call fsc_fname%kill
            return
        endif
        fsc = file2rarr(fsc_fname)
        call fsc_fname%kill
        if( fsc_box < 1 ) fsc_box = params%box
        res = get_resarr(fsc_box, params%smpd)
        if( size(res) < size(fsc) ) THROW_HARD('startup FSC/box size mismatch; '//workflow_label//' registration pass')
        call get_resolution_at_fsc(fsc, res(:size(fsc)), params%regpass_fsc, res_pass)
        deallocate(fsc, res)
        if( res_pass <= TINY )then
            write(logfhandle,'(A,F5.2)') '>>> '//workflow_label// &
                &' REGISTRATION PASS SKIPPED: startup FSC never reaches ', params%regpass_fsc
            return
        endif
        write(logfhandle,'(A,F5.3,A,F7.3,A,I0,A)') '>>> '//workflow_label// &
            &' REGISTRATION PASS: refine=greedy (exhaustive, incumbent included), all particles, band-limited at FSC=', &
            &params%regpass_fsc, ' -> lpstop=', res_pass, ' A, nspace=', REGPASS_NSPACE
        ! Save the poses before the pass for the reassignment diagnostic.
        call pass_proj%read_segment(params%oritype, params%projfile)
        call os_before%copy(pass_proj%os_ptcl3D)
        call pass_proj%kill
        ! Greedy includes the incumbent: a particle moves only when the current
        ! objective has a better pose among the searched candidates.
        cline_pass = cline
        call cline_pass%set('prg',        'refine3D')
        call cline_pass%set('refine',     'greedy')
        call cline_pass%set('nspace',     REGPASS_NSPACE)
        call cline_pass%delete('nspace_sub')
        call cline_pass%delete('prob_athres')
        if( l_user_lpstop ) res_pass = max(res_pass, user_lpstop)
        ! Use lpstop, never lp: l_lpset can merge the half references.
        call cline_pass%set('lpstop',     res_pass)
        call cline_pass%set('maxits',     1)
        call cline_pass%set('minits',     1)
        call cline_pass%set('startit',    1)
        call cline_pass%set('which_iter', 1)
        call cline_pass%delete('endit')
        call cline_pass%delete('continue')
        ! Use every particle, regardless of the main run's sampling policy.
        call cline_pass%delete('update_frac')
        call cline_pass%delete('nsample')
        call xrefine3D%execute(cline_pass)
        call cline_pass%kill
        ! Reassignment is measured against the basin width at the pass band.
        call pass_proj%read_segment(params%oritype, params%projfile)
        call os_after%copy(pass_proj%os_ptcl3D)
        call pass_proj%kill
        n        = os_after%get_noris()
        delta    = res_pass / (0.5 * params%mskdiam)   ! radians
        n_active = 0
        n_moved1 = 0
        n_moved2 = 0
        n_shift  = 0
        do i = 1, n
            if( os_after%get_state(i) < 1 ) cycle
            n_active = n_active + 1
            call os_before%get_ori(i, o_before)
            dist = os_after%euldist(i, o_before)
            if( dist > delta      ) n_moved1 = n_moved1 + 1
            if( dist > 2. * delta ) n_moved2 = n_moved2 + 1
            sh_before = os_before%get_2Dshift(i)
            sh_after  = os_after%get_2Dshift(i)
            if( sqrt(sum((sh_after - sh_before)**2)) > 1. ) n_shift = n_shift + 1
        enddo
        call o_before%kill
        call os_before%kill
        call os_after%kill
        if( n_active > 0 ) write(logfhandle,'(A,F6.2,A,F6.2,A,F6.2,A,F6.2,A)') &
            &'>>> '//workflow_label//' REGISTRATION PASS REASSIGNED: ', 100.*real(n_moved1)/real(n_active), &
            &' % of directions moved > ', rad2deg(delta), ' deg (basin width at the pass band), ', &
            &100.*real(n_moved2)/real(n_active), ' % moved > twice that, ', &
            &100.*real(n_shift)/real(n_active), ' % of shifts moved > 1 px'
        ! Continue from the pass's assembled reference and saved poses as iteration 2.
        call cline%set('vol1',       refine3D_state_vol_fname(1))
        call cline%set('startit',    2)
        call cline%set('which_iter', 2)
        call cline%delete('endit')
        if( l_user_lpstop )then
            call cline%set('lpstop', user_lpstop)
        else
            call cline%delete('lpstop')
        endif
    end subroutine run_refine3D_registration_pass

    subroutine exec_refine3D_states( self, cline )
        use simple_final_rec,                  only: calc_final_rec
        use simple_commanders_rec,             only: commander_rec3D
        use simple_refine3D_states_controller, only: refine3D_states_controller, refine3D_states_policy
        class(commander_refine3D_states), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        type(gui_communicator)           :: gui_comm
        type(refine3D_states_controller) :: states_ctrl
        type(refine3D_states_policy)     :: policy
        type(string)                     :: m_estimator_arg
        real,             parameter :: STATE_OVERLAP_NEIGH_REFINE3D_STATES = 0.99
        real,             parameter :: LPSTART_REFINE3D_STATES             = 10.0
        real,             parameter :: LPSTOP_REFINE3D_STATES              = 6.0
        character(len=*), parameter :: WORKFLOW_LABEL                      = 'REFINE3D_STATES'
        logical :: l_gui_comm_active
        ! commanders
        type(commander_rec3D)           :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        type(commander_refine3D)        :: xrefine3D
        policy%label            = WORKFLOW_LABEL
        ! the gridding M-step for the initialization: no PCG kernels and no solvent support needed
        policy%flex_rec_backend = 'gridding'
        call states_ctrl%new(policy)
        call cline%set('prg', 'refine3D_states')
        ! no weighted mode: the refinement moves the labels the frozen flex_pca weights were computed for;
        ! one state is refined with its weights by refine3D_auto state=X m_estimator=flex
        if( cline%defined('m_estimator') )then
            m_estimator_arg = cline%get_carg('m_estimator')
            if( m_estimator_arg%to_char() /= 'no' ) THROW_HARD(WORKFLOW_LABEL//' has no weighted mode; use refine3D_auto state=X m_estimator=flex')
            call m_estimator_arg%kill
        endif
        call states_ctrl%normalize_inputs(cline)
        if( .not. cline%defined('filt_mode')     ) call cline%set('filt_mode', 'nonuniform_lpset')
        call states_ctrl%resolve_states(cline)
        ! hard defaults (balance falls back to none without selected class averages)
        call cline%set('balance',        'class')
        call cline%set('greedy_sampling', 'no')
        call cline%set('frac_best',       1.0)
        call cline%set('trail_rec',     'yes')
        call cline%set('objfun',      'euclid')
        call cline%set('lplim_crit',       0.5)
        call cline%set('incrreslim',      'no')
        ! overridable defaults
        if( .not. cline%defined('combine_eo')      ) call cline%set('combine_eo',        'no')
        if( .not. cline%defined('envfsc')          ) call cline%set('envfsc',            'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('center')          ) call cline%set('center',            'no')
        if( .not. cline%defined('prob_inpl')       ) call cline%set('prob_inpl',        'yes')
        call states_ctrl%resolve_sampling(cline)
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',               'yes')
        if( .not. cline%defined('lpstart')     ) call cline%set('lpstart', LPSTART_REFINE3D_STATES)
        if( .not. cline%defined('lpstop')      ) call cline%set('lpstop',  LPSTOP_REFINE3D_STATES)
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',               'no')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap', STATE_OVERLAP_NEIGH_REFINE3D_STATES)
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',               'no')
        call states_ctrl%read_maxits(cline)
        call params%new(cline)
        l_gui_comm_active = params%niceserver .ne. "" .and. params%niceprocid > 0
        if( l_gui_comm_active ) call gui_comm%new(params)
        call states_ctrl%validate(params)
        call cline%set('mkdir', 'no')
        call states_ctrl%initialize_states(params, cline)
        call states_ctrl%plan(params, cline)
        call states_ctrl%startup_maps(params, cline)
        call states_ctrl%run_march(params, cline, xrefine3D, gui_comm)
        call states_ctrl%ensure_all_updated(params, cline, xrefine3D)
        ! the shared ending (simple_final_rec): final all-particle
        ! reconstruction at native sampling with reused or bootstrapped sigmas,
        ! project registration, final products and reprojections
        call calc_final_rec(params, spproj, params%projfile, cline, xrec3D, xbootstrap_rec3D, &
            &l_postprocess=.true., lp_snapshot=params%lpstop)
        call gui_comm%add_metadata(spproj, oritype='cls3D', stage=0, selection=.true.)
        call spproj%kill
        call states_ctrl%kill
        call gui_comm%kill()
        call simple_end('**** SIMPLE_REFINE3D_STATES NORMAL STOP ****')
    end subroutine exec_refine3D_states

    !> refine3D_states on the sub-class averages of cls_expansion, in a temporary work project built as
    !! solve3D_cavgs builds its own
    subroutine exec_refine3D_states_cavgs( self, cline )
        use simple_final_rec,                  only: calc_final_rec
        use simple_commanders_rec,             only: commander_rec3D
        use simple_cavgs_work_project,         only: cavgs_work_project
        use simple_refine3D_states_controller, only: refine3D_states_controller, refine3D_states_policy
        class(commander_refine3D_states_cavgs), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters)                 :: params
        type(sp_project)                 :: spproj, work_proj
        type(gui_communicator)           :: gui_comm
        type(cmdline)                    :: work_cline
        type(refine3D_states_controller) :: states_ctrl
        type(refine3D_states_policy)     :: policy
        type(cavgs_work_project)         :: wp
        type(string)                     :: input_projfile, mean_vol, rec_backend_arg, ref_vols(1), checkpoint_vols(1)
        type(string),        allocatable :: final_vols(:)
        logical,             allocatable :: ptcl_active(:)
        integer,             allocatable :: cls_label(:), cls_state(:), new_state(:)
        real,                allocatable :: cls_euls(:,:), cls_shifts(:,:), cls_corr(:), cls_proj(:)
        real,             parameter :: STATE_OVERLAP_TARGET = 0.99
        real,             parameter :: LPSTART_DEFAULT      = 10.0
        real,             parameter :: LPSTOP_DEFAULT       = 6.0
        real,             parameter :: MIN_STATE_FRAC       = 0.1
        real,             parameter :: SMPD_TOL             = 0.001
        integer,          parameter :: MIN_SUBAVGS          = 100
        integer,          parameter :: FLEX_MIN_NEFF        = 20
        integer,          parameter :: FLEX_MIN_FIT_ROWS    = 100
        character(len=*), parameter :: WORKFLOW_LABEL       = 'REFINE3D_STATES_CAVGS'
        character(len=*), parameter :: WORK_PROJFILE        = 'refine3D_states_cavgs_tmpproj.simple'
        character(len=*), parameter :: WORK_SIGMA2_STATE    = 'refine3D_states_cavgs_sigma2_state.bin'
        character(len=*), parameter :: MEAN_VOL_FNAME       = 'refine3D_states_cavgs_mean.mrc'
        integer :: ncavgs, nsel
        logical :: l_gui_comm_active, l_cls3D_labels
        ! commanders
        type(commander_rec3D)           :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        type(commander_refine3D)        :: xrefine3D
        call cline%set('prg', 'refine3D_states_cavgs')
        call refuse_inputs()
        if( .not. cline%defined('rec_backend') ) call cline%set('rec_backend', 'pcg')
        rec_backend_arg = cline%get_carg('rec_backend')
        ! the class-average lineage: one backend everywhere (PCG unless overridden), flex at the native box on
        ! the march's lpstart band and in shared memory, every row in every iteration, no stage maps for the GUI
        policy%label                   = WORKFLOW_LABEL
        policy%flex_rec_backend        = rec_backend_arg%to_char()
        policy%flex_rec_states_backend = rec_backend_arg%to_char()
        policy%flex_min_neff           = FLEX_MIN_NEFF
        policy%l_flex_native_box       = .true.
        policy%l_flex_lp_lpstart       = .true.
        policy%l_flex_shmem            = .true.
        policy%l_full_batch            = .true.
        policy%l_gui_stages            = .false.
        call states_ctrl%new(policy)
        call states_ctrl%normalize_inputs(cline)
        ! fixed settings of the class-average route
        call cline%set('filt_mode',       'none')
        call cline%set('automsk',         'no')
        call cline%set('envfsc',          'no')
        call cline%set('combine_eo',      'no')
        call cline%set('objfun',          'euclid')
        call cline%set('trail_rec',       'no')
        call cline%set('balance',         'none')
        call cline%set('greedy_sampling', 'no')
        call cline%set('frac_best',       1.0)
        call cline%set('lplim_crit',      0.5)
        call cline%set('incrreslim',      'no')
        ! overridable defaults
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('center')         ) call cline%set('center',         'no')
        if( .not. cline%defined('prob_inpl')      ) call cline%set('prob_inpl',      'yes')
        if( .not. cline%defined('ml_reg')         ) call cline%set('ml_reg',         'yes')
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',        LPSTART_DEFAULT)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',         LPSTOP_DEFAULT)
        if( .not. cline%defined('overlap')        ) call cline%set('overlap',        STATE_OVERLAP_TARGET)
        if( .not. cline%defined('keepvol')        ) call cline%set('keepvol',        'no')
        if( .not. cline%defined('min_state_frac') ) call cline%set('min_state_frac', MIN_STATE_FRAC)
        call states_ctrl%resolve_sampling(cline)
        call states_ctrl%read_maxits(cline)
        call params%new(cline)
        input_projfile    = params%projfile
        l_gui_comm_active = params%niceserver .ne. "" .and. params%niceprocid > 0
        if( l_gui_comm_active ) call gui_comm%new(params)
        call states_ctrl%validate(params)
        call cline%set('mkdir', 'no')
        ! the input project is read once here and written once at the end
        call spproj%read(input_projfile)
        call preflight()
        ! the hand-back covers the ptcl2D selection, which a rerun on the output finds unchanged
        ptcl_active = spproj%os_ptcl2D%get_all_asint('state') > 0
        ! the work project, the command line of every child, and the parameters that follow it
        call wp%new(spproj, string(WORK_PROJFILE), string(WORK_SIGMA2_STATE))
        ncavgs     = wp%get_ncavgs()
        work_cline = cline
        call work_cline%set('projfile', WORK_PROJFILE)
        call work_cline%set('oritype',  'ptcl3D')
        params%projfile = wp%get_projfile()
        params%oritype  = 'ptcl3D'
        params%nptcls   = 2 * ncavgs
        call wp%cap_nparts(params, work_cline)
        if( .not. l_cls3D_labels ) call check_flex_fits()
        call seed_poses()
        if( l_cls3D_labels ) call copy_cls3D_labels()
        call states_ctrl%resolve_states(work_cline)
        params%nstates = work_cline%get_iarg('nstates')
        if( l_cls3D_labels )then
            call states_ctrl%initialize_states(params, work_cline)
        else
            call states_ctrl%initialize_states(params, work_cline, flex_mean=mean_vol)
        endif
        call states_ctrl%plan(params, work_cline)
        call states_ctrl%startup_maps(params, work_cline)
        ! the pose seeding ran iteration 1
        call states_ctrl%run_march(params, work_cline, xrefine3D, gui_comm, iters_done=1)
        call states_ctrl%ensure_all_updated(params, work_cline, xrefine3D)
        call pick_better_rows()
        call renumber_states()
        ! the final maps, from every active row as assigned, at native sampling and without postprocessing
        call calc_final_rec(params, work_proj, params%projfile, work_cline, xrec3D, xbootstrap_rec3D, &
            &l_postprocess=.false., lp_snapshot=params%lpstop)
        call write_final_maps()
        call commit()
        if( l_gui_comm_active ) call gui_comm%add_metadata(spproj, oritype='cls3D', stage=0, selection=.true.)
        call wp%delete_files
        call simple_rmdir(string(STKPARTSDIR))
        call wp%kill
        call states_ctrl%kill
        call spproj%kill
        call work_proj%kill
        call work_cline%kill
        call gui_comm%kill()
        call simple_end('**** SIMPLE_REFINE3D_STATES_CAVGS NORMAL STOP ****')

    contains

        subroutine refuse_inputs()
            type(string) :: val
            if( cline%defined('m_estimator') )then
                val = cline%get_carg('m_estimator')
                if( val%to_char() /= 'no' ) THROW_HARD(WORKFLOW_LABEL//' has no weighted mode')
            endif
            call refuse_key('nsample')
            call refuse_key('update_frac')
            call refuse_key('cohort_sampling')
            call refuse_key('trail_seed')
            call refuse_key('pcg_mskfile')
            call require_value('filt_mode',   'none')
            call require_value('automsk',     'no')
            call require_value('envfsc',      'no')
            call require_value('combine_eo',  'no')
            call require_value('trail_rec',   'no')
            call require_value('pcg_solvent', 'no')
            call require_value('objfun',      'euclid')
            if( cline%defined('rec_backend') )then
                val = cline%get_carg('rec_backend')
                select case(val%to_char())
                    case('pcg', 'gridding')
                        ! supported
                    case default
                        THROW_HARD(WORKFLOW_LABEL//' supports rec_backend=pcg|gridding')
                end select
            endif
            if( .not. cline%defined('mskdiam') ) THROW_HARD(WORKFLOW_LABEL//' needs mskdiam, the support of the PCG solves')
            call val%kill
        end subroutine refuse_inputs

        subroutine refuse_key( key )
            character(len=*), intent(in) :: key
            character(len=:), allocatable :: errmsg
            if( cline%defined(key) )then
                errmsg = WORKFLOW_LABEL//' sets '//key//' itself; remove it from the command line'
                THROW_HARD(errmsg)
            endif
        end subroutine refuse_key

        subroutine require_value( key, val )
            character(len=*), intent(in) :: key, val
            character(len=:), allocatable :: errmsg
            type(string)                  :: given
            if( .not. cline%defined(key) ) return
            given = cline%get_carg(key)
            if( given%to_char() /= val )then
                errmsg = WORKFLOW_LABEL//' requires '//key//'='//val
                THROW_HARD(errmsg)
            endif
            call given%kill
        end subroutine require_value

        !> everything checked before the first write: stacks, consensus map, particle fields, selection, the
        !! populated cls3D labels and, without multi-state labels, flex's state count and population floor
        subroutine preflight()
            use simple_fileio,  only: fname2ext
            use simple_imghead, only: find_img_smpd
            type(string)                  :: stk, stk_eo(2), ext, consensus_vol
            integer,          allocatable :: cls2D_states(:), cls3D_states(:), ptcl2D_states(:), ptcl3D_states(:)
            character(len=:), allocatable :: errmsg
            integer :: ncls_stk, ldim(3), nimgs, vol_box, nrows, nmin, ieo, state, nlabels, icls
            real    :: smpd_stk, smpd_hdr, vol_smpd
            call spproj%get_cavgs_stk(stk, ncls_stk, smpd_stk)
            if( spproj%os_cls2D%get_noris() /= ncls_stk .or. spproj%os_cls3D%get_noris() /= ncls_stk )then
                THROW_HARD(WORKFLOW_LABEL//' needs a cls2D and a cls3D row for every class average')
            endif
            if( .not. file_exists(stk) ) THROW_HARD('class-average stack does not exist: '//stk%to_char())
            call find_ldim_nptcls(stk, ldim, nimgs)
            if( nimgs /= ncls_stk )then
                THROW_HARD(WORKFLOW_LABEL//' needs one image in the class-average stack per cls2D row')
            endif
            if( ldim(1) /= params%box .or. abs(smpd_stk - params%smpd) > SMPD_TOL )then
                THROW_HARD(WORKFLOW_LABEL//' needs the class averages at the box and sampling of the particles')
            endif
            ! the even and odd companions match the merged stack in box, count and header sampling
            smpd_hdr  = find_img_smpd(stk)
            ext       = string('.')//fname2ext(stk)
            stk_eo(1) = add2fbody(stk, ext, '_even')
            stk_eo(2) = add2fbody(stk, ext, '_odd')
            do ieo = 1,2
                if( .not. file_exists(stk_eo(ieo)) )then
                    THROW_HARD('class-average companion does not exist: '//stk_eo(ieo)%to_char())
                endif
                call find_ldim_nptcls(stk_eo(ieo), ldim, nimgs)
                if( ldim(1) /= params%box .or. nimgs /= ncls_stk .or. &
                    &abs(find_img_smpd(stk_eo(ieo)) - smpd_hdr) > SMPD_TOL )then
                    THROW_HARD('class-average companion differs from the merged stack: '//stk_eo(ieo)%to_char())
                endif
            enddo
            if( .not. spproj%isthere_in_osout('vol', 1) )then
                THROW_HARD(WORKFLOW_LABEL//' needs the consensus map (vol, state 1) in the project')
            endif
            call spproj%get_vol('vol', 1, consensus_vol, vol_smpd, vol_box)
            if( .not. file_exists(consensus_vol) ) THROW_HARD('consensus map does not exist: '//consensus_vol%to_char())
            if( vol_box /= params%box .or. abs(vol_smpd - params%smpd) > SMPD_TOL )then
                THROW_HARD(WORKFLOW_LABEL//' needs the consensus map at the box and sampling of the particles')
            endif
            ref_vols(1) = simple_abspath(consensus_vol)
            ptcl2D_states = spproj%os_ptcl2D%get_all_asint('state')
            ptcl3D_states = spproj%os_ptcl3D%get_all_asint('state')
            if( size(ptcl2D_states) /= size(ptcl3D_states) )then
                THROW_HARD(WORKFLOW_LABEL//' needs congruent ptcl2D and ptcl3D fields')
            endif
            ! a rerun on the output finds the particles of uncommitted sub-classes inactive in ptcl3D only
            if( any(ptcl3D_states > 0 .and. ptcl2D_states <= 0) )then
                THROW_HARD(WORKFLOW_LABEL//' needs every particle active in ptcl3D to be active in ptcl2D')
            endif
            cls2D_states = nint(spproj%os_cls2D%get_all('state'))
            nsel = count(cls2D_states > 0)
            if( nsel < MIN_SUBAVGS )then
                errmsg = WORKFLOW_LABEL//' needs at least '//int2str(MIN_SUBAVGS)//&
                    &' selected sub-class averages, found '//int2str(nsel)
                THROW_HARD(errmsg)
            endif
            ! the populated cls3D labels of the selected classes, numbered consecutively (an empty label is
            ! dropped); two or more are refined as they are, otherwise flex_pca initializes the states
            cls3D_states = nint(spproj%os_cls3D%get_all('state'))
            where( cls2D_states <= 0 .or. cls3D_states < 0 ) cls3D_states = 0
            allocate(cls_label(ncls_stk), source=0)
            nlabels = 0
            do state = 1,maxval(cls3D_states)
                if( count(cls3D_states == state) == 0 ) cycle
                nlabels = nlabels + 1
                where( cls3D_states == state ) cls_label = nlabels
            enddo
            l_cls3D_labels = nlabels > 1
            if( l_cls3D_labels )then
                write(logfhandle,'(A,I0,A)') '>>> '//WORKFLOW_LABEL//' REFINES THE ', nlabels, &
                    &' POPULATED STATE LABELS OF cls3D'
                if( cline%defined('nstates') .and. params%nstates /= nlabels )then
                    errmsg = WORKFLOW_LABEL//' found '//int2str(nlabels)//' populated cls3D state labels, not nstates='//&
                        &int2str(params%nstates)
                    THROW_HARD(errmsg)
                endif
            else
                if( .not. cline%defined('nstates') .or. params%nstates < 3 )then
                    THROW_HARD(WORKFLOW_LABEL//' needs nstates >= 3 for the flex_pca initialization')
                endif
                ! flex's own population-floor check, on the active rows: two per selected class, less
                ! the halves cls2D flags with state 0
                nrows = 0
                do icls = 1,ncls_stk
                    if( cls2D_states(icls) <= 0 ) cycle
                    nrows = nrows + half_active(icls, 'state_even') + half_active(icls, 'state_odd')
                enddo
                nmin  = max(1, nint(params%min_state_frac * real(nrows)))
                if( params%nstates * nmin > nrows )then
                    THROW_HARD(WORKFLOW_LABEL//' has too few selected sub-class averages for nstates at min_state_frac')
                endif
            endif
        end subroutine preflight

        !> 1 when the half of class icls is active, 0 when cls2D flags it with state 0 (flag absent: active)
        integer function half_active( icls, flag )
            integer,          intent(in) :: icls
            character(len=*), intent(in) :: flag
            half_active = 1
            if( spproj%os_cls2D%isthere(icls, flag) )then
                if( spproj%os_cls2D%get_int(icls, flag) == 0 ) half_active = 0
            endif
        end function half_active

        !> flex's two fits split the rows by flex_pca_half_of and each needs FLEX_MIN_FIT_ROWS active rows,
        !! counted on the work project's rows as flex will see them
        subroutine check_flex_fits()
            use simple_flex_pca_stages, only: flex_pca_half_of
            character(len=:), allocatable :: errmsg
            integer :: irow, nfit(2)
            call work_proj%read_segment('ptcl3D', params%projfile)
            nfit = 0
            do irow = 1,work_proj%os_ptcl3D%get_noris()
                if( work_proj%os_ptcl3D%get_state(irow) == 0 ) cycle
                nfit(flex_pca_half_of(irow)) = nfit(flex_pca_half_of(irow)) + 1
            enddo
            call work_proj%kill
            if( any(nfit < FLEX_MIN_FIT_ROWS) )then
                errmsg = WORKFLOW_LABEL//' gives a flex_pca fit fewer than '//int2str(FLEX_MIN_FIT_ROWS)//' rows: '//&
                    &int2str(nfit(1))//' and '//int2str(nfit(2))
                THROW_HARD(errmsg)
            endif
        end subroutine check_flex_fits

        !> every row against the consensus map; the checkpoint map, from the class averages, is flex's mean.
        !! The rows are pose-less (the work project gives the class-average halves no orientation or shift)
        subroutine seed_poses()
            integer :: nactive_rows
            call work_proj%read_segment('ptcl3D', params%projfile)
            nactive_rows = work_proj%os_ptcl3D%count_state_gt_zero()
            call work_proj%kill
            call initialize_poses_against_external_references(params, work_cline, xrefine3D, xrec3D, &
                &nactive_rows, ref_vols, checkpoint_vols, 1, lp_init=params%lpstart, l_all_rows=.true., l_poseless=.true.)
            if( l_cls3D_labels ) return
            ! flex's delivered maps reuse the checkpoint's file name
            call simple_copy_file(checkpoint_vols(1), string(MEAN_VOL_FNAME))
            mean_vol = simple_abspath(string(MEAN_VOL_FNAME))
        end subroutine seed_poses

        !> the consecutive cls3D labels of the selected classes onto the active rows of each class (an inactive
        !! half stays inactive)
        subroutine copy_cls3D_labels()
            call work_proj%read_segment('ptcl3D', params%projfile)
            call wp%label_rows(work_proj%os_ptcl3D, cls_label)
            call work_proj%write_segment_inside('ptcl3D', params%projfile)
            call work_proj%kill
        end subroutine copy_cls3D_labels

        !> per class, the pose, state and score of the better-scoring active row, as solve3D_cavgs picks them
        !! (a tie goes to the odd row); a class without an active row stays at state 0
        subroutine pick_better_rows()
            integer :: icls, ie, io, ibest, state_even, state_odd
            call work_proj%read_segment('ptcl3D', params%projfile)
            allocate(cls_state(ncavgs), source=0)
            allocate(cls_euls(3,ncavgs), cls_shifts(2,ncavgs), cls_corr(ncavgs), cls_proj(ncavgs), source=0.)
            do icls = 1,ncavgs
                ie         = wp%even_row(icls)
                io         = wp%odd_row(icls)
                state_even = work_proj%os_ptcl3D%get_state(ie)
                state_odd  = work_proj%os_ptcl3D%get_state(io)
                if( state_even == 0 .and. state_odd == 0 ) cycle
                if( state_odd == 0 )then
                    ibest = ie
                else if( state_even == 0 )then
                    ibest = io
                else if( work_proj%os_ptcl3D%get(ie, 'corr') > work_proj%os_ptcl3D%get(io, 'corr') )then
                    ibest = ie
                else
                    ibest = io
                endif
                cls_state(icls)    = work_proj%os_ptcl3D%get_state(ibest)
                cls_euls(:,icls)   = work_proj%os_ptcl3D%get_euler(ibest)
                cls_shifts(:,icls) = work_proj%os_ptcl3D%get_2Dshift(ibest)
                cls_corr(icls)     = work_proj%os_ptcl3D%get(ibest, 'corr')
                cls_proj(icls)     = work_proj%os_ptcl3D%get(ibest, 'proj')
            enddo
            write(logfhandle,'(A,F6.1,A)') '>>> '//WORKFLOW_LABEL//' EVEN/ODD STATE OVERLAP: ', &
                &wp%eo_state_overlap(work_proj%os_ptcl3D), '%'
            write(logfhandle,'(A,F6.1)') '>>> '//WORKFLOW_LABEL//' EVEN/ODD AVERAGE ANGULAR DISTANCE: ', &
                &wp%eo_angular_distance(work_proj%os_ptcl3D, trim(params%pgrp))
            call work_proj%kill
        end subroutine pick_better_rows

        !> rec_final_stateNN.mrc, NN the renumbered state, for every populated state with a map, and its
        !! low-pass snapshot at FSC 0.143, never finer than lpstop because the even/odd FSC of class averages
        !! reads optimistic
        subroutine write_final_maps()
            use simple_halfmap_diagnostics, only: copy_support_provenance
            use simple_solve3D_utils,       only: solve3D_rec_fbody, solve3D_state_fsc_lowpass, &
                &write_solve3D_lowpass_snapshot
            type(string) :: vol_name, vol_final_lp
            real    :: lp_snapshot
            integer :: state
            allocate(final_vols(params%nstates))
            call work_proj%read_segment('out', params%projfile)
            do state = 1,params%nstates
                if( new_state(state) == 0 ) cycle
                if( .not. work_proj%isthere_in_osout('vol', state) ) cycle
                vol_name = refine3D_state_vol_fname(state)
                if( .not. file_exists(vol_name) ) cycle
                final_vols(state) = solve3D_rec_fbody()//int2str_pad(new_state(state),2)//MRC_EXT
                call simple_copy_file(vol_name, final_vols(state))
                call copy_support_provenance(vol_name, final_vols(state))
                lp_snapshot  = max(params%lpstop, solve3D_state_fsc_lowpass(state, params%box, params%smpd, params%lpstop))
                vol_final_lp = add2fbody(final_vols(state), MRC_EXT, LP_SUFFIX)
                call write_solve3D_lowpass_snapshot(final_vols(state), lp_snapshot, vol_final_lp, params%smpd, box=params%box)
            enddo
            call work_proj%kill
            call vol_name%kill
            call vol_final_lp%kill
        end subroutine write_final_maps

        !> the states that hold committed sub-class averages, numbered consecutively; an empty state is reported
        !! and dropped, never an error
        subroutine renumber_states()
            integer :: state, icls, nnew
            allocate(new_state(params%nstates), source=0)
            nnew = 0
            do state = 1,params%nstates
                if( count(cls_state == state) == 0 )then
                    write(logfhandle,'(A,I0,A)') '>>> '//WORKFLOW_LABEL//' STATE ', state, ' HOLDS NO SUB-CLASS AVERAGE; DROPPED'
                    cycle
                endif
                nnew = nnew + 1
                new_state(state) = nnew
                write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> '//WORKFLOW_LABEL//' STATE ', state, ' -> ', nnew, &
                    &' (', count(cls_state == state), ' SUB-CLASS AVERAGES)'
            enddo
            if( nnew < 2 ) write(logfhandle,'(A,I0,A)') '>>> '//WORKFLOW_LABEL//' ONLY ', nnew, &
                &' STATE HOLDS SUB-CLASS AVERAGES'
            do icls = 1,ncavgs
                if( cls_state(icls) > 0 ) cls_state(icls) = new_state(cls_state(icls))
            enddo
        end subroutine renumber_states

        !> cls3D in place (cluster and accept kept), particle states (poses kept), vol_cavg replaced as a set;
        !! one write of the input project
        subroutine commit()
            integer :: icls, state
            do icls = 1,ncavgs
                call spproj%os_cls3D%set_state(icls, cls_state(icls))
                if( cls_state(icls) == 0 ) cycle
                call spproj%os_cls3D%set_euler(icls, cls_euls(:,icls))
                call spproj%os_cls3D%set_shift(icls, cls_shifts(:,icls))
                call spproj%os_cls3D%set(icls, 'corr', cls_corr(icls))
                call spproj%os_cls3D%set(icls, 'proj', cls_proj(icls))
            enddo
            call spproj%map_cls3D_states2ptcls(ptcl_active)
            do state = 1,MAXS
                call spproj%remove_entry_from_osout('vol_cavg', state)
            enddo
            do state = 1,params%nstates
                if( new_state(state) == 0 ) cycle
                if( .not. final_vols(state)%is_allocated() )then
                    THROW_HARD(WORKFLOW_LABEL//' found no final map for a state that holds sub-class averages')
                endif
                call spproj%add_vol2os_out(simple_abspath(final_vols(state)), params%smpd, new_state(state), 'vol_cavg')
            enddo
            call spproj%write(input_projfile)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> '//WORKFLOW_LABEL//' COMMITTED ', count(cls_state > 0), &
                &' OF ', ncavgs, ' SUB-CLASS AVERAGES'
            write(logfhandle,'(A,I0,A,I0,A)') '>>> '//WORKFLOW_LABEL//' HANDED BACK STATES TO ', &
                &spproj%os_ptcl3D%count_state_gt_zero(), ' OF ', count(ptcl_active), ' PARTICLES SELECTED IN ptcl2D'
        end subroutine commit

    end subroutine exec_refine3D_states_cavgs

    subroutine exec_classify3D_refs( self, cline )
        use simple_final_rec,      only: calc_final_rec
        use simple_commanders_rec, only: commander_rec3D
        use simple_estimate_ssnr,  only: lpstages_setlims
        class(commander_classify3D_refs), intent(inout) :: self
        class(cmdline),                intent(inout)    :: cline
        integer, parameter :: NSAMPLE_PER_STATE_CLASSIFY3D_REFS = 10000
        integer, parameter :: NSAMPLE_CLASSIFY3D_REFS_CAP = 100000
        integer, parameter :: NSPACE                      = 5000
        integer, parameter :: NSPACE_SUB                  = 500
        integer, parameter :: MINITS_CLASSIFY3D_REFS      = 10
        integer, parameter :: MAXITS_CLASSIFY3D_REFS      = 50
        real,    parameter :: TARGET_UPDATES_PER_PARTICLE = 4.0
        real,    parameter :: LPSTART_CLASSIFY3D_REFS     = 10.0
        real,    parameter :: LPSTOP_CLASSIFY3D_REFS      = 6.0
        character(len=*), parameter :: WORKFLOW_LABEL     = 'CLASSIFY3D_REFS'
        type(commander_rec3D)           :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        type(commander_refine3D)  :: xrefine3D
        type(parameters)          :: params
        type(sp_project)          :: spproj
        type(lp_crop_inf)         :: lpinfo_master(1)
        type(string), allocatable :: init_vols(:)
        integer,      allocatable :: state_pops(:)
        integer :: nptcls_eff, iter_glob, state
        call cline%set('prg', 'classify3D_refs')
        ! hard defaults (balance falls back to none without selected class averages)
        call cline%set('balance',         'class')
        call cline%set('greedy_sampling', 'no')
        call cline%set('frac_best',       1.0)
        call cline%set('trail_rec',       'yes')
        call cline%set('objfun',          'euclid')
        call cline%set('lplim_crit',      0.5)
        call cline%set('incrreslim',      'no')
        call cline%set('combine_eo',      'no')
        call cline%set('multivol_mode',   'independent')
        ! staged frequency marching always uses planner-driven downscaling
        if( cline%defined('autoscale') ) write(logfhandle,'(A)') &
            &'>>> '//WORKFLOW_LABEL//' IGNORES autoscale: frequency marching owns downscaling'
        call cline%set('autoscale', 'yes')
        ! overridable defaults
        if( .not. cline%defined('filt_mode')       ) call cline%set('filt_mode',       'nonuniform_lpset')
        if( .not. cline%defined('envfsc')          ) call cline%set('envfsc',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('center')          ) call cline%set('center',          'no')
        if( .not. cline%defined('prob_inpl')       ) call cline%set('prob_inpl',       'yes')
        if( .not. cline%defined('refine')          ) call cline%set('refine',          'prob_neigh')
        if( .not. cline%defined('prob_neigh_mode') ) call cline%set('prob_neigh_mode', 'state')
        if( .not. cline%defined('ml_reg')          ) call cline%set('ml_reg',          'yes')
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart', LPSTART_CLASSIFY3D_REFS)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',  LPSTOP_CLASSIFY3D_REFS)
        if( .not. cline%defined('automsk')         ) call cline%set('automsk',         'no')
        if( .not. cline%defined('nsample')         ) call cline%set('nsample', NSAMPLE_PER_STATE_CLASSIFY3D_REFS)
        if( .not. cline%defined('keepvol')         ) call cline%set('keepvol',         'no')
        if( .not. cline%defined('nspace')          ) call cline%set('nspace',          NSPACE)
        if( .not. cline%defined('nspace_sub')      ) call cline%set('nspace_sub',      NSPACE_SUB)
        call params%new(cline)
        if( params%nstates < 2 ) THROW_HARD('nstates must be >= 2 for '//WORKFLOW_LABEL)
        call cline%set('mkdir', 'no')
        call spproj%read( params%projfile )
        ! Search planning
        iter_glob = 0
        call set_classify3D_refs_nstates(nptcls_eff, state_pops)
        call set_classify3D_refs_nsample
        call validate_classify3D_refs_filtering
        call validate_classify3D_refs_search_mode
        call set_classify3D_refs_sampling
        call prepare_classify3D_refs_class_sampling
        call spproj%os_ptcl2D%kill
        call set_classify3D_refs_downscaling
        if( spproj%os_ptcl3D%has_been_sampled() )then
            call spproj%os_ptcl3D%clean_entry('sampled', 'updatecnt')
        endif
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        call initialize_state_volumes
        ! Search
        call cline%set('prg',         'refine3D')
        call cline%set('maxits_glob', max(1, params%maxits_glob))
        call cline%set('balance',     params%balance)
        call run_classify3D_refs
        ! Final particle mapping pass
        call ensure_all_active_particles_updated
        ! Final reconstruction
        call reconstruct_all_particles_volumes
        if( allocated(state_pops) ) deallocate(state_pops)
        if( allocated(init_vols) )then
            do state = 1,size(init_vols)
                call init_vols(state)%kill
            enddo
            deallocate(init_vols)
        endif
        call spproj%kill
        call simple_end('**** SIMPLE_CLASSIFY3D_REFS NORMAL STOP ****')
    contains

        subroutine set_classify3D_refs_nstates( nptcls_eff, state_pops )
            integer,              intent(inout) :: nptcls_eff
            integer, allocatable, intent(inout) :: state_pops(:)
            integer :: state, nstates_labels
            nptcls_eff = spproj%os_ptcl3D%count_state_gt_zero()
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nstates_labels = spproj%os_ptcl3D%get_n('state')
            if( nstates_labels == 1 )then
                allocate(state_pops(params%nstates), source=0)
                state_pops(1) = nptcls_eff
            elseif( nstates_labels == params%nstates )then
                call spproj%os_ptcl3D%get_pops(state_pops, 'state', maxn=params%nstates)
                write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' NSTATES FROM PROJECT: ', nstates_labels
            else
                THROW_HARD('command-line nstates does not match project state labels for '//WORKFLOW_LABEL)
            endif
            do state = 1, size(state_pops)
                write(logfhandle,'(A,I3,A,I0)') '>>> '//WORKFLOW_LABEL//' STATE POPULATION: ', state, '/', state_pops(state)
            enddo
        end subroutine set_classify3D_refs_nstates

        subroutine set_classify3D_refs_nsample()
            params%nsample = min(NSAMPLE_CLASSIFY3D_REFS_CAP, params%nsample * params%nstates)
            params%nsample = min(nptcls_eff, params%nsample)
            if( params%nsample < 1 ) THROW_HARD('nsample must be >= 1 for '//WORKFLOW_LABEL)
            write(logfhandle,'(A,I0,A,I0)')'>>> '//WORKFLOW_LABEL//' NSAMPLE: ', params%nsample
        end subroutine set_classify3D_refs_nsample

        subroutine validate_classify3D_refs_filtering()
            select case(trim(params%filt_mode))
                case('nonuniform_lpset', 'none', 'nonuniform')
                    ! supported
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports filt_mode=nonuniform_lpset|none')
            end select
        end subroutine validate_classify3D_refs_filtering

        subroutine validate_classify3D_refs_search_mode()
            select case(trim(params%refine))
                case('prob_neigh', 'prob', 'shc')
                    call cline%set('refine', params%refine)
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports refine=prob_neigh|prob|shc')
            end select
            call cline%set('refine', params%refine)
            if( trim(params%refine).eq.'prob_neigh' )then
                select case(trim(params%prob_neigh_mode))
                case('geom', 'state')
                    call cline%set('prob_neigh_mode', params%prob_neigh_mode)
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports prob_neigh_mode=geom|state')
                end select
            endif
        end subroutine validate_classify3D_refs_search_mode

        subroutine set_classify3D_refs_sampling()
            real    :: update_frac_auto
            integer :: maxits_auto, nptcls_per_iter, stage_cap
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nptcls_per_iter = min(nptcls_eff, params%nsample)
            if( nptcls_eff <= params%nsample )then
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)')'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &nptcls_eff, '/', params%nsample, ' -> FULL UPDATE'
            else
                update_frac_auto = real(params%nsample) / real(nptcls_eff)
                if( update_frac_auto <= 0.99 )then
                    params%update_frac   = update_frac_auto
                    params%l_update_frac = .true.
                    params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                    call cline%set('update_frac', update_frac_auto)
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)')'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                        &nptcls_eff, '/', params%nsample, '/', update_frac_auto
                else
                    params%update_frac   = 1.0
                    params%l_update_frac = .false.
                    params%l_trail_rec   = .false.
                    call cline%delete('update_frac')
                    write(logfhandle,'(A,I0,A,I0,A)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                        &nptcls_eff, '/', params%nsample, ' -> FULL UPDATE'
                endif
            endif
            if( cline%defined('maxits') )then
                write(logfhandle,'(A,I0)') &
                    &'>>> '//WORKFLOW_LABEL//' STAGE MAXITS COMMAND-LINE OVERRIDE: ', params%maxits
            else
                maxits_auto   = ceiling((TARGET_UPDATES_PER_PARTICLE * real(nptcls_eff)) / real(nptcls_per_iter))
                stage_cap     = max(MINITS_CLASSIFY3D_REFS, min(MAXITS_CLASSIFY3D_REFS, maxits_auto))
                params%maxits = stage_cap
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,F5.1,A)') '>>> '//WORKFLOW_LABEL//' STAGE MAXITS: ', &
                    &stage_cap, ' FOR ~', TARGET_UPDATES_PER_PARTICLE, ' UPDATES/PARTICLE'
            endif
        end subroutine set_classify3D_refs_sampling

        subroutine prepare_classify3D_refs_class_sampling()
            use simple_view_partition_sampling, only: make_class_samples, report_class_sample_coverage, has_selected_cavgs
            type(class_sample), allocatable :: clssmp(:)
            if( params%update_frac > 0.99 )then
                write(logfhandle,'(A)') '>>> FORCING FULL ACTIVE SAMPLING (NO FRACTIONAL OR TRAILING UPDATE)'
                params%balance = 'none'
            else if( spproj%is_virgin_field('ptcl2D') .or. .not. has_selected_cavgs(spproj) )then
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' NO SELECTED CLASS AVERAGES: balance=none'
                params%balance = 'none'
            endif
            if( trim(params%balance) == 'none' ) return
            call make_class_samples(params, spproj, params%balance, params%nthr, .true., clssmp)
            call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
            call report_class_sample_coverage(clssmp, params%balance, min(nptcls_eff, params%nsample), params%maxits, 'iteration')
            call deallocate_class_samples(clssmp)
        end subroutine prepare_classify3D_refs_class_sampling

        subroutine set_classify3D_refs_downscaling()
            call lpstages_setlims(params%box, 1, params%smpd, params%lpstart, params%lpstop, lpinfo_master(1))
            params%trs = lpinfo_master(1)%trslim
            call cline%set('trs', params%trs)
            if( lpinfo_master(1)%l_autoscale )then
                params%box_crop  = lpinfo_master(1)%box_crop
                params%smpd_crop = lpinfo_master(1)%smpd_crop
                call cline%set('box_crop', params%box_crop)
                call cline%set('smpd_crop', params%smpd_crop)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//WORKFLOW_LABEL//' DOWNSCALING BOX/SMPD_CROP: ', &
                    &params%box, '/', params%box_crop, '/', params%smpd_crop
            else
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' DOWNSCALING: native sampling retained'
            endif
        end subroutine set_classify3D_refs_downscaling

        subroutine initialize_state_volumes()
            type(string), allocatable :: checkpoint_vols(:)
            logical :: vols_defined(params%nstates)
            integer :: state
            do state = 1,params%nstates
                vols_defined(state) = cline%defined('vol'//int2str(state))
            enddo
            if( .not. all(vols_defined) )then
                THROW_HARD(WORKFLOW_LABEL//' requires a complete external vol1..volN reference set')
            endif
            allocate(init_vols(params%nstates), checkpoint_vols(params%nstates))
            call validate_input_volumes()
            init_vols(1:params%nstates) = params%vols(1:params%nstates)
            write(logfhandle,'(A)')'>>> '//WORKFLOW_LABEL//&
                &' EXTERNAL REFERENCES ARE UNTRUSTED UNTIL CC POSE INITIALIZATION COMPLETES'
            write(logfhandle,'(A,F7.2)') '>>> '//WORKFLOW_LABEL//' INITIALIZING POSES WITH LP:', params%lpstart
            call initialize_poses_against_external_references(params, cline, xrefine3D, xrec3D, &
                &nptcls_eff, init_vols, checkpoint_vols, 1, lp_init=params%lpstart)
            params%vols(1:params%nstates) = checkpoint_vols
            iter_glob = 1
            ! update command line
            do state = 1,params%nstates
                call cline%set('vol'//int2str(state), params%vols(state))
            enddo
            deallocate(checkpoint_vols)
        end subroutine initialize_state_volumes

        subroutine validate_input_volumes()
            type(string) :: vol
            integer :: state, ldim(3), nptcls_dummy
            real    :: vol_smpd
            do state = 1,params%nstates
                vol = params%vols(state)
                if( .not. file_exists(vol) ) THROW_HARD('Input volume does not exist: '//vol%to_char())
                call find_ldim_nptcls(vol, ldim, nptcls_dummy)
                vol_smpd = find_img_smpd(vol)
                if( any(ldim /= [params%box,params%box,params%box]) .or. abs(vol_smpd - params%smpd) > 1.e-6 )then
                    THROW_HARD('Input state volumes must have same dimensions/sampling as the project particles')
                endif
            enddo
            call vol%kill
        end subroutine validate_input_volumes

        subroutine run_classify3D_refs()
            integer, parameter  :: STEP_IT = 3
            type(refine3D_stage_plan_entry), allocatable :: stages(:)
            integer :: startit, lastit, stage
            call plan_refine3D_frequency_stages(params%box, params%smpd, params%lpstart, params%lpstop, &
                &params%maxits, STEP_IT, iter_glob + 1, lpinfo_master(1), stages)
            write(logfhandle,'(A,I0,A,F7.2,A,F7.2)') &
                &'>>> '//WORKFLOW_LABEL//' FREQUENCY MARCHING NSTAGES/LPSTART/LPSTOP: ',&
                &size(stages), '/', real(params%lpstart), '/', real(params%lpstop)
            do stage = 1,size(stages)
                write(logfhandle,'(A,I3,A,F7.2)')'>>> '//WORKFLOW_LABEL//' PLANNED STAGE/LP: ',&
                    &stage, ' / ', stages(stage)%lpinfo%lp
            enddo
            lastit = 0
            do stage = 1,size(stages)
                startit = iter_glob + 1
                call cline%set('minits',     1)
                call cline%set('maxits',     stages(stage)%nits)
                call cline%set('startit',    startit)
                call cline%set('which_iter', startit)
                call cline%set('extr_iter',  startit)
                call cline%set('trs',        stages(stage)%lpinfo%trslim)
                call cline%set('lp',         stages(stage)%lpinfo%lp)
                call cline%set('lpstop',     stages(stage)%lpinfo%lp)
                if( stages(stage)%lpinfo%l_autoscale )then
                    call cline%set('box_crop',  stages(stage)%lpinfo%box_crop)
                    call cline%set('smpd_crop', stages(stage)%lpinfo%smpd_crop)
                else
                    call cline%delete('box_crop')
                    call cline%delete('smpd_crop')
                endif
                write(logfhandle,'(A,I0,A,F7.2)')'>>> '//WORKFLOW_LABEL//' ENTERING STAGE ', &
                    &stage, ' LP: ', stages(stage)%lpinfo%lp
                call xrefine3D%execute(cline)
                lastit    = cline%get_iarg('endit')
                iter_glob = lastit
                call cline%delete('endit')
            enddo
            write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' EXITING FREQUENCY MARCHING AT ITERATION ', iter_glob
            deallocate(stages)
        end subroutine run_classify3D_refs

        subroutine ensure_all_active_particles_updated()
            integer :: nactive, nupdated, nmissing
            call read_update_coverage(nactive, nupdated, nmissing)
            if( nactive < 1 )then
                THROW_HARD(WORKFLOW_LABEL//' has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                call run_classify3D_refs_missing_update(nmissing, nactive)
                call read_update_coverage(nactive, nupdated, nmissing)
                if( nmissing > 0 )then
                    THROW_HARD(WORKFLOW_LABEL//' final missing-update pass failed to update every active particle')
                endif
            endif
        end subroutine ensure_all_active_particles_updated

        subroutine read_update_coverage( nactive, nupdated, nmissing )
            type(sp_project) :: update_proj
            integer, allocatable :: states(:), updatecnts(:)
            integer, intent(out) :: nactive, nupdated, nmissing
            call update_proj%read_segment('ptcl3D', params%projfile)
            if( .not. update_proj%os_ptcl3D%isthere('updatecnt') )then
                call update_proj%kill
                THROW_HARD(WORKFLOW_LABEL//' cannot finish before active particles are updated')
            endif
            states     = update_proj%os_ptcl3D%get_all_asint('state')
            updatecnts = update_proj%os_ptcl3D%get_all_asint('updatecnt')
            nactive    = count(states > 0)
            nupdated   = count(states > 0 .and. updatecnts > 0)
            nmissing   = nactive - nupdated
            write(logfhandle,'(A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLE UPDATE COVERAGE UPDATED/ACTIVE/MISSING: ', &
                &nupdated, '/', nactive, '/', nmissing
            if( allocated(states)     ) deallocate(states)
            if( allocated(updatecnts) ) deallocate(updatecnts)
            call update_proj%kill
        end subroutine read_update_coverage

        subroutine run_classify3D_refs_missing_update( nmissing, nactive )
            integer, intent(in) :: nmissing, nactive
            type(cmdline) :: cline_missing
            iter_glob = iter_glob + 1
            write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' FINAL MISSING-UPDATE ASSIGNMENT', &
                &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_glob
            call flush(logfhandle)
            cline_missing = cline
            call cline_missing%set('prg',           'refine3D')
            call cline_missing%set('mkdir',              'no')
            call cline_missing%set('balance',          'none')
            call cline_missing%set('frac_best',           1.0)
            call cline_missing%set('fillin',             'no')
            call cline_missing%set('update_frac',         1.0)
            call cline_missing%set('trail_rec',          'no')
            call cline_missing%set('volrec',             'no')
            call cline_missing%set('maxits',                1)
            call cline_missing%set('startit',       iter_glob)
            call cline_missing%set('which_iter',    iter_glob)
            call cline_missing%set('extr_iter',     iter_glob)
            call cline_missing%set('refine',         'greedy')
            call cline_missing%set('greedy_sampling',   'yes')
            call cline_missing%set('update_missing',    'yes')
            call cline_missing%delete('endit')
            call xrefine3D%execute(cline_missing)
            call cline_missing%kill
        end subroutine run_classify3D_refs_missing_update

        subroutine reconstruct_all_particles_volumes
            ! the shared ending (simple_final_rec): final all-particle
            ! reconstruction at native sampling with reused or bootstrapped
            ! sigmas, project registration, final products and reprojections
            call calc_final_rec(params, spproj, params%projfile, cline, xrec3D, xbootstrap_rec3D, &
                &l_postprocess=.true., lp_snapshot=params%lpstop)
        end subroutine reconstruct_all_particles_volumes

    end subroutine exec_classify3D_refs

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_refine3D( self, cline )
        use simple_core_module_api
        use simple_refine3D_strategy
        use simple_eul_prob_tab_utils, only: prob_beta_at
        use simple_defs,               only: PROB_BETA_TAIL
        class(commander_refine3D), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(refine3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        type(string)     :: filt_mode_arg, refine_arg
        character(len=STDLEN) :: errmsg
        logical          :: converged, l_anneal
        integer          :: niters, it_tail
        real             :: beta_start
        ! the modes the 3D matcher implements, refused here so that a distributed master never
        ! waits for workers that stop on an unknown mode
        if( cline%defined('refine') )then
            refine_arg = cline%get_carg('refine')
            select case(trim(refine_arg%to_char()))
                case('shc','neigh','greedy','greedy_inpl','eval','prob','prob_neigh','sigma','cont')
                case DEFAULT
                    errmsg = 'refine3D refine='//trim(refine_arg%to_char())//&
                        &' unsupported; one of shc|neigh|greedy|greedy_inpl|eval|prob|prob_neigh|sigma|cont'
                    THROW_HARD(trim(errmsg))
            end select
            call refine_arg%kill
        endif
        ! sanity check: multiple input volumes require nstates > 1
        if( cline%defined('vol2') )then
            if( .not. cline%defined('nstates') )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES is not set; set NSTATES to the number of states')
            else if( cline%get_iarg('nstates') <= 1 )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES <= 1; set NSTATES to the number of states')
            endif
        endif
        if( cline%defined('filt_mode') .and. cline%defined('nstates') )then
            filt_mode_arg = cline%get_carg('filt_mode')
            if( cline%get_iarg('nstates') > 1 .and. trim(filt_mode_arg%to_char()).eq.'uniform' )then
                THROW_HARD('filt_mode=uniform is disabled for multi-state refine3D search')
            endif
            call filt_mode_arg%kill
        endif
        ! local defaults (kept consistent with previous distributed master)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%set('prg', 'refine3D')
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_refine3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        if( params%nstates > 1 .and. trim(params%filt_mode).eq.'uniform' )then
            THROW_HARD('filt_mode=uniform is disabled for multi-state refine3D search')
        endif
        ! per-state M-estimation refines one state on its work project; a multi-state search moves the
        ! labels the frozen weights were computed for
        if( params%l_m_estimator_flex .and. params%nstates > 1 )then
            THROW_HARD('m_estimator=flex refines one state; use refine3D_auto state=X m_estimator=flex')
        endif
        ! Main loop counter semantics:
        !   - params%maxits is the *number of iterations to run* in this invocation.
        !   - params%which_iter starts at params%startit.
        niters            = 0
        params%which_iter = params%startit - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        beta_start = params%prob_beta
        l_anneal   = params%l_prob_align_mode .and. abs(params%prob_beta_end - beta_start) > TINY
        do
            niters            = niters + 1
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            ! the probability tables' posterior draws: the iteration's tempering factor, on the command line for the
            ! workers (check_conv3D reports it in its settings box). A run with prob_beta_end searches flat at
            ! prob_beta until its resolution stops gaining or its budget nears its end (check_conv3D starts the
            ! tail and keeps its first iteration), then ramps geometrically to prob_beta_end over the PROB_BETA_TAIL
            ! iterations of the tail and ends
            if( params%l_prob_align_mode )then
                if( l_anneal )then
                    it_tail = strategy%anneal_tail_start()
                    if( it_tail > 0 )then
                        params%prob_beta = prob_beta_at(params%which_iter - it_tail + 1, PROB_BETA_TAIL,&
                            &beta_start, params%prob_beta_end)
                    else
                        params%prob_beta = beta_start
                    endif
                endif
                call cline%set('prob_beta', params%prob_beta)
            endif
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged .or. niters >= params%maxits ) exit
        end do
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        ! Global teardown (strategies may have built different toolboxes)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call build%pftc%kill
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    !> Distributed worker (single-iteration execution). This should be the command
    !> invoked by the scheduler for each partition.
    subroutine exec_refine3D_distr_worker( self, cline )
        use simple_core_module_api
        use simple_strategy3D_matcher, only: refine3D_exec
        class(commander_refine3D_distr_worker), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        logical          :: l_write_partial_recs
        ! Flags required for worker execution
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        ! Worker needs the alignment toolboxes
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        if( params%which_iter < 1 ) params%which_iter = max(1, params%startit)
        if( .not. cline%defined('extr_iter') ) params%extr_iter = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        l_write_partial_recs = trim(params%volrec) .eq. 'yes'
        call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
    end subroutine exec_refine3D_distr_worker

    !> Sigma2 bootstrap + reconstruction for 3D poses without a consumable sigma2 state (calc_final_rec at a new
    !! sampling, or standalone): image-power seed, euclid ML gridding bootstrap map, one residual pass (refine=sigma)
    !! committing the next generation, then the shipped euclid ML map. Returns vol1..N and which_iter+1 on cline.
    subroutine exec_bootstrap_rec3D( self, cline )
        use simple_commanders_rec,    only: commander_rec3D
        use simple_commanders_euclid, only: commander_calc_pspec
        use simple_sigma2_bootstrap,  only: prepare_pspec_cline, prepare_residual_sigma2_pass_cline
        use simple_solve3D_utils,     only: configure_final_pcg_solve_budget, strip_pcg_backend_keys
        class(commander_bootstrap_rec3D), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(commander_rec3D)      :: xrec3D
        type(commander_calc_pspec) :: xcalc_pspec
        type(commander_refine3D)   :: xrefine3D
        type(cmdline)              :: cline_rec, cline_pspec, cline_sigma
        type(parameters)           :: params
        type(string), allocatable  :: seed_vols(:)
        integer                    :: state, which_iter
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')     ) call cline%set('nstates',          1)
        call warn_for_forced_bootstrap_overrides(cline)
        if( .not. cline%defined('which_iter')  ) call cline%set('which_iter',       1)
        if( .not. cline%defined('postprocess') ) call cline%set('postprocess',  'yes')
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',    'no')
        if( .not. cline%defined('envfsc')      ) call cline%set('envfsc',        'no')
        call cline%delete('objfun')
        call cline%delete('ml_reg')
        call params%new(cline)
        which_iter = max(1, params%which_iter)
        call cline%set('which_iter', which_iter)
        call cline%set('mkdir', 'no') ! child calls must not create nested run directories
        ! 1. One sigma2 basis for every bootstrap: the particle power spectra, what a fresh refinement
        ! seeds from, committed to the canonical state; which_iter only numbers the residual pass.
        call prepare_pspec_cline(cline, params%projfile, 'ptcl3D', which_iter, cline_pspec)
        write(logfhandle,'(A,I0)') '>>> BOOTSTRAP_REC3D SIGMA2 FROM PARTICLE POWER SPECTRA, ITERATION ', which_iter
        call xcalc_pspec%execute(cline_pspec)
        call cline_pspec%kill
        ! 2. Bootstrap map: one euclid ML-regularized gridding assembly on the seed, the residual pass's
        ! reference; it keeps the caller's filt_mode/automsk, as residual sigmas depend on its regularization.
        cline_rec = cline
        call prepare_bootstrap_rec_cline(cline_rec, which_iter, l_final=.false.)
        write(logfhandle,'(A,I0)') '>>> BOOTSTRAP_REC3D: GRIDDING ML-REGULARIZED BOOTSTRAP MAP ON THE IMAGE-POWER SEED, SIGMA ITERATION ', &
            &which_iter
        call xrec3D%execute(cline_rec)
        call cline_rec%kill
        ! 3. No refinement iteration follows, so the image-power seed is
        ! upgraded here: one residual sigma2 pass (refine=sigma: no search,
        ! no volume assembly, no orientation output) against the bootstrap
        ! map at this sampling (simple_sigma2_bootstrap).
        allocate(seed_vols(params%nstates))
        do state = 1, params%nstates
            seed_vols(state) = refine3D_state_vol_fname(state)
        enddo
        call prepare_residual_sigma2_pass_cline(cline, which_iter, params%nstates, seed_vols, cline_sigma)
        write(logfhandle,'(A,I0)') '>>> BOOTSTRAP_REC3D: RESIDUAL SIGMA2 PASS AGAINST THE BOOTSTRAP MAP, ITERATION ', which_iter
        call xrefine3D%execute(cline_sigma)
        call cline_sigma%kill
        do state = 1, params%nstates
            call seed_vols(state)%kill
        enddo
        deallocate(seed_vols)
        ! 4. The residual pass has committed the next canonical generation.
        ! 5. The shipped map: euclid ML-regularized reconstruction on the
        ! residual sigmas, on the caller's backend. A PCG solve at the native
        ! box starts from nothing (no warm start exists at this sampling), so
        ! it gets the cold-solve iteration budget (at least
        ! FINAL_PCG_MAXITS_FLOOR) whoever the caller is.
        cline_rec = cline
        call prepare_bootstrap_rec_cline(cline_rec, which_iter + 1, l_final=.true.)
        write(logfhandle,'(A,I0)') '>>> BOOTSTRAP_REC3D: EUCLID ML-REGULARIZED RECONSTRUCTION ON RESIDUAL SIGMAS, SIGMA ITERATION ', &
            &which_iter + 1
        if( cline_rec%defined('rec_backend') )then
            if( cline_rec%get_carg('rec_backend') == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> BOOTSTRAP_REC3D: FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_rec%get_iarg('maxits_pcg')
        endif
        call xrec3D%execute(cline_rec)
        call cline_rec%kill
        call register_bootstrap_rec_outputs()
        do state = 1, params%nstates
            call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
        enddo
        call cline%set('which_iter', which_iter + 1)
        call simple_end('**** SIMPLE_BOOTSTRAP_REC3D NORMAL STOP ****', print_simple=.false.)

    contains

        !> euclid ML-regularized reconstruct3D on the sigma2 of iter. Bootstrap map (l_final=.false.): gridding,
        !! no postprocess, caller's filt_mode/automsk. Shipped map: caller's backend, postprocess and automsk,
        !! filt_mode=none unless automsk=nu; on PCG the cold-solve budget and the bootstrap map as lag-one vol.
        subroutine prepare_bootstrap_rec_cline( cline_rec, iter, l_final )
            class(cmdline), intent(inout) :: cline_rec
            integer,        intent(in)    :: iter
            logical,        intent(in)    :: l_final
            integer :: istate
            logical :: l_automsk_nu
            l_automsk_nu = .false.
            if( cline_rec%defined('automsk') ) l_automsk_nu = cline_rec%get_carg('automsk') == 'nu'
            if( l_final )then
                ! the shipped map is classical (no nonuniform filtering, a
                ! matching-reference feature) except under automsk=nu, which
                ! is valid only in an NU filt_mode: there the caller's
                ! filt_mode is kept, as on the direct route (simple_final_rec)
                if( .not. l_automsk_nu ) call cline_rec%set('filt_mode', 'none')
                if( cline_rec%defined('rec_backend') )then
                    if( cline_rec%get_carg('rec_backend') == 'pcg' ) &
                        &call configure_final_pcg_solve_budget(cline, cline_rec)
                endif
            else
                call strip_pcg_backend_keys(cline_rec)
                call cline_rec%delete('pcg_mskfile')
                call cline_rec%set('rec_backend', 'gridding')
                call cline_rec%set('postprocess',       'no')
            endif
            call cline_rec%set('prg',       'reconstruct3D')
            call cline_rec%set('mkdir',              'no')
            call cline_rec%set('oritype', params%oritype)
            call cline_rec%set('nstates', params%nstates)
            call cline_rec%set('which_iter',       iter)
            call cline_rec%set('trail_rec',        'no')
            call cline_rec%set('combine_eo',       'no')
            call cline_rec%set('objfun',       'euclid')
            call cline_rec%set('ml_reg',          'yes')
            call cline_rec%delete('refine')
            call cline_rec%delete('update_frac')
            call cline_rec%delete('fillin')
            call cline_rec%delete('ufrac_trec')
            call cline_rec%delete('endit')
            call cline_rec%delete('vol_even')
            call cline_rec%delete('vol_odd')
            call cline_rec%delete('refs')
            call cline_rec%delete('refs_even')
            call cline_rec%delete('refs_odd')
            do istate = 1, params%nstates
                call cline_rec%delete('vol'//int2str(istate))
            enddo
            ! Under automsk the shipped PCG map is estimated on the density envelope like every refinement iteration
            ! (build_pcg_state_support derives it from vol<state>); the step-2 bootstrap map is that reference.
            if( l_final )then
                if( cline_rec%defined('rec_backend') )then
                    if( cline_rec%get_carg('rec_backend') == 'pcg' )then
                        do istate = 1, params%nstates
                            call cline_rec%set('vol'//int2str(istate), refine3D_state_vol_fname(istate))
                        enddo
                    endif
                endif
            endif
        end subroutine prepare_bootstrap_rec_cline

        subroutine register_bootstrap_rec_outputs()
            type(sp_project) :: spproj
            type(string)     :: volname, fscname
            integer          :: istate, pop
            character(len=16) :: imgkind
            call spproj%read_segment('out', params%projfile)
            call spproj%read_segment(params%oritype, params%projfile)
            select case(trim(params%oritype))
                case('cls3D')
                    imgkind = 'vol_cavg'
                case DEFAULT
                    imgkind = 'vol'
            end select
            do istate = 1, params%nstates
                select case(trim(params%oritype))
                    case('cls3D')
                        pop = spproj%os_cls3D%get_pop(istate, 'state')
                    case DEFAULT
                        pop = spproj%os_ptcl3D%get_pop(istate, 'state')
                end select
                if( pop == 0 )cycle
                volname = refine3D_state_vol_fname(istate)
                if( .not. file_exists(volname) )then
                    call volname%kill
                    cycle
                endif
                fscname = refine3D_fsc_fname(istate)
                ! params%box_crop/smpd_crop are the effective reconstruction
                ! sampling resolved by params%new or explicitly pinned by
                ! callers that must avoid staged downsampling leakage.
                call spproj%add_vol2os_out(volname, params%smpd_crop, istate, trim(imgkind), pop=pop)
                if( file_exists(fscname) ) call spproj%add_fsc2os_out(fscname, istate, params%box_crop)
                call volname%kill
                call fscname%kill
            enddo
            call spproj%write_segment_inside('out', params%projfile)
            call spproj%kill
        end subroutine register_bootstrap_rec_outputs

        subroutine warn_for_forced_bootstrap_overrides( cline_in )
            class(cmdline), intent(inout) :: cline_in
            if( cline_in%defined('objfun') ) THROW_WARN('bootstrap_rec3D controls objfun internally; ignoring input objfun')
            if( cline_in%defined('ml_reg') ) THROW_WARN('bootstrap_rec3D controls ml_reg internally; ignoring input ml_reg')
        end subroutine warn_for_forced_bootstrap_overrides

    end subroutine exec_bootstrap_rec3D

    subroutine exec_state_tree( self, cline )
        use simple_state_tree,      only: state_tree
        use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
        class(commander_state_tree), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)          :: params
        type(sp_project)          :: spproj
        type(state_tree)          :: tree
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: labels(:), names(:)
        real(dp),     allocatable :: counts(:)
        type(string) :: mapdir, vol
        real         :: smpd_vol
        integer      :: box_vol, s, nlab, n
        logical      :: l_unfil
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('hp')      ) call cline%set('hp',      40.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',      10.)
        if( .not. cline%defined('linkage') ) call cline%set('linkage', 'complete')
        if( .not. cline%defined('minpop')  ) call cline%set('minpop',  0)
        call cline%set('oritype', 'ptcl3D')
        call params%new(cline)
        call spproj%read_segment('ptcl3D', params%projfile)
        labels = spproj%os_ptcl3D%get_all_asint('state')
        if( cline%defined('dir') )then
            mapdir = params%dir//'/'
        else
            call spproj%read_segment('out', params%projfile)
            call spproj%get_vol('vol', 1, vol, smpd_vol, box_vol)
            mapdir = get_fpath(vol)
        endif
        nlab = maxval(labels)
        if( nlab < 2 ) THROW_HARD('state_tree needs a project with at least two states')
        l_unfil = .true.
        do s = 1, nlab
            if( count(labels == s) == 0 ) cycle
            if( .not. file_exists(mapdir//refine3D_state_halfvol_fname(s, 'even', .false.)) ) cycle
            l_unfil = l_unfil .and. file_exists(mapdir//refine3D_state_halfvol_fname(s, 'even', .true.)) .and. &
                &file_exists(mapdir//refine3D_state_halfvol_fname(s, 'odd', .true.))
        end do
        allocate(evens(nlab), odds(nlab), names(nlab), counts(nlab))
        n = 0
        do s = 1, nlab
            if( count(labels == s) == 0 ) cycle
            if( .not. file_exists(mapdir//refine3D_state_halfvol_fname(s, 'even', l_unfil)) ) cycle
            if( .not. file_exists(mapdir//refine3D_state_halfvol_fname(s, 'odd',  l_unfil)) ) cycle
            n = n + 1
            evens(n)  = mapdir//refine3D_state_halfvol_fname(s, 'even', l_unfil)
            odds(n)   = mapdir//refine3D_state_halfvol_fname(s, 'odd',  l_unfil)
            names(n)  = s
            counts(n) = real(count(labels == s), dp)
        end do
        write(logfhandle,'(A)') '>>> STATE_TREE half-maps '//merge('unfiltered','filtered  ',l_unfil)//' from '//mapdir%to_char()
        call tree%new(evens(1:n), odds(1:n), names(1:n), counts(1:n), find_img_smpd(evens(1)), params%hp, params%lp)
        call tree%build(params%linkage, real(params%minpop))
        call tree%write(labels)
        call tree%kill
        call spproj%kill
        call simple_end('**** SIMPLE_STATE_TREE NORMAL STOP ****')
    end subroutine exec_state_tree

    subroutine exec_state_mergeup( self, cline )
        use simple_state_tree,    only: state_tree
        use simple_state_mergeup, only: state_mergeup, EM_LIM
        use simple_eo_partition,  only: ensure_field_halves
        class(commander_state_mergeup), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(sp_project)     :: spproj
        type(state_tree)     :: tree
        type(state_mergeup)  :: mergeup
        integer, allocatable :: labels(:), units(:)
        real    :: smpd
        integer :: box, i
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp')    ) call cline%set('hp',    40.)
        if( .not. cline%defined('lp')    ) call cline%set('lp',    10.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('box_crop') )then
            if( .not. cline%defined('projfile') ) THROW_HARD('state_mergeup needs projfile (or box_crop)')
            call spproj%read_segment('stk', cline%get_carg('projfile'))
            box  = spproj%get_box()
            smpd = spproj%get_smpd()
            call spproj%kill
            call cline%set('box_crop', min(box, max(26, 2*ceiling(real(box)*smpd/cline%get_rarg('lp')))))
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call tree%read(params%infile)
        labels = build%spproj_field%get_all_asint('state')
        call mergeup%new(params, build, EM_LIM)
        if( cline%defined('infile2') ) call mergeup%leaf_focus(params%infile2)
        call mergeup%run(tree, labels, units)
        do i = 1, size(units)
            call build%spproj_field%set_state(i, units(i))
        end do
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call mergeup%kill
        call tree%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_STATE_MERGEUP NORMAL STOP ****')
    end subroutine exec_state_mergeup

    subroutine exec_state_split( self, cline )
        use simple_state_mergeup, only: state_mergeup, EM_LIM, locality_guard, effect_floor, write_split_stats
        use simple_eo_partition,  only: ensure_field_halves
        class(commander_state_split), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(sp_project)     :: spproj
        type(state_mergeup)  :: split
        integer, allocatable :: labels(:), rows(:), child(:), junk(:), out(:)
        type(string), allocatable :: qlines(:), plines(:)
        real,         allocatable :: eul_in(:,:)
        real    :: smpd, st(6)
        integer :: box, i, kstar, kpre, funit, ios
        logical :: l_q, l_ex
        character(len=8) :: envv
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp')    ) call cline%set('hp',    40.)
        if( .not. cline%defined('lp')    ) call cline%set('lp',    10.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('box_crop') )then
            if( .not. cline%defined('projfile') ) THROW_HARD('state_split needs projfile (or box_crop)')
            call spproj%read_segment('stk', cline%get_carg('projfile'))
            box  = spproj%get_box()
            smpd = spproj%get_smpd()
            call spproj%kill
            call cline%set('box_crop', min(box, max(26, 2*ceiling(real(box)*smpd/cline%get_rarg('lp')))))
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        labels = build%spproj_field%get_all_asint('state')
        call get_environment_variable('SIMPLE_SPLIT_CONTRAST', envv, status=ios)
        if( ios /= 0 ) envv = ''
        call split%new(params, build, EM_LIM, contrast=any(trim(envv) == ['yes  ', 'maps ']), &
            &amplitude=any(trim(envv) == ['yes  ', 'estep']))
        if( cline%defined('pcg_mskfile') ) call split%focus(params%pcg_mskfile)
        if( params%merge_insert_audit == 'yes' .and. params%merge_audit_request%strlen() > 0 )then
            call split%audit_em_wrapper()
            call split%kill()
            call build%kill_general_tbox()
            call simple_end('**** SIMPLE_EM_WRAPPER_CHECK NORMAL STOP ****')
            return
        endif
        rows = split%get_rows()
        call split%test_once(labels(rows), maxval(labels), child, kstar, junk, params%split_ntests, stats=st)
        ! the quality flag describes the split the test proposed, before the guard; it gates nothing
        call split%quality_flag(child, kstar, qlines, l_q)
        kpre = kstar
        allocate(out(size(labels)), source=0)
        out(rows) = child
        call fopen(funit, string('split_children_pre.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('state_split: cannot write split_children_pre.txt', ios)
        do i = 1, size(out)
            write(funit,'(I0)') out(i)
        end do
        call fclose(funit)
        deallocate(out)
        call locality_guard(labels(rows), child, kstar, params%mskdiam, split)
        call effect_floor(labels(rows), child, kstar, params%mskdiam)
        call write_split_stats('split_stats.txt', st, kstar, params%split_ntests, kpre)
        call split%write_weights('ptcl_weights.txt')
        inquire(file='../n.simple', exist=l_ex)
        if( l_ex )then
            call spproj%read_segment('ptcl3D', string('../n.simple'))
            if( spproj%os_ptcl3D%get_noris() == size(labels) )then
                allocate(eul_in(3,size(labels)))
                do i = 1, size(labels)
                    eul_in(:,i) = spproj%os_ptcl3D%get_euler(i)
                end do
            endif
            call spproj%kill
        endif
        call split%pose_lines(labels(rows), child, kstar, eul_in, trim(params%pgrp), plines)
        call fopen(funit, string('split_stats.txt'), status='old', action='write', position='append', iostat=ios)
        call fileiochk('state_split: cannot append to split_stats.txt', ios)
        do i = 1, size(qlines)
            write(funit,'(A)') '# '//qlines(i)%to_char()
        end do
        do i = 1, size(plines)
            write(funit,'(A)') plines(i)%to_char()
            write(logfhandle,'(A)') plines(i)%to_char()
        end do
        call fclose(funit)
        allocate(out(size(labels)), source=0)
        out(rows) = child
        call fopen(funit, string('split_children.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('state_split: cannot write split_children.txt', ios)
        do i = 1, size(out)
            write(funit,'(I0)') out(i)
        end do
        call fclose(funit)
        call fopen(funit, string('split_junk.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('state_split: cannot write split_junk.txt', ios)
        do i = 1, size(junk)
            write(funit,'(I0)') junk(i)
        end do
        call fclose(funit)
        call split%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_STATE_SPLIT NORMAL STOP ****')
    end subroutine exec_state_split

    subroutine exec_site_rank( self, cline )
        use simple_state_mergeup,     only: state_mergeup, merge_poles, merge_poles_subspace, gram_corr, gram_resid, fisher_z, &
                                          &ivec, HALVES_KSIG, mask_neff, seed_poles
        use simple_eo_partition,      only: ensure_field_halves
        use simple_srch_sort_loc,     only: hpsort
        use simple_string_utils,      only: list_of_ints2arr
        use simple_refine3D_strategy, only: PROB_MULTISTATE_MIN_POP
        use simple_rnd,               only: ran3
        class(commander_site_rank), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        character(len=2), parameter :: TAILS(2) = ['lo', 'hi']
        type(parameters)     :: params
        type(builder)        :: build
        type(sp_project)     :: spproj
        type(state_mergeup)  :: site
        type(string), allocatable :: names(:)
        real,    allocatable :: z(:,:), post(:), curve(:,:), posts(:,:), key(:), sp(:,:), posth(:), curveh(:,:)
        real,    allocatable :: cls2(:), hpk(:), zself(:), zdist(:), zred(:), rho(:,:), zres(:,:)
        real(dp),allocatable :: gee(:,:), goo(:,:), geo(:,:)
        integer, allocatable :: rows(:), ord(:), top(:), eo(:), ordh(:), sel(:), dirs(:), grp(:), seeds(:), nsx(:)
        integer, allocatable :: fidx(:), ntopx(:), epmx(:)
        logical, allocatable :: seed(:), fit(:), starts(:,:), took(:)
        complex, allocatable :: devs(:,:,:,:), dfull(:,:,:,:,:)
        character(len=:), allocatable :: tag, members
        character(len=2048) :: line
        real    :: smpd
        integer :: box, i, j, itail, nstar, ntail, n, funit, ios, nrows, ih, idir, ie, h, ne, nstarq, nfit, ntop, band(2), minseed
        real    :: neff
        logical :: l_halves
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    10.)
        call cline%set('hp',      40.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('box_crop') )then
            call spproj%read_segment('stk', cline%get_carg('projfile'))
            box  = spproj%get_box()
            smpd = spproj%get_smpd()
            call spproj%kill
            call cline%set('box_crop', min(box, max(26, 2*ceiling(real(box)*smpd/cline%get_rarg('lp')))))
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        nrows = build%spproj_field%get_noris()
        call site%new(params, build, em_lim=15.)
        call site%fit_contrast
        rows  = site%get_rows()
        n     = size(rows)
        if( cline%defined('site_dirs') )then
            dirs = list_of_ints2arr(trim(params%site_dirs))
            if( size(dirs) == 0 ) THROW_HARD('site_rank: no direction in site_dirs')
            call read_latent
        else if( cline%defined('infile') )then
            allocate(dirs(1), source=0)
            allocate(z(nrows,1))
            call fopen(funit, params%infile, status='old', action='read', iostat=ios)
            call fileiochk('site_rank: cannot read '//params%infile%to_char(), ios)
            do i = 1, nrows
                read(funit,*) z(i,1)
            end do
            call fclose(funit)
        else
            allocate(dirs(0), z(nrows,0))
        endif
        allocate(fit(n), source=.true.)
        if( params%nsample > 0 .and. params%nsample < n )then
            fidx = [(i, i=1,n)]
            do i = 1, params%nsample
                j       = min(n, i + int(ran3() * real(n - i + 1)))
                h       = fidx(i)
                fidx(i) = fidx(j)
                fidx(j) = h
            end do
            fit = .false.
            fit(fidx(1:params%nsample)) = .true.
        endif
        nfit = count(fit)
        if( nfit < n ) write(logfhandle,'(A)') '>>> SITE_RANK the EMs fit on '//int2str(nfit)//' of the '//int2str(n)//&
            &' particles, drawn at random; every particle is scored with the converged maps'
        eo    = [(build%spproj_field%get_eo(rows(i)), i=1,n)]
        if( params%site_size .eq. 'yes' )then
            call size_rankings
            call site%kill
            call build%kill_general_tbox
            call simple_end('**** SIMPLE_SITE_RANK NORMAL STOP ****')
            return
        endif
        ntail = min(15000, nfit / 6)
        ! the floor refine3D itself needs, applied to every state site_rank writes (poles, the rest and the pole rule's rest bucket)
        minseed = PROB_MULTISTATE_MIN_POP + 1
        call starts_of
        ne = size(names)
        if( ne == 0 ) THROW_HARD('site_rank: no start; give infile (with site_dirs) or infile2')
        allocate(posts(n,ne), nsx(ne), ntopx(ne), cls2(ne), hpk(ne), epmx(ne))
        l_halves = params%pole_merge <= 0.
        do ie = 1, ne
            tag  = names(ie)%to_char()
            seed = starts(:,ie)
            write(logfhandle,'(A)') '>>> SITE_RANK '//tag//' start: '//int2str(count(seed))//' of '//int2str(n)//' particles'
            call site%site_em(seed, fit, params%pcg_mskfile, post, nstar, curve, devs, [real ::], nfit < n)
            if( .not. allocated(dfull) ) allocate(dfull(size(devs,1),size(devs,2),size(devs,3),2,ne))
            dfull(:,:,:,:,ie) = devs
            posts(:,ie) = post
            nsx(ie)     = nstar
            ntop        = nint(real(nstar) * real(n) / real(nfit))
            ntopx(ie)   = ntop
            cls2(ie)    = sum(post) / real(n)
            hpk(ie)     = maxval(curve(2,:))
            epmx(ie)    = nint(maxval(curve(4,:)) * real(n) / real(nfit))
            allocate(top(nrows), source=0)
            do ih = 0, 1
                sel  = pack([(i, i=1,n)], eo == ih)
                ordh = [(i, i=1,size(sel))]
                key  = -post(sel)
                call hpsort(key, ordh)
                top(rows(sel(ordh(1:min(ntop/2, size(sel)))))) = 1
            end do
            call write_rows('site_post_'//tag//'.txt', post)
            call write_ints('site_top_'//tag//'.txt', top)
            call fopen(funit, string('site_curve_'//tag//'.txt'), status='replace', action='write', iostat=ios)
            call fileiochk('site_rank: cannot write site_curve_'//tag//'.txt', ios)
            write(funit,'(A)') '# N holdout_corr holdout_amplitude equivalent_pure_members'
            do i = 1, size(curve,2)
                write(funit,'(I0,1X,F8.4,1X,ES12.5,1X,F10.1)') nint(curve(1,i)), curve(2,i), curve(3,i), curve(4,i)
            end do
            call fclose(funit)
            if( nfit < n )then
                write(logfhandle,'(A,A,A,I0,A,I0,A,F6.3)') '>>> SITE_RANK ', tag, ' pole: N* ', nstar, ' of the fit, ', ntop, &
                    &' over the node; site class share ', cls2(ie)
            else
                write(logfhandle,'(A,A,A,I0,A,F7.1)') '>>> SITE_RANK ', tag, ' pole: N* ', nstar, ', site class weight ', sum(post)
            endif
            deallocate(top)
        end do
        call merge_poles(posts, merge(2., params%pole_merge, l_halves), PROB_MULTISTATE_MIN_POP + 1, sp, grp, seeds)
        band = site%holdout_band()
        neff = mask_neff(params%pcg_mskfile, params%smpd, params%box_crop, params%smpd_crop, params%msk_crop, band(1), band(2))
        if( l_halves )then
            call site%dev_gram(dfull, gee, goo, geo)
            ! a class that grew past half the node has found no minority state: it is no pole and claims no seeds
            took = cls2 > 0.5
            call merge_poles_subspace(gee, goo, geo, neff, posts, minseed, grp, seeds, zself, zdist, rho, zres, zred, l_ok=.not. took)
            call write_subspace
            if( maxval([0, seeds]) > maxval([0, grp]) ) call test_rest
            call fclose(funit)
            if( maxval([0, grp]) > 0 ) call pole_rule
        endif
        do i = 1, size(names)
            if( grp(i) /= 0 ) cycle
            if( l_halves )then
                if( i <= ne )then
                    if( took(i) )then
                        write(logfhandle,'(A)') '>>> SITE_RANK endpoint '//names(i)%to_char()//' rejected: its site class grew '//&
                            &'past half the node (share '//fstr(cls2(i))//'), no minority state'
                        cycle
                    endif
                endif
                if( zself(i) < HALVES_KSIG )then
                    write(logfhandle,'(A)') '>>> SITE_RANK endpoint '//names(i)%to_char()//' rejected: its cross-half deviation '//&
                        &'does not reproduce (z '//fstr(zself(i))//' below '//fstr(HALVES_KSIG)//')'
                    cycle
                else if( zred(i) < HALVES_KSIG )then
                    write(logfhandle,'(A)') '>>> SITE_RANK endpoint '//names(i)%to_char()//' redundant: in the span of the other '//&
                        &'poles (residual z '//fstr(zred(i))//')'
                    cycle
                endif
            endif
            write(logfhandle,'(A)') '>>> SITE_RANK endpoint '//names(i)%to_char()//' dropped (technical): '//&
                &'its pole drew fewer than '//int2str(minseed)//' seeds, the fewest refine3D refines as a state'
        end do
        call fopen(funit, string('site_neff.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('site_rank: cannot write site_neff.txt', ios)
        write(funit,'(A)') '# box_crop smpd_crop msk_crop hold-out shells, effective independent coefficients of the masked band'
        write(funit,'(I0,2(1X,F9.5),2(1X,I0),1X,F12.3)') params%box_crop, params%smpd_crop, params%msk_crop, band, neff
        call fclose(funit)
        call fopen(funit, string('site_starts.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('site_rank: cannot write site_starts.txt', ios)
        write(funit,'(A)') '# start, seed particles, N*, pole it reached (0 = rejected or dropped), N* over the node, site class '//&
            &'share, hold-out peak and its SE ((1 - r^2) / sqrt(n_eff - 3), n_eff in site_neff.txt), equivalent pure members '//&
            &'over the node'
        do i = 1, ne
            write(funit,'(A,4(1X,I0),3(1X,F7.4),1X,I0)') names(i)%to_char(), count(starts(:,i)), nsx(i), grp(i), ntopx(i), &
                &cls2(i), hpk(i), (1. - hpk(i)**2) / sqrt(max(neff - 3., 1.)), epmx(i)
        end do
        call fclose(funit)
        call fopen(funit, string('site_spearman.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('site_rank: cannot write site_spearman.txt', ios)
        do i = 1, ne
            do j = i + 1, ne
                write(logfhandle,'(A,F7.3)') '>>> SITE_RANK Spearman '//names(i)%to_char()//' vs '//names(j)%to_char()//' ', sp(i,j)
                write(funit,'(A,F7.3)') names(i)%to_char()//' '//names(j)%to_char()//' ', sp(i,j)
            end do
        end do
        call fclose(funit)
        if( params%junk_screen .eq. 'yes' ) call site%junk_screen(seeds)
        allocate(top(nrows), source=0)
        top(rows) = seeds
        call write_ints('site_seeds.txt', top)
        call fopen(funit, string('site_poles.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('site_rank: cannot write site_poles.txt', ios)
        do j = 1, maxval(grp)
            members = ''
            do i = 1, size(names)
                if( grp(i) == j ) members = members//' '//names(i)%to_char()
            end do
            write(funit,'(I0,1X,I0,A)') j, count(seeds == j), members
            write(logfhandle,'(A)') '>>> SITE_RANK pole '//int2str(j)//':'//members//', '//int2str(count(seeds == j))//' seeds'
        end do
        j = maxval([0, seeds])
        if( j > maxval([0, grp]) )then
            write(funit,'(I0,1X,I0,A)') j, count(seeds == j), ' rest'
            write(logfhandle,'(A)') '>>> SITE_RANK rest state '//int2str(j)//': '//int2str(count(seeds == j))//&
                &' seeds (best pole posterior < 0.5)'
        endif
        call fclose(funit)
        call site%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_SITE_RANK NORMAL STOP ****')

    contains

        ! site_halves.txt in the controller's format: endpoint rows (name, N* twice, Spearman-Brown reliability twice | its even
        ! half's correlation with every endpoint's odd half | self z | pole), merge lines and accept lines (row = reproducibility,
        ! column = distinctness); the endpoint Grams go to site_gram.txt
        subroutine write_subspace
            integer :: a, b, fu, io
            real    :: r, rel
            call fopen(funit, string('site_halves.txt'), status='replace', action='write', iostat=ios)
            call fileiochk('site_rank: cannot write site_halves.txt', ios)
            write(funit,'(A)') '# endpoint, N* (twice), Spearman-Brown reliability of its deviation (twice), correlation of its even '//&
                &'half with every endpoint''s odd half, self z (Fisher, n_eff '//fstr(neff)//'), pole (0 = rejected or dropped)'
            do a = 1, ne
                r   = gram_corr(gee, goo, geo, a, a)
                rel = 2. * r / max(1. + r, 1.e-6)
                write(line,'(1X,I0,1X,I0,2(1X,F6.3))') nsx(a), nsx(a), rel, rel
                members = names(a)%to_char()//trim(line)//' |'
                do b = 1, ne
                    write(line,'(1X,F6.3)') gram_corr(gee, goo, geo, a, b)
                    members = members//trim(line)
                end do
                write(line,'(F8.3)') min(zself(a), 999.)
                members = members//' | '//trim(adjustl(line))//' | '//int2str(grp(a))
                write(funit,'(A)') members
                write(logfhandle,'(A)') '>>> SITE_RANK subspace '//members
            end do
            write(funit,'(A)') '# merge tests: endpoint pair, disattenuated correlation, residual z of the first after the second and '//&
                &'of the second after the first, their max (one pole when the correlation is positive and both are below '//&
                &int2str(nint(HALVES_KSIG))//'; groups merge closest first when every cross pair is such a pair), merged'
            do a = 1, ne
                do b = a + 1, ne
                    write(line,'(4(1X,F0.3))') rho(a,b), min(zres(a,b), 999.), min(zres(b,a), 999.), min(max(zres(a,b), zres(b,a)), 999.)
                    members = 'merge '//names(a)%to_char()//' '//names(b)%to_char()//trim(line)//' '//&
                        &trim(merge('yes', 'no ', grp(a) > 0 .and. grp(a) == grp(b)))
                    write(funit,'(A)') members
                    write(logfhandle,'(A)') '>>> SITE_RANK subspace '//members
                end do
            end do
            write(funit,'(A)') '# accept tests: endpoint, row = its cross-half correlation, Fisher SE and z (accepted at '//&
                &int2str(nint(HALVES_KSIG))//'); column = the smallest residual z against a positively correlated other pole '//&
                &'(twice, SE 1; its own z when there is none), accepted (reproducible, not in the span of the other poles and '//&
                &'its site class under half the node)'
            do a = 1, ne
                call accept_lines(names(a)%to_char(), gram_corr(gee, goo, geo, a, a), zself(a), zdist(a), &
                    &zself(a) >= HALVES_KSIG .and. zred(a) >= HALVES_KSIG .and. .not. took(a))
            end do
            call fopen(fu, string('site_gram.txt'), status='replace', action='write', iostat=io)
            call fileiochk('site_rank: cannot write site_gram.txt', io)
            write(fu,'(A)') '# hold-out band Gram matrices of the endpoint deviation maps: n_eff, shells, names, then even.even, '//&
                &'odd.odd, even.odd (rows x columns)'
            write(fu,'(F12.3,2(1X,I0))') neff, band
            members = ''
            do a = 1, ne
                members = members//' '//names(a)%to_char()
            end do
            write(fu,'(A)') trim(adjustl(members))
            do a = 1, ne
                write(fu,'(*(1X,ES14.6))') gee(a,:)
            end do
            do a = 1, ne
                write(fu,'(*(1X,ES14.6))') goo(a,:)
            end do
            do a = 1, ne
                write(fu,'(*(1X,ES14.6))') geo(a,:)
            end do
            call fclose(fu)
        end subroutine write_subspace

        subroutine accept_lines( nm, r, z, zd, l_ok )
            character(len=*), intent(in) :: nm
            real,             intent(in) :: r, z, zd
            logical,          intent(in) :: l_ok
            write(line,'(3(1X,F0.3))') r, 1. / sqrt(max(neff - 3., 1.)), min(z, 999.)
            write(funit,'(A)') 'accept '//nm//' row'//trim(line)//' '//trim(merge('yes', 'no ', l_ok))
            write(logfhandle,'(A)') '>>> SITE_RANK subspace accept '//nm//' row'//trim(line)//' '//trim(merge('yes', 'no ', l_ok))
            write(line,'(3(1X,F0.3))') min(zd, 999.), 1., min(zd, 999.)
            write(funit,'(A)') 'accept '//nm//' column'//trim(line)//' '//trim(merge('yes', 'no ', l_ok))
            write(logfhandle,'(A)') '>>> SITE_RANK subspace accept '//nm//' column'//trim(line)//' '//trim(merge('yes', 'no ', l_ok))
        end subroutine accept_lines

        function fstr( x ) result( s )
            real, intent(in) :: x
            character(len=:), allocatable :: s
            character(len=32) :: b
            write(b,'(F0.2)') x
            s = trim(adjustl(b))
        end function fstr

        ! mean disattenuated correlation of map a with the maps in mem, from the Grams g2e / g2o / g2x
        real function dis_corr( g2e, g2o, g2x, a, mem )
            real(dp), intent(in) :: g2e(:,:), g2o(:,:), g2x(:,:)
            integer,  intent(in) :: a, mem(:)
            integer :: k, f
            dis_corr = 0.
            do k = 1, size(mem)
                f = mem(k)
                dis_corr = dis_corr + 0.5 * (gram_corr(g2e, g2o, g2x, a, f) + gram_corr(g2e, g2o, g2x, f, a)) / &
                    &sqrt(max(gram_corr(g2e, g2o, g2x, a, a) * gram_corr(g2e, g2o, g2x, f, f), 1.e-12)) / real(size(mem))
            end do
        end function dis_corr

        ! the distinctness z of map a against the direction of the maps in mem: the larger of a's residual after them and the
        ! smallest residual of one of them after a (both must fail for a to be the same direction)
        real function dist_z( g2e, g2o, g2x, a, mem )
            real(dp), intent(in) :: g2e(:,:), g2o(:,:), g2x(:,:)
            integer,  intent(in) :: a, mem(:)
            type(ivec) :: dm(1), da(1)
            real    :: z2
            integer :: k
            dm(1)%v = mem
            da(1)%v = [a]
            dist_z  = fisher_z(gram_resid(g2e, g2o, g2x, a, dm), neff - 1.)
            z2      = huge(1.)
            do k = 1, size(mem)
                z2 = min(z2, fisher_z(gram_resid(g2e, g2o, g2x, mem(k), da), neff - 1.))
            end do
            dist_z = max(dist_z, z2)
        end function dist_z

        ! the rest (no accepted pole claims it at 0.5), ranked by 1 - its best pole posterior and scored once: a state when its
        ! deviation reproduces and is distinct from every pole it correlates with positively; else every particle joins its best pole
        subroutine test_rest
            type(string), allocatable :: nm2(:)
            complex,      allocatable :: d2(:,:,:,:,:)
            real(dp),     allocatable :: g2e(:,:), g2o(:,:), g2x(:,:), gp(:,:), gp2(:,:)
            real,         allocatable :: rscore(:), p2(:,:)
            integer,      allocatable :: plab(:), mem(:)
            real    :: zs, zd, r
            integer :: m, np, g, f
            logical :: l_ok
            m  = size(names)
            np = maxval(grp)
            allocate(gp(n,np), source=0._dp)
            do f = 1, m
                if( grp(f) > 0 ) gp(:,grp(f)) = gp(:,grp(f)) + real(posts(:,f),dp) / real(count(grp == grp(f)),dp)
            end do
            rscore = real(1._dp - maxval(gp, dim=2))
            write(logfhandle,'(A)') '>>> SITE_RANK rest: '//int2str(count(seeds == np + 1))//' particles at best pole posterior '//&
                &'< 0.5, ranked by 1 - that posterior and scored once'
            call site%site_em(fit, fit, params%pcg_mskfile, posth, nstarq, curveh, devs, rscore, .false.)
            allocate(d2(size(dfull,1),size(dfull,2),size(dfull,3),2,ne+1))
            d2(:,:,:,:,1:ne) = dfull
            d2(:,:,:,:,ne+1) = devs
            call site%dev_gram(d2, g2e, g2o, g2x)
            r  = gram_corr(g2e, g2o, g2x, ne+1, ne+1)
            zs = fisher_z(r, neff)
            zd = huge(1.)
            do g = 1, np
                mem = pack([(f, f=1,ne)], grp(1:ne) == g)
                if( size(mem) == 0 ) cycle
                if( dis_corr(g2e, g2o, g2x, ne+1, mem) <= 0. ) cycle
                zd = min(zd, dist_z(g2e, g2o, g2x, ne+1, mem))
            end do
            if( zd == huge(1.) ) zd = zs
            l_ok = zs >= HALVES_KSIG .and. zd >= HALVES_KSIG
            write(line,'(1X,I0,1X,I0,2(1X,F6.3))') nstarq, nstarq, 2. * r / max(1. + r, 1.e-6), 2. * r / max(1. + r, 1.e-6)
            members = 'rest'//trim(line)//' |'
            do f = 1, ne
                write(line,'(1X,F6.3)') gram_corr(g2e, g2o, g2x, ne+1, f)
                members = members//trim(line)
            end do
            write(line,'(F8.3)') min(zs, 999.)
            members = members//' | '//trim(adjustl(line))//' | '//int2str(merge(np + 1, 0, l_ok))
            write(funit,'(A)') members
            write(logfhandle,'(A)') '>>> SITE_RANK subspace '//members
            call accept_lines('rest', r, zs, zd, l_ok)
            write(logfhandle,'(A)') '>>> SITE_RANK rest state (best pole posterior < 0.5): '//trim(merge('accepted as a state       ', &
                &'rejected, to its best pole', l_ok))
            nm2 = [names, string('rest')]
            call move_alloc(nm2, names)
            allocate(p2(n,m+1))
            p2(:,1:m) = posts
            p2(:,m+1) = rscore
            call move_alloc(p2, posts)
            zself = [zself, zs]
            zdist = [zdist, zd]
            zred  = [zred, huge(1.)]
            if( l_ok )then
                allocate(gp2(n,np+1))
                gp2(:,1:np) = gp
                gp2(:,np+1) = real(rscore,dp)
                call seed_poles(gp2, minseed, .false., seeds, plab)
                grp = [grp, np + 1]
            else
                ! a rejected minority rest joins its best pole; a majority rest stays the complement of the poles
                call seed_poles(gp, minseed, real(count(seeds == np + 1)) >= 0.5 * real(n), seeds, plab)
                grp = [grp, 0]
            endif
            do f = 1, m + 1
                if( grp(f) > 0 ) grp(f) = plab(grp(f))
            end do
        end subroutine test_rest

        ! a pole keeps its seeds when they hold under half of the node and their own deviation (scored by the pole's posterior on
        ! them, no EM) reproduces and is distinct from every other pole's direction it correlates with positively; its strength
        ! (the margin strongest_pole ranks sites by) is its weakest endpoint's hold-out reproducibility and distinctness
        subroutine pole_rule
            complex,  allocatable :: d2(:,:,:,:,:)
            real(dp), allocatable :: gp(:,:), gpp(:,:), g2e(:,:), g2o(:,:), g2x(:,:)
            real,     allocatable :: sc(:), shs(:), zrs(:), zcs(:), zse(:), zde(:)
            integer,  allocatable :: plab(:), sd(:), kid(:), cnt0(:), fin(:), pos(:), mem(:), newidx(:)
            logical,  allocatable :: keep(:), lpole(:)
            real    :: share, zr, zc
            integer :: g, h, np, m, k, f, fu, io, nmap
            np = maxval(grp)
            m  = size(names)
            allocate(keep(np), source=.true.)
            allocate(shs(np), zrs(np), zcs(np), zse(np), zde(np), source=0.)
            allocate(cnt0(np), fin(np), pos(np), source=0)
            lpole = [(any([(grp(f) == g .and. names(f)%to_char() /= 'rest', f=1,m)]), g=1,np)]
            allocate(gpp(n,np), source=0._dp)
            do f = 1, m
                if( grp(f) > 0 ) gpp(:,grp(f)) = gpp(:,grp(f)) + real(posts(:,f),dp) / real(count(grp == grp(f)),dp)
            end do
            ! one scoring pass per minority pole; its map joins the endpoint maps in the Grams
            allocate(d2(size(dfull,1),size(dfull,2),size(dfull,3),2,ne+np))
            d2(:,:,:,:,1:ne) = dfull
            nmap = ne
            do g = 1, np
                if( .not. lpole(g) ) cycle
                if( real(count(seeds == g)) / real(n) >= 0.5 ) cycle
                sc = real(gpp(:,g))
                where( seeds /= g ) sc = 0.
                write(logfhandle,'(A)') '>>> SITE_RANK pole rule: pole '//int2str(g)//' scored on its '//int2str(count(seeds == g))//&
                    &' seeds'
                call site%site_em(fit, fit, params%pcg_mskfile, posth, nstarq, curveh, devs, sc, .false.)
                nmap = nmap + 1
                d2(:,:,:,:,nmap) = devs
                pos(g) = nmap
            end do
            call site%dev_gram(d2(:,:,:,:,1:nmap), g2e, g2o, g2x)
            do g = 1, np
                if( .not. lpole(g) ) cycle
                share = real(count(seeds == g)) / real(n)
                zr    = -huge(zr)
                zc    = -huge(zc)
                if( pos(g) > 0 )then
                    zr = fisher_z(gram_corr(g2e, g2o, g2x, pos(g), pos(g)), neff)
                    zc = huge(1.)
                    do h = 1, np
                        if( h == g .or. .not. lpole(h) ) cycle
                        if( pos(h) > 0 )then
                            mem = [pos(h)]
                        else
                            mem = pack([(f, f=1,ne)], grp(1:ne) == h)
                        endif
                        if( size(mem) == 0 ) cycle
                        if( dis_corr(g2e, g2o, g2x, pos(g), mem) <= 0. ) cycle
                        zc = min(zc, dist_z(g2e, g2o, g2x, pos(g), mem))
                    end do
                    if( zc == huge(1.) ) zc = zr
                endif
                keep(g) = share < 0.5 .and. zr >= HALVES_KSIG .and. zc >= HALVES_KSIG
                shs(g)  = share
                zse(g)  = max(min(zr, 999.), -999.)
                zde(g)  = max(min(zc, 999.), -999.)
                zrs(g)  = min(minval(zself(1:ne), mask=grp(1:ne) == g), 999.)
                zcs(g)  = min(minval(zdist(1:ne), mask=grp(1:ne) == g), 999.)
                cnt0(g) = count(seeds == g)
                write(line,'(A,I0,A,I0,A,F6.3,A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)') '>>> SITE_RANK pole rule: pole ', g, ', ', &
                    &count(seeds == g), ' seeds (', share, ' of the node), seed deviation z ', zse(g), ' / distinct z ', zde(g), &
                    &'; strength (weakest endpoint hold-out z / distinct z) ', zrs(g), ' / ', zcs(g), ': '//&
                    &trim(merge('kept    ', 'rejected', keep(g)))
                write(logfhandle,'(A)') trim(line)
            end do
            if( all(keep) )then
                fin = [(g, g=1,np)]
            else
                kid = pack([(g, g=1,np)], keep)
                if( .not. any(keep .and. lpole) )then
                    seeds = 1
                    grp   = 0
                    write(logfhandle,'(A)') '>>> SITE_RANK pole rule: no pole kept, one state'
                else
                    ! a kept pole keeps the seed rows it was tested on; a rejected pole's rows join their best kept pole at
                    ! posterior >= 0.5, else the rest
                    kid = pack([(g, g=1,np)], keep .and. lpole)
                    allocate(gp(n,size(kid)), source=0._dp)
                    do k = 1, size(kid)
                        do f = 1, m
                            if( grp(f) == kid(k) ) gp(:,k) = gp(:,k) + real(posts(:,f),dp) / real(count(grp == kid(k)),dp)
                        end do
                    end do
                    sd   = [(0, f=1,n)]
                    plab = [(0, g=1,np)]
                    do k = 1, size(kid)
                        plab(kid(k)) = k
                    end do
                    do f = 1, n
                        if( seeds(f) <= 0 ) cycle
                        if( plab(seeds(f)) > 0 )then
                            sd(f) = plab(seeds(f))
                        else if( lpole(seeds(f)) )then
                            k = maxloc(gp(f,:), dim=1)
                            sd(f) = size(kid) + 1
                            if( gp(f,k) >= 0.5_dp ) sd(f) = k
                        else
                            sd(f) = size(kid) + 1
                        endif
                    end do
                    do f = 1, m
                        if( grp(f) == 0 ) cycle
                        if( plab(grp(f)) > 0 )then
                            grp(f) = plab(grp(f))
                        else if( lpole(grp(f)) )then
                            grp(f) = 0
                        else
                            grp(f) = merge(size(kid) + 1, 0, any(sd == size(kid) + 1))
                        endif
                    end do
                    seeds = sd
                    fin   = plab
                endif
            endif
            call fold_small_states(newidx)
            do g = 1, np
                if( fin(g) > 0 )then
                    if( newidx(fin(g)) == 0 ) keep(g) = .false.
                    fin(g) = newidx(fin(g))
                endif
            end do
            call fopen(fu, string('site_pole_rule.txt'), status='replace', action='write', iostat=io)
            call fileiochk('site_rank: cannot write site_pole_rule.txt', io)
            write(fu,'(A)') '# pole, seeds, share of the node, strength = its weakest endpoint''s hold-out reproducibility z and '//&
                &'distinctness z (strongest_pole ranks by their min), its state in site_seeds.txt after the rule and that '//&
                &'state''s seeds (0 = none; a rejection re-seeds the kept poles at posterior >= 0.5), Fisher z of the seeds'' own '//&
                &'deviation (reproducibility, distinctness), kept (share < 0.5 and both seed z >= '//int2str(nint(HALVES_KSIG))//')'
            do g = 1, np
                if( .not. lpole(g) ) cycle
                write(fu,'(I0,1X,I0,1X,F6.3,2(1X,F0.3),2(1X,I0),2(1X,F0.3),1X,A)') g, cnt0(g), shs(g), zrs(g), zcs(g), fin(g), &
                    &merge(count(seeds == fin(g)), 0, fin(g) > 0), zse(g), zde(g), trim(merge('kept    ', 'rejected', keep(g)))
            end do
            call fclose(fu)
        end subroutine pole_rule

        ! states below the seed floor fold into the remaining state with the highest mean posterior for their rows (a rest
        ! bucket, which has no endpoints, is the complement of the others); newidx maps the old state numbers (0 = folded)
        subroutine fold_small_states( newidx )
            integer, allocatable, intent(out) :: newidx(:)
            real(dp), allocatable :: gp(:,:)
            integer,  allocatable :: cnt(:), oldidx(:)
            logical,  allocatable :: alive(:)
            integer :: ns, s, t, f, i, m, k
            m  = size(names)
            ns = maxval([0, seeds])
            allocate(newidx(max(ns,1)), source=0)
            if( ns == 0 ) return
            allocate(alive(ns), source=.true.)
            do
                cnt = [(count(seeds == s), s=1,ns)]
                s   = 0
                do t = 1, ns
                    if( .not. alive(t) .or. cnt(t) == 0 .or. cnt(t) >= minseed ) cycle
                    if( s == 0 ) s = t
                    if( cnt(t) < cnt(s) ) s = t
                end do
                if( s == 0 .or. count(alive .and. cnt > 0) < 2 ) exit
                allocate(gp(n,ns), source=0._dp)
                do f = 1, m
                    if( grp(f) > 0 ) gp(:,grp(f)) = gp(:,grp(f)) + real(posts(:,f),dp) / real(count(grp == grp(f)),dp)
                end do
                do t = 1, ns
                    if( .not. any(grp == t) ) gp(:,t) = 1._dp - maxval(gp, dim=2)
                end do
                alive(s) = .false.
                do i = 1, n
                    if( seeds(i) /= s ) cycle
                    seeds(i) = maxloc(gp(i,:), dim=1, mask=alive .and. cnt > 0)
                end do
                where( grp == s ) grp = 0
                write(logfhandle,'(A)') '>>> SITE_RANK state '//int2str(s)//' folded: '//int2str(cnt(s))//' seeds, fewer than the '//&
                    &int2str(minseed)//' refine3D refines as a state; its rows join their best remaining state'
                deallocate(gp)
            end do
            k = 0
            do s = 1, ns
                if( alive(s) .and. count(seeds == s) > 0 )then
                    k = k + 1
                    newidx(s) = k
                endif
            end do
            oldidx = seeds
            do i = 1, n
                if( oldidx(i) > 0 ) seeds(i) = newidx(oldidx(i))
            end do
            do f = 1, m
                if( grp(f) > 0 ) grp(f) = newidx(grp(f))
            end do
        end subroutine fold_small_states

        subroutine size_rankings
            type(string), allocatable :: fn(:), ff(:)
            real,         allocatable :: v(:), sc(:)
            character(len=4096) :: ln
            integer :: fu, io, k0, j0, fu2
            if( .not. cline%defined('infile2') ) THROW_HARD('site_rank: site_size needs the rankings in infile2')
            allocate(fn(0), ff(0))
            call fopen(fu, params%infile2, status='old', action='read', iostat=io)
            call fileiochk('site_rank: cannot read '//params%infile2%to_char(), io)
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                ln = adjustl(ln)
                if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
                j0 = index(trim(ln), ' ')
                if( j0 == 0 ) THROW_HARD('site_rank: a ranking needs a name and a file: '//trim(ln))
                fn = [fn, string(ln(1:j0-1))]
                ff = [ff, string(trim(adjustl(ln(j0+1:))))]
            end do
            call fclose(fu)
            call fopen(fu2, string('site_size.txt'), status='replace', action='write', iostat=io)
            call fileiochk('site_rank: cannot write site_size.txt', io)
            write(fu2,'(A)') '# ranking, N* (hold-out correlation peak) over the fit and over the node, equivalent pure members '//&
                &'(hold-out amplitude x N / amplitude at the smallest N, its maximum) over the fit and over the node, hold-out peak'
            allocate(v(nrows))
            do k0 = 1, size(fn)
                call fopen(fu, ff(k0), status='old', action='read', iostat=io)
                call fileiochk('site_rank: cannot read '//ff(k0)%to_char(), io)
                read(fu,*) v
                call fclose(fu)
                sc  = max(v(rows), 0.)
                tag = fn(k0)%to_char()
                write(logfhandle,'(A)') '>>> SITE_RANK sizing the ranking '//tag
                call site%site_em(fit, fit, params%pcg_mskfile, post, nstar, curve, devs, sc, .false.)
                call write_top('site_top_'//tag//'.txt', nint(real(nstar) * real(n) / real(nfit)), sc)
                call write_top('site_epm_'//tag//'.txt', nint(maxval(curve(4,:)) * real(n) / real(nfit)), sc)
                call fopen(funit, string('site_curve_'//tag//'.txt'), status='replace', action='write', iostat=ios)
                call fileiochk('site_rank: cannot write site_curve_'//tag//'.txt', ios)
                write(funit,'(A)') '# N holdout_corr holdout_amplitude equivalent_pure_members'
                do i = 1, size(curve,2)
                    write(funit,'(I0,1X,F8.4,1X,ES12.5,1X,F10.1)') nint(curve(1,i)), curve(2,i), curve(3,i), curve(4,i)
                end do
                call fclose(funit)
                write(fu2,'(A,4(1X,I0),1X,F7.4)') tag, nstar, nint(real(nstar) * real(n) / real(nfit)), nint(maxval(curve(4,:))), &
                    &nint(maxval(curve(4,:)) * real(n) / real(nfit)), maxval(curve(2,:))
                write(logfhandle,'(A,A,A,I0,A,I0,A)') '>>> SITE_RANK ', tag, ': N* ', nint(real(nstar) * real(n) / real(nfit)), &
                    &', equivalent pure members ', nint(maxval(curve(4,:)) * real(n) / real(nfit)), ' over the node'
            end do
            call fclose(fu2)
        end subroutine size_rankings

        subroutine write_top( fname, nt, sc )
            character(len=*), intent(in) :: fname
            integer,          intent(in) :: nt
            real,             intent(in) :: sc(:)
            integer, allocatable :: tp(:), s2(:), o2(:)
            real,    allocatable :: k2(:)
            integer :: hh, ii
            allocate(tp(nrows), source=0)
            do hh = 0, 1
                s2 = pack([(ii, ii=1,n)], eo == hh)
                o2 = [(ii, ii=1,size(s2))]
                k2 = -sc(s2)
                call hpsort(k2, o2)
                tp(rows(s2(o2(1:min(max(nt,0)/2, size(s2)))))) = 1
            end do
            call write_ints(fname, tp)
        end subroutine write_top

        subroutine starts_of
            type(string), allocatable :: fn(:), ff(:)
            real,         allocatable :: v(:)
            character(len=4096) :: ln
            integer :: fu, io, k0, j0, nd2
            allocate(fn(0), ff(0))
            if( cline%defined('infile2') )then
                call fopen(fu, params%infile2, status='old', action='read', iostat=io)
                call fileiochk('site_rank: cannot read '//params%infile2%to_char(), io)
                do
                    read(fu,'(A)',iostat=io) ln
                    if( io /= 0 ) exit
                    ln = adjustl(ln)
                    if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
                    j0 = index(trim(ln), ' ')
                    if( j0 == 0 ) THROW_HARD('site_rank: a start needs a name and a seed file: '//trim(ln))
                    fn = [fn, string(ln(1:j0-1))]
                    ff = [ff, string(trim(adjustl(ln(j0+1:))))]
                end do
                call fclose(fu)
            endif
            nd2 = 2 * size(dirs)
            allocate(names(nd2 + size(fn)))
            allocate(starts(n, nd2 + size(fn)), source=.false.)
            sel = pack([(i, i=1,n)], fit)
            do idir = 1, size(dirs)
                ord = [(i, i=1,nfit)]
                key = z(rows(sel),idir)
                call hpsort(key, ord)
                do itail = 1, 2
                    k0  = 2 * (idir - 1) + itail
                    tag = TAILS(itail)
                    if( dirs(idir) > 0 ) tag = 'z'//int2str(dirs(idir))//tag
                    names(k0) = tag
                    if( itail == 1 )then
                        starts(sel(ord(1:ntail)),k0) = .true.
                    else
                        starts(sel(ord(nfit-ntail+1:nfit)),k0) = .true.
                    endif
                end do
            end do
            do k0 = 1, size(fn)
                allocate(v(nrows))
                call fopen(fu, ff(k0), status='old', action='read', iostat=io)
                call fileiochk('site_rank: cannot read '//ff(k0)%to_char(), io)
                do j0 = 1, nrows
                    read(fu,*) v(j0)
                end do
                call fclose(fu)
                starts(:,nd2+k0) = v(rows) > 0. .and. fit
                names(nd2+k0)    = fn(k0)
                deallocate(v)
            end do
        end subroutine starts_of

        subroutine read_latent
            real,    allocatable :: zq(:)
            character(len=8192)  :: line
            integer :: pind, ieo, ilab, nq, io
            real    :: res1, res2
            allocate(z(nrows,size(dirs)), source=0.)
            call fopen(funit, params%infile, status='old', action='read', iostat=ios)
            call fileiochk('site_rank: cannot read '//params%infile%to_char(), ios)
            nq = 0
            do
                read(funit,'(A)',iostat=io) line
                if( io /= 0 ) exit
                if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
                if( nq == 0 )then
                    nq = count_words(line) - 5
                    if( any(dirs < 1 .or. dirs > nq) ) THROW_HARD('site_rank: site_dirs beyond the '//int2str(nq)//' latent directions')
                    allocate(zq(nq))
                endif
                read(line,*) pind, ieo, ilab, res1, res2, zq
                if( pind >= 1 .and. pind <= nrows ) z(pind,:) = zq(dirs)
            end do
            call fclose(funit)
        end subroutine read_latent

        integer function count_words( line )
            character(len=*), intent(in) :: line
            integer :: k
            count_words = 0
            do k = 1, len_trim(line)
                if( line(k:k) /= ' ' .and. (k == 1 .or. line(max(1,k-1):max(1,k-1)) == ' ') ) count_words = count_words + 1
            end do
        end function count_words

        subroutine write_rows( fname, vals )
            character(len=*), intent(in) :: fname
            real,             intent(in) :: vals(:)
            real, allocatable :: full(:)
            integer :: fu, io, k
            allocate(full(nrows), source=-1.)
            full(rows) = vals
            call fopen(fu, string(fname), status='replace', action='write', iostat=io)
            call fileiochk('site_rank: cannot write '//fname, io)
            do k = 1, nrows
                write(fu,'(F8.5)') full(k)
            end do
            call fclose(fu)
        end subroutine write_rows

        subroutine write_ints( fname, vals )
            character(len=*), intent(in) :: fname
            integer,          intent(in) :: vals(:)
            integer :: fu, io, k
            call fopen(fu, string(fname), status='replace', action='write', iostat=io)
            call fileiochk('site_rank: cannot write '//fname, io)
            do k = 1, size(vals)
                write(fu,'(I0)') vals(k)
            end do
            call fclose(fu)
        end subroutine write_ints

    end subroutine exec_site_rank

    subroutine exec_landscape_compare( self, cline )
        use simple_landscape_compare, only: landscape_compare
        class(commander_landscape_compare), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    15.)
        call params%new(cline)
        if( cline%defined('infile2') )then
            call landscape_compare(params%infile, params%pcg_mskfile, params%lp, 'landscape_compare.txt', params%infile2)
        else
            call landscape_compare(params%infile, params%pcg_mskfile, params%lp, 'landscape_compare.txt')
        endif
        call simple_end('**** SIMPLE_LANDSCAPE_COMPARE NORMAL STOP ****')
    end subroutine exec_landscape_compare

    !> the add-only per-site split of a would-be leaf (labels and statistics only)
    subroutine exec_site_split( self, cline )
        use simple_site_split, only: run_site_split
        class(commander_site_split), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call run_site_split(cline)
        call simple_end('**** SIMPLE_SITE_SPLIT NORMAL STOP ****')
    end subroutine exec_site_split

    subroutine exec_site_discover( self, cline )
        use simple_hetero_tree,     only: discover_sites, discover_islands
        use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
        class(commander_site_discover), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)          :: params
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: labels(:), names(:)
        real(dp),     allocatable :: counts(:)
        integer :: s, i, funit, ios, nsite
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        allocate(labels(nlines(params%infile)))
        call fopen(funit, params%infile, status='old', action='read', iostat=ios)
        call fileiochk('site_discover: cannot read '//params%infile%to_char(), ios)
        do i = 1, size(labels)
            read(funit,*) labels(i)
        end do
        call fclose(funit)
        allocate(evens(0), odds(0), names(0), counts(0))
        do s = 1, maxval(labels)
            if( count(labels == s) == 0 ) cycle
            evens  = [evens, params%dir//'/'//refine3D_state_halfvol_fname(s, 'even', .true.)]
            odds   = [odds,  params%dir//'/'//refine3D_state_halfvol_fname(s, 'odd',  .true.)]
            names  = [names, s]
            counts = [counts, real(count(labels == s),dp)]
        end do
        if( size(names) < 2 ) THROW_HARD('site_discover needs at least two states')
        call discover_sites(evens, odds, counts, names, params%mskdiam, '.', nsite, trim(params%pgrp))
        write(logfhandle,'(A)') '>>> SITE_DISCOVER '//int2str(nsite)//' sites over '//int2str(size(names))//' states'
        if( params%site_split .eq. 'yes' ) call discover_islands(evens, odds, counts, params%mskdiam, '.')
        call simple_end('**** SIMPLE_SITE_DISCOVER NORMAL STOP ****')
    end subroutine exec_site_discover

    subroutine exec_landscape_sample( self, cline )
        use simple_landscape_sample, only: landscape_sample
        use simple_eo_partition,     only: ensure_field_halves
        class(commander_landscape_sample), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)       :: params
        type(builder)          :: build
        type(sp_project)       :: spproj
        type(landscape_sample) :: lsamp
        real    :: smpd
        integer :: box
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    10.)
        call cline%set('hp',      40.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('box_crop') )then
            call spproj%read_segment('stk', cline%get_carg('projfile'))
            box  = spproj%get_box()
            smpd = spproj%get_smpd()
            call spproj%kill
            call cline%set('box_crop', min(box, max(26, 2*ceiling(real(box)*smpd/cline%get_rarg('lp')))))
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call lsamp%new(params, build)
        if( cline%defined('infile2') )then
            call lsamp%run(params%infile, params%infile2, params%dir)
        else
            call lsamp%run(params%infile, string(''), params%dir)
        endif
        call lsamp%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_LANDSCAPE_SAMPLE NORMAL STOP ****')
    end subroutine exec_landscape_sample

    subroutine exec_hetero_tree( self, cline )
        use simple_hetero_tree, only: hetero_tree
        class(commander_hetero_tree), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters), target :: params
        type(hetero_tree) :: tree
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('kfloor')   ) call cline%set('kfloor',   5000)
        if( .not. cline%defined('minpop')   ) call cline%set('minpop',   15000)
        if( .not. cline%defined('maxdepth') ) call cline%set('maxdepth', 4)
        if( .not. cline%defined('nodepar')  ) call cline%set('nodepar',  2)
        if( .not. cline%defined('nparts')   ) call cline%set('nparts',   8)
        if( .not. cline%defined('nthr')     ) call cline%set('nthr',     8)
        if( .not. cline%defined('pgrp')     ) call cline%set('pgrp',     'c1')
        if( .not. cline%defined('site_dirs')) call cline%set('site_dirs', '1,3')
        if( .not. cline%defined('site_disc')) call cline%set('site_disc', 'yes')
        ! node-sized projects by default, with or without the site layer
        if( .not. cline%defined('node_proj') ) call cline%set('node_proj', 'yes')
        call cline%set('oritype', 'ptcl3D')
        call params%new(cline)
        call tree%new(params, cline%defined('infile'), cline)
        call tree%run
        call tree%kill
        call simple_end('**** SIMPLE_HETERO_TREE NORMAL STOP ****')
    end subroutine exec_hetero_tree

end module simple_commanders_refine3D
