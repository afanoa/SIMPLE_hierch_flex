!@descr: recursive heterogeneity tree: a resumable node queue whose mean, flex, focused refinement and split-test steps run as queue jobs, with cheap children, falloff-prior focus masks, node decisions and the leaves
module simple_hetero_tree
use simple_defs,            only: dp, logfhandle, PI, CWD_GLOB
use simple_defs_fname,      only: JOB_INFO_EXT
use simple_error,           only: simple_exception
use simple_string,          only: string
use simple_string_utils,    only: int2str
use simple_fileio,          only: fopen, fclose, fileiochk, nlines, read_exit_code, simple_copy_file
use simple_syslib,          only: file_exists, dir_exists, simple_mkdir, simple_chdir, simple_getcwd, simple_rename, simple_abspath, del_file, &
                                & exec_cmdline, simple_list_files_regexp, simple_sync_dir
use simple_imghead,         only: find_ldim_nptcls, find_img_smpd
use simple_rnd,             only: seed_rnd_fixed, ran3
use simple_linalg,          only: eigh, jacobi
use simple_srch_sort_loc,   only: hpsort
use simple_parameters,      only: parameters
use simple_cmdline,         only: cmdline
use simple_sp_project,      only: sp_project
use simple_oris,            only: oris
use simple_qsys_env,        only: qsys_env
use simple_image,           only: image
use simple_state_tree,      only: state_tree
use simple_tree_wiener,     only: tree_wiener_maps
use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
use simple_stat,            only: avg_ranks
use simple_state_mergeup,   only: mask_neff, split_zcrit
use simple_refine3D_strategy, only: PROB_MULTISTATE_MIN_POP
implicit none

public :: hetero_tree, discover_sites
private
#include "simple_local_flags.inc"

integer,  parameter :: KMAX = 8, NCHILD = 10, PERBIN = 1300, NBIN_MIN = 20, NKM_ITS = 40, POLL_SECS = 5, RND_SEED = 0
integer,  parameter :: DMINVOX = 200, DSEEDMIN = 1000, DFACT_ITS = 20
real,     parameter :: DMINSHARE = 0.05, DFACT_COS = 0.8, DFACT_SHARE = 0.1, DDICE = 0.9, DMAXMOL = 0.3, DMATE = 0.99
integer,  parameter :: ITA0 = 1, ITA = 3, ITB = 3, BOXA = 64, BOXB = 88, NSPA = 2500, NSPB = 5000, NPROBE = 2
integer,  parameter :: NPARTS_MIN = 4, PTCLS_PER_PART = 20000, SITE_NTHR = 32, SBATCH_PTCLS = 200000, PTCLS_PER_THR = 2500
! the merge-up job: a thread per 10k particles, up to a 48-core node (32 -> 48 threads cut a 692k merge-up by about 20 %)
integer,  parameter :: MERGEUP_NTHR = 48, MERGEUP_PTCLS = 10000
real,     parameter :: LPF = 15., LPA = 15., LPB = 10., MLP = 15., MFRAC = 0.2, MSIG = 1.5, MDIL = 3., MFALL = 6., MSUP = 0.1
real,     parameter :: ATHR = 0.01, LPD = 10.
integer,  parameter :: S_QUEUED = 0, S_MEAN = 1, S_FLEX = 2, S_CHEAP = 3, S_FFLEX = 4, S_STA0 = 5, S_RECB0 = 6, S_SITE = 7, &
                     & S_SEEDS = 8, S_STAGEA = 9, S_RECB = 10, S_STAGEB = 11, S_SPLIT = 12, S_DONE = 13, S_FAILED = -1
integer,  parameter :: S_SUBTREE = 14, SUB_NTHR = 48, S_MOTION = 15, S_MBATCH = 16, S_GATEW = 17, S_SBATCH = 18, S_OFFP = 19
integer,  parameter :: NMBIN = 5, NVSHUF = 20, GATE_BOX = 100, GATE_MINGRP = 50, INH0 = 100, NULL_PERM = 4
real,     parameter :: MRES = 15., NB_SW = 0.030, ZRES_SW = 5., ZLOC_MIN = 5., GRADED = 0.40, VIEW_EXC = 0.40, DFX_EXC = 0.20
real,     parameter :: GATE_LP = 15., GATE_EDGE = 6., GATE_NSD = 3.
character(len=*), parameter :: FLOOR_FILE = 'motion_floor.txt'
real(dp), parameter :: SLOW_FACT = 3._dp, SLOW_FIRST = 14400._dp
character(len=6), parameter :: STEP_DIRS(S_SPLIT) = [character(len=6) :: 'mean', 'flex', 'cheap', 'fflex', 'stA0', &
                                                     & 'recB0', 'site', 'seeds', 'stageA', 'recB', 'stageB', 'split']
character(len=*), parameter :: NODE_TABLE = 'hetero_tree_nodes.txt', LEAF_TABLE = 'hetero_tree_leaves.txt'
character(len=*), parameter :: STEP_EXIT = 'hetero_step.exit', NODE_PROJ = 'n.simple'
character(len=*), parameter :: NR_PROJ = 'node_rows.simple', NR_ROWS = 'node_rows.txt', NR_NAME = 'node_rows_proj.txt'
character(len=*), parameter :: POOL_PENDING = 'hetero_pool.pending'
character(len=*), parameter :: NODE_LOCAL = 'node_local.txt'
integer,  parameter :: SITE_NFIT = 100000, F_FAIL = -1, F_RUN = 1, F_TEST = 2, F_UP = 3, F_LEVEL = 4, F_COPY = 5, F_NOPOLE = 6, &
                     & F_DUP = 7, F_KEEP = 8, F_MAPS = 9, F_DONE = 10, F_REC = 11, F_REF = 12, F_READY = 13, F_PREC = 14, &
                     & F_PREF = 15, F_SIZE = 16
real,     parameter :: FGATE = 0.5
character(len=*), parameter :: SITE_TABLE = 'hetero_tree_sites', FDIR = 'factorial'
character(len=8), parameter :: FST_NAMES(F_FAIL:F_SIZE) = [character(len=8) :: 'failed', '-', 'ranking', 'ranked', 'pooling', &
                                                           & 'level', 'copy', 'no_pole', 'dup', 'kept', 'maps', 'done', 'refclass', &
                                                           & 'refrank', 'ready', 'prerec', 'preref', 'sizing']

type :: tnode
    character(len=:), allocatable :: id, pproj, prows
    integer  :: depth = 0, step = S_QUEUED, sseen = -1, nr = 0
    real(dp) :: tseen = -1._dp
    logical  :: lslow = .false.
end type tnode

type :: fsite
    character(len=:), allocatable :: name, home, level, target
    integer :: k = 0, st = F_RUN, ref = 0
    real    :: repro = 0., rho = 0., dice = 0.
end type fsite

type :: hetero_tree
    private
    class(parameters), pointer :: p => null()
    type(tnode), allocatable :: nodes(:)
    type(fsite), allocatable :: fs(:)
    type(string), allocatable :: subd(:)
    logical,      allocatable :: tjunk(:)
    type(qsys_env) :: qenv, qsite
    type(qsys_env), allocatable :: qpool
    character(len=:), allocatable :: root, sitemask
    real    :: smpd = 0.
    integer :: box = 0, nrows = 0, npcur = 0, ntests = 1
    logical :: l_subtree = .false., l_lsdone = .false., l_mudone = .false., l_fsread = .false., l_oneflex = .false.
    logical :: l_nodeproj = .false., l_nodeseed = .false., l_motion = .false., l_strong = .false., l_twopole = .false.
    logical :: l_gate = .false., l_inherit = .false., l_multi = .true.
    character(len=:), allocatable :: subroot, ntab, ltab
    type(qsys_env) :: qsub
    type(cmdline), allocatable :: cl0
    real(dp) :: srate(S_SBATCH) = 0._dp
  contains
    procedure :: new
    procedure :: run
    procedure :: kill
    procedure, private :: start
    procedure, private :: submit_subtree
    procedure, private :: merge_subtree
    procedure, private :: watch
    procedure, private :: prev_step
    procedure, private :: advance
    procedure, private :: submit
    procedure, private :: submit_in
    procedure, private :: submit_fflex
    procedure, private :: set_parts
    procedure, private :: final
    procedure, private :: plain
    procedure, private :: landscape
    procedure, private :: mergeup
    procedure, private :: step_status
    procedure, private :: pending
    procedure, private :: decide
    procedure, private :: finish
    procedure, private :: read_table
    procedure, private :: write_table
    procedure, private :: ndir
    procedure, private :: fs_sync
    procedure, private :: factorial
    procedure, private :: fs_rank
    procedure, private :: fs_climb
    procedure, private :: fs_maps
    procedure, private :: fs_output
    procedure, private :: fs_hd
    procedure, private :: fs_rd
    procedure, private :: fs_rho
    procedure, private :: fs_eval_rows
    procedure, private :: fs_labels
    procedure, private :: dir_status
    procedure, private :: fs_nrun
    procedure, private :: fs_cd
    procedure, private :: fs_md
    procedure, private :: fs_fd
    procedure, private :: fs_refine
    procedure, private :: fs_guard
    procedure, private :: fs_prerefine
    procedure, private :: fs_size
    procedure, private :: fs_sets
    procedure, private :: fs_mask
    procedure, private :: fs_zd
    procedure, private :: tree_junk
    procedure, private :: fs_read_table
    procedure, private :: fs_write_table
    procedure, private :: step_dir
    procedure, private :: root_axes
    procedure, private :: batch_status
    procedure, private :: off_work
    procedure, private :: anc_wait
    procedure, private :: floor_ready
    procedure, private :: floor_value
end type hetero_tree

contains

    subroutine new( self, params, l_rows, cline )
        class(hetero_tree),        intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        logical,                   intent(in)    :: l_rows
        class(cmdline), optional,  intent(in)    :: cline
        type(sp_project)     :: spproj
        type(string)         :: cwd, absproj
        integer, allocatable :: states(:), rows(:), sel(:)
        character(len=8)     :: zs
        character(len=:), allocatable :: how
        integer :: i, nroot
        call self%kill
        self%p         => params
        self%l_subtree = l_rows
        self%subroot   = ''
        if( params%subtree_node%strlen() > 0 ) self%subroot = params%subtree_node%to_char()
        self%l_nodeseed = params%node_seed .eq. 'yes' .or. params%subtree_max > 0 .or. len(self%subroot) > 0
        if( len(self%subroot) > 0 ) call simple_chdir(string('../..'))
        self%ntab = NODE_TABLE
        self%ltab = LEAF_TABLE
        if( len(self%subroot) > 0 )then
            self%ntab = 'node_'//self%subroot//'/subtree/subtree_nodes.txt'
            self%ltab = 'node_'//self%subroot//'/subtree/subtree_leaves.txt'
        else if( params%subtree_max > 0 )then
            if( params%site_layer .eq. 'yes' ) THROW_HARD('subtree_max > 0 needs site_layer=no')
            if( .not. present(cline) ) THROW_HARD('subtree jobs need the command line; new')
            allocate(self%cl0)
            self%cl0 = cline
        endif
        self%l_oneflex = params%oneflex .eq. 'yes'
        self%l_motion  = params%node_motion .eq. 'yes'
        self%l_strong  = params%strongest_pole .eq. 'yes' .and. .not. self%l_motion
        self%l_twopole = params%two_poles .eq. 'yes' .and. .not. self%l_motion
        self%l_gate    = self%l_motion .and. params%motion_gate .eq. 'yes'
        self%l_inherit = self%l_motion .and. params%motion_inherit .eq. 'yes'
        self%l_multi   = params%motion_multi .ne. 'no'
        self%l_nodeproj = params%node_proj .eq. 'yes' .and. params%site_layer .ne. 'yes'
        if( params%node_proj .eq. 'yes' .and. .not. self%l_nodeproj ) write(logfhandle,'(A)') '>>> HETERO_TREE node_proj=yes '//&
            &'needs site_layer=no: the node projects stay the size of the input project'
        call simple_getcwd(cwd)
        self%root     = cwd%to_char()
        self%sitemask = ''
        if( params%pcg_mskfile%strlen() > 0 )then
            absproj       = simple_abspath(params%pcg_mskfile)
            self%sitemask = absproj%to_char()
        endif
        call spproj%read_segment('stk', params%projfile)
        self%box  = spproj%get_box()
        self%smpd = spproj%get_smpd()
        call spproj%read_segment('ptcl3D', params%projfile)
        states     = spproj%os_ptcl3D%get_all_asint('state')
        self%nrows = size(states)
        call spproj%kill
        if( file_exists(string(self%ntab)) )then
            call self%read_table
            write(logfhandle,'(A)') '>>> HETERO_TREE resumed from '//self%ntab//': '//int2str(size(self%nodes))//' nodes'
        else if( len(self%subroot) > 0 )then
            call subtree_root_node
            call self%write_table
            write(logfhandle,'(A)') '>>> HETERO_TREE subtree job of node '//self%subroot//', depth '//&
                &int2str(self%nodes(1)%depth)
        else
            if( l_rows )then
                sel = read_ints(params%infile%to_char())
                if( size(sel) /= self%nrows ) THROW_HARD('infile must hold one value per particle of the project')
                rows = pack([(i, i=1,self%nrows)], sel > 0 .and. states > 0)
            else
                rows = pack([(i, i=1,self%nrows)], states > 0)
            endif
            if( size(rows) < 2 ) THROW_HARD('no particles for the root node')
            absproj = simple_abspath(params%projfile)
            allocate(self%nodes(1))
            self%nodes(1)%id    = 'r'
            self%nodes(1)%depth = 0
            self%nodes(1)%step  = S_QUEUED
            self%nodes(1)%pproj = absproj%to_char()
            self%nodes(1)%prows = '-'
            call simple_mkdir(string(self%ndir(1)))
            call write_ints(self%ndir(1)//'/rows.txt', rows)
            call self%write_table
            write(logfhandle,'(A)') '>>> HETERO_TREE root node: '//int2str(size(rows))//' particles of '//int2str(self%nrows)
        endif
        ! one Bonferroni count for the whole tree, from the root's particles (a subtree job runs in the main tree's root)
        if( .not. file_exists(string(self%root//'/node_r/rows.txt')) ) THROW_HARD('node_r/rows.txt is missing: the root''s &
            &particles set the split re-test''s Bonferroni count')
        nroot       = nlines(string(self%root//'/node_r/rows.txt'))
        self%ntests = max(1, (2 * nroot) / max(1, params%minpop) - 1)
        how         = '2 x '//int2str(nroot)//' root particles / minpop '//int2str(params%minpop)//' - 1'
        if( present(cline) )then
            ! a developer's split_ntests replaces the count; 1 is the uncorrected test
            if( cline%defined('split_ntests') )then
                self%ntests = max(1, params%split_ntests)
                how         = 'split_ntests on the command line'
            endif
        endif
        write(zs,'(F5.2)') split_zcrit(self%ntests)
        write(logfhandle,'(A)') '>>> HETERO_TREE split re-test: Bonferroni over '//int2str(self%ntests)//' split tests ('//how//&
            &'), z_crit '//trim(adjustl(zs))
        if( len(self%subroot) > 0 )then
            call self%qenv%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr, qsys_name=string('local'))
            call self%qsite%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=max(params%nthr, SITE_NTHR), &
                &qsys_name=string('local'))
        else
            call self%qenv%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr, &
                &base_only=params%node_worker_pool == 'yes')
            call self%qsite%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=max(params%nthr, SITE_NTHR), &
                &base_only=params%node_worker_pool == 'yes')
        endif

    contains

        subroutine subtree_root_node
            character(len=4096) :: line
            integer :: funit, ios, il
            call fopen(funit, string(NODE_TABLE), status='old', action='read', iostat=ios)
            call fileiochk('hetero_tree: cannot read '//NODE_TABLE, ios)
            do il = 1, nlines(string(NODE_TABLE))
                read(funit,'(A)') line
                if( word(line, 1) /= self%subroot ) cycle
                allocate(self%nodes(1))
                self%nodes(1)%id    = self%subroot
                self%nodes(1)%depth = str2int(word(line, 2))
                self%nodes(1)%step  = S_QUEUED
                self%nodes(1)%pproj = word(line, 4)
                self%nodes(1)%prows = word(line, 5)
                exit
            end do
            call fclose(funit)
            if( .not. allocated(self%nodes) ) THROW_HARD('node '//self%subroot//' is not in '//NODE_TABLE)
        end subroutine subtree_root_node

    end subroutine new

    subroutine run( self )
        class(hetero_tree), intent(inout) :: self
        integer(kind=8) :: t0, t1, rate
        integer :: i, nact, code, nsub, noff
        character(len=:), allocatable :: sdone
        logical :: l_fin, l_pl
        call seed_rnd_fixed(RND_SEED)
        do
            if( allocated(self%qpool) ) call self%qpool%service_persistent_worker_warmup()
            nact = 0
            nsub = 0
            i    = 0
            do while( i < size(self%nodes) )
                i = i + 1
                if( self%nodes(i)%step == S_SUBTREE )then
                    if( .not. self%merge_subtree(i) ) nsub = nsub + 1
                    cycle
                endif
                if( .not. live_step(self%nodes(i)%step) ) cycle
                code = self%step_status(i, self%nodes(i)%step)
                call self%watch(i, code)
                if( code == -1 )then
                    ! a node waiting for the gate floor holds no slot: the floor may wait on a node not yet started
                    if( self%nodes(i)%step /= S_GATEW ) nact = nact + 1
                    cycle
                endif
                if( code == -2 )then
                    call self%start(i, self%nodes(i)%step)
                    nact = nact + 1
                    call self%write_table
                    cycle
                endif
                if( code /= 0 )then
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': step '//&
                        &step_label(self%nodes(i)%step)//' FAILED (exit code '//int2str(code)//'), node stopped'
                    self%nodes(i)%step = S_FAILED
                    call self%write_table
                    cycle
                endif
                call system_clock(t0, rate)
                call expand_proj(self%step_dir(i, self%nodes(i)%step))
                sdone = self%step_dir(i, self%nodes(i)%step)
                call self%advance(i)
                if( .not. keep_maps() )then
                    call prune_step(sdone)
                    if( self%nodes(i)%step == S_DONE ) call prune_node(self, i)
                endif
                call system_clock(t1)
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': controller '//&
                    &rstr(real(t1 - t0) / real(max(rate, 1_8)))//' s after its step'
                if( live_step(self%nodes(i)%step) ) nact = nact + 1
                call self%write_table
            end do
            do i = 1, size(self%nodes)
                if( self%nodes(i)%step /= S_QUEUED ) cycle
                if( self%p%subtree_max > 0 .and. len(self%subroot) == 0 .and. self%nodes(i)%depth > 0 )then
                    if( nlines(string(self%ndir(i)//'/rows.txt')) <= self%p%subtree_max )then
                        call self%submit_subtree(i)
                        nsub = nsub + 1
                        call self%write_table
                        cycle
                    endif
                endif
                if( nact >= max(1, self%p%nodepar) ) cycle
                call system_clock(t0, rate)
                call self%start(i, S_MEAN)
                call system_clock(t1)
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': controller '//&
                    &rstr(real(t1 - t0) / real(max(rate, 1_8)))//' s to start'
                nact = nact + 1
                call self%write_table
            end do
            if( self%p%site_layer .eq. 'yes' )then
                call self%fs_sync
            else if( len(self%subroot) == 0 )then
                l_pl = self%plain()
            endif
            call self%off_work(noff)
            if( nact == 0 .and. nsub == 0 .and. .not. any(self%nodes(:)%step == S_QUEUED .or. self%nodes(:)%step == S_GATEW) )then
                l_fin = self%final()
                if( l_fin .and. noff == 0 ) exit
            endif
            call sleep(POLL_SECS)
        end do
        if( len(self%subroot) == 0 .and. .not. keep_maps() ) call prune_maps(self)
    end subroutine run

    subroutine start( self, i, istep )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, istep
        type(cmdline)        :: cl
        integer, allocatable :: rows(:)
        character(len=:), allocatable :: d
        integer :: k
        if( istep > S_MEAN )then
            k = self%prev_step(i, istep)
            if( k >= S_MEAN )then
                call self%advance(i, k)
                return
            endif
        endif
        d = self%ndir(i)
        if( .not. file_exists(string(d//'/'//NODE_PROJ)) )then
            call node_project(self, i)
            if( len(self%subroot) > 0 ) call stamp_local(d//'/'//NODE_PROJ)
        endif
        rows = read_ints(d//'/rows.txt')
        k    = min(KMAX, max(2, size(rows) / max(1, self%p%kfloor)))
        call self%set_parts(size(rows))
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' (depth '//int2str(self%nodes(i)%depth)//'): '//&
            &int2str(size(rows))//' particles, '//int2str(k)//' flex states, '//int2str(self%npcur)//' parts per step'
        if( self%l_oneflex .and. self%nodes(i)%depth > 0 )then
            call oneflex_inputs(self, i)
            call self%submit_fflex(i, d//'/pmean.mrc')
            return
        endif
        call simple_mkdir(string(d//'/mean'))
        call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/mean/m.simple'))
        call cl%set('prg',      'reconstruct3D')
        call cl%set('projfile', 'm.simple')
        call cl%set('nstates',  1)
        call common_keys(self, cl)
        call cl%set('postprocess', 'no')
        call self%submit(i, S_MEAN, cl)
    end subroutine start

    subroutine submit_subtree( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        type(cmdline) :: cl
        type(string)  :: absproj
        character(len=:), allocatable :: sd
        integer :: nr, np, nc
        sd = self%ndir(i)//'/subtree'
        call simple_mkdir(string(sd))
        if( .not. self%qsub%exists() ) call self%qsub%new(self%p, 1, exec_bin=string('simple_exec'), qsys_nthr=SUB_NTHR, &
            &base_only=self%p%node_worker_pool == 'yes')
        nr = nlines(string(self%ndir(i)//'/rows.txt'))
        np = min(self%p%nparts, max(NPARTS_MIN, ceiling(real(nr) / real(PTCLS_PER_PART))))
        nc = min(SUB_NTHR, max(self%p%nthr, max(1, self%p%nodepar) * np * self%p%nthr))
        cl = self%cl0
        if( self%p%node_worker_pool == 'yes' )then
            call cl%set('node_worker_pool', 'no')
            call cl%set('qsys_name', self%qenv%get_qsys())
            call cl%delete('worker_server')
            call cl%delete('worker_server_nthr')
        endif
        absproj = simple_abspath(self%p%projfile)
        call cl%set('projfile',     absproj%to_char())
        call cl%set('subtree_node', self%nodes(i)%id)
        call cl%set('mkdir',        'no')
        call cl%delete('infile')
        if( len(self%sitemask) > 0 ) call cl%set('pcg_mskfile', self%sitemask)
        call simple_chdir(string(sd))
        CWD_GLOB = sd
        call self%qsub%qdescr%set('job_memory_per_task', int2str(max(64000, 8000 + nint(0.1 * real(nr)))))
        call self%qsub%qdescr%set('job_cpus_per_task',   int2str(nc))
        call self%qsub%exec_simple_prg_in_queue_async(cl, string('./subtree.sh'), string('subtree.log'), &
            &exit_code_fname=string(STEP_EXIT))
        call simple_chdir(string(self%root))
        CWD_GLOB = self%root
        if( .not. allocated(self%subd) ) allocate(self%subd(0))
        self%subd = [self%subd, string(sd)]
        self%nodes(i)%step = S_SUBTREE
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' (depth '//int2str(self%nodes(i)%depth)//', '//&
            &int2str(nr)//' particles): subtree job submitted ('//int2str(nc)//' cores)'
    end subroutine submit_subtree

    ! a finished subtree job's nodes and leaves join the tree (merged.txt: its leaves were appended)
    logical function merge_subtree( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        type(tnode), allocatable :: tmp(:)
        character(len=4096) :: line
        character(len=:), allocatable :: sd, id
        integer :: code, funit, fout, ios, il, nleaf, nnew, jj
        merge_subtree = .false.
        sd   = self%ndir(i)//'/subtree'
        code = self%dir_status(sd)
        if( code == -1 ) return
        if( code == -2 )then
            call self%submit_subtree(i)
            return
        endif
        merge_subtree = .true.
        if( code /= 0 )then
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': subtree job FAILED (exit code '//&
                &int2str(code)//'), node stopped'
            self%nodes(i)%step = S_FAILED
            call self%write_table
            return
        endif
        nleaf = 0
        if( .not. file_exists(string(sd//'/merged.txt')) )then
            if( file_exists(string(sd//'/subtree_leaves.txt')) )then
                call fopen(funit, string(sd//'/subtree_leaves.txt'), status='old', action='read', iostat=ios)
                call fileiochk('hetero_tree: cannot read '//sd//'/subtree_leaves.txt', ios)
                call fopen(fout, string(self%root//'/'//LEAF_TABLE), status='unknown', position='append', action='write', &
                    &iostat=ios)
                call fileiochk('hetero_tree: cannot append to '//LEAF_TABLE, ios)
                do il = 1, nlines(string(sd//'/subtree_leaves.txt'))
                    read(funit,'(A)') line
                    write(fout,'(A)') trim(line)
                    nleaf = nleaf + 1
                end do
                call fclose(fout)
                call fclose(funit)
            endif
            call write_ints(sd//'/merged.txt', [nleaf])
        endif
        nnew = 0
        call fopen(funit, string(sd//'/subtree_nodes.txt'), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//sd//'/subtree_nodes.txt', ios)
        do il = 1, nlines(string(sd//'/subtree_nodes.txt'))
            read(funit,'(A)') line
            id = word(line, 1)
            if( id == self%nodes(i)%id )then
                self%nodes(i)%step = str2int(word(line, 3))
                cycle
            endif
            if( any([(self%nodes(jj)%id == id, jj=1,size(self%nodes))]) ) cycle
            allocate(tmp(size(self%nodes) + 1))
            tmp(1:size(self%nodes)) = self%nodes
            tmp(size(tmp))%id    = id
            tmp(size(tmp))%depth = str2int(word(line, 2))
            tmp(size(tmp))%step  = str2int(word(line, 3))
            tmp(size(tmp))%pproj = word(line, 4)
            tmp(size(tmp))%prows = word(line, 5)
            call move_alloc(tmp, self%nodes)
            nnew = nnew + 1
        end do
        call fclose(funit)
        call self%write_table
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': subtree job done, '//int2str(nnew)//&
            &' descendant nodes and its leaves joined the tree'
    end function merge_subtree

    ! report-only: a step running over SLOW_FACT times this run's slowest rate of its kind (SLOW_FIRST s without one)
    subroutine watch( self, i, code )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, code
        integer(kind=8) :: cnt, rate
        character(len=4096) :: rec(3)
        character(len=:), allocatable :: sd
        real(dp) :: now, el, lim
        integer  :: s, funit, ios, j
        s = self%nodes(i)%step
        if( self%nodes(i)%sseen /= s )then
            self%nodes(i)%sseen = s
            self%nodes(i)%tseen = -1._dp
            self%nodes(i)%lslow = .false.
            self%nodes(i)%nr    = max(1, nlines(string(self%ndir(i)//'/rows.txt')))
        endif
        sd = self%step_dir(i, s)
        call system_clock(cnt, rate)
        now = real(cnt,dp) / real(max(rate, 1_8),dp)
        if( self%nodes(i)%tseen < 0._dp )then
            if( code /= -1 .or. .not. file_exists(string(sd//'/'//STEP_EXIT//JOB_INFO_EXT)) ) return
            self%nodes(i)%tseen = now
        endif
        el = now - self%nodes(i)%tseen
        if( code == 0 ) self%srate(s) = max(self%srate(s), el / real(self%nodes(i)%nr,dp))
        if( code /= -1 .or. self%nodes(i)%lslow ) return
        lim = SLOW_FIRST
        if( self%srate(s) > 0._dp ) lim = SLOW_FACT * self%srate(s) * real(self%nodes(i)%nr,dp)
        if( el <= lim ) return
        self%nodes(i)%lslow = .true.
        rec = '?'
        call fopen(funit, string(sd//'/'//STEP_EXIT//JOB_INFO_EXT), status='old', action='read', iostat=ios)
        if( ios == 0 )then
            do j = 1, 3
                read(funit,'(A)',iostat=ios) rec(j)
                if( ios /= 0 ) exit
            end do
            call fclose(funit)
        endif
        call fopen(funit, string(sd//'/hetero_step.slow'), status='replace', action='write', iostat=ios)
        if( ios == 0 )then
            write(funit,'(A)') 'running '//int2str(nint(el))//' s, limit '//int2str(nint(lim))//' s, job '//trim(rec(3))//&
                &' on '//trim(rec(2))
            call fclose(funit)
        endif
        write(logfhandle,'(A)') '>>> HETERO_TREE WARNING node '//self%nodes(i)%id//': step '//step_label(s)//' running '//&
            &int2str(nint(el))//' s, over its limit '//int2str(nint(lim))//' s ('//trim(rec(3))//' on '//trim(rec(2))//')'
    end subroutine watch

    subroutine advance( self, i, done_step )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        integer, optional,  intent(in)    :: done_step
        type(cmdline)             :: cl
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: rows(:), st(:), pres(:), ch(:), seeds(:), jk(:)
        real(dp),     allocatable :: w(:)
        character(len=:), allocatable :: d, mask, sitemask, sd, tmsg
        type(string), allocatable :: bdirs(:)
        integer,      allocatable :: hrow(:)
        real    :: fl
        integer :: s, k, j, ks
        logical :: l_tp, l_fl, l_wait
        s = self%nodes(i)%step
        if( present(done_step) ) s = done_step
        d        = self%ndir(i)
        mask     = d//'/mask.mrc'
        sitemask = self%sitemask
        if( len(sitemask) == 0 ) sitemask = mask
        rows = read_ints(d//'/rows.txt')
        k    = min(KMAX, max(2, size(rows) / max(1, self%p%kfloor)))
        call self%set_parts(size(rows))
        select case(s)
            case(S_MEAN)
                call simple_mkdir(string(d//'/flex'))
                call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/flex/f.simple'))
                call flex_keys(self, cl, d//'/mean/recvol_state01.mrc', k)
                if( local_rows(d) ) call cl%set('rowmap', d//'/rows.txt')
                call self%submit(i, S_FLEX, cl)
            case(S_FLEX)
                if( self%nodes(i)%depth == 0 .and. .not. self%l_subtree )then
                    st   = project_states(d//'/flex/f.simple')
                    pres = present_labels(st)
                    allocate(evens(size(pres)), odds(size(pres)), w(size(pres)))
                    do j = 1, size(pres)
                        evens(j) = d//'/flex/'//hv(pres(j), 'even')
                        odds(j)  = d//'/flex/'//hv(pres(j), 'odd')
                        w(j)     = real(count(st == pres(j)), dp)
                    end do
                    call falloff_mask(self, evens, odds, w, mask)
                    call self%submit_fflex(i, d//'/mean/recvol_state01.mrc')
                else
                    call cheap_children(self, i, ch)
                    call simple_mkdir(string(d//'/cheap'))
                    call write_states(d//'/'//NODE_PROJ, d//'/cheap/c.simple', ch)
                    call cl%set('prg',      'reconstruct3D')
                    call cl%set('projfile', 'c.simple')
                    call cl%set('nstates',  NCHILD)
                    call common_keys(self, cl)
                    call halves_only_keys(cl)
                    call self%submit(i, S_CHEAP, cl)
                endif
            case(S_CHEAP)
                st = project_states(d//'/cheap/c.simple')
                allocate(evens(0), odds(0), w(0))
                do j = 1, NCHILD
                    if( count(st == j) == 0 ) cycle
                    if( .not. file_exists(string(d//'/cheap/'//hv(j, 'even'))) ) cycle
                    evens = [evens, string(d//'/cheap/'//hv(j, 'even'))]
                    odds  = [odds,  string(d//'/cheap/'//hv(j, 'odd'))]
                    w     = [w, real(count(st == j), dp)]
                end do
                if( size(evens) < 4 ) THROW_HARD('fewer than four cheap children maps for the node mask')
                call falloff_mask(self, evens, odds, w, mask)
                call self%submit_fflex(i, d//'/mean/recvol_state01.mrc')
            case(S_FFLEX)
                if( self%p%site_disc .ne. 'yes' )then
                    call site_step(sitemask)
                    return
                endif
                st   = project_states(d//'/fflex/f.simple')
                pres = present_labels(st)
                ! refine3D's multi-state startup needs more than PROB_MULTISTATE_MIN_POP particles per state; the rows stay in the node
                do j = 1, size(pres)
                    if( count(st == pres(j)) > PROB_MULTISTATE_MIN_POP .or. count(st > 0 .and. st /= pres(j)) == 0 ) cycle
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': flex state '//int2str(pres(j))//' ('//&
                        &int2str(count(st == pres(j)))//' particles) is under refine3D''s minimum, unassigned for stA0'
                    where( st == pres(j) ) st = 0
                end do
                pres = present_labels(st)
                ks   = size(pres)
                allocate(ch(size(st)), source=0)
                do j = 1, ks
                    where( st == pres(j) ) ch = j
                end do
                call simple_mkdir(string(d//'/stA0'))
                call write_states(d//'/'//NODE_PROJ, d//'/stA0/r.simple', ch)
                call cl%set('prg',      'refine3D')
                call cl%set('projfile', 'r.simple')
                call cl%set('nstates',  ks)
                do j = 1, ks
                    call cl%set('vol'//int2str(j), d//'/fflex/recvol_state'//pad2(pres(j))//'.mrc')
                end do
                call refine_keys(self, cl, LPA, ITA0, NSPA, BOXA, mask)
                if( self%p%sta0_npeaks > 0 ) call cl%set('npeaks_inpl', self%p%sta0_npeaks)
                call self%submit(i, S_STA0, cl)
            case(S_STA0)
                ks = maxval(project_states(d//'/stA0/r.simple'))
                call simple_mkdir(string(d//'/recB0'))
                call simple_copy_file(string(d//'/stA0/r.simple'), string(d//'/recB0/b.simple'))
                call cl%set('prg',      'reconstruct3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                call common_keys(self, cl)
                call halves_only_keys(cl)
                call self%submit(i, S_RECB0, cl)
            case(S_RECB0)
                st   = project_states(d//'/recB0/b.simple')
                pres = with_maps(d//'/recB0', present_labels(st))
                allocate(evens(size(pres)), odds(size(pres)), w(size(pres)))
                do j = 1, size(pres)
                    evens(j) = d//'/recB0/'//hv(pres(j), 'even')
                    odds(j)  = d//'/recB0/'//hv(pres(j), 'odd')
                    w(j)     = real(count(st == pres(j)), dp)
                end do
                call simple_mkdir(string(d//'/sites'))
                call discover_sites(evens, odds, w, pres, self%p%mskdiam, d//'/sites', ks, trim(self%p%pgrp))
                if( ks > 1 .and. (self%l_motion .or. self%l_strong) )then
                    call site_batch(ks)
                    return
                endif
                if( ks > 1 .and. size(rows) > SBATCH_PTCLS )then
                    call site_batch(ks)
                    call write_ints(d//'/sites/rest_first.txt', [1])
                    return
                endif
                if( ks > 0 )then
                    ! a site step restarted after losing its record ranks the site its directory belongs to, not site 1
                    j = seed_site(d)
                    if( j > ks ) THROW_HARD('node '//self%nodes(i)%id//': seed site '//int2str(j)//' beyond the '//int2str(ks)//' sites found again')
                    sitemask = d//'/sites/site_'//pad2(j)//'.mrc'
                else
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': no local site found, the node mask ranks the poles'
                    sitemask = mask
                endif
                call site_step(sitemask)
            case(S_SITE, S_SBATCH)
                sd    = seed_dir(d)
                seeds = read_ints(sd//'/site_seeds.txt')
                ks    = maxval([0, seeds])
                if( s == S_SBATCH )then
                    do j = 1, nsites(d)
                        call log_site(site_dir_of(d, j), 'site '//int2str(j)//' ')
                    end do
                else
                    call log_site(sd, '')
                endif
                j = str2int(trim(adjustl(first_line(d//'/sites/seed_site.txt', '1'))))
                if( s == S_SBATCH .and. .not. (self%l_motion .or. self%l_strong) )then
                    ! the remaining sites were ranked at once: the first that gives a split is the one the chain stops at
                    j = str2int(trim(adjustl(first_line(d//'/sites/rest_first.txt', '2'))))
                    do while( j < nsites(d) )
                        if( maxval([0, read_ints(site_dir_of(d, j)//'/site_seeds.txt')]) >= 2 ) exit
                        j = j + 1
                    end do
                    call write_ints(d//'/sites/seed_site.txt', [j])
                    sd    = seed_dir(d)
                    seeds = read_ints(sd//'/site_seeds.txt')
                    ks    = maxval([0, seeds])
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': site '//int2str(j)//' of the batch '//&
                        &merge('gives the split   ', 'is the last site  ', ks >= 2)
                endif
                if( (self%l_motion .or. self%l_strong) .and. j < nsites(d) )then
                    if( ks >= 2 )then
                        call next_site(j, 'keeps a pole; every site is ranked')
                    else
                        call next_site(j, 'gives no split')
                    endif
                    return
                endif
                if( self%l_motion )then
                    jk = kept_sites(d)
                    do j = 1, size(jk)
                        call add_cand(d, jk(j))
                    end do
                    call motion_next
                    return
                endif
                if( self%l_strong )then
                    jk = kept_sites(d)
                    ks = 0
                    if( size(jk) > 0 )then
                        j     = strongest_site(d, jk)
                        sd    = site_dir_of(d, j)
                        seeds = read_ints(sd//'/site_seeds.txt')
                        ks    = maxval(seeds)
                        call write_ints(d//'/sites/chosen_site.txt', [j])
                        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': site '//int2str(j)//' has the '//&
                            &'strongest kept pole (smaller pole-rule margin per kept site: '//kept_margins(d, jk)//'), its seeds '//&
                            &'seed the refinement'
                    endif
                endif
                l_tp = .false.
                if( self%l_twopole .and. ks >= 2 )then
                    if( two_pole_seeds(d, j, self%p%minpop, ch, jk, tmsg) )then
                        seeds = ch
                        ks    = 2
                        l_tp  = .true.
                        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': two-pole seeds from site '//&
                            &int2str(j)//': '//tmsg
                    endif
                endif
                if( ks < 2 )then
                    if( .not. self%l_strong .and. j < nsites(d) )then
                        call rest_batch(j)
                        return
                    endif
                    ch = two_groups(d)
                    jk = junk_rows(d)
                    if( size(jk) > 0 .and. size(ch) > 0 ) ch(jk) = 0
                    if( count(ch == 1) > PROB_MULTISTATE_MIN_POP .and. count(ch == 2) > PROB_MULTISTATE_MIN_POP )then
                        seeds = ch
                        ks    = 2
                        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': no site gives a split, the two '//&
                            &'groups of its fflex states seed the refinement'
                    endif
                endif
                call seed_node(seeds, ks, sd//'/site_contrast.txt')
                if( l_tp ) call write_ints(d//'/seeds/start.txt', jk)
            case(S_MOTION, S_MBATCH, S_GATEW)
                call motion_next
            case(S_OFFP)
                call off_next
            case(S_SEEDS)
                seeds = read_ints(d//'/seeds/seeds.txt')
                ks    = maxval(seeds)
                if( file_exists(string(d//'/seeds/start.txt')) ) seeds = read_ints(d//'/seeds/start.txt')
                call simple_mkdir(string(d//'/stageA'))
                call write_states(d//'/'//NODE_PROJ, d//'/stageA/r.simple', seeds)
                call cl%set('prg',      'refine3D')
                call cl%set('projfile', 'r.simple')
                call cl%set('nstates',  ks)
                do j = 1, ks
                    call cl%set('vol'//int2str(j), d//'/seeds/recvol_state'//pad2(j)//'.mrc')
                end do
                call refine_keys(self, cl, LPA, ITA, NSPA, BOXA, mask)
                call self%submit(i, S_STAGEA, cl)
            case(S_STAGEA)
                st = project_states(d//'/stageA/r.simple')
                ch = compact(fold_small(st))
                if( maxval(ch) < maxval(compact(st)) ) write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//&
                    &': stage-A states under refine3D''s minimum fold into the largest, '//int2str(maxval(ch))//' remain'
                ks = maxval(ch)
                call simple_mkdir(string(d//'/recB'))
                call write_states(d//'/stageA/r.simple', d//'/recB/b.simple', ch)
                call cl%set('prg',      'reconstruct3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                call common_keys(self, cl)
                call cl%set('postprocess', 'no')
                call crop_keys(self, cl, BOXB)
                call self%submit(i, S_RECB, cl)
            case(S_RECB)
                ch = compact(fold_small(project_states(d//'/stageA/r.simple')))
                ks = maxval(ch)
                call simple_mkdir(string(d//'/stageB'))
                call write_states(d//'/stageA/r.simple', d//'/stageB/b.simple', ch)
                call cl%set('prg',      'refine3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                do j = 1, ks
                    call cl%set('vol'//int2str(j), d//'/recB/recvol_state'//pad2(j)//'.mrc')
                end do
                call refine_keys(self, cl, LPB, ITB, NSPB, BOXB, mask)
                call self%submit(i, S_STAGEB, cl)
            case(S_STAGEB)
                st   = project_states(d//'/stageB/b.simple')
                pres = present_labels(st)
                if( .not. all([(file_exists(string(d//'/stageB/'//hv(pres(j), 'even'))) .and. &
                    &file_exists(string(d//'/stageB/'//hv(pres(j), 'odd'))) .and. &
                    &file_exists(string(d//'/stageB/recvol_state'//pad2(pres(j))//'.mrc')), j=1,size(pres))]) )then
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': stage B maps missing, node stopped'
                    self%nodes(i)%step = S_FAILED
                    return
                endif
                call simple_mkdir(string(d//'/split'))
                if( maxval([0, st]) < 2 )then
                    call write_ints(d//'/split/split_children.txt', merge(1, 0, st > 0))
                    call write_ints(d//'/split/split_junk.txt', [integer ::])
                    call self%decide(i)
                    self%nodes(i)%step = S_DONE
                    return
                endif
                call simple_copy_file(string(d//'/stageB/b.simple'), string(d//'/split/s.simple'))
                call cl%set('prg',         'state_split')
                call cl%set('projfile',    's.simple')
                call cl%set('mskdiam',     self%p%mskdiam)
                call cl%set('pcg_mskfile', mask)
                call cl%set('pgrp',        trim(self%p%pgrp))
                call cl%set('split_ntests', self%ntests)
                call cl%set('junk_screen', trim(self%p%junk_screen))
                if( self%p%merge_geom_reuse == 'yes' ) call cl%set('merge_geom_reuse', 'yes')
                if( self%p%merge_insert_reuse == 'yes' ) call cl%set('merge_insert_reuse', 'yes')
                call cl%set('nthr',        max(self%p%nthr, min(SITE_NTHR, ceiling(real(count(st > 0)) / real(PTCLS_PER_PART)))))
                call self%submit(i, S_SPLIT, cl, count(st > 0))
            case(S_SPLIT)
                call self%decide(i)
                self%nodes(i)%step = S_DONE
        end select

    contains

        subroutine site_step( smask )
            character(len=*), intent(in) :: smask
            call site_submit(smask, seed_dir(d))
            self%nodes(i)%step = S_SITE
        end subroutine site_step

        subroutine site_submit( smask, sdir )
            character(len=*), intent(in) :: smask, sdir
            type(cmdline) :: cj
            call simple_mkdir(string(sdir))
            call simple_copy_file(string(d//'/'//NODE_PROJ), string(sdir//'/s.simple'))
            call cj%set('prg',         'site_rank')
            call cj%set('projfile',    's.simple')
            call cj%set('infile',      d//'/fflex/flex_pca_coordinates.txt')
            call cj%set('site_dirs',   trim(self%p%site_dirs))
            call cj%set('pole_merge',  self%p%pole_merge)
            if( self%p%pooled_poles .eq. 'yes' ) call cj%set('pooled_poles', 'yes')
            call cj%set('junk_screen', trim(self%p%junk_screen))
            call cj%set('pcg_mskfile', smask)
            call cj%set('pgrp',        trim(self%p%pgrp))
            call cj%set('mskdiam',     self%p%mskdiam)
            ! a site rank gets a thread per PTCLS_PER_THR particles (the root keeps SITE_NTHR)
            call cj%set('nthr',        max(self%p%nthr, min(SITE_NTHR, ceiling(real(size(rows)) / real(PTCLS_PER_THR)))))
            if( local_rows(d) ) call cj%set('rowmap', d//'/rows.txt')
            call self%submit_in(sdir, cj, size(rows))
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': '//sdir(index(sdir, '/', back=.true.)+1:)//&
                &' submitted'
        end subroutine site_submit

        ! every site of the node ranked at once, one job each (the rankings are independent); the batch ends as the chain did
        subroutine site_batch( kss )
            integer, intent(in) :: kss
            type(string), allocatable :: sds(:)
            character(len=:), allocatable :: sdir
            integer :: kk, nsub
            allocate(sds(0))
            nsub = 0
            do kk = 1, kss
                if( kk > 1 ) call simple_mkdir(string(d//'/nsite_'//self%nodes(i)%id//'_s'//int2str(kk)))
                sdir = site_dir_of(d, kk)
                sds  = [sds, string(sdir)]
                if( self%dir_status(sdir) /= -2 ) cycle
                call site_submit(d//'/sites/site_'//pad2(kk)//'.mrc', sdir)
                nsub = nsub + 1
            end do
            call write_ints(d//'/sites/seed_site.txt', [kss])
            call write_lines(d//'/sites/batch.txt', sds)
            self%nodes(i)%step = S_SBATCH
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': its '//int2str(kss)//' sites are ranked at '//&
                &'once, '//int2str(nsub)//' jobs submitted'
        end subroutine site_batch

        subroutine rest_batch( jf )
            integer, intent(in) :: jf
            type(string), allocatable :: sds(:)
            character(len=:), allocatable :: sdir
            integer :: kk, nsub, ns
            ns = nsites(d)
            allocate(sds(0))
            nsub = 0
            do kk = jf + 1, ns
                call simple_mkdir(string(d//'/nsite_'//self%nodes(i)%id//'_s'//int2str(kk)))
                sdir = site_dir_of(d, kk)
                sds  = [sds, string(sdir)]
                if( self%dir_status(sdir) /= -2 ) cycle
                call site_submit(d//'/sites/site_'//pad2(kk)//'.mrc', sdir)
                nsub = nsub + 1
            end do
            call write_ints(d//'/sites/rest_first.txt', [jf + 1])
            call write_lines(d//'/sites/batch.txt', sds)
            call write_ints(d//'/sites/seed_site.txt', [ns])
            self%nodes(i)%step = S_SBATCH
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': site '//int2str(jf)//' gives no split, sites '//&
                &int2str(jf + 1)//'..'//int2str(ns)//' are ranked at once, '//int2str(nsub)//' jobs submitted'
        end subroutine rest_batch

        subroutine log_site( sdl, pre )
            character(len=*), intent(in) :: sdl, pre
            call log_lines(sdl//'/site_poles.txt', pre//'pole (index, seeds, endpoints)')
            call log_lines(sdl//'/site_halves.txt', pre//'halves')
            call log_lines(sdl//'/site_pole_rule.txt', pre//'pole rule')
            call log_lines(sdl//'/site_spearman.txt', pre//'endpoint Spearman (merged at >= '//&
                &trim(adjustl(rfmt2(self%p%pole_merge)))//')')
            call log_lines(sdl//'/site_junk_report.txt', pre//'junk screen')
        end subroutine log_site

        subroutine next_site( j0, what )
            integer,          intent(in) :: j0
            character(len=*), intent(in) :: what
            call simple_mkdir(string(d//'/nsite_'//self%nodes(i)%id//'_s'//int2str(j0 + 1)))
            call write_ints(d//'/sites/seed_site.txt', [j0 + 1])
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': site '//int2str(j0)//' '//what//', '//&
                &'ranking site '//int2str(j0 + 1)
            call site_step(d//'/sites/site_'//pad2(j0 + 1)//'.mrc')
        end subroutine next_site

        subroutine seed_node( sds, kss, cfile )
            integer, allocatable, intent(inout) :: sds(:)
            integer,              intent(inout) :: kss
            character(len=*),     intent(in)    :: cfile
            integer, allocatable :: jj(:)
            integer :: m
            if( kss < 2 )then
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': fewer than two seed states (accepted '//&
                    &'poles and the rest), the node is refined as one state and becomes a leaf'
                if( kss < 1 )then
                    sds = merge(1, 0, project_states(d//'/'//NODE_PROJ) > 0)
                    jj  = junk_rows(d)
                    if( size(jj) > 0 ) sds(jj) = 0
                endif
                kss = 1
            endif
            if( any([(count(sds == m) == 0, m=1,kss)]) )then
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': ERROR, seed states 1..'//int2str(kss)//&
                    &' hold '//joined([(count(sds == m), m=1,kss)])//' particles; a state without particles cannot be '//&
                    &'reconstructed, node stopped'
                self%nodes(i)%step = S_FAILED
                return
            endif
            call write_contrast(d//'/'//NODE_PROJ, cfile)
            call simple_mkdir(string(d//'/seeds'))
            call write_ints(d//'/seeds/seeds.txt', sds)
            call write_states(d//'/'//NODE_PROJ, d//'/seeds/s.simple', sds)
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', 's.simple')
            call cl%set('nstates',  kss)
            call common_keys(self, cl)
            call cl%set('postprocess', 'no')
            call crop_keys(self, cl, BOXA)
            call self%submit(i, S_SEEDS, cl)
        end subroutine seed_node

        ! node_motion: every candidate gets its 5-bin reconstruction and call (the fflex fallback only when no site is a switch or
        ! swing); a switch is two seeded poles, a swing its poles and middle, and their product cells seed the split
        subroutine motion_next
            integer, allocatable :: cands(:), fac(:), sv(:,:), cv(:,:), g(:), jj(:), sds(:), start(:), lv(:), cv1(:), ord(:), off(:)
            character(len=:), allocatable :: how, w, one
            integer :: ic, kss, m, nr, nb
            logical :: l_anc
            allocate(bdirs(0), off(0))
            l_wait = .false.
            l_fl   = self%floor_ready()
            fl     = -huge(fl)
            if( l_fl .and. self%l_gate .and. self%nodes(i)%depth > 0 ) fl = self%floor_value()
            if( self%l_gate )then
                if( .not. file_exists(string(d//'/motion/gate_halves.txt')) ) call gate_halves
                hrow = read_ints(d//'/motion/gate_halves.txt')
                if( self%p%motion_floor <= 0. .and. floor_node(self%nodes(i)%id) ) call null_work
            endif
            l_anc = self%anc_wait(i)
            if( self%l_inherit .and. .not. l_anc .and. .not. file_exists(string(d//'/motion/inherited.txt')) ) call inherit_cands
            cands = read_cands(d)
            if( self%l_multi .or. size(cands) < 2 )then
                do ic = 1, size(cands)
                    call cand_work(cands(ic))
                end do
            else
                ! one factor is the strongest passing candidate: the strongest goes first, and when it passes the rest run off the path
                ord = strength_order(cands)
                nb  = size(bdirs)
                call cand_work(ord(1))
                if( size(bdirs) == nb .and. .not. l_wait )then
                    if( size(passing(ord(1:1))) == 1 )then
                        off = ord(2:)
                    else
                        do ic = 2, size(ord)
                            call cand_work(ord(ic))
                        end do
                    endif
                endif
            endif
            if( size(bdirs) == 0 .and. .not. l_wait .and. .not. l_anc .and. .not. any(cands == 0) .and. &
                &size(passing(cands)) == 0 )then
                g  = two_groups(d)
                jj = junk_rows(d)
                if( size(jj) > 0 .and. size(g) > 0 ) g(jj) = 0
                if( count(g == 1) > PROB_MULTISTATE_MIN_POP .and. count(g == 2) > PROB_MULTISTATE_MIN_POP )then
                    call add_cand(d, 0)
                    cands = read_cands(d)
                    call cand_work(0)
                endif
            endif
            if( size(bdirs) > 0 )then
                call write_lines(d//'/motion/batch.txt', bdirs)
                self%nodes(i)%step = S_MBATCH
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: '//int2str(size(bdirs))//&
                    &' jobs submitted or running'
                return
            endif
            if( l_wait )then
                if( self%nodes(i)%step /= S_GATEW ) write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//&
                    &': gate scores done, waiting for the tree''s floor (the permutation splits of r, r1 and r2)'
                self%nodes(i)%step = S_GATEW
                return
            endif
            if( l_anc )then
                if( self%nodes(i)%step /= S_GATEW ) write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//&
                    &': waiting for its ancestors'' off-path candidates, whose axes it inherits'
                self%nodes(i)%step = S_GATEW
                return
            endif
            if( size(off) > 0 ) cands = ord(1:1)
            do ic = 1, size(cands)
                call write_line(motion_dir(d, cands(ic))//'/verdict.txt', verdict(cands(ic)))
            end do
            fac = passing(cands)
            kss = 0
            how = 'no candidate calls switch or swing'
            one = ''
            if( .not. self%l_multi .and. size(fac) > 1 )then
                m   = strongest_of(fac)
                one = 'one factor (motion_multi=no), the strongest by pole-rule margin of '//int2str(size(fac))//' passing: '
                fac = [fac(m)]
            endif
            if( size(off) > 0 )then
                one = 'one factor (motion_multi=no), the strongest by pole-rule margin passes ('//int2str(size(off))//&
                    &' others run off the path): '
                call write_ints(d//'/motion/offpath.txt', off)
            endif
            if( size(fac) > 0 )then
                nr = size(project_states(d//'/'//NODE_PROJ))
                allocate(sv(nr,size(fac)), cv(nr,size(fac)), source=0)
                how = one
                do m = 1, size(fac)
                    call motion_factor(d, fac(m), motion_word(d, fac(m)) == 'swing', self%p%minpop, lv, cv1, w)
                    if( size(lv) == nr ) sv(:,m) = lv
                    if( size(cv1) == nr ) cv(:,m) = cv1
                    how = how//merge('; ', '  ', m > 1)//cand_label(d, fac(m))//' '//w
                end do
                start = product_cells(sv, 0)
                how   = trim(adjustl(how))//' -> product cells before folding (levels per motion: rows)'
                do m = 1, maxval([0, start])
                    jj  = sv(findloc(start, m, dim=1),:)
                    how = how//' ('//joined(jj, ',')//'): '//int2str(count(start == m))
                end do
                start = product_cells(sv, self%p%minpop)
                kss   = maxval([0, start])
                allocate(sds(nr), source=0)
                where( all(cv > 0, dim=2) ) sds = start
                do m = 1, kss
                    if( count(sds == m) < self%p%minpop ) where( start == m ) sds = m
                end do
                how = how//' -> K '//int2str(kss)//' after folding cells under minpop into their nearest: '//&
                    &joined([(count(start == m), m=1,kss)])//' rows, seeded from '//joined([(count(sds == m), m=1,kss)])//&
                    &' (its core when that holds minpop, else the whole cell)'
                do m = 1, size(fac)
                    if( fac(m) > 0 .and. fac(m) <= INH0 )then
                        call write_ints(d//'/sites/chosen_site.txt', [fac(m)])
                        exit
                    endif
                end do
            else
                allocate(sds(0))
            endif
            call node_motion_summary(d, read_cands(d), how)
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion hand-back: '//how
            call seed_node(sds, kss, d//'/motion/site_contrast.txt')
            if( kss >= 2 ) call write_ints(d//'/seeds/start.txt', start)
            if( size(off) > 0 ) write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: '//&
                &int2str(size(off))//' candidates run off the path: '//off_labels(off, .false.)
        end subroutine motion_next

        subroutine off_next
            integer, allocatable :: oc(:)
            integer :: ic, fu, io
            allocate(bdirs(0))
            l_wait = .false.
            l_fl   = self%floor_ready()
            fl     = -huge(fl)
            if( l_fl .and. self%l_gate .and. self%nodes(i)%depth > 0 ) fl = self%floor_value()
            if( self%l_gate ) hrow = read_ints(d//'/motion/gate_halves.txt')
            oc = read_ints(d//'/motion/offpath.txt')
            do ic = 1, size(oc)
                call cand_work(oc(ic))
            end do
            call write_lines(d//'/motion/offbatch.txt', bdirs)
            if( size(bdirs) > 0 .or. l_wait ) return
            call fopen(fu, string(d//'/motion_call.txt'), status='unknown', position='append', action='write', iostat=io)
            do ic = 1, size(oc)
                call write_line(motion_dir(d, oc(ic))//'/verdict.txt', verdict(oc(ic)))
                if( io == 0 ) write(fu,'(A)') 'off the path, '//cand_label(d, oc(ic))//': '//&
                    &first_line(motion_dir(d, oc(ic))//'/motion_call.txt', '-')//' | gate '//&
                    &first_line(motion_dir(d, oc(ic))//'/verdict.txt', '-')
            end do
            if( io == 0 ) call fclose(fu)
            call write_ints(d//'/motion/offpath_done.txt', oc)
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: the '//int2str(size(oc))//&
                &' off-path candidates are called: '//off_labels(oc, .true.)
        end subroutine off_next

        function off_labels( oc, lw ) result( str )
            integer, intent(in) :: oc(:)
            logical, intent(in) :: lw
            character(len=:), allocatable :: str
            integer :: ic
            str = ''
            do ic = 1, size(oc)
                str = str//merge(', ', '  ', ic > 1)//cand_label(d, oc(ic))
                if( lw ) str = str//' '//trim(motion_word(d, oc(ic)))//' ('//first_line(motion_dir(d, oc(ic))//'/verdict.txt', '-')//')'
            end do
            str = trim(adjustl(str))
        end function off_labels

        logical function gated( c )
            integer, intent(in) :: c
            gated = self%l_gate .and. self%nodes(i)%depth > 0 .and. c <= INH0
        end function gated

        real function cand_strength( c )
            integer, intent(in) :: c
            character(len=4096) :: ln
            integer :: fu, io
            cand_strength = -3.e6
            if( c > 0 .and. c <= INH0 )then
                cand_strength = kept_margin(site_dir_of(d, c))
            else if( c > INH0 )then
                cand_strength = -2.e6
                call fopen(fu, string(d//'/motion/inherited.txt'), status='old', action='read', iostat=io)
                if( io == 0 )then
                    do
                        read(fu,'(A)',iostat=io) ln
                        if( io /= 0 ) exit
                        if( word(ln, 1) /= int2str(c) .or. nwords(ln) < 3 ) cycle
                        if( str2int(word(ln, 3)) > 0 ) cand_strength = kept_margin(site_dir_of(self%root//'/node_'//&
                            &word(ln, 2), str2int(word(ln, 3)))) - 1.e6
                    end do
                    call fclose(fu)
                endif
            endif
        end function cand_strength

        integer function strongest_of( fac )
            integer, intent(in) :: fac(:)
            real    :: s, best
            integer :: m
            strongest_of = 1
            best = -huge(best)
            do m = 1, size(fac)
                s = cand_strength(fac(m))
                if( s > best )then
                    best = s
                    strongest_of = m
                endif
            end do
        end function strongest_of

        function strength_order( cs ) result( ord )
            integer, intent(in)  :: cs(:)
            integer, allocatable :: ord(:)
            real,    allocatable :: sv(:)
            real    :: x
            integer :: a, b, t
            ord = cs
            allocate(sv(size(cs)))
            do a = 1, size(cs)
                sv(a) = cand_strength(cs(a))
            end do
            do a = 2, size(ord)
                b = a
                do while( b > 1 )
                    if( sv(b) <= sv(b - 1) ) exit
                    t = ord(b)
                    ord(b) = ord(b - 1)
                    ord(b - 1) = t
                    x = sv(b)
                    sv(b) = sv(b - 1)
                    sv(b - 1) = x
                    b = b - 1
                end do
            end do
        end function strength_order

        function passing( cands ) result( fac )
            integer, intent(in)  :: cands(:)
            integer, allocatable :: fac(:)
            integer :: ic
            allocate(fac(0))
            do ic = 1, size(cands)
                if( .not. any(motion_word(d, cands(ic)) == [character(len=12) :: 'switch', 'insufficient', 'swing']) ) cycle
                if( index(verdict(cands(ic)), 'fail') == 1 ) cycle
                fac = [fac, cands(ic)]
            end do
        end function passing

        function verdict( c ) result( v )
            integer, intent(in) :: c
            character(len=:), allocatable :: v, gs
            real :: r
            integer :: ios
            if( .not. gated(c) )then
                v = 'ungated'
                return
            endif
            gs = first_line(motion_dir(d, c)//'/gate/gate.txt', 'fail: no gate score')
            read(gs,*,iostat=ios) r
            if( ios /= 0 )then
                v = gs
                if( index(v, 'fail') /= 1 ) v = 'fail: '//v
            else if( .not. l_fl )then
                v = 'pending: split-half corr '//gs//', the floor is not known yet'
            else if( r >= fl )then
                v = 'pass: split-half corr '//gs//' >= floor '//rstr(fl)
            else
                v = 'fail: split-half corr '//gs//' < floor '//rstr(fl)
            endif
        end function verdict

        subroutine cand_work( c )
            integer, intent(in) :: c
            character(len=1), parameter :: HN(2) = ['A', 'B']
            type(string), allocatable :: ends(:)
            integer,      allocatable :: ga(:), gb(:), grp(:)
            character(len=:), allocatable :: md, gd, ln
            real    :: r
            integer :: h
            logical :: ldone
            md = motion_dir(d, c)
            call simple_mkdir(string(md))
            if( .not. file_exists(string(md//'/motion_call.txt')) )then
                select case(self%dir_status(md))
                    case(-2)
                        call motion_submit(c)
                    case(-1)
                        bdirs = [bdirs, string(md)]
                    case default
                        call motion_stats(md, mask, nsym_of(trim(self%p%pgrp)), cand_label(d, c), root_axes(self, i, c), &
                            &cand_origin(d, c), ln)
                        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: '//ln
                end select
            endif
            if( .not. gated(c) ) return
            gd = md//'/gate'
            if( file_exists(string(gd//'/gate.txt')) )then
                if( .not. l_fl ) l_wait = .true.
                return
            endif
            call simple_mkdir(string(gd))
            ldone = .true.
            do h = 1, 2
                select case(self%dir_status(gd//'/'//HN(h)))
                    case(-2)
                        if( c == 0 )then
                            call gate_flex(gd//'/'//HN(h), h)
                        else
                            call gate_refit(gd//'/'//HN(h), c, h, d//'/motion/coords_'//HN(h)//'.txt')
                        endif
                        ldone = .false.
                    case(-1)
                        bdirs = [bdirs, string(gd//'/'//HN(h))]
                        ldone = .false.
                end select
            end do
            if( .not. ldone ) return
            select case(self%dir_status(gd//'/rec'))
                case(-2)
                    if( c == 0 )then
                        ga = two_groups_in(gd//'/A')
                        gb = two_groups_in(gd//'/B')
                    else
                        ends = kept_ends(site_dir_of(d, c))
                        ga   = post_groups(mean_post(gd//'/A', ends))
                        gb   = post_groups(mean_post(gd//'/B', ends))
                    endif
                    grp = groups4(ga, gb)
                    if( size(grp) > 0 )then
                        call gate_rec(gd//'/rec', grp)
                    else
                        call write_line(gd//'/gate.txt', 'fail: a half group holds under '//int2str(GATE_MINGRP)//' particles')
                        if( .not. l_fl ) l_wait = .true.
                    endif
                case(-1)
                    bdirs = [bdirs, string(gd//'/rec')]
                case default
                    r = gate_corr(gd//'/rec', self%p%mskdiam)
                    if( r < -1.5 )then
                        call write_line(gd//'/gate.txt', 'fail: the half reconstructions are missing')
                    else
                        if( c == 0 ) r = abs(r)
                        call write_line(gd//'/gate.txt', rstr(r))
                    endif
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' gate: '//cand_label(d, c)//&
                        &', split-half corr(D_A, D_B) at 15 A: '//first_line(gd//'/gate.txt', '-')
                    if( .not. l_fl ) l_wait = .true.
            end select
        end subroutine cand_work

        ! the gate's noise: random fifths of each particle half (permutation splits), reconstructed and scored as an axis
        subroutine null_work
            integer, allocatable :: grp(:)
            character(len=:), allocatable :: nd0, rd
            real    :: r
            integer :: j, fu, io, nperm
            logical :: ldone
            nd0 = d//'/motion/null'
            if( file_exists(string(nd0//'/gate.txt')) ) return
            call simple_mkdir(string(nd0))
            ldone = .true.
            do j = 1, NULL_PERM
                rd = nd0//'/perm/rec_p'//int2str(j)
                if( file_exists(string(rd//'/score.txt')) ) cycle
                select case(self%dir_status(rd))
                    case(-2)
                        call simple_mkdir(string(nd0//'/perm'))
                        grp = perm_groups(hrow, gate_seed(self%nodes(i)%id) + 9000_8 + int(j, 8))
                        if( size(grp) > 0 )then
                            call gate_rec(rd, grp)
                            ldone = .false.
                        else
                            call simple_mkdir(string(rd))
                            call write_line(rd//'/score.txt', 'fail')
                        endif
                    case(-1)
                        bdirs = [bdirs, string(rd)]
                        ldone = .false.
                    case default
                        r = gate_corr(rd, self%p%mskdiam)
                        if( r < -1.5 )then
                            call write_line(rd//'/score.txt', 'fail')
                        else
                            call write_line(rd//'/score.txt', rstr(r))
                        endif
                end select
            end do
            if( .not. ldone ) return
            call fopen(fu, string(nd0//'/gate.tmp'), status='replace', action='write', iostat=io)
            call fileiochk('hetero_tree: cannot write '//nd0//'/gate.tmp', io)
            nperm = 0
            do j = 1, NULL_PERM
                rd = first_line(nd0//'/perm/rec_p'//int2str(j)//'/score.txt', 'fail')
                if( rd == 'fail' ) cycle
                write(fu,'(A)') 'perm:p'//int2str(j)//' '//rd//' - 1'
                nperm = nperm + 1
            end do
            call fclose(fu)
            call simple_rename(string(nd0//'/gate.tmp'), string(nd0//'/gate.txt'))
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' gate null: '//int2str(nperm)//' of '//&
                &int2str(NULL_PERM)//' permutation splits scored'
        end subroutine null_work

        ! one random 50/50 split of the node's flex rows for every gate test of the node
        subroutine gate_halves
            type(string), allocatable :: hdr(:), body(:)
            integer,      allocatable :: perm(:), hl(:), hr(:)
            character(len=8192) :: ln
            character(len=:), allocatable :: cf
            integer(kind=8) :: seed
            integer :: fu, io, n, nh, j, p, nrp
            cf = d//'/fflex/flex_pca_coordinates.txt'
            call fopen(fu, string(cf), status='old', action='read', iostat=io)
            call fileiochk('hetero_tree: cannot read '//cf, io)
            n  = 0
            nh = 0
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                if( len_trim(ln) == 0 ) cycle
                if( ln(1:1) == '#' )then
                    nh = nh + 1
                else
                    n = n + 1
                endif
            end do
            allocate(hdr(nh), body(n))
            rewind(fu)
            n  = 0
            nh = 0
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                if( len_trim(ln) == 0 ) cycle
                if( ln(1:1) == '#' )then
                    nh      = nh + 1
                    hdr(nh) = string(trim(ln))
                else
                    n       = n + 1
                    body(n) = string(trim(ln))
                endif
            end do
            call fclose(fu)
            seed = gate_seed(self%nodes(i)%id)
            perm = lcg_perm(n, seed)
            allocate(hl(n), source=2)
            hl(perm(1:n/2)) = 1
            nrp = size(project_states(d//'/'//NODE_PROJ))
            allocate(hr(nrp), source=0)
            do j = 1, n
                p = str2int(word(body(j)%to_char(), 1))
                if( p >= 1 .and. p <= nrp ) hr(p) = hl(j)
            end do
            call simple_mkdir(string(d//'/motion'))
            call write_coords(d//'/motion/coords_A.txt', hdr, pack(body, hl == 1))
            call write_coords(d//'/motion/coords_B.txt', hdr, pack(body, hl == 2))
            call write_ints(d//'/motion/gate_halves.txt', hr)
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' gate: random halves of its '//int2str(n)//&
                &' flex rows, '//int2str(count(hl == 1))//' / '//int2str(count(hl == 2))
        end subroutine gate_halves

        subroutine inherit_cands
            integer, allocatable :: ac(:), mx(:), ab(:), s0(:), abn(:), s0n(:)
            real,    allocatable :: sc(:), scn(:)
            character(len=:), allocatable :: id, ad, md, wd
            integer :: lp, ia, c, j, fu, io, ninh
            id   = self%nodes(i)%id
            c    = INH0
            ninh = 0
            call simple_mkdir(string(d//'/motion'))
            call fopen(fu, string(d//'/motion/inherited.tmp'), status='replace', action='write', iostat=io)
            call fileiochk('hetero_tree: cannot write '//d//'/motion/inherited.tmp', io)
            do lp = 1, len(id) - 1
                ad = self%root//'/node_'//id(1:lp)
                ac = read_cands(ad)
                if( size(ac) == 0 ) cycle
                mx = anc_index(ad, d)
                do ia = 1, size(ac)
                    if( ac(ia) > INH0 ) cycle
                    wd = trim(motion_word(ad, ac(ia)))
                    if( wd /= 'switch' .and. wd /= 'swing' ) cycle
                    if( index(first_line(motion_dir(ad, ac(ia))//'/verdict.txt', 'ungated'), 'fail') == 1 ) cycle
                    sc = read_reals(motion_dir(ad, ac(ia))//'/score.txt')
                    ab = read_ints(motion_dir(ad, ac(ia))//'/ab.txt')
                    s0 = read_ints(motion_dir(ad, ac(ia))//'/st0.txt')
                    if( size(mx) == 0 .or. size(sc) == 0 .or. size(ab) /= size(sc) .or. size(s0) /= size(sc) ) cycle
                    if( maxval(mx) > size(sc) ) cycle
                    allocate(scn(size(mx)), source=-1.)
                    allocate(abn(size(mx)), s0n(size(mx)), source=0)
                    do j = 1, size(mx)
                        if( mx(j) < 1 ) cycle
                        scn(j) = sc(mx(j))
                        abn(j) = ab(mx(j))
                        s0n(j) = s0(mx(j))
                    end do
                    if( count(scn >= 0.) >= NMBIN * GATE_MINGRP )then
                        c  = c + 1
                        md = motion_dir(d, c)
                        call simple_mkdir(string(md))
                        call write_reals(md//'/inherit_score.txt', scn)
                        call write_ints(md//'/inherit_ab.txt', abn)
                        call write_ints(md//'/inherit_st0.txt', s0n)
                        call write_line(md//'/inherit.txt', id(1:lp)//' '//cand_name(ac(ia))//' ('//wd//' there)')
                        call add_cand(d, c)
                        write(fu,'(A)') int2str(c)//' '//id(1:lp)//' '//int2str(ac(ia))//' '//wd//' '//int2str(count(scn >= 0.))
                        ninh = ninh + 1
                    endif
                    deallocate(scn, abn, s0n)
                end do
            end do
            call fclose(fu)
            call simple_rename(string(d//'/motion/inherited.tmp'), string(d//'/motion/inherited.txt'))
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: '//int2str(ninh)//' ancestor axes '//&
                &'(switch or swing there) join the candidates'
        end subroutine inherit_cands

        subroutine motion_submit( c )
            integer, intent(in) :: c
            type(cmdline) :: cj
            real,    allocatable :: sc(:)
            integer, allocatable :: bins(:)
            character(len=:), allocatable :: md
            integer :: m
            md = motion_dir(d, c)
            call simple_mkdir(string(md))
            sc   = cand_score(d, c)
            bins = quint_bins(sc)
            call write_reals(md//'/score.txt', sc)
            call write_ints(md//'/bins.txt', bins)
            call write_states(d//'/recB0/b.simple', md//'/m.simple', bins)
            call cj%set('prg',         'reconstruct3D')
            call cj%set('projfile',    'm.simple')
            call cj%set('nstates',     NMBIN)
            call cj%set('postprocess', 'no')
            call common_keys(self, cj)
            call self%submit_in(md, cj, count(bins > 0))
            bdirs = [bdirs, string(md)]
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' motion: '//cand_label(d, c)//', '//&
                &int2str(count(bins > 0))//' rows in '//int2str(NMBIN)//' score bins of '//&
                &joined([(count(bins == m), m=1,NMBIN)])//', their reconstruction submitted'
        end subroutine motion_submit

        subroutine gate_refit( sdir, kk, h, cfile )
            character(len=*), intent(in) :: sdir, cfile
            integer,          intent(in) :: kk, h
            type(cmdline) :: cj
            call simple_mkdir(string(sdir))
            call write_states(d//'/'//NODE_PROJ, sdir//'/s.simple', merge(1, 0, hrow == h))
            call cj%set('prg',         'site_rank')
            call cj%set('projfile',    's.simple')
            call cj%set('infile',      cfile)
            call cj%set('site_dirs',   trim(self%p%site_dirs))
            call cj%set('pole_merge',  self%p%pole_merge)
            if( self%p%pooled_poles .eq. 'yes' ) call cj%set('pooled_poles', 'yes')
            call cj%set('junk_screen', trim(self%p%junk_screen))
            call cj%set('pcg_mskfile', cand_mask(d, kk))
            call cj%set('pgrp',        trim(self%p%pgrp))
            call cj%set('mskdiam',     self%p%mskdiam)
            call cj%set('nthr',        max(self%p%nthr, SITE_NTHR))
            if( local_rows(d) ) call cj%set('rowmap', d//'/rows.txt')
            call self%submit_in(sdir, cj, count(hrow == h))
            bdirs = [bdirs, string(sdir)]
        end subroutine gate_refit

        subroutine gate_flex( sdir, h )
            character(len=*), intent(in) :: sdir
            integer,          intent(in) :: h
            type(cmdline) :: cj
            character(len=:), allocatable :: vm
            integer :: nh
            call simple_mkdir(string(sdir))
            nh = count(hrow == h)
            call write_states(d//'/'//NODE_PROJ, sdir//'/f.simple', merge(1, 0, hrow == h))
            if( self%l_oneflex .and. self%nodes(i)%depth > 0 )then
                vm = d//'/pmean.mrc'
            else
                vm = d//'/mean/recvol_state01.mrc'
            endif
            call self%set_parts(nh)
            call flex_keys(self, cj, vm, min(KMAX, max(2, nh / max(1, self%p%kfloor))))
            call cj%set('pcg_mskfile', d//'/mask.mrc')
            if( self%p%flex_min_neff > 0 ) call cj%set('min_neff', self%p%flex_min_neff)
            if( local_rows(d) ) call cj%set('rowmap', d//'/rows.txt')
            call self%submit_in(sdir, cj)
            call self%set_parts(size(rows))
            bdirs = [bdirs, string(sdir)]
        end subroutine gate_flex

        subroutine gate_rec( sdir, grp )
            character(len=*), intent(in) :: sdir
            integer,          intent(in) :: grp(:)
            type(cmdline) :: cj
            call simple_mkdir(string(sdir))
            call write_states(d//'/'//NODE_PROJ, sdir//'/r.simple', grp)
            call cj%set('prg',         'reconstruct3D')
            call cj%set('projfile',    'r.simple')
            call cj%set('nstates',     4)
            call cj%set('pgrp',        trim(self%p%pgrp))
            call cj%set('mskdiam',     self%p%mskdiam)
            call cj%set('postprocess', 'no')
            call cj%set('nthr',        max(self%p%nthr, SITE_NTHR))
            if( self%box > GATE_BOX )then
                call cj%set('box_crop',  GATE_BOX)
                call cj%set('smpd_crop', real(self%box) * self%smpd / real(GATE_BOX))
            endif
            call self%submit_in(sdir, cj, count(grp > 0))
            bdirs = [bdirs, string(sdir)]
        end subroutine gate_rec

        subroutine log_lines( fname, what )
            character(len=*), intent(in) :: fname, what
            character(len=4096) :: ln
            integer :: fu, io
            if( .not. file_exists(string(fname)) ) return
            call fopen(fu, string(fname), status='old', action='read', iostat=io)
            if( io /= 0 ) return
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' '//what//': '//trim(ln)
            end do
            call fclose(fu)
        end subroutine log_lines

        function rfmt2( x ) result( s )
            real, intent(in) :: x
            character(len=8) :: s
            write(s,'(F5.2)') x
        end function rfmt2

    end subroutine advance

    integer function prev_step( self, i, istep )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i, istep
        prev_step = istep - 1
        select case(istep)
            case(S_FFLEX)
                if( self%l_oneflex .and. self%nodes(i)%depth > 0 )then
                    prev_step = S_QUEUED
                else if( self%nodes(i)%depth == 0 .and. .not. self%l_subtree )then
                    prev_step = S_FLEX
                endif
            case(S_SITE)
                if( self%p%site_disc .ne. 'yes' ) prev_step = S_FFLEX
            case(S_SBATCH)
                prev_step = merge(S_RECB0, S_SITE, self%l_motion .or. self%l_strong .or. &
                    &first_line(self%ndir(i)//'/sites/rest_first.txt', '2') == '1')
            case(S_MOTION, S_MBATCH, S_GATEW)
                prev_step = S_SITE
        end select
    end function prev_step

    subroutine submit_fflex( self, i, vmean )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        character(len=*),   intent(in)    :: vmean
        type(cmdline)        :: cl
        integer, allocatable :: rows(:)
        character(len=:), allocatable :: d
        d    = self%ndir(i)
        rows = read_ints(d//'/rows.txt')
        call simple_mkdir(string(d//'/fflex'))
        call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/fflex/f.simple'))
        call flex_keys(self, cl, vmean, min(KMAX, max(2, size(rows) / max(1, self%p%kfloor))))
        call cl%set('pcg_mskfile', d//'/mask.mrc')
        if( self%p%flex_min_neff > 0 ) call cl%set('min_neff', self%p%flex_min_neff)
        if( local_rows(d) ) call cl%set('rowmap', d//'/rows.txt')
        ! its state maps feed only stA0's references and the two-group fallback's 40-10 A tree: stage A's box suffices
        if( self%box > BOXA ) call cl%set('rec_states_box', BOXA)
        call self%submit(i, S_FFLEX, cl)
    end subroutine submit_fflex

    subroutine submit( self, i, istep, cl, nmem )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, istep
        class(cmdline),     intent(inout) :: cl
        integer, optional,  intent(in)    :: nmem
        type(string) :: pf
        if( self%p%node_rows .eq. 'yes' .and. any(istep == [S_MEAN, S_CHEAP, S_STA0, S_RECB0, S_SEEDS, S_STAGEA, S_RECB, &
            &S_STAGEB]) )then
            pf = cl%get_carg('projfile')
            call compact_proj(self%ndir(i)//'/'//trim(STEP_DIRS(istep)), pf%to_char(), self%p%node_subset_io .eq. 'yes')
        endif
        call self%submit_in(self%ndir(i)//'/'//trim(STEP_DIRS(istep)), cl, nmem)
        self%nodes(i)%step = istep
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': '//trim(STEP_DIRS(istep))//' submitted'
    end subroutine submit

    subroutine submit_in( self, sd, cl, nmem )
        class(hetero_tree), intent(inout) :: self
        character(len=*),   intent(in)    :: sd
        class(cmdline),     intent(inout) :: cl
        integer, optional,  intent(in)    :: nmem
        call cl%set('mkdir', 'no')
        call simple_chdir(string(sd))
        CWD_GLOB = sd
        if( cl%get_iarg('nthr') > self%p%nthr )then
            call queue_in(self%qsite)
        else
            call queue_in(self%qenv)
        endif
        call simple_chdir(string(self%root))
        CWD_GLOB = self%root
        if( .not. allocated(self%subd) ) allocate(self%subd(0))
        self%subd = [self%subd, string(sd)]

    contains

        ! a one-job step over nmem particles asks 8 GB + 2.5 GB per 100k (measured about 2 GB per 100k at 692k)
        subroutine queue_in( q )
            use simple_persistent_worker_server, only: persistent_worker
            type(qsys_env), intent(inout) :: q
            type(string) :: mem0, memclaim, prg, cpu0
            integer :: task_mem, task_nthr, fu, io, ncpu
            logical :: l_pool
            if( present(nmem) )then
                mem0 = q%qdescr%get('job_memory_per_task')
                call q%qdescr%set('job_memory_per_task', int2str(8000 + nint(0.025 * real(nmem))))
            endif
            ncpu = 0
            if( cl%defined('nthr') ) ncpu = cl%get_iarg('nthr')
            if( ncpu > 0 )then
                cpu0 = q%qdescr%get('job_cpus_per_task')
                call q%qdescr%set('job_cpus_per_task', int2str(ncpu))
            endif
            l_pool = .false.
            if( self%p%node_worker_pool == 'yes' )then
                call cl%set('qsys_name', self%qenv%get_qsys())
                prg = cl%get_carg('prg')
                select case(prg%to_char())
                    case('state_split', 'state_mergeup', 'site_rank')
                        task_nthr = cl%get_iarg('nthr')
                        memclaim = q%qdescr%get('job_memory_per_task')
                        task_mem = str2int(memclaim%to_char())
                        l_pool = task_nthr > 0 .and. task_nthr <= self%p%worker_nthr .and. &
                            &task_mem > 0 .and. task_mem <= self%p%node_worker_mem
                end select
            endif
            if( l_pool )then
                if( .not. allocated(self%qpool) )then
                    if( associated(persistent_worker%server) ) THROW_HARD('node_worker_pool requires exclusive server ownership')
                    ! Launch from the root so worker scripts have stable paths across task directories.
                    call simple_chdir(string(self%root))
                    CWD_GLOB = self%root
                    allocate(self%qpool)
                    call self%qpool%new(self%p, 1, exec_bin=string('simple_exec'), &
                        &qsys_name=string(self%qenv%get_qsys()//'_worker'), qsys_nthr=self%p%worker_nthr, &
                        &qsys_memory=self%p%node_worker_mem, strict_dispatch=.true.)
                    write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> HETERO_TREE worker pool: ', self%p%workers, &
                        &' slots, ', self%p%worker_nthr, ' threads and ', self%p%node_worker_mem, ' MB per slot'
                    call simple_chdir(string(sd))
                    CWD_GLOB = sd
                endif
                ! A complete slot is reserved, so concurrent tasks cannot overcommit its memory.
                call cl%set('qsys_name', 'local')
                call cl%delete('nparts')
                call cl%delete('worker_server')
                call cl%delete('worker_server_nthr')
                if( file_exists(string(POOL_PENDING)) )then
                    if( .not. file_exists(string(STEP_EXIT)) ) THROW_HARD('worker task already submitted in '//sd)
                    call del_file(string(POOL_PENDING))
                endif
                if( file_exists(string(STEP_EXIT)) ) call del_file(string(STEP_EXIT))
                if( file_exists(string(STEP_EXIT//JOB_INFO_EXT)) ) call del_file(string(STEP_EXIT//JOB_INFO_EXT))
                call fopen(fu, string(POOL_PENDING), status='new', action='write', iostat=io)
                call fileiochk('unfinished or duplicate worker submission in '//sd, io)
                write(fu,'(A)') prg%to_char()//' '//int2str(task_nthr)//' '//int2str(task_mem)
                call fclose(fu)
                call self%qpool%exec_simple_prg_in_queue_async(cl, string('./hetero_step.sh'), string('hetero_step.log'), &
                    &exit_code_fname=string(STEP_EXIT))
                write(logfhandle,'(A)') '>>> HETERO_TREE pooled '//prg%to_char()//' in '//sd
            else
                call q%exec_simple_prg_in_queue_async(cl, string('./hetero_step.sh'), string('hetero_step.log'), &
                    &exit_code_fname=string(STEP_EXIT))
            endif
            if( present(nmem) ) call q%qdescr%set('job_memory_per_task', mem0)
            if( ncpu > 0 ) call q%qdescr%set('job_cpus_per_task', cpu0)
        end subroutine queue_in

    end subroutine submit_in

    ! the distributed steps of a node get one part per PTCLS_PER_PART particles, at least NPARTS_MIN, at most nparts
    subroutine set_parts( self, n )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: n
        self%npcur = min(self%p%nparts, max(NPARTS_MIN, ceiling(real(n) / real(PTCLS_PER_PART))))
    end subroutine set_parts

    ! a job submitted by this run writes its record only once the queue starts it
    logical function pending( self, sd )
        class(hetero_tree), intent(in) :: self
        character(len=*),   intent(in) :: sd
        integer :: j
        pending = .false.
        if( .not. allocated(self%subd) ) return
        do j = 1, size(self%subd)
            if( self%subd(j)%to_char() == sd ) pending = .true.
        end do
    end function pending

    logical function final( self )
        class(hetero_tree), intent(inout) :: self
        logical :: l_fac
        final = .false.
        if( .not. self%plain() ) return
        if( len(self%subroot) > 0 )then
            if( .not. self%l_lsdone )then
                if( .not. self%landscape() ) return
                self%l_lsdone = .true.
            endif
            final = .true.
            return
        endif
        l_fac = self%factorial()
        if( .not. self%l_lsdone )then
            if( .not. self%landscape() ) return
            self%l_lsdone = .true.
        endif
        ! the merge-up job goes first and the final maps are made while it runs (neither reads the other's files)
        if( .not. self%l_mudone ) self%l_mudone = self%mergeup()
        if( .not. file_exists(string(self%root//'/final_maps/choice.txt')) ) call final_maps(self)
        if( .not. self%l_mudone ) return
        final = l_fac
    end function final

    ! the landscape between the accepted poles of every finished node (landscape_sample in node_<id>/landscape)
    logical function landscape( self )
        class(hetero_tree), intent(inout) :: self
        type(cmdline) :: cl
        character(len=:), allocatable :: nd, ld
        integer :: i, nrun, np, code
        logical :: err
        landscape = .true.
        nrun      = 0
        do i = 1, size(self%nodes)
            if( self%nodes(i)%step /= S_DONE ) cycle
            nd = self%ndir(i)
            ld = nd//'/landscape'
            if( .not. file_exists(string(nd//'/split/split_children.txt')) ) cycle
            if( file_exists(string(ld//'/'//STEP_EXIT)) ) cycle
            np = count_poles(seed_dir(nd)//'/site_poles.txt')
            if( np < 2 )then
                call simple_mkdir(string(ld))
                call write_ints(ld//'/landscape_pairs.txt', [integer ::])
                call write_ints(ld//'/'//STEP_EXIT, [0])
                cycle
            endif
            landscape = .false.
            if( file_exists(string(ld//'/'//STEP_EXIT//JOB_INFO_EXT)) .or. self%pending(ld) )then
                nrun = nrun + 1
                cycle
            endif
            if( nrun >= max(1, self%p%nodepar) ) cycle
            call simple_mkdir(string(ld))
            call cl%set('prg',      'landscape_sample')
            call cl%set('projfile', nd//'/stageB/b.simple')
            call cl%set('dir',      nd//'/stageB')
            call cl%set('infile',   seed_dir(nd)//'/site_poles.txt')
            call cl%set('infile2',  nd//'/split/split_children.txt')
            call cl%set('mskdiam',  self%p%mskdiam)
            call cl%set('pgrp',     trim(self%p%pgrp))
            call cl%set('nthr',     self%p%nthr)
            call self%submit_in(ld, cl, nlines(string(nd//'/rows.txt')))
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': landscape submitted, '//int2str(np)//&
                &' accepted poles'
            nrun = nrun + 1
        end do
        if( .not. landscape ) return
        do i = 1, size(self%nodes)
            if( self%nodes(i)%step /= S_DONE ) cycle
            ld = self%ndir(i)//'/landscape'
            if( .not. file_exists(string(ld//'/'//STEP_EXIT)) ) cycle
            np = count_poles(seed_dir(self%ndir(i))//'/site_poles.txt')
            call read_exit_code(string(ld//'/'//STEP_EXIT), code, err)
            if( np < 2 )then
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': landscape skipped, '//int2str(np)//&
                    &' accepted pole(s)'
            else if( code /= 0 .or. err )then
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': landscape FAILED (exit code '//&
                    &int2str(code)//')'
            else
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': landscape of '//int2str(np)//&
                    &' accepted poles, '//int2str(nlines(string(ld//'/landscape_pairs.txt')))//' pairs'
            endif
        end do
    end function landscape

    logical function plain( self )
        class(hetero_tree), intent(inout) :: self
        type(cmdline)        :: cl
        type(string), allocatable :: lines(:), nds(:)
        integer, allocatable :: lab(:)
        character(len=4096)  :: line
        character(len=:), allocatable :: pd
        integer :: nl, il, inode, nrun, nleaf
        plain = .true.
        if( .not. file_exists(string(self%ltab)) ) return
        call leaf_lines(lines, nds, self%ltab)
        nl   = size(lines)
        nrun = 0
        do inode = 1, size(nds)
            if( any([(nds(il) == nds(inode), il=1,inode-1)]) ) cycle
            pd = nds(inode)%to_char()//'/plain'
            if( file_exists(string(pd//'/'//STEP_EXIT)) ) cycle
            plain = .false.
            if( file_exists(string(pd//'/'//STEP_EXIT//JOB_INFO_EXT)) .or. self%pending(pd) )then
                nrun = nrun + 1
                cycle
            endif
            if( nrun >= max(1, self%p%nodepar) ) cycle
            allocate(lab(self%nrows), source=0)
            nleaf = 0
            do il = 1, nl
                if( nds(il) /= nds(inode) ) cycle
                nleaf = nleaf + 1
                line  = lines(il)%to_char()
                lab(read_ints(word(line, 8))) = nleaf
            end do
            call simple_mkdir(string(pd))
            if( local_rows(nds(inode)%to_char()) ) lab = lab(read_ints(nds(inode)%to_char()//'/rows.txt'))
            call write_states(nds(inode)%to_char()//'/stageB/b.simple', pd//'/p.simple', lab)
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', 'p.simple')
            call cl%set('nstates',  nleaf)
            call self%set_parts(count(lab > 0))
            call common_keys(self, cl)
            call self%submit_in(pd, cl)
            write(logfhandle,'(A)') '>>> HETERO_TREE plain leaf maps of '//pd//': '//int2str(nleaf)//' leaves submitted'
            nrun = nrun + 1
            deallocate(lab)
        end do
    end function plain

    logical function mergeup( self )
        class(hetero_tree), intent(inout) :: self
        type(cmdline)        :: cl
        integer, allocatable :: units(:), leaves(:), ju(:)
        type(string), allocatable :: lines(:), nds(:)
        character(len=:), allocatable :: md, f
        logical :: err
        integer :: code, nn, j, np
        md      = self%root//'/mergeup'
        f       = md//'/'//STEP_EXIT
        mergeup = .true.
        if( file_exists(string(f)) )then
            call read_exit_code(string(f), code, err)
            if( code == 0 .and. file_exists(string(md//'/merged_units.txt')) )then
                units  = read_ints(md//'/merged_units.txt')
                leaves = read_ints(self%root//'/hetero_tree_labels.txt')
                ju     = pack(leaves, leaves > 0 .and. units == 0)
                do j = 1, size(ju)
                    if( count(ju(1:j-1) == ju(j)) > 0 ) cycle
                    where( leaves == ju(j) .and. units == 0 ) units = maxval(units) + 1
                end do
                call write_ints(self%root//'/hetero_tree_merged.txt', units)
                write(logfhandle,'(A)') '>>> HETERO_TREE merge-up of the leaves: '//int2str(size(present_labels(units)))//&
                    &' units over '//int2str(count(units > 0))//' particles: hetero_tree_merged.txt'
            else
                write(logfhandle,'(A)') '>>> HETERO_TREE merge-up of the leaves FAILED (exit code '//int2str(code)//')'
            endif
            return
        endif
        mergeup = .false.
        if( file_exists(string(f//JOB_INFO_EXT)) .or. self%pending(md) ) return
        call self%finish(nn)
        if( nn < 2 )then
            mergeup = .true.
            return
        endif
        call simple_mkdir(string(md))
        call final_project(self, md//'/final.simple')
        call leaf_masks(md//'/leaf_masks.txt')
        call leaf_lines(lines, nds)
        np = sum([(merge(str2int(word(lines(j)%to_char(), 3)), 0, word(lines(j)%to_char(), 2) /= 'junk'), j=1,size(lines))])
        call cl%set('prg',      'state_mergeup')
        call cl%set('projfile', 'final.simple')
        call cl%set('infile',   self%root//'/leaf_tree/merge_tree.txt')
        call cl%set('infile2',  md//'/leaf_masks.txt')
        call cl%set('mskdiam',  self%p%mskdiam)
        call cl%set('pgrp',     trim(self%p%pgrp))
        call cl%set('junk_screen', trim(self%p%junk_screen))
        if( self%p%merge_geom_reuse == 'yes' ) call cl%set('merge_geom_reuse', 'yes')
        if( self%p%merge_insert_reuse == 'yes' ) call cl%set('merge_insert_reuse', 'yes')
        call cl%set('nthr',     max(self%p%nthr, min(MERGEUP_NTHR, ceiling(real(np) / real(MERGEUP_PTCLS)))))
        call self%submit_in(md, cl, np)
        write(logfhandle,'(A)') '>>> HETERO_TREE merge-up of the '//int2str(nn)//' leaves submitted'
    end function mergeup

    subroutine fs_sync( self )
        class(hetero_tree), intent(inout) :: self
        type(fsite), allocatable :: tmp(:)
        character(len=4096) :: ln
        character(len=:), allocatable :: f, sd
        integer :: i, j, funit, ios, nrun, code
        logical :: l_reg, l_chg
        if( .not. self%l_fsread )then
            call self%fs_read_table
            self%l_fsread = .true.
        endif
        l_chg = .false.
        do i = 1, size(self%nodes)
            f = self%ndir(i)//'/sites/sites.txt'
            if( .not. file_exists(string(f)) ) cycle
            l_reg = .false.
            do j = 1, size(self%fs)
                if( self%fs(j)%home == self%nodes(i)%id ) l_reg = .true.
            end do
            if( l_reg ) cycle
            call fopen(funit, string(f), status='old', action='read', iostat=ios)
            if( ios /= 0 ) cycle
            do
                read(funit,'(A)',iostat=ios) ln
                if( ios /= 0 ) exit
                if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
                allocate(tmp(size(self%fs) + 1))
                tmp(1:size(self%fs)) = self%fs
                j = size(tmp)
                tmp(j)%k      = str2int(word(ln, 1))
                tmp(j)%name   = self%nodes(i)%id//'_s'//int2str(tmp(j)%k)
                tmp(j)%home   = self%nodes(i)%id
                tmp(j)%level  = self%nodes(i)%id
                tmp(j)%target = '-'
                tmp(j)%repro  = str2real(word(ln, 3))
                call move_alloc(tmp, self%fs)
                l_chg = .true.
            end do
            call fclose(funit)
        end do
        nrun = self%fs_nrun()
        do j = 1, size(self%fs)
            if( self%fs(j)%st /= F_RUN ) cycle
            sd   = self%fs_hd(j)
            code = self%dir_status(sd)
            if( code == -1 ) cycle
            if( code == -2 )then
                if( self%fs(j)%k == 1 .or. nrun >= max(1, self%p%nodepar) ) cycle
                call self%fs_rank(sd, self%fs(j)%home, site_mask(self, j), '', '', SITE_NFIT, .false.)
                write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': ranking at its node submitted'
                nrun = nrun + 1
                cycle
            endif
            self%fs(j)%st = merge(F_TEST, F_FAIL, code == 0)
            if( code /= 0 ) write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': ranking FAILED (exit code '//&
                &int2str(code)//')'
            l_chg = .true.
        end do
        if( l_chg ) call self%fs_write_table
    end subroutine fs_sync

    logical function factorial( self )
        class(hetero_tree), intent(inout) :: self
        integer, allocatable :: ev(:), ord(:)
        real,    allocatable :: key(:)
        character(len=:), allocatable :: sd
        integer :: j, a, b, dmin, nrun, code, best, np
        real    :: r, rbest
        logical :: l_chg, l_dup, l_ok
        factorial = .true.
        if( self%p%site_layer .ne. 'yes' ) return
        call self%fs_sync
        factorial = .false.
        l_chg     = .false.
        nrun      = self%fs_nrun()
        do j = 1, size(self%fs)
            if( self%fs(j)%st /= F_UP ) cycle
            sd   = self%fs_cd(j)
            code = self%dir_status(sd)
            if( code == -1 ) cycle
            if( code == -2 )then
                if( nrun >= max(1, self%p%nodepar) ) cycle
                call self%fs_climb(j)
                nrun = nrun + 1
                cycle
            endif
            l_chg = .true.
            np    = count_poles(sd//'/site_poles.txt')
            r     = 0.
            if( code == 0 ) r = self%fs_rho(sd, self%fs_hd(j), self%fs_eval_rows(self%fs(j)%home), .true.)
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': pooling test '//self%fs(j)%level//' -> '//&
                &parent_of(self%fs(j)%level)//', first pole Spearman on the '//self%fs(j)%home//' particles '//rstr(r)//', '//&
                &int2str(np)//' pole(s) pass the halves test there'//trim(merge(' (run FAILED)', '             ', code /= 0))
            if( code /= 0 .or. r < FGATE )then
                self%fs(j)%st = F_LEVEL
                cycle
            endif
            self%fs(j)%level = parent_of(self%fs(j)%level)
            self%fs(j)%rho   = r
            if( len(self%fs(j)%level) == 1 )then
                self%fs(j)%st = F_LEVEL
            else if( nrun < max(1, self%p%nodepar) )then
                call self%fs_climb(j)
                nrun = nrun + 1
            endif
        end do
        do j = 1, size(self%fs)
            if( any(self%fs(j)%st == [F_PREC, F_PREF]) ) call self%fs_prerefine(j, nrun, l_chg)
        end do
        dmin = huge(dmin)
        do j = 1, size(self%fs)
            if( any(self%fs(j)%st == [F_RUN, F_TEST, F_UP, F_PREC, F_PREF]) ) dmin = min(dmin, len(self%fs(j)%home) - 1)
        end do
        if( dmin < huge(dmin) )then
            key = [(-self%fs(j)%repro, j=1,size(self%fs))]
            ord = [(j, j=1,size(self%fs))]
            call hpsort(key, ord)
            do b = 1, size(ord)
                j = ord(b)
                if( self%fs(j)%st /= F_TEST .or. len(self%fs(j)%home) - 1 /= dmin ) cycle
                if( npoles_of(self%fs_hd(j)) == 0 .and. self%fs(j)%ref == 0 )then
                    l_ok = best_minz(self%fs_hd(j)) >= 2. .and. best_minz(self%fs_hd(j)) < 3.
                    if( l_ok )then
                        if( nrun < max(1, self%p%nodepar) ) call self%fs_prerefine(j, nrun, l_chg)
                        cycle
                    endif
                endif
                if( npoles_of(self%fs_hd(j)) == 0 )then
                    self%fs(j)%st = F_NOPOLE
                    l_chg = .true.
                    write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': no accepted pole, not a site'
                    cycle
                endif
                ev    = self%fs_eval_rows(self%fs(j)%home)
                rbest = 0.
                best  = 0
                do a = 1, size(self%fs)
                    if( self%fs(a)%st /= F_LEVEL .or. len(self%fs(a)%home) - 1 >= dmin ) cycle
                    if( .not. is_anc(self%fs(a)%level, self%fs(j)%home) ) cycle
                    r = self%fs_rho(self%fs_rd(a), self%fs_hd(j), ev, .false.)
                    if( r > rbest )then
                        rbest = r
                        best  = a
                    endif
                end do
                if( best > 0 .and. rbest >= FGATE )then
                    self%fs(j)%st     = F_COPY
                    self%fs(j)%target = self%fs(best)%name
                    self%fs(j)%rho    = rbest
                    self%fs(j)%dice   = mask_dice(site_mask(self, j), site_mask(self, best))
                    write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//' = site '//self%fs(best)%name//&
                        &' (pooled at '//self%fs(best)%level//'): Spearman on the '//self%fs(j)%home//' particles '//rstr(rbest)//&
                        &', mask Dice '//rstr(self%fs(j)%dice)
                    l_chg = .true.
                    cycle
                endif
                rbest = 0.
                best  = 0
                do a = 1, size(self%fs)
                    if( a == j .or. self%fs(a)%home /= self%fs(j)%home ) cycle
                    if( .not. any(self%fs(a)%st == [F_UP, F_LEVEL, F_COPY]) ) cycle
                    r = self%fs_rho(self%fs_hd(a), self%fs_hd(j), ev, .false.)
                    if( r > rbest )then
                        rbest = r
                        best  = a
                    endif
                end do
                if( best > 0 .and. rbest >= FGATE )then
                    self%fs(j)%st     = F_DUP
                    self%fs(j)%target = self%fs(best)%name
                    self%fs(j)%rho    = rbest
                    self%fs(j)%dice   = mask_dice(site_mask(self, j), site_mask(self, best))
                    write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//' duplicates site '//self%fs(best)%name//&
                        &' at its node '//self%fs(j)%home//', not pooled: Spearman '//rstr(rbest)//', mask Dice '//&
                        &rstr(self%fs(j)%dice)
                else if( len(self%fs(j)%home) == 1 )then
                    self%fs(j)%st = F_LEVEL
                else
                    if( nrun >= max(1, self%p%nodepar) ) cycle
                    self%fs(j)%st  = F_UP
                    self%fs(j)%rho = rbest
                    call self%fs_climb(j)
                    nrun = nrun + 1
                endif
                l_chg = .true.
            end do
            if( l_chg ) call self%fs_write_table
            return
        endif
        if( any([(self%fs(j)%st == F_LEVEL, j=1,size(self%fs))]) )then
            key = [(-self%fs(j)%repro, j=1,size(self%fs))]
            ord = [(j, j=1,size(self%fs))]
            call hpsort(key, ord)
            do b = 1, size(ord)
                j = ord(b)
                if( self%fs(j)%st /= F_LEVEL ) cycle
                ev    = self%fs_eval_rows(self%fs(j)%level)
                l_dup = .false.
                do a = 1, size(self%fs)
                    if( self%fs(a)%st /= F_KEEP .or. self%fs(a)%level /= self%fs(j)%level ) cycle
                    r = self%fs_rho(self%fs_rd(a), self%fs_rd(j), ev, .false.)
                    if( r < FGATE ) cycle
                    self%fs(j)%st     = F_DUP
                    self%fs(j)%target = self%fs(a)%name
                    self%fs(j)%rho    = r
                    self%fs(j)%dice   = mask_dice(site_mask(self, j), site_mask(self, a))
                    write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//' duplicates site '//self%fs(a)%name//&
                        &' at '//self%fs(j)%level//': Spearman '//rstr(r)//', mask Dice '//rstr(self%fs(j)%dice)
                    l_dup = .true.
                    exit
                end do
                if( .not. l_dup ) self%fs(j)%st = F_KEEP
            end do
            l_chg = .true.
        endif
        do j = 1, size(self%fs)
            if( any(self%fs(j)%st == [F_KEEP, F_REC, F_REF]) )then
                call self%fs_refine(j, nrun, l_chg)
            else if( self%fs(j)%st == F_READY .or. (self%fs(j)%st == F_DONE .and. &
                &self%dir_status(self%fs_zd(j)) == -2) )then
                if( nrun >= max(1, self%p%nodepar) ) cycle
                call self%fs_size(j)
                nrun  = nrun + 1
                l_chg = .true.
            else if( self%fs(j)%st == F_SIZE )then
                code = self%dir_status(self%fs_zd(j))
                if( code == -1 ) cycle
                if( code == -2 )then
                    if( nrun >= max(1, self%p%nodepar) ) cycle
                    call self%fs_size(j)
                    nrun = nrun + 1
                    cycle
                endif
                if( code /= 0 ) write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': sizing FAILED, its '//&
                    &'labels and maps keep the top-N* sets of its ranking'
                l_chg = .true.
                if( self%dir_status(self%fs_md(j)) == 0 )then
                    self%fs(j)%st = F_DONE
                else
                    if( nrun >= max(1, self%p%nodepar) ) cycle
                    call self%fs_maps(j)
                    nrun = nrun + 1
                endif
            else if( self%fs(j)%st == F_MAPS )then
                sd   = self%fs_md(j)
                code = self%dir_status(sd)
                if( code == -1 ) cycle
                if( code == -2 )then
                    if( nrun >= max(1, self%p%nodepar) ) cycle
                    call self%fs_maps(j)
                    nrun = nrun + 1
                    cycle
                endif
                if( code /= 0 ) write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': pole maps FAILED'
                self%fs(j)%st = F_DONE
                l_chg = .true.
            endif
        end do
        if( l_chg ) call self%fs_write_table
        if( any([(any(self%fs(j)%st == [F_KEEP, F_REC, F_REF, F_READY, F_MAPS, F_SIZE]), j=1,size(self%fs))]) ) return
        if( .not. file_exists(string(self%root//'/'//FDIR//'/sites.txt')) .or. l_chg ) call self%fs_output
        factorial = .true.
    end function factorial

    subroutine fs_rank( self, sd, id, mask, startf, srcproj, nfit, l_size )
        class(hetero_tree), intent(inout) :: self
        character(len=*),   intent(in)    :: sd, id, mask, startf, srcproj
        integer,            intent(in)    :: nfit
        logical,            intent(in)    :: l_size
        type(cmdline)        :: cl
        integer, allocatable :: rows(:), jk(:), st(:)
        character(len=:), allocatable :: nd
        nd = self%root//'/node_'//id
        call simple_mkdir(string(sd))
        if( len(srcproj) > 0 )then
            call simple_copy_file(string(srcproj), string(sd//'/s.simple'))
            st = project_states(sd//'/s.simple')
        else
            rows = read_ints(nd//'/rows.txt')
            allocate(st(self%nrows), source=0)
            st(rows) = 1
            jk = junk_rows(nd)
            if( size(jk) > 0 ) st(jk) = 0
            call write_states(nd//'/'//NODE_PROJ, sd//'/s.simple', st)
        endif
        call cl%set('prg',         'site_rank')
        call cl%set('projfile',    's.simple')
        if( len(startf) > 0 )then
            call cl%set('infile2',   startf)
        else
            call cl%set('infile',    nd//'/fflex/flex_pca_coordinates.txt')
            call cl%set('site_dirs', trim(self%p%site_dirs))
        endif
        call cl%set('nsample',     nfit)
        if( l_size ) call cl%set('site_size', 'yes')
        call cl%set('pole_merge',  self%p%pole_merge)
        if( self%p%pooled_poles .eq. 'yes' ) call cl%set('pooled_poles', 'yes')
        call cl%set('junk_screen', trim(self%p%junk_screen))
        call cl%set('pcg_mskfile', mask)
        call cl%set('pgrp',        trim(self%p%pgrp))
        call cl%set('mskdiam',     self%p%mskdiam)
        call cl%set('nthr',        max(self%p%nthr, min(SITE_NTHR, ceiling(real(count(st > 0)) / real(PTCLS_PER_PART)))))
        call self%submit_in(sd, cl)
    end subroutine fs_rank

    subroutine fs_climb( self, j )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        integer, allocatable :: rows(:), v(:)
        real,    allocatable :: pp(:,:)
        logical, allocatable :: tp(:,:)
        character(len=:), allocatable :: sd, p
        integer :: a, funit, ios
        p  = parent_of(self%fs(j)%level)
        sd = self%fs_cd(j)
        call site_poles(self%fs_rd(j), self%nrows, rows, pp, tp)
        call simple_mkdir(string(sd))
        call fopen(funit, string(sd//'/starts.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//sd//'/starts.txt', ios)
        do a = 1, size(pp,2)
            allocate(v(self%nrows), source=0)
            v(pack(rows, tp(:,a))) = 1
            call write_ints(sd//'/start_p'//int2str(a)//'.txt', v)
            write(funit,'(A)') 'p'//int2str(a)//' '//sd//'/start_p'//int2str(a)//'.txt'
            deallocate(v)
        end do
        call fclose(funit)
        call self%fs_rank(sd, p, site_mask(self, j), sd//'/starts.txt', '', SITE_NFIT, .false.)
        write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': pooling run at '//p//' submitted, started from its '//&
            &int2str(size(pp,2))//' pole(s) at '//self%fs(j)%level
    end subroutine fs_climb

    subroutine fs_maps( self, j )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        type(cmdline)        :: cl
        integer, allocatable :: lab(:), pres(:), st(:)
        character(len=:), allocatable :: wd
        integer :: a
        wd   = self%fs_md(j)
        lab  = self%fs_sets(j, 1)
        pres = present_labels(lab)
        if( size(pres) == 0 )then
            self%fs(j)%st = F_DONE
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': no pole particles, no maps'
            return
        endif
        allocate(st(size(lab)), source=0)
        do a = 1, size(pres)
            where( lab == pres(a) ) st = a
        end do
        call simple_mkdir(string(wd))
        call write_states(self%root//'/node_'//self%fs(j)%level//'/stageB/b.simple', wd//'/p.simple', st)
        call cl%set('prg',      'reconstruct3D')
        call cl%set('projfile', 'p.simple')
        call cl%set('nstates',  size(pres))
        call self%set_parts(count(st > 0))
        call common_keys(self, cl)
        call cl%set('postprocess', 'no')
        call self%submit_in(wd, cl)
        self%fs(j)%st = F_MAPS
        write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': maps of its '//int2str(size(pres))//' pole(s) at '//&
            &self%fs(j)%level//' submitted'
    end subroutine fs_maps

    subroutine fs_output( self )
        class(hetero_tree), intent(inout) :: self
        integer, allocatable :: lab(:), rows(:), pres(:), epm(:)
        real,    allocatable :: pp(:,:), z(:)
        logical, allocatable :: tp(:,:)
        character(len=:), allocatable :: fd, sd, md, line, nm
        character(len=64) :: nums
        integer :: j, a, i, fsit, ftab, fpost, ios, np, nkeep
        fd = self%root//'/'//FDIR
        call simple_mkdir(string(fd))
        call simple_mkdir(string(fd//'/labels'))
        call simple_mkdir(string(fd//'/post'))
        call simple_mkdir(string(fd//'/maps'))
        call simple_mkdir(string(fd//'/maps/halves'))
        call fopen(fsit, string(fd//'/sites.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//fd//'/sites.txt', ios)
        call fopen(ftab, string(fd//'/table.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//fd//'/table.txt', ios)
        write(fsit,'(A)') '# site home candidate discovery_repro state level joins Spearman mask_Dice mask refined | kept sites: '//&
            &'poles | particles per pole | halves z per pole (the weakest accept test of its endpoints) | hard share per pole '//&
            &'(posterior >= 0.5 over the level''s non-junk particles) | equivalent pure members per pole (its own top-EPM set; '//&
            &'far fewer than its particles marks a pole labelled by the complement)'
        write(ftab,'(A)') '# site level labels (0 neither, j pole j) posteriors (one column per pole, -1 outside the level) rows'
        call self%tree_junk
        nkeep = 0
        do j = 1, size(self%fs)
            nm = self%fs(j)%name
            write(nums,'(3(1X,F6.3))') self%fs(j)%repro, self%fs(j)%rho, self%fs(j)%dice
            md   = site_mask(self, j)
            if( self%fs(j)%ref == 1 ) md = self%root//'/node_'//self%fs(j)%level//'/fsite_'//nm//'/refine_mask.mrc'
            if( self%fs(j)%ref == 2 ) md = self%root//'/node_'//self%fs(j)%home//'/fsite_'//nm//'/pre_mask.mrc'
            line = nm//' '//self%fs(j)%home//' '//int2str(self%fs(j)%k)//' '//word(nums, 1)//' '//&
                &trim(FST_NAMES(self%fs(j)%st))//' '//self%fs(j)%level//' '//self%fs(j)%target//' '//word(nums, 2)//' '//&
                &word(nums, 3)//' '//md//' '//trim(merge('yes  ', merge('first', 'no   ', self%fs(j)%ref == 2), self%fs(j)%ref == 1))
            if( self%fs(j)%st == F_DONE )then
                nkeep = nkeep + 1
                sd    = self%fs_fd(j)
                lab   = self%fs_labels(j)
                call site_poles(sd, self%nrows, rows, pp, tp)
                np = size(pp,2)
                z  = pole_z(sd, np)
                call write_ints(fd//'/labels/'//nm//'.txt', lab)
                call fopen(fpost, string(fd//'/post/'//nm//'.txt'), status='replace', action='write', iostat=ios)
                call fileiochk('hetero_tree: cannot write '//fd//'/post/'//nm//'.txt', ios)
                if( allocated(pres) ) deallocate(pres)
                allocate(pres(self%nrows), source=0)
                pres(rows) = [(i, i=1,size(rows))]
                do i = 1, self%nrows
                    if( pres(i) == 0 .or. self%tjunk(i) )then
                        write(fpost,'(*(F8.4,:,1X))') [(-1., a=1,np)]
                    else
                        write(fpost,'(*(F8.5,:,1X))') pp(pres(i),:)
                    endif
                end do
                call fclose(fpost)
                deallocate(pres)
                write(ftab,'(A)') nm//' '//self%fs(j)%level//' '//fd//'/labels/'//nm//'.txt '//fd//'/post/'//nm//'.txt '//&
                    &self%root//'/node_'//self%fs(j)%level//'/rows.txt'
                line = line//' | '//int2str(np)//' |'
                do a = 1, np
                    line = line//' '//int2str(count(lab == a))
                end do
                line = line//' |'
                do a = 1, np
                    line = line//' '//rstr(z(a))
                end do
                line = line//' |'
                do a = 1, np
                    line = line//' '//rstr(real(count(pp(:,a) >= 0.5 .and. .not. self%tjunk(rows))) / &
                        &real(max(1, count(.not. self%tjunk(rows)))))
                end do
                line = line//' |'
                do a = 1, np
                    epm = read_ints(self%fs_zd(j)//'/site_epm_p'//int2str(a)//'.txt')
                    if( size(epm) == self%nrows )then
                        line = line//' '//int2str(count(epm(rows) > 0 .and. .not. self%tjunk(rows)))
                    else
                        line = line//' -'
                    endif
                end do
                md   = self%fs_md(j)
                pres = present_labels(lab)
                do a = 1, size(pres)
                    if( .not. (file_exists(string(md//'/'//hv(a, 'even'))) .and. file_exists(string(md//'/'//hv(a, 'odd')))) ) cycle
                    call pole_map(md//'/'//hv(a, 'even'), md//'/'//hv(a, 'odd'), fd//'/maps/'//nm//'_pole'//int2str(pres(a))//&
                        &'.mrc', fd//'/maps/halves/'//nm//'_pole'//int2str(pres(a)), self%p%mskdiam)
                end do
            endif
            write(fsit,'(A)') line
        end do
        call fclose(fsit)
        call fclose(ftab)
        write(logfhandle,'(A)') '>>> HETERO_TREE per-site layer: '//int2str(nkeep)//' sites kept of '//int2str(size(self%fs))//&
            &' candidates ('//int2str(count_st(F_COPY))//' joined an ancestor''s site, '//int2str(count_st(F_DUP))//&
            &' duplicated a kept site, '//int2str(count_st(F_NOPOLE))//' without a pole, '//int2str(count_st(F_FAIL))//&
            &' failed): '//fd//'/sites.txt, table.txt, maps'

    contains

        integer function count_st( st )
            integer, intent(in) :: st
            integer :: jj
            count_st = 0
            do jj = 1, size(self%fs)
                if( self%fs(jj)%st == st ) count_st = count_st + 1
            end do
        end function count_st

    end subroutine fs_output

    function fs_hd( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        if( self%fs(j)%ref == 2 )then
            d = self%root//'/node_'//self%fs(j)%home//'/fsite_'//self%fs(j)%name//'/pre'
        else if( self%fs(j)%k == 1 )then
            d = self%root//'/node_'//self%fs(j)%home//'/site'
        else
            d = self%root//'/node_'//self%fs(j)%home//'/fsite_'//self%fs(j)%name
        endif
    end function fs_hd

    function fs_rd( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        if( self%fs(j)%level == self%fs(j)%home )then
            d = self%fs_hd(j)
        else
            d = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name
        endif
    end function fs_rd

    function fs_cd( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        d = self%root//'/node_'//parent_of(self%fs(j)%level)//'/fsite_'//self%fs(j)%name
    end function fs_cd

    function fs_fd( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        if( self%fs(j)%ref == 1 )then
            d = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name//'/refine'
        else
            d = self%fs_rd(j)
        endif
    end function fs_fd

    subroutine fs_refine( self, j, nrun, l_chg )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        integer,            intent(inout) :: nrun
        logical,            intent(inout) :: l_chg
        type(cmdline)             :: cl
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: lab(:), st(:), rows(:), pres(:)
        real(dp),     allocatable :: w(:)
        character(len=:), allocatable :: wd, sf
        integer :: code, s
        logical :: l_ok
        wd = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name
        select case(self%fs(j)%st)
            case(F_KEEP)
                if( self%fs(j)%ref == 2 )then
                    self%fs(j)%st = F_READY
                    l_chg = .true.
                    return
                endif
                if( count_poles(self%fs_rd(j)//'/site_poles.txt') == 0 )then
                    self%fs(j)%st = F_READY
                    l_chg = .true.
                    write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': no pole passes the halves test at '//&
                        &self%fs(j)%level//', not refined'
                    return
                endif
                if( nrun >= max(1, self%p%nodepar) ) return
                call submit_classes
            case(F_REC)
                code = self%dir_status(wd//'/refine_rec')
                if( code == -1 ) return
                if( code == -2 )then
                    if( nrun < max(1, self%p%nodepar) ) call submit_classes
                    return
                endif
                if( code /= 0 )then
                    call finish_ref(.false., 'its class maps FAILED')
                    return
                endif
                if( nrun >= max(1, self%p%nodepar) ) return
                st   = project_states(wd//'/refine_rec/c.simple')
                pres = with_maps(wd//'/refine_rec', present_labels(st))
                if( size(pres) < 2 )then
                    call finish_ref(.false., 'its class maps FAILED')
                    return
                endif
                allocate(evens(size(pres)), odds(size(pres)), w(size(pres)))
                do s = 1, size(pres)
                    evens(s) = wd//'/refine_rec/'//hv(pres(s), 'even')
                    odds(s)  = wd//'/refine_rec/'//hv(pres(s), 'odd')
                    w(s)     = real(count(st == pres(s)), dp)
                end do
                call falloff_mask(self, evens, odds, w, wd//'/refine_mask.mrc')
                sf = ''
                if( self%fs(j)%level /= self%fs(j)%home )then
                    sf = wd//'/refine/starts.txt'
                    call simple_mkdir(string(wd//'/refine'))
                    call simple_copy_file(string(self%fs_rd(j)//'/starts.txt'), string(sf))
                endif
                call self%fs_rank(wd//'/refine', self%fs(j)%level, wd//'/refine_mask.mrc', sf, self%fs_rd(j)//'/s.simple', &
                    &merge(0, SITE_NFIT, self%fs(j)%k == 1 .and. self%fs(j)%level == self%fs(j)%home), .false.)
                nrun = nrun + 1
                self%fs(j)%st = F_REF
                l_chg = .true.
                write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': ranking with the mask from its '//&
                    &int2str(size(pres))//' class maps submitted'
            case(F_REF)
                code = self%dir_status(wd//'/refine')
                if( code == -1 ) return
                if( code == -2 )then
                    self%fs(j)%st = F_REC
                    return
                endif
                if( code /= 0 )then
                    call finish_ref(.false., 'its ranking FAILED')
                    return
                endif
                l_ok = self%fs_guard(j, self%fs_rd(j), wd//'/refine', .false.)
                call finish_ref(l_ok, trim(merge('kept    ', 'rejected', l_ok)))
        end select

    contains

        subroutine submit_classes
            lab  = self%fs_sets(j, 0)
            rows = read_ints(self%root//'/node_'//self%fs(j)%level//'/rows.txt')
            allocate(st(self%nrows), source=0)
            st(rows) = maxval([0, lab]) + 1
            where( lab > 0 ) st = lab
            where( self%tjunk ) st = 0
            pres = present_labels(st)
            lab  = st
            do s = 1, size(pres)
                where( lab == pres(s) ) st = s
            end do
            call simple_mkdir(string(wd//'/refine_rec'))
            call write_states(self%root//'/node_'//self%fs(j)%level//'/'//NODE_PROJ, wd//'/refine_rec/c.simple', st)
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', 'c.simple')
            call cl%set('nstates',  size(pres))
            call self%set_parts(count(st > 0))
            call common_keys(self, cl)
            call halves_only_keys(cl)
            call self%submit_in(wd//'/refine_rec', cl)
            nrun = nrun + 1
            self%fs(j)%st = F_REC
            l_chg = .true.
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': maps of its '//int2str(size(pres))//&
                &' classes for the mask refinement submitted'
        end subroutine submit_classes

        subroutine finish_ref( l_keep, what )
            logical,          intent(in) :: l_keep
            character(len=*), intent(in) :: what
            self%fs(j)%ref = merge(1, -1, l_keep)
            self%fs(j)%st  = F_READY
            l_chg = .true.
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': refined mask '//what
        end subroutine finish_ref

    end subroutine fs_refine

    logical function fs_guard( self, j, od, nd, l_halves )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        character(len=*),   intent(in)    :: od, nd
        logical,            intent(in)    :: l_halves
        type(string), allocatable :: ends(:)
        integer,      allocatable :: pole(:)
        real,         allocatable :: v0(:), v1(:)
        logical,      allocatable :: ok(:)
        character(len=4096) :: ln
        character(len=:), allocatable :: e
        real    :: sc, mc, n0, n1, r0, r1, se, sh0, sh1
        integer :: k, fu, io, b(2), bc
        logical :: l_pass
        fs_guard = npoles_of(nd) > 0
        if( .not. fs_guard )then
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': the refined ranking accepts no pole'
            return
        endif
        call fopen(fu, string(nd//'/site_neff.txt'), status='old', action='read', iostat=io)
        call fileiochk('hetero_tree: cannot read '//nd//'/site_neff.txt', io)
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( ln(1:1) /= '#' ) exit
        end do
        call fclose(fu)
        read(ln,*) bc, sc, mc, b(1), b(2), n1
        n0 = mask_neff(string(site_mask(self, j)), self%smpd, bc, sc, mc, b(1), b(2))
        call pole_ends(od, ends, pole, l_halves)
        call self%tree_junk
        do k = 1, size(ends)
            e   = ends(k)%to_char()
            r0  = curve_peak(od//'/site_curve_'//e//'.txt')
            r1  = curve_peak(nd//'/site_curve_'//e//'.txt')
            se  = sqrt(((1. - r0**2) / sqrt(max(n0 - 3., 1.)))**2 + ((1. - r1**2) / sqrt(max(n1 - 3., 1.)))**2)
            v0  = read_reals(od//'/site_post_'//e//'.txt')
            v1  = read_reals(nd//'/site_post_'//e//'.txt')
            if( size(v0) /= self%nrows .or. size(v1) /= self%nrows )then
                fs_guard = .false.
                cycle
            endif
            ok  = v0 >= 0. .and. v1 >= 0. .and. .not. self%tjunk
            sh0 = real(count(ok .and. v0 >= 0.5)) / real(max(1, count(ok)))
            sh1 = real(count(ok .and. v1 >= 0.5)) / real(max(1, count(ok)))
            l_pass = r0 - r1 <= se .and. ((sh0 < 0.5) .eqv. (sh1 < 0.5))
            if( .not. l_pass ) fs_guard = .false.
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//' endpoint '//e//': hold-out peak '//rstr(r0)//&
                &' -> '//rstr(r1)//' (combined SE '//rstr(se)//', n_eff '//rstr(n0)//' -> '//rstr(n1)//'), class share '//&
                &rstr(sh0)//' -> '//rstr(sh1)//trim(merge(' passes', ' fails ', l_pass))
        end do
    end function fs_guard

    subroutine fs_prerefine( self, j, nrun, l_chg )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        integer,            intent(inout) :: nrun
        logical,            intent(inout) :: l_chg
        type(cmdline)             :: cl
        type(string), allocatable :: evens(:), odds(:), ends(:)
        integer,      allocatable :: st(:), pres(:), pole(:), top(:), lab(:)
        real,         allocatable :: v(:), pm(:), best(:)
        real(dp),     allocatable :: w(:)
        character(len=:), allocatable :: wd, hd
        integer :: code, s, e, np
        wd = self%root//'/node_'//self%fs(j)%home//'/fsite_'//self%fs(j)%name
        hd = self%fs_hd(j)
        select case(self%fs(j)%st)
            case(F_TEST)
                call submit_classes
            case(F_PREC)
                code = self%dir_status(wd//'/pre_rec')
                if( code == -1 ) return
                if( code == -2 )then
                    if( nrun < max(1, self%p%nodepar) ) call submit_classes
                    return
                endif
                if( code /= 0 )then
                    call reject('its class maps FAILED')
                    return
                endif
                if( nrun >= max(1, self%p%nodepar) ) return
                st   = project_states(wd//'/pre_rec/c.simple')
                pres = with_maps(wd//'/pre_rec', present_labels(st))
                if( size(pres) < 2 )then
                    call reject('its class maps FAILED')
                    return
                endif
                allocate(evens(size(pres)), odds(size(pres)), w(size(pres)))
                do s = 1, size(pres)
                    evens(s) = wd//'/pre_rec/'//hv(pres(s), 'even')
                    odds(s)  = wd//'/pre_rec/'//hv(pres(s), 'odd')
                    w(s)     = real(count(st == pres(s)), dp)
                end do
                call falloff_mask(self, evens, odds, w, wd//'/pre_mask.mrc')
                call self%fs_rank(wd//'/pre', self%fs(j)%home, wd//'/pre_mask.mrc', '', hd//'/s.simple', &
                    &merge(0, SITE_NFIT, self%fs(j)%k == 1), .false.)
                nrun = nrun + 1
                self%fs(j)%st = F_PREF
                l_chg = .true.
                write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': fails the halves test narrowly, ranking with '//&
                    &'the mask from its '//int2str(size(pres))//' class maps submitted'
            case(F_PREF)
                code = self%dir_status(wd//'/pre')
                if( code == -1 ) return
                if( code == -2 )then
                    self%fs(j)%st = F_PREC
                    return
                endif
                if( code /= 0 )then
                    call reject('its ranking FAILED')
                    return
                endif
                if( .not. self%fs_guard(j, hd, wd//'/pre', .true.) )then
                    call reject('its refined ranking fails the halves test or the guard')
                    return
                endif
                self%fs(j)%ref = 2
                self%fs(j)%st  = F_TEST
                l_chg = .true.
                write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': refined first, its refined ranking passes '//&
                    &'the halves test and the guard and becomes its ranking'
        end select

    contains

        subroutine submit_classes
            call pole_ends(hd, ends, pole, .true.)
            np = maxval([0, pole])
            if( np == 0 )then
                call reject('no endpoint')
                return
            endif
            v = read_reals(hd//'/site_post_'//ends(1)%to_char()//'.txt')
            allocate(lab(self%nrows), source=0)
            allocate(best(self%nrows), source=-1.)
            where( v >= 0. ) lab = np + 1
            do s = 1, np
                allocate(pm(self%nrows), source=0.)
                allocate(top(self%nrows), source=0)
                do e = 1, size(ends)
                    if( pole(e) /= s ) cycle
                    pm = pm + read_reals(hd//'/site_post_'//ends(e)%to_char()//'.txt') / real(count(pole == s))
                    top = max(top, read_ints(hd//'/site_top_'//ends(e)%to_char()//'.txt'))
                end do
                where( top > 0 .and. v >= 0. .and. pm > best )
                    lab  = s
                    best = pm
                end where
                deallocate(pm, top)
            end do
            call self%tree_junk
            where( self%tjunk ) lab = 0
            pres = present_labels(lab)
            allocate(st(self%nrows), source=0)
            do s = 1, size(pres)
                where( lab == pres(s) ) st = s
            end do
            call simple_mkdir(string(wd//'/pre_rec'))
            call write_states(self%root//'/node_'//self%fs(j)%home//'/'//NODE_PROJ, wd//'/pre_rec/c.simple', st)
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', 'c.simple')
            call cl%set('nstates',  size(pres))
            call self%set_parts(count(st > 0))
            call common_keys(self, cl)
            call halves_only_keys(cl)
            call self%submit_in(wd//'/pre_rec', cl)
            nrun = nrun + 1
            self%fs(j)%st = F_PREC
            l_chg = .true.
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': fails the halves test narrowly (best endpoint '//&
                &rstr(best_minz(hd))//' SE), maps of its '//int2str(size(pres))//' classes for refining first submitted'
        end subroutine submit_classes

        subroutine reject( why )
            character(len=*), intent(in) :: why
            self%fs(j)%st = F_NOPOLE
            l_chg = .true.
            write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': not a site after refining first, '//why
        end subroutine reject

    end subroutine fs_prerefine

    function fs_zd( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        d = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name//'/size'
    end function fs_zd

    function fs_mask( self, j ) result( f )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: f
        select case(self%fs(j)%ref)
            case(1)
                f = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name//'/refine_mask.mrc'
            case(2)
                f = self%root//'/node_'//self%fs(j)%home//'/fsite_'//self%fs(j)%name//'/pre_mask.mrc'
            case default
                f = site_mask(self, j)
        end select
    end function fs_mask

    subroutine fs_size( self, j )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        integer, allocatable :: rows(:)
        real,    allocatable :: pp(:,:), v(:)
        logical, allocatable :: tp(:,:)
        character(len=:), allocatable :: sd
        integer :: a, funit, ios
        sd = self%fs_zd(j)
        call site_poles(self%fs_fd(j), self%nrows, rows, pp, tp)
        call simple_mkdir(string(sd))
        call fopen(funit, string(sd//'/rankings.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//sd//'/rankings.txt', ios)
        allocate(v(self%nrows))
        do a = 1, size(pp,2)
            v       = -1.
            v(rows) = pp(:,a)
            call write_reals(sd//'/pole_p'//int2str(a)//'.txt', v)
            write(funit,'(A)') 'p'//int2str(a)//' '//sd//'/pole_p'//int2str(a)//'.txt'
        end do
        call fclose(funit)
        call self%fs_rank(sd, self%fs(j)%level, self%fs_mask(j), sd//'/rankings.txt', self%fs_fd(j)//'/s.simple', SITE_NFIT, .true.)
        self%fs(j)%st = F_SIZE
        write(logfhandle,'(A)') '>>> HETERO_TREE site '//self%fs(j)%name//': sizing of its '//int2str(size(pp,2))//&
            &' pole ranking(s) at '//self%fs(j)%level//' submitted'
    end subroutine fs_size

    function fs_md( self, j ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: d
        d = self%root//'/node_'//self%fs(j)%level//'/fsite_'//self%fs(j)%name//'/maps'
    end function fs_md

    integer function dir_status( self, sd )
        class(hetero_tree), intent(in) :: self
        character(len=*),   intent(in) :: sd
        logical :: err
        integer :: ios
        call simple_sync_dir(sd, status=ios)
        if( file_exists(string(sd//'/'//STEP_EXIT)) )then
            call read_exit_code(string(sd//'/'//STEP_EXIT), dir_status, err)
            if( err ) dir_status = -1
        else if( file_exists(string(sd//'/'//POOL_PENDING)) .and. .not. self%pending(sd) )then
            THROW_HARD('unfinished worker task from an earlier controller; resolve before restart: '//sd)
        else if( file_exists(string(sd//'/'//STEP_EXIT//JOB_INFO_EXT)) .or. self%pending(sd) )then
            dir_status = -1
        else
            dir_status = -2
        endif
    end function dir_status

    integer function batch_status( self, bf, lfail )
        class(hetero_tree), intent(in) :: self
        character(len=*),   intent(in) :: bf
        logical, optional,  intent(in) :: lfail
        character(len=4096) :: ln
        integer :: fu, io, c, cf
        batch_status = 0
        cf = 0
        if( .not. file_exists(string(bf)) )then
            batch_status = -2
            return
        endif
        call fopen(fu, string(bf), status='old', action='read', iostat=io)
        if( io /= 0 )then
            batch_status = -1
            return
        endif
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) == 0 ) cycle
            c = self%dir_status(trim(ln))
            if( c == -2 )then
                batch_status = -2
                exit
            endif
            if( c == -1 ) batch_status = -1
            if( c > 0 .and. cf == 0 ) cf = c
        end do
        call fclose(fu)
        if( batch_status /= 0 .or. .not. present(lfail) ) return
        if( lfail ) batch_status = cf
    end function batch_status

    subroutine off_work( self, noff )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(out)   :: noff
        integer :: i
        noff = 0
        if( self%l_multi .or. .not. self%l_motion ) return
        do i = 1, size(self%nodes)
            if( self%nodes(i)%step == S_SUBTREE ) cycle
            if( .not. off_open(self%ndir(i)) ) cycle
            noff = noff + 1
            if( self%batch_status(self%ndir(i)//'/motion/offbatch.txt') == -1 ) cycle
            call self%advance(i, S_OFFP)
        end do
    end subroutine off_work

    logical function anc_wait( self, i )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i
        integer :: lp
        anc_wait = .false.
        if( .not. self%l_inherit ) return
        if( file_exists(string(self%ndir(i)//'/motion/inherited.txt')) ) return
        do lp = 1, len(self%nodes(i)%id) - 1
            if( off_open(self%root//'/node_'//self%nodes(i)%id(1:lp)) ) anc_wait = .true.
        end do
    end function anc_wait

    logical function floor_ready( self )
        class(hetero_tree), intent(in) :: self
        character(len=2), parameter :: FN(3) = ['r ', 'r1', 'r2']
        integer :: f, j, jj
        floor_ready = .true.
        if( .not. self%l_gate .or. self%p%motion_floor > 0. ) return
        if( file_exists(string(self%root//'/'//FLOOR_FILE)) ) return
        do f = 1, 3
            j = 0
            do jj = 1, size(self%nodes)
                if( self%nodes(jj)%id == trim(FN(f)) ) j = jj
            end do
            if( j == 0 ) cycle
            if( file_exists(string(self%root//'/node_'//trim(FN(f))//'/motion/null/gate.txt')) ) cycle
            if( self%nodes(j)%step == S_DONE .or. self%nodes(j)%step == S_FAILED ) cycle
            if( self%nodes(j)%step >= S_SEEDS .and. self%nodes(j)%step <= S_SPLIT ) cycle
            floor_ready = .false.
            return
        end do
    end function floor_ready

    real function floor_value( self )
        class(hetero_tree), intent(in) :: self
        character(len=2), parameter :: FN(3) = ['r ', 'r1', 'r2']
        real, allocatable :: v(:)
        character(len=4096) :: ln
        character(len=256)  :: buf
        character(len=:), allocatable :: ff, gf, nodes, vals
        real    :: x, med, mads
        integer :: f, fu, io, nkept, ndrop, nempty, nperm
        floor_value = self%p%motion_floor
        if( floor_value > 0. ) return
        ff = self%root//'/'//FLOOR_FILE
        if( file_exists(string(ff)) )then
            buf = first_line(ff, '-1.0e30')
            read(buf,*,iostat=io) floor_value
            if( io /= 0 ) floor_value = -huge(1.)
            if( floor_value < -1.e29 ) floor_value = -huge(1.)
            return
        endif
        allocate(v(0))
        nodes = ''
        vals  = ''
        do f = 1, 3
            gf = self%root//'/node_'//trim(FN(f))//'/motion/null/gate.txt'
            if( .not. file_exists(string(gf)) ) cycle
            nkept  = 0
            ndrop  = 0
            nempty = 0
            nperm  = 0
            call fopen(fu, string(gf), status='old', action='read', iostat=io)
            if( io /= 0 ) cycle
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                if( nwords(ln) < 4 ) cycle
                if( word(ln, 4) == '0' )then
                    ndrop = ndrop + 1
                    cycle
                endif
                if( word(ln, 4) /= '1' .and. word(ln, 4) /= 'e' ) cycle
                buf = word(ln, 2)
                read(buf,*,iostat=io) x
                if( io /= 0 ) cycle
                v     = [v, abs(x)]
                vals  = vals//' '//word(ln, 2)
                if( word(ln, 4) == 'e' )then
                    nempty = nempty + 1
                else
                    nkept = nkept + 1
                endif
                if( index(ln, 'perm:') == 1 ) nperm = nperm + 1
            end do
            call fclose(fu)
            if( nperm == nkept + nempty )then
                nodes = nodes//' '//trim(FN(f))//' ('//int2str(nperm)//' permutation splits)'
            else
                nodes = nodes//' '//trim(FN(f))//' ('//int2str(nkept - nperm)//' fake-direction refits kept, '//int2str(ndrop)//&
                    &' dropped, '//int2str(nempty)//' empty, '//int2str(nperm)//' permutation splits)'
            endif
        end do
        call fopen(fu, string(ff), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write '//ff, io)
        if( size(v) == 0 )then
            floor_value = -huge(1.)
            write(fu,'(A)') '-1.0e30'
            write(fu,'(A)') 'no null split scored: every gated axis passes'
            write(logfhandle,'(A)') '>>> HETERO_TREE WARNING gate floor: no null split scored, every gated axis passes'
        else
            med  = quant_of(v, 0.5)
            mads = 1.4826 * quant_of(abs(v - med), 0.5)
            floor_value = med + GATE_NSD * mads
            write(fu,'(A)') rstr(floor_value)
            write(fu,'(A)') 'median |corr| '//rstr(med)//' + '//rstr(GATE_NSD)//' MAD-sd '//rstr(mads)//' of '//&
                &int2str(size(v))//' null splits on'//nodes//':'//vals
            write(logfhandle,'(A)') '>>> HETERO_TREE gate floor '//rstr(floor_value)//': median |corr| '//rstr(med)//' + '//&
                &rstr(GATE_NSD)//' MAD-sd '//rstr(mads)//' of '//int2str(size(v))//' null splits on'//nodes
        endif
        call fclose(fu)
    end function floor_value

    integer function fs_nrun( self )
        class(hetero_tree), intent(in) :: self
        integer :: j
        fs_nrun = 0
        do j = 1, size(self%fs)
            select case(self%fs(j)%st)
                case(F_RUN)
                    if( self%fs(j)%k > 1 .and. self%dir_status(self%fs_hd(j)) == -1 ) fs_nrun = fs_nrun + 1
                case(F_UP)
                    if( self%dir_status(self%fs_cd(j)) == -1 ) fs_nrun = fs_nrun + 1
                case(F_MAPS)
                    if( self%dir_status(self%fs_md(j)) == -1 ) fs_nrun = fs_nrun + 1
            end select
        end do
    end function fs_nrun

    real function fs_rho( self, sda, sdb, ev, l_first )
        class(hetero_tree), intent(in) :: self
        character(len=*),   intent(in) :: sda, sdb
        integer,            intent(in) :: ev(:)
        logical,            intent(in) :: l_first
        integer,  allocatable :: ra(:), rb(:), ia(:), ib(:), sel(:)
        real,     allocatable :: pa(:,:), pb(:,:)
        logical,  allocatable :: ta(:,:), tb(:,:)
        real(dp), allocatable :: ka(:,:), kb(:,:)
        integer :: a, b, i, n
        fs_rho = 0.
        call site_poles(sda, self%nrows, ra, pa, ta)
        call site_poles(sdb, self%nrows, rb, pb, tb)
        if( size(pa,2) == 0 .or. size(pb,2) == 0 ) return
        allocate(ia(self%nrows), ib(self%nrows), source=0)
        ia(ra) = [(i, i=1,size(ra))]
        ib(rb) = [(i, i=1,size(rb))]
        sel = pack(ev, ia(ev) > 0 .and. ib(ev) > 0)
        n   = size(sel)
        if( n < 3 ) return
        allocate(ka(n,size(pa,2)), kb(n,size(pb,2)))
        do a = 1, size(pa,2)
            ka(:,a) = avg_ranks(pa(ia(sel),a))
            ka(:,a) = ka(:,a) - sum(ka(:,a)) / real(n,dp)
        end do
        do b = 1, size(pb,2)
            kb(:,b) = avg_ranks(pb(ib(sel),b))
            kb(:,b) = kb(:,b) - sum(kb(:,b)) / real(n,dp)
        end do
        do a = 1, size(pa,2)
            do b = 1, size(pb,2)
                if( l_first .and. (a > 1 .or. b > 1) ) cycle
                fs_rho = max(fs_rho, real(abs(sum(ka(:,a) * kb(:,b))) / max(sqrt(sum(ka(:,a)**2) * sum(kb(:,b)**2)), 1.e-30_dp)))
            end do
        end do
    end function fs_rho

    function fs_eval_rows( self, id ) result( ev )
        class(hetero_tree), intent(in) :: self
        character(len=*),   intent(in) :: id
        integer, allocatable :: ev(:), jk(:)
        logical, allocatable :: m(:)
        ev = read_ints(self%root//'/node_'//id//'/rows.txt')
        jk = junk_rows(self%root//'/node_'//id)
        if( size(jk) == 0 ) return
        allocate(m(self%nrows), source=.false.)
        m(jk) = .true.
        ev = pack(ev, .not. m(ev))
    end function fs_eval_rows

    function fs_labels( self, j ) result( lab )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j
        integer, allocatable :: lab(:)
        lab = self%fs_sets(j, 2)
    end function fs_labels

    function fs_sets( self, j, kind ) result( lab )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: j, kind
        integer, allocatable :: lab(:), rows(:), top(:), mi(:)
        real,    allocatable :: pp(:,:)
        logical, allocatable :: tp(:,:), minor(:)
        character(len=:), allocatable :: sz
        integer :: i, a, best
        real    :: pb, fmin
        call site_poles(self%fs_fd(j), self%nrows, rows, pp, tp)
        call self%tree_junk
        sz = self%fs_zd(j)
        if( kind > 0 .and. self%dir_status(sz) == 0 )then
            do a = 1, size(pp,2)
                if( .not. file_exists(string(sz//'/site_'//merge('top', 'epm', kind == 1)//'_p'//int2str(a)//'.txt')) ) cycle
                top     = read_ints(sz//'/site_'//merge('top', 'epm', kind == 1)//'_p'//int2str(a)//'.txt')
                tp(:,a) = top(rows) > 0
            end do
            if( kind == 2 )then
                minor = [(count(pp(:,a) >= 0.5 .and. .not. self%tjunk(rows)) < count(.not. self%tjunk(rows)) - &
                    &count(pp(:,a) >= 0.5 .and. .not. self%tjunk(rows)), a=1,size(pp,2))]
                if( any(minor) .and. .not. all(minor) )then
                    mi   = pack([(i, i=1,size(pp,2))], minor)
                    fmin = real(count(any(tp(:,mi), dim=2) .and. .not. self%tjunk(rows))) / real(max(1, count(.not. self%tjunk(rows))))
                    do a = 1, size(pp,2)
                        if( minor(a) ) cycle
                        if( real(count(pp(:,a) >= 0.5 .and. .not. self%tjunk(rows))) / real(max(1, count(.not. self%tjunk(rows)))) &
                            &+ fmin > 1. ) cycle
                        tp(:,a) = .not. any(tp(:,mi), dim=2) .and. pp(:,a) > maxval(pp(:,mi), dim=2)
                    end do
                endif
            endif
        endif
        allocate(lab(self%nrows), source=0)
        do i = 1, size(rows)
            best = 0
            pb   = -1.
            do a = 1, size(pp,2)
                if( tp(i,a) .and. pp(i,a) > pb )then
                    best = a
                    pb   = pp(i,a)
                endif
            end do
            lab(rows(i)) = best
        end do
        call self%tree_junk
        where( self%tjunk ) lab = 0
    end function fs_sets

    subroutine tree_junk( self )
        class(hetero_tree), intent(inout) :: self
        type(string), allocatable :: lines(:), nds(:)
        integer,      allocatable :: rows(:)
        logical,      allocatable :: inleaf(:)
        character(len=4096) :: line
        integer :: il
        if( allocated(self%tjunk) ) return
        allocate(inleaf(self%nrows), source=.false.)
        if( file_exists(string(self%root//'/'//LEAF_TABLE)) )then
            call leaf_lines(lines, nds)
            do il = 1, size(lines)
                line = lines(il)%to_char()
                if( word(line, 2) == 'junk' ) cycle
                rows = read_ints(word(line, 8))
                inleaf(rows) = .true.
            end do
        endif
        allocate(self%tjunk(self%nrows), source=.false.)
        rows = read_ints(self%ndir(1)//'/rows.txt')
        self%tjunk(rows) = .not. inleaf(rows)
    end subroutine tree_junk

    subroutine fs_read_table( self )
        class(hetero_tree), intent(inout) :: self
        character(len=4096) :: line
        character(len=:), allocatable :: f
        integer :: funit, ios, n, j
        if( allocated(self%fs) ) deallocate(self%fs)
        f = self%root//'/'//SITE_TABLE//'.txt'
        if( .not. file_exists(string(f)) )then
            allocate(self%fs(0))
            return
        endif
        n = nlines(string(f))
        allocate(self%fs(n))
        call fopen(funit, string(f), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//f, ios)
        do j = 1, n
            read(funit,'(A)') line
            self%fs(j)%name   = word(line, 1)
            self%fs(j)%home   = word(line, 2)
            self%fs(j)%k      = str2int(word(line, 3))
            self%fs(j)%st     = str2int(word(line, 4))
            self%fs(j)%level  = word(line, 5)
            self%fs(j)%target = word(line, 6)
            self%fs(j)%repro = str2real(word(line, 7))
            self%fs(j)%rho   = str2real(word(line, 8))
            self%fs(j)%dice  = str2real(word(line, 9))
            if( nwords(line) >= 10 ) self%fs(j)%ref = str2int(word(line, 10))
        end do
        call fclose(funit)
        write(logfhandle,'(A)') '>>> HETERO_TREE per-site layer resumed from '//f//': '//int2str(n)//' sites'
    end subroutine fs_read_table

    subroutine fs_write_table( self )
        class(hetero_tree), intent(in) :: self
        character(len=64) :: nums
        integer :: funit, ios, j
        call fopen(funit, string(self%root//'/'//SITE_TABLE//'.txt.tmp'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//SITE_TABLE//'.txt', ios)
        do j = 1, size(self%fs)
            write(nums,'(3(1X,F7.4))') self%fs(j)%repro, self%fs(j)%rho, self%fs(j)%dice
            write(funit,'(A)') self%fs(j)%name//' '//self%fs(j)%home//' '//int2str(self%fs(j)%k)//' '//&
                &int2str(self%fs(j)%st)//' '//self%fs(j)%level//' '//self%fs(j)%target//trim(nums)//' '//int2str(self%fs(j)%ref)
        end do
        call fclose(funit)
        call simple_rename(string(self%root//'/'//SITE_TABLE//'.txt.tmp'), string(self%root//'/'//SITE_TABLE//'.txt'))
    end subroutine fs_write_table

    integer function step_status( self, i, istep )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i, istep
        character(len=:), allocatable :: f, sdir
        logical :: err
        integer :: ios
        if( istep == S_MBATCH )then
            step_status = self%batch_status(self%ndir(i)//'/motion/batch.txt')
            return
        else if( istep == S_SBATCH )then
            step_status = self%batch_status(self%ndir(i)//'/sites/batch.txt', .true.)
            return
        else if( istep == S_GATEW )then
            step_status = merge(0, -1, self%floor_ready() .and. .not. self%anc_wait(i))
            return
        endif
        sdir = self%step_dir(i, istep)
        f = sdir//'/'//STEP_EXIT
        ! opening the step directory revalidates NFS's cached "absent" answer, so a finished step is seen at once
        call simple_sync_dir(sdir, status=ios)
        if( file_exists(string(f)) )then
            call read_exit_code(string(f), step_status, err)
            if( err ) step_status = -1
        else if( file_exists(string(sdir//'/'//POOL_PENDING)) .and. .not. self%pending(sdir) )then
            THROW_HARD('unfinished worker task from an earlier controller; resolve before restart: '//sdir)
        else if( file_exists(string(f//JOB_INFO_EXT)) .or. self%pending(sdir) )then
            step_status = -1
        else
            step_status = -2
        endif
    end function step_status

    subroutine decide( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        integer, allocatable :: rows(:), lab(:), ch(:), junk(:), states(:), to(:), pop(:), mem(:)
        type(tnode), allocatable :: tmp(:)
        character(len=:), allocatable :: d, name, comp
        character(len=8) :: ev
        logical :: l_split, l_re
        integer :: s, g, ng, n, fg, big, j, ios
        call get_environment_variable('SIMPLE_TREE_REENTER', ev, status=ios)
        l_re   = ios == 0 .and. ev == 'yes'
        d      = self%ndir(i)
        rows   = read_ints(d//'/rows.txt')
        lab    = project_states(d//'/stageB/b.simple')
        ch     = read_ints(d//'/split/split_children.txt')
        junk   = read_ints(d//'/split/split_junk.txt')
        if( size(ch) /= size(lab) ) THROW_HARD('split children do not match the node; decide')
        states = present_labels(lab)
        allocate(pop(max(1,maxval([0,ch]))), source=0)
        do g = 1, size(pop)
            pop(g) = count(ch == g)
        end do
        fg = 0
        if( count(pop > 0) >= 2 .and. all(pack(pop, pop > 0) >= self%p%minpop) .and. &
            &self%nodes(i)%depth < self%p%maxdepth )then
            do g = 1, size(pop)
                if( pop(g) == 0 ) cycle
                comp = ''
                do j = 1, size(states)
                    n = count(ch == g .and. lab == states(j))
                    if( n > 0 ) comp = comp//' '//int2str(states(j))//':'//int2str(n)
                end do
                call add_child([(ch(j) == g, j=1,size(lab))], 'split group '//int2str(g)//' ('//int2str(pop(g))//&
                    &' particles; stage-B states'//comp//')')
            end do
        else
            allocate(to(size(states)), source=0)
            do j = 1, size(states)
                s = states(j)
                if( any(junk == s) )then
                    to(j) = -1
                    cycle
                endif
                if( count(lab == s .and. ch > 0) > 0 ) to(j) = mode_of(pack(ch, lab == s .and. ch > 0))
            end do
            big = 1
            if( any(to > 0) ) big = maxloc(pop, dim=1, mask=[(any(to == g), g=1,size(pop))])
            where( to == 0 ) to = big
            ng      = count([(any(to == g), g=1,size(pop))])
            l_split = ng >= 2
            do g = 1, size(pop)
                if( .not. any(to == g) ) cycle
                mem = pack(states, to == g)
                n   = count([(any(mem == lab(j)), j=1,size(lab))])
                if( l_split .and. n >= self%p%minpop .and. self%nodes(i)%depth < self%p%maxdepth )then
                    call add_child([(any(mem == lab(j)), j=1,size(lab))], 'states '//joined(mem)//' ('//int2str(n)//&
                        &' particles)')
                else if( l_re .and. ng == 1 .and. any([(count(lab == junk(j)) > 0, j=1,size(junk))]) .and. &
                    &n >= self%p%minpop .and. self%nodes(i)%depth < self%p%maxdepth .and. &
                    &.not. reentered_line(self%root, self%nodes(i)%id) )then
                    ! one group beside the junk set aside re-enters with the rest leaf's rows, once per lineage (test switch)
                    call add_child(rest_sel(), 'states '//joined(mem)//' ('//int2str(n)//' particles), the one group '//&
                        &'beside the junk set aside, re-enters')
                    call write_ints(self%ndir(size(self%nodes))//'/reentered.txt', [1])
                else
                    call add_leaf(self, i, mem, .false., ng == 1)
                endif
            end do
        endif
        do j = 1, size(junk)
            if( count(lab == junk(j)) > 0 ) call add_leaf(self, i, [junk(j)], .true.)
        end do

    contains

        subroutine add_child( sel, what )
            logical,          intent(in) :: sel(:)
            character(len=*), intent(in) :: what
            integer, allocatable :: cr(:)
            integer :: jj
            fg   = fg + 1
            name = self%nodes(i)%id//int2str(fg)
            allocate(tmp(size(self%nodes) + 1))
            tmp(1:size(self%nodes)) = self%nodes
            tmp(size(tmp))%id    = name
            tmp(size(tmp))%depth = self%nodes(i)%depth + 1
            tmp(size(tmp))%step  = S_QUEUED
            tmp(size(tmp))%pproj = d//'/stageB/b.simple'
            tmp(size(tmp))%prows = '-'
            call move_alloc(tmp, self%nodes)
            call simple_mkdir(string(self%ndir(size(self%nodes))))
            cr = pack([(jj, jj=1,size(sel))], sel)
            if( local_rows(d) ) cr = rows(cr)
            call write_ints(self%ndir(size(self%nodes))//'/rows.txt', cr)
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': '//what//' -> node '//name
        end subroutine add_child

        function rest_sel() result( sel )
            logical, allocatable :: sel(:)
            integer, allocatable :: jk(:)
            integer :: jj
            allocate(sel(size(lab)), source=local_rows(d))
            if( .not. local_rows(d) ) sel(rows) = .true.
            jk = junk_rows(d)
            if( size(jk) > 0 ) sel(jk) = .false.
            sel = [(any(mem == lab(jj)) .or. (lab(jj) == 0 .and. sel(jj)), jj=1,size(lab))]
        end function rest_sel

    end subroutine decide

    subroutine finish( self, nn )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(out)   :: nn
        type(state_tree)          :: tree
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: labels(:), rows(:), names(:)
        real(dp),     allocatable :: counts(:)
        character(len=4096) :: line
        integer :: funit, ios, nl, il, ndup
        if( any(self%nodes(:)%step == S_FAILED) ) write(logfhandle,'(A)') '>>> HETERO_TREE WARNING: failed nodes'//&
            &failed_list(self)
        allocate(labels(self%nrows), source=0)
        ndup = 0
        nl   = 0
        if( file_exists(string(LEAF_TABLE)) ) nl = nlines(string(LEAF_TABLE))
        allocate(evens(nl), odds(nl), names(nl), counts(nl))
        nn = 0
        if( nl > 0 )then
            call fopen(funit, string(LEAF_TABLE), status='old', action='read', iostat=ios)
            call fileiochk('hetero_tree: cannot read '//LEAF_TABLE, ios)
            do il = 1, nl
                read(funit,'(A)') line
                rows = read_ints(word(line, 8))
                ndup = ndup + count(labels(rows) > 0)
                labels(rows) = il
                if( word(line, 2) == 'junk' ) cycle
                nn = nn + 1
                evens(nn)  = word(line, 6)
                odds(nn)   = word(line, 7)
                names(nn)  = il
                counts(nn) = real(str2int(word(line, 3)), dp)
            end do
            call fclose(funit)
        endif
        call write_ints('hetero_tree_labels.txt', labels)
        call unplaced_rows(self, read_ints(self%ndir(1)//'/rows.txt'), labels, ndup)
        write(logfhandle,'(A)') '>>> HETERO_TREE '//int2str(nl)//' leaves ('//int2str(nl - nn)//' junk) over '//&
            &int2str(count(labels > 0))//' particles: hetero_tree_labels.txt, '//LEAF_TABLE
        if( nn >= 2 )then
            call simple_mkdir(string('leaf_tree'))
            call simple_chdir(string('leaf_tree'))
            call tree%new(evens(1:nn), odds(1:nn), names(1:nn), counts(1:nn), find_img_smpd(evens(1)), 40., 10.)
            call tree%build('complete', 0.)
            call tree%write(labels)
            call tree%kill
            call simple_chdir(string(self%root))
        endif
    end subroutine finish

    subroutine unplaced_rows( self, root, labels, ndup )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: root(:), labels(:), ndup
        integer, allocatable :: seed(:), nr(:), sv(:), jk(:), miss(:)
        logical, allocatable :: junk(:)
        character(len=:), allocatable :: nd
        integer :: i, k, nj, nu, no, fu, io
        allocate(seed(self%nrows), source=-1)
        allocate(junk(self%nrows), source=.false.)
        do i = 1, size(self%nodes)
            nd = self%ndir(i)
            if( .not. file_exists(string(nd//'/seeds/seeds.txt')) ) cycle
            nr = read_ints(nd//'/rows.txt')
            sv = read_ints(nd//'/seeds/seeds.txt')
            if( local_rows(nd) )then
                seed(nr) = sv
            else
                seed(nr) = sv(nr)
            endif
            jk = junk_rows(nd)
            if( size(jk) == 0 ) cycle
            if( local_rows(nd) ) jk = nr(jk)
            junk(jk) = .true.
        end do
        miss = pack(root, labels(root) == 0)
        nj   = count(junk(miss))
        nu   = count(.not. junk(miss) .and. seed(miss) == 0)
        no   = size(miss) - nj - nu
        write(logfhandle,'(A)') '>>> HETERO_TREE leaf rows: '//int2str(count(labels > 0))//' particles in leaves, '//&
            &int2str(ndup)//' in two leaves, '//int2str(count(labels > 0) - count(labels(root) > 0))//' outside the root; '//&
            &int2str(size(miss))//' of the root''s in none: '//int2str(nj)//' junk, '//int2str(nu)//' unseeded, '//&
            &int2str(no)//' other'
        if( nu + no > 0 ) write(logfhandle,'(A)') '>>> HETERO_TREE WARNING: '//int2str(nu + no)//' particles are in no leaf '//&
            &'and not junk; rows and causes in hetero_tree_unplaced.txt'
        if( size(miss) == 0 ) return
        call fopen(fu, string('hetero_tree_unplaced.txt'), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write hetero_tree_unplaced.txt', io)
        write(fu,'(A)') '# input-project row, cause at the deepest node holding it (junk, unseeded, other)'
        do k = 1, size(miss)
            write(fu,'(I0,1X,A)') miss(k), trim(merge('junk    ', merge('unseeded', 'other   ', seed(miss(k)) == 0), &
                &junk(miss(k))))
        end do
        call fclose(fu)
    end subroutine unplaced_rows

    subroutine read_table( self )
        class(hetero_tree), intent(inout) :: self
        character(len=4096) :: line
        integer :: funit, ios, n, i
        n = nlines(string(self%ntab))
        allocate(self%nodes(n))
        call fopen(funit, string(self%ntab), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//self%ntab, ios)
        do i = 1, n
            read(funit,'(A)') line
            self%nodes(i)%id    = word(line, 1)
            self%nodes(i)%depth = str2int(word(line, 2))
            self%nodes(i)%step  = str2int(word(line, 3))
            self%nodes(i)%pproj = word(line, 4)
            self%nodes(i)%prows = word(line, 5)
        end do
        call fclose(funit)
    end subroutine read_table

    subroutine write_table( self )
        class(hetero_tree), intent(in) :: self
        integer :: funit, ios, i
        call fopen(funit, string(self%ntab//'.tmp'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//self%ntab, ios)
        do i = 1, size(self%nodes)
            write(funit,'(A)') self%nodes(i)%id//' '//int2str(self%nodes(i)%depth)//' '//int2str(self%nodes(i)%step)//&
                &' '//self%nodes(i)%pproj//' '//self%nodes(i)%prows
        end do
        call fclose(funit)
        call simple_rename(string(self%ntab//'.tmp'), string(self%ntab))
    end subroutine write_table

    function ndir( self, i ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i
        character(len=:), allocatable :: d
        d = self%root//'/node_'//self%nodes(i)%id
    end function ndir

    subroutine kill( self )
        class(hetero_tree), intent(inout) :: self
        if( allocated(self%qpool) )then
            call self%qpool%stop_exclusive_workers()
            call self%qpool%kill()
            deallocate(self%qpool)
        endif
        if( allocated(self%nodes) ) deallocate(self%nodes)
        if( allocated(self%fs)    ) deallocate(self%fs)
        if( allocated(self%tjunk) ) deallocate(self%tjunk)
        if( allocated(self%subd)  ) deallocate(self%subd)
        if( allocated(self%cl0)   ) deallocate(self%cl0)
        self%l_fsread = .false.
        call self%qenv%kill
        call self%qsite%kill
        call self%qsub%kill
        self%p => null()
    end subroutine kill

    subroutine final_project( self, fname )
        class(hetero_tree), intent(inout) :: self
        character(len=*),   intent(in)    :: fname
        type(sp_project)     :: master, nproj
        integer, allocatable :: rows(:), ix(:)
        character(len=4096)  :: line
        character(len=:), allocatable :: lrows, nd, nd_read
        integer :: funit, ios, nl, il, j, np
        call master%read(string(self%nodes(1)%pproj))
        np = master%os_ptcl3D%get_noris()
        do j = 1, np
            call master%os_ptcl3D%set_state(j, 0)
            call master%os_ptcl2D%set_state(j, 0)
        end do
        nl = nlines(string(LEAF_TABLE))
        nd_read = ''
        call fopen(funit, string(LEAF_TABLE), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//LEAF_TABLE, ios)
        do il = 1, nl
            read(funit,'(A)') line
            if( word(line, 2) == 'junk' ) cycle
            lrows = word(line, 8)
            nd    = lrows(1:index(lrows, '/leaf_') - 1)
            if( nd /= nd_read )then
                call nproj%kill
                call nproj%read_segment('ptcl3D', string(nd//'/stageB/b.simple'))
                ix      = node_index(self%nrows, nd)
                nd_read = nd
            endif
            rows = read_ints(lrows)
            do j = 1, size(rows)
                call master%os_ptcl3D%transfer_ori(rows(j), nproj%os_ptcl3D, ix(rows(j)))
                call master%os_ptcl3D%set_state(rows(j), il)
                call master%os_ptcl2D%set_state(rows(j), il)
            end do
        end do
        call fclose(funit)
        call master%write(string(fname))
        call master%kill
        call nproj%kill
    end subroutine final_project

    subroutine leaf_masks( fname )
        character(len=*), intent(in) :: fname
        character(len=4096)  :: line
        character(len=:), allocatable :: lrows, nd
        integer :: funit, fout, ios, il
        call fopen(funit, string(LEAF_TABLE), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//LEAF_TABLE, ios)
        call fopen(fout, string(fname), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//fname, ios)
        do il = 1, nlines(string(LEAF_TABLE))
            read(funit,'(A)') line
            if( word(line, 2) == 'junk' ) cycle
            lrows = word(line, 8)
            nd    = lrows(1:index(lrows, '/leaf_') - 1)
            if( file_exists(string(nd//'/mask.mrc')) ) write(fout,'(A)') int2str(il)//' '//nd//'/mask.mrc'
        end do
        call fclose(fout)
        call fclose(funit)
    end subroutine leaf_masks

    subroutine leaf_lines( lines, nds, tab )
        type(string), allocatable,  intent(out) :: lines(:), nds(:)
        character(len=*), optional, intent(in)  :: tab
        character(len=4096) :: line
        character(len=:), allocatable :: lrows, t
        integer :: funit, ios, nl, il
        t = LEAF_TABLE
        if( present(tab) ) t = tab
        nl = nlines(string(t))
        allocate(lines(nl), nds(nl))
        call fopen(funit, string(t), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//t, ios)
        do il = 1, nl
            read(funit,'(A)') line
            lines(il) = trim(line)
            lrows     = word(line, 8)
            nds(il)   = lrows(1:index(lrows, '/leaf_') - 1)
        end do
        call fclose(funit)
    end subroutine leaf_lines

    logical function keep_maps()
        character(len=8) :: ev
        integer :: ios
        call get_environment_variable('SIMPLE_TREE_KEEPMAPS', ev, status=ios)
        keep_maps = ios == 0 .and. ev == 'yes'
    end function keep_maps

    ! a step's partials, neighbour tables, reprojection models, embedding and sigma2 state (but stage B's) are read only inside it
    subroutine prune_step( dir )
        character(len=*), intent(in) :: dir
        type(string), allocatable :: list(:)
        character(len=:), allocatable :: rx
        integer :: j
        rx = '(_part[0-9]+_(even|odd)[.]mrc|^dist_part_neigh_[0-9]+[.]dat|^reprojection_model_(even|odd)[.]bin|^flex_pca_embedding[.]bin'
        if( index(dir, '/stageB', back=.true.) /= len(dir) - 6 ) rx = rx//'|^sigma2_state[.]bin'
        call simple_list_files_regexp(string(dir), rx//')$', list)
        do j = 1, size(list)
            call del_file(list(j))
        end do
    end subroutine prune_step

    subroutine prune_node( self, i )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i
        type(string), allocatable :: subs(:), list(:)
        character(len=:), allocatable :: sd, nm
        integer :: j, k
        logical :: l_kids
        if( self%p%site_layer .eq. 'yes' .or. self%l_motion ) return
        l_kids = any([(self%nodes(j)%pproj == self%ndir(i)//'/stageB/b.simple', j=1,size(self%nodes))])
        call del_file(string(self%ndir(i)//'/'//NODE_PROJ))
        call del_file(string(self%ndir(i)//'/pmean.mrc'))
        call simple_list_files_regexp(string(self%ndir(i)), &
            &'^(mean|flex|cheap|fflex|stA0|recB0|site|nsite_.+|sites|seeds|stageA|recB|split)$', subs)
        do j = 1, size(subs)
            sd = subs(j)%to_char()
            nm = sd(index(sd, '/', back=.true.)+1:)
            if( nm == 'recB0' .and. l_kids )then
                call simple_list_files_regexp(string(sd), '([.]bin|^recvol_state[0-9]+(_even|_odd)?[.]mrc)$', list)
            else if( nm == 'recB' .and. l_kids )then
                call simple_list_files_regexp(string(sd), '([.](simple|bin)|_(even|odd)(_unfil)?[.]mrc)$', list)
            else if( nm == 'sites' )then
                call simple_list_files_regexp(string(sd), '[.]mrc$', list)
            else
                call simple_list_files_regexp(string(sd), '[.](simple|mrc|bin)$', list)
            endif
            do k = 1, size(list)
                call del_file(list(k))
            end do
        end do
    end subroutine prune_node

    ! once the tree is done, the maps under its node directories that no leaf uses are removed
    subroutine prune_maps( self )
        class(hetero_tree), intent(inout) :: self
        type(string), allocatable :: keep(:), maps(:)
        character(len=4096) :: line
        character(len=:), allocatable :: lst, f
        integer(kind=8) :: sz, freed
        integer :: funit, ios, il, j, ndel, nfound
        allocate(keep(0), maps(0))
        if( .not. file_exists(string(self%root//'/'//self%ltab)) ) return
        call fopen(funit, string(self%root//'/'//self%ltab), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do il = 1, nlines(string(self%root//'/'//self%ltab))
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            keep = [keep, string(word(line, 5)), string(word(line, 6)), string(word(line, 7))]
        end do
        call fclose(funit)
        lst = self%root//'/prune_maps_list.txt'
        call exec_cmdline('find '//self%root//' -mindepth 2 -path "*/node_*" -type f -name "*.mrc" -o -path "*/node_*" -type f -name '//&
            &'"reprojection_model_*.bin" -o -path "*/node_*" -type f -name "dist_part_neigh_*.dat" > '//lst, suppress_errors=.true.)
        if( .not. file_exists(string(lst)) ) return
        call fopen(funit, string(lst), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do il = 1, nlines(string(lst))
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            maps = [maps, string(trim(line))]
        end do
        call fclose(funit)
        call del_file(string(lst))
        nfound = count([(any(maps == keep(j)), j=1,size(keep))])
        if( nfound < count([(file_exists(keep(j)), j=1,size(keep))]) )then
            write(logfhandle,'(A)') '>>> HETERO_TREE WARNING: leaf maps not found by path, intermediate maps kept'
            return
        endif
        freed = 0
        ndel  = 0
        do j = 1, size(maps)
            f = maps(j)%to_char()
            if( any(keep == maps(j)) ) cycle
            if( index(f, '/leafmaps/') > 0 .or. index(f, '/landscape/') > 0 ) cycle
            if( index(f, '/mask.mrc', back=.true.) == len(f) - 8 ) cycle
            if( index(f, '/plain/') > 0 .and. index(f, '_part') == 0 ) cycle
            inquire(file=f, size=sz)
            call del_file(maps(j))
            freed = freed + max(0_8, sz)
            ndel  = ndel + 1
        end do
        write(logfhandle,'(A)') '>>> HETERO_TREE removed '//int2str(ndel)//' intermediate maps and reprojection models ('//&
            &int2str(nint(real(freed,dp) / 1.e9_dp))//' GB); the leaves'' maps, plain leaf maps and final maps are kept'
    end subroutine prune_maps

    subroutine final_maps( self )
        class(hetero_tree), intent(inout) :: self
        type(state_tree)          :: tree
        type(string), allocatable :: lines(:), nds(:), evens(:), odds(:)
        real(dp),     allocatable :: counts(:)
        integer,      allocatable :: names(:), steps(:,:)
        logical,      allocatable :: ljunk(:)
        character(len=4096) :: line
        character(len=:), allocatable :: md, ev, od
        integer :: il, j, k, funit, ios
        md = self%root//'/final_maps'
        call leaf_lines(lines, nds)
        allocate(evens(0), odds(0), counts(0), names(0), ljunk(0))
        do il = 1, size(lines)
            j  = count([(nds(k) == nds(il), k=1,il)])
            ev = nds(il)%to_char()//'/plain/'//hv(j, 'even')
            od = nds(il)%to_char()//'/plain/'//hv(j, 'odd')
            line = lines(il)%to_char()
            if( .not. (file_exists(string(ev)) .and. file_exists(string(od))) )then
                write(logfhandle,'(A)') '>>> HETERO_TREE no plain leaf maps for '//word(line, 1)//': left out of the final maps'
                cycle
            endif
            evens  = [evens, string(ev)]
            odds   = [odds,  string(od)]
            counts = [counts, real(str2int(word(line, 3)), dp)]
            names  = [names, il]
            ljunk  = [ljunk, word(line, 2) == 'junk']
        end do
        call simple_mkdir(string(md))
        call fopen(funit, string(md//'/units_index.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write units_index.txt', ios)
        write(funit,'(A)') '# unit leaf particles kind'
        do j = 1, size(names)
            line = lines(names(j))%to_char()
            write(funit,'(A)') pad2(names(j))//' '//word(line, 1)//' '//word(line, 3)//' '//word(line, 2)
        end do
        call fclose(funit)
        if( size(evens) < 2 )then
            call fopen(funit, string(md//'/choice.txt'), status='replace', action='write', iostat=ios)
            call fileiochk('hetero_tree: cannot write choice.txt', ios)
            write(funit,'(A)') 'none'
            call fclose(funit)
            return
        endif
        call tree%new(evens, odds, names, counts, find_img_smpd(evens(1)), 40., 10.)
        call tree%build('complete', 0.)
        call tree%get_steps(steps)
        call tree%kill
        call tree_wiener_maps(evens, odds, counts, names, ljunk, steps, self%p%mskdiam, LPD, md)
    end subroutine final_maps

    integer function seed_site( nd )
        character(len=*), intent(in) :: nd
        if( file_exists(string(nd//'/sites/chosen_site.txt')) )then
            seed_site = str2int(trim(adjustl(first_line(nd//'/sites/chosen_site.txt', '1'))))
        else
            seed_site = str2int(trim(adjustl(first_line(nd//'/sites/seed_site.txt', '1'))))
        endif
    end function seed_site

    function seed_dir( nd ) result( sd )
        character(len=*), intent(in) :: nd
        character(len=:), allocatable :: sd
        sd = site_dir_of(nd, seed_site(nd))
    end function seed_dir

    function site_dir_of( nd, k ) result( sd )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: k
        character(len=:), allocatable :: sd
        sd = nd//'/site'
        if( k > 1 )then
            sd = nd//'/nsite_'//nd(index(nd, '/node_', back=.true.)+6:)//'_s'//int2str(k)
            ! fsite_ is the per-site layer's ranking of the same site; trees begun before nsite_ ranked there too
            if( .not. dir_exists(sd) ) sd = nd//'/fsite_'//nd(index(nd, '/node_', back=.true.)+6:)//'_s'//int2str(k)
        endif
    end function site_dir_of

    function first_line( fname, dflt ) result( ln )
        character(len=*), intent(in) :: fname, dflt
        character(len=:), allocatable :: ln
        character(len=256) :: buf
        integer :: fu, io
        ln = dflt
        if( .not. file_exists(string(fname)) ) return
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        read(fu,'(A)',iostat=io) buf
        if( io == 0 .and. len_trim(buf) > 0 ) ln = trim(buf)
        call fclose(fu)
    end function first_line

    integer function nsites( nd )
        character(len=*), intent(in) :: nd
        character(len=4096) :: ln
        integer :: fu, io
        nsites = 0
        if( .not. file_exists(string(nd//'/sites/sites.txt')) ) return
        call fopen(fu, string(nd//'/sites/sites.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) > 0 .and. ln(1:1) /= '#' ) nsites = nsites + 1
        end do
        call fclose(fu)
    end function nsites

    function two_groups( nd ) result( lab )
        character(len=*), intent(in) :: nd
        integer, allocatable :: lab(:)
        lab = two_groups_in(nd//'/fflex')
    end function two_groups

    function two_groups_in( fd ) result( lab )
        character(len=*), intent(in) :: fd
        integer, allocatable :: lab(:), st(:), pres(:), grp(:)
        type(string), allocatable :: evens(:), odds(:)
        real(dp),     allocatable :: cnt(:)
        type(state_tree) :: tree
        integer :: j
        allocate(lab(0))
        if( .not. file_exists(string(fd//'/f.simple')) ) return
        st   = project_states(fd//'/f.simple')
        pres = with_maps(fd, present_labels(st))
        deallocate(lab)
        allocate(lab(size(st)), source=0)
        if( size(pres) < 2 ) return
        evens = [(string(fd//'/'//hv(pres(j), 'even')), j=1,size(pres))]
        odds  = [(string(fd//'/'//hv(pres(j), 'odd')), j=1,size(pres))]
        cnt   = [(real(count(st == pres(j)), dp), j=1,size(pres))]
        call tree%new(evens, odds, pres, cnt, find_img_smpd(evens(1)), 40., 10.)
        call tree%build('complete', 0.)
        call tree%cut(2, grp)
        call tree%kill
        do j = 1, size(pres)
            where( st == pres(j) ) lab = grp(j)
        end do
        call nearest_group(fd//'/flex_pca_coordinates.txt', lab)
    end function two_groups_in

    ! every coordinate row still at 0 joins the nearer group centre in the flex latent, per axis scaled by the pooled spread
    subroutine nearest_group( fname, lab )
        character(len=*), intent(in)    :: fname
        integer,          intent(inout) :: lab(:)
        real(dp), allocatable :: z(:,:), mu(:,:), sd(:), d(:)
        integer,  allocatable :: pinds(:), g(:)
        integer :: n, q, k, i, a
        call flex_coords(fname, pinds, z)
        n = size(pinds)
        q = size(z, 2)
        if( n == 0 ) return
        g = lab(pinds)
        if( q < 1 .or. count(g == 1) < 2 .or. count(g == 2) < 2 .or. count(g == 0) == 0 ) return
        allocate(mu(q,2), sd(q), d(2))
        do a = 1, 2
            mu(:,a) = [(sum(z(:,k), mask=g == a) / real(count(g == a), dp), k=1,q)]
        end do
        do k = 1, q
            sd(k) = sqrt(max((sum((z(:,k) - mu(k,1))**2, mask=g == 1) + sum((z(:,k) - mu(k,2))**2, mask=g == 2)) / &
                &real(count(g > 0) - 2, dp), 1.e-12_dp))
        end do
        do i = 1, n
            if( g(i) /= 0 ) cycle
            d = [(sum(((z(i,:) - mu(:,a)) / sd)**2), a=1,2)]
            lab(pinds(i)) = minloc(d, dim=1)
        end do
        write(logfhandle,'(A)') '>>> HETERO_TREE '//fname(1:index(fname, '/', back=.true.)-1)//': the two-group fallback seeds '//&
            &'its '//int2str(count(g == 0))//' rows flex left unassigned by the nearer group centre'
    end subroutine nearest_group

    subroutine flex_coords( fname, pinds, z )
        character(len=*),      intent(in)  :: fname
        integer,  allocatable, intent(out) :: pinds(:)
        real(dp), allocatable, intent(out) :: z(:,:)
        character(len=8192) :: line
        real    :: e(3), u(1)
        integer :: funit, ios, n, q, k
        allocate(pinds(0), z(0,0))
        if( .not. file_exists(string(fname)) ) return
        call fopen(funit, string(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        n = 0
        q = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            n = n + 1
            if( q == 0 ) q = nwords(line) - 5
        end do
        rewind(funit)
        deallocate(pinds, z)
        allocate(pinds(n))
        allocate(z(n,max(q,0)), source=0._dp)
        k = 0
        do while( k < n )
            read(funit,'(A)') line
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            k = k + 1
            read(line,*) pinds(k), e(1), e(2), e(3), u(1), z(k,1:q)
        end do
        call fclose(funit)
    end subroutine flex_coords

    function with_maps( dir, pres ) result( kept )
        character(len=*), intent(in) :: dir
        integer,          intent(in) :: pres(:)
        integer, allocatable :: kept(:)
        integer :: j
        kept = pack(pres, [(file_exists(string(dir//'/'//hv(pres(j), 'even'))) .and. &
            &file_exists(string(dir//'/'//hv(pres(j), 'odd'))), j=1,size(pres))])
    end function with_maps

    function compact( st ) result( ch )
        integer, intent(in)  :: st(:)
        integer, allocatable :: ch(:), pres(:)
        integer :: j
        pres = present_labels(st)
        allocate(ch(size(st)), source=0)
        do j = 1, size(pres)
            where( st == pres(j) ) ch = j
        end do
    end function compact

    function fold_small( st ) result( fs )
        integer, intent(in)  :: st(:)
        integer, allocatable :: fs(:), pres(:)
        integer :: j, big
        fs   = st
        pres = present_labels(st)
        if( size(pres) < 2 ) return
        big  = pres(maxloc([(count(st == pres(j)), j=1,size(pres))], dim=1))
        do j = 1, size(pres)
            if( count(st == pres(j)) <= PROB_MULTISTATE_MIN_POP ) where( st == pres(j) ) fs = big
        end do
    end function fold_small

    subroutine node_project( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        integer, allocatable :: rows(:), st(:), ix(:)
        integer :: j
        rows = read_ints(self%ndir(i)//'/rows.txt')
        if( self%l_nodeproj )then
            ix = [(j, j=1,self%nrows)]
            if( self%nodes(i)%depth > 0 ) ix = node_index(self%nrows, parent_dir(self%nodes(i)%pproj))
            if( any(ix(rows) == 0) ) THROW_HARD('a node row is missing from its parent''s project; node_project')
            call write_cut(self%nodes(i)%pproj, self%ndir(i)//'/n_cut.simple', ix(rows), self%p%node_subset_io .eq. 'yes')
            call write_ints(self%ndir(i)//'/'//NODE_LOCAL, [size(rows)])
            call simple_rename(string(self%ndir(i)//'/n_cut.simple'), string(self%ndir(i)//'/'//NODE_PROJ))
            return
        endif
        allocate(st(self%nrows), source=0)
        st(rows) = 1
        if( size(project_states(self%nodes(i)%pproj)) /= self%nrows )then
            THROW_HARD('the parent project does not hold every particle; node_project')
        endif
        call write_states(self%nodes(i)%pproj, self%ndir(i)//'/'//NODE_PROJ, st)
    end subroutine node_project

    subroutine oneflex_inputs( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        type(string), allocatable :: evens(:), odds(:)
        type(image)           :: vol
        integer,  allocatable :: rows(:), lab(:), pres(:)
        real(dp), allocatable :: w(:)
        real,     allocatable :: acc(:,:,:)
        character(len=:), allocatable :: d, pd, fn
        integer :: j, ldim(3), ifoo
        real    :: smpd, wsum
        d    = self%ndir(i)
        pd   = parent_dir(self%nodes(i)%pproj)
        rows = read_ints(d//'/rows.txt')
        lab  = node_index(self%nrows, pd)
        rows = lab(rows)
        allocate(evens(0), odds(0), w(0))
        if( file_exists(string(pd//'/recB0/b.simple')) )then
            lab  = project_states(pd//'/recB0/b.simple')
            pres = present_labels(lab(rows))
            do j = 1, size(pres)
                if( .not. file_exists(string(pd//'/recB0/'//hv(pres(j), 'even'))) ) cycle
                evens = [evens, string(pd//'/recB0/'//hv(pres(j), 'even'))]
                odds  = [odds,  string(pd//'/recB0/'//hv(pres(j), 'odd'))]
                w     = [w, real(count(lab(rows) == pres(j)), dp)]
            end do
        endif
        if( size(evens) >= 2 )then
            call falloff_mask(self, evens, odds, w, d//'/mask.mrc')
        else
            call simple_copy_file(string(pd//'/mask.mrc'), string(d//'/mask.mrc'))
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': its parent''s mask is its focus mask'
        endif
        lab  = project_states(self%nodes(i)%pproj)
        pres = present_labels(lab(rows))
        wsum = 0.
        do j = 1, size(pres)
            fn = pd//'/recB/recvol_state'//pad2(pres(j))//'.mrc'
            call find_ldim_nptcls(string(fn), ldim, ifoo)
            smpd = find_img_smpd(string(fn))
            call vol%new(ldim, smpd)
            call vol%read(string(fn))
            if( j == 1 )then
                acc = real(count(lab(rows) == pres(j))) * vol%get_rmat()
            else
                acc = acc + real(count(lab(rows) == pres(j))) * vol%get_rmat()
            endif
            wsum = wsum + real(count(lab(rows) == pres(j)))
            call vol%kill
        end do
        call vol%new(ldim, smpd)
        call vol%set_rmat(acc / max(wsum, 1.), .false.)
        call vol%write(string(d//'/pmean.mrc'))
        call vol%kill
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': flex mean from the parent''s maps of states '//&
            &joined(pres)//', mask from '//int2str(size(evens))//' parent state maps'
    end subroutine oneflex_inputs

    function parent_dir( pproj ) result( pd )
        character(len=*), intent(in) :: pproj
        character(len=:), allocatable :: pd
        pd = pproj(1:index(pproj, '/', back=.true.)-1)
        pd = pd(1:index(pd, '/', back=.true.)-1)
    end function parent_dir

    subroutine cheap_children( self, i, child )
        class(hetero_tree),   intent(inout) :: self
        integer,              intent(in)    :: i
        integer, allocatable, intent(out)   :: child(:)
        type(sp_project)     :: spproj
        real,    allocatable :: z(:,:), zc(:,:), y(:), cen(:,:), v(:,:), bm(:,:), cov(:,:), vals(:), evec(:,:)
        integer, allocatable :: pinds(:), lab(:), bins(:), nbin(:), inds(:), sel(:)
        character(len=8192)  :: line
        character(len=:), allocatable :: fname
        real    :: e(3), u(3), ev(1)
        integer :: n, q, nb, it, j, b, funit, ios, k, nsel, r, np, t
        if( self%l_nodeseed ) call seed_rnd_fixed(RND_SEED + id_seed(self%nodes(i)%id))
        fname = self%ndir(i)//'/flex/flex_pca_coordinates.txt'
        call fopen(funit, string(fname), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//fname, ios)
        n = 0
        q = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            n = n + 1
            if( q == 0 ) q = nwords(line) - 5
        end do
        rewind(funit)
        allocate(pinds(n), z(n,q))
        k = 0
        do while( k < n )
            read(funit,'(A)') line
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            k = k + 1
            read(line,*) pinds(k), e(1), e(2), e(3), u(1), z(k,:)
        end do
        call fclose(funit)
        call spproj%read_segment('ptcl3D', string(self%ndir(i)//'/'//NODE_PROJ))
        np = spproj%os_ptcl3D%get_noris()
        allocate(v(3,n))
        do j = 1, n
            e      = spproj%os_ptcl3D%get_euler(pinds(j)) * real(PI) / 180.
            v(:,j) = [sin(e(2))*cos(e(1)), sin(e(2))*sin(e(1)), cos(e(2))]
        end do
        call spproj%kill
        nb  = max(NBIN_MIN, nint(real(n) / real(PERBIN)))
        sel = [(j, j=1,n)]
        do j = 1, nb
            r      = j + int(ran3() * real(n - j + 1))
            r      = min(r, n)
            t      = sel(j)
            sel(j) = sel(r)
            sel(r) = t
        end do
        cen = v(:, sel(1:nb))
        allocate(lab(n))
        do it = 1, NKM_ITS + 1
            do j = 1, n
                lab(j) = maxloc(matmul(v(:,j), cen), dim=1)
            end do
            if( it > NKM_ITS ) exit
            do b = 1, nb
                if( count(lab == b) > 0 )then
                    u = [sum(v(1,:), mask=lab == b), sum(v(2,:), mask=lab == b), sum(v(3,:), mask=lab == b)]
                    cen(:,b) = u / max(norm2(u), 1.e-12)
                else
                    cen(:,b) = v(:, min(n, 1 + int(ran3() * real(n))))
                endif
            end do
        end do
        allocate(bins(n), source=0)
        nsel = 0
        do b = 1, nb
            if( count(lab == b) == 0 ) cycle
            nsel = nsel + 1
            where( lab == b ) bins = nsel
        end do
        nb   = nsel
        nbin = [(count(bins == b), b=1,nb)]
        allocate(bm(nb,q), source=0.)
        do j = 1, n
            bm(bins(j),:) = bm(bins(j),:) + z(j,:)
        end do
        do b = 1, nb
            bm(b,:) = bm(b,:) / real(nbin(b))
        end do
        allocate(zc(n,q))
        do j = 1, n
            zc(j,:) = z(j,:) - bm(bins(j),:)
        end do
        do k = 1, q
            zc(:,k) = zc(:,k) - sum(zc(:,k)) / real(n)
        end do
        cov = matmul(transpose(zc), zc) / real(max(n - 1, 1))
        allocate(evec(q,1))
        call eigh(q, cov, 1, ev, evec)
        y = matmul(zc, evec(:,1))
        allocate(child(np), source=0)
        do b = 1, nb
            sel  = pack([(j, j=1,n)], bins == b)
            vals = y(sel)
            inds = [(j, j=1,size(sel))]
            call hpsort(vals, inds)
            do r = 1, size(sel)
                child(pinds(sel(inds(r)))) = min((r - 1) * NCHILD / size(sel), NCHILD - 1) + 1
            end do
        end do
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' cheap children: '//int2str(n)//&
            &' particles, '//int2str(nb)//' direction bins, '//int2str(q)//'-dim latent, populations '//&
            &joined([(count(child == b), b=1,NCHILD)])
    end subroutine cheap_children

    subroutine falloff_mask( self, evens, odds, w, fname )
        class(hetero_tree), intent(inout) :: self
        type(string),       intent(in)    :: evens(:), odds(:)
        real(dp),           intent(in)    :: w(:)
        character(len=*),   intent(in)    :: fname
        type(image)          :: img
        real,    allocatable :: vs(:,:,:), m(:,:,:)
        logical, allocatable :: sph(:,:,:), sup(:,:,:), core(:,:,:)
        integer :: n
        real    :: smpd
        call half_variance(evens, odds, w, self%p%mskdiam, n, smpd, sph, sup, vs)
        core = vs > MFRAC * maxval(vs)
        if( .not. any(core) )then
            m = merge(1., 0., sph)
        else
            m    = dist_px(core, n, smpd)
            core = m <= MDIL
            m    = dist_px(core, n, smpd)
            m    = exp(-m**2 / (2. * MFALL**2))
            where( .not. sph ) m = 0.
        endif
        call img%new([n,n,n], smpd)
        call img%set_rmat(m, .false.)
        call img%write(string(fname))
        call img%set_rmat(vs, .false.)
        call img%write(string(fname(1:len(fname)-4)//'_variance.mrc'))
        call img%kill
        write(logfhandle,'(A)') '>>> HETERO_TREE focus mask '//fname//' from '//int2str(size(evens))//' maps: '//&
            &int2str(count(m > 0.5))//' voxels above 0.5 ('//int2str(nint(100. * real(count(m > 0.5 .and. sup)) / &
            &real(max(1,count(sup)))))//' % of the molecule support)'
    end subroutine falloff_mask

    subroutine half_variance( evens, odds, w, mskdiam, n, smpd, sph, sup, vs, de, dv )
        type(string),                   intent(in)  :: evens(:), odds(:)
        real(dp),                       intent(in)  :: w(:)
        real,                           intent(in)  :: mskdiam
        integer,                        intent(out) :: n
        real,                           intent(out) :: smpd
        logical, allocatable,           intent(out) :: sph(:,:,:), sup(:,:,:)
        real,    allocatable,           intent(out) :: vs(:,:,:)
        real,    allocatable, optional, intent(out) :: de(:,:,:,:), dv(:,:,:,:)
        real,     allocatable :: e(:,:,:,:), o(:,:,:,:), em(:,:,:), om(:,:,:), v(:,:,:)
        real(dp), allocatable :: wn(:)
        integer :: ldim(3), nm, j, ifoo
        real    :: sc
        call find_ldim_nptcls(evens(1), ldim, ifoo)
        n    = ldim(1)
        nm   = size(evens)
        smpd = find_img_smpd(evens(1))
        wn   = w / sum(w)
        allocate(e(n,n,n,nm), o(n,n,n,nm))
        do j = 1, nm
            call read_lp(evens(j), e(:,:,:,j))
            call read_lp(odds(j),  o(:,:,:,j))
        end do
        sph = sphere(n, mskdiam / 2. / smpd)
        em  = wsum(e)
        om  = wsum(o)
        sup = (0.5 * (em + om) > MSUP * maxval(0.5 * (em + om), mask=sph)) .and. sph
        do j = 1, nm
            sc = sum(e(:,:,:,j) * em, mask=sup) / max(sum(e(:,:,:,j)**2, mask=sup), 1.e-30)
            e(:,:,:,j) = e(:,:,:,j) * sc
            sc = sum(o(:,:,:,j) * om, mask=sup) / max(sum(o(:,:,:,j)**2, mask=sup), 1.e-30)
            o(:,:,:,j) = o(:,:,:,j) * sc
        end do
        em = wsum(e)
        om = wsum(o)
        allocate(v(n,n,n), source=0.)
        do j = 1, nm
            v = v + real(wn(j)) * (e(:,:,:,j) - em) * (o(:,:,:,j) - om)
        end do
        vs = filtered_vol(v, n, smpd, 2)
        where( .not. sph ) vs = 0.
        if( present(de) .and. present(dv) )then
            do j = 1, nm
                e(:,:,:,j) = e(:,:,:,j) - em
                o(:,:,:,j) = o(:,:,:,j) - om
            end do
            call move_alloc(e, de)
            call move_alloc(o, dv)
        endif

    contains

        function wsum( x ) result( s )
            real, intent(in) :: x(:,:,:,:)
            real, allocatable :: s(:,:,:)
            integer :: jj
            allocate(s(n,n,n), source=0.)
            do jj = 1, nm
                s = s + real(wn(jj)) * x(:,:,:,jj)
            end do
        end function wsum

        subroutine read_lp( f, x )
            type(string), intent(in)  :: f
            real,         intent(out) :: x(:,:,:)
            type(image) :: vol
            call vol%new([n,n,n], smpd)
            call vol%read(f)
            call fourier_filter(vol, n, smpd, 1)
            x = vol%get_rmat()
            call vol%kill
        end subroutine read_lp

    end subroutine half_variance

    subroutine discover_sites( evens, odds, counts, names, mskdiam, outdir, nsite, pgrp )
        type(string),     intent(in)  :: evens(:), odds(:)
        real(dp),         intent(in)  :: counts(:)
        integer,          intent(in)  :: names(:)
        real,             intent(in)  :: mskdiam
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: nsite
        character(len=*), optional, intent(in) :: pgrp
        type :: dcand
            real, allocatable :: m(:,:,:)
            real    :: r = 0., share = 0., mol = 0., rad = 0.
            integer :: comp = 0, fac = 0, ncore = 0, seed = 0, vox = 0
        end type dcand
        type :: dreg
            integer, allocatable :: v(:)
            real    :: cut = 0.
            integer :: comp = 0
            logical :: top = .true.
        end type dreg
        type(dcand), allocatable :: cands(:), tmp(:)
        type(dreg),  allocatable :: regs(:), rtmp(:)
        type(image)              :: img
        real,     allocatable :: vs(:,:,:), de(:,:,:,:), dv(:,:,:,:), pn(:,:), cc(:,:), vv(:)
        logical,  allocatable :: sph(:,:,:), sup(:,:,:), core(:,:,:), cm(:,:,:), b(:,:,:), okk(:)
        integer,  allocatable :: lab(:,:,:), lb(:,:,:), asg(:), vidx(:,:), ord(:), rep(:), memb(:)
        real(dp), allocatable :: wn(:)
        real,     allocatable :: rots(:,:,:)
        character(len=:), allocatable :: mlist
        real    :: smpd, tot, cosv, sh(2), cut2, dm
        integer :: n, nk, nl, ic, nv, iv, it, j, i, k, x, y, z, nb, funit, ios, ir, nc0, np2
        logical :: l_sib, l_fac
        wn  = counts / sum(counts)
        nk  = size(evens)
        okk = counts >= real(DSEEDMIN,dp)
        call half_variance(evens, odds, counts, mskdiam, n, smpd, sph, sup, vs, de, dv)
        core = vs > MFRAC * maxval(vs)
        call label6(core, lab, nl)
        tot  = sum(vs, mask=core)
        allocate(cands(0), regs(0))
        do ic = 1, nl
            call push(lab == ic, MFRAC * maxval(vs), .true., ic)
        end do
        ir  = 0
        np2 = nl
        do while( ir < size(regs) )
            ir = ir + 1
            if( allocated(cm) ) deallocate(cm)
            allocate(cm(n,n,n), source=.false.)
            call set_voxels(cm, regs(ir)%v)
            if( count(cm) < DMINVOX ) cycle
            if( regs(ir)%top .and. sum(vs, mask=cm) / tot < DMINSHARE ) cycle
            ic    = regs(ir)%comp
            nc0   = size(cands)
            l_fac = .false.
            if( count(cm) >= 2 * DMINVOX )then
                nv = count(cm)
                allocate(pn(nv,nk), vidx(3,nv), vv(nv))
                iv = 0
                do z = 1, n
                    do y = 1, n
                        do x = 1, n
                            if( .not. cm(x,y,z) ) cycle
                            iv = iv + 1
                            vidx(:,iv) = [x, y, z]
                            vv(iv)     = vs(x,y,z)
                            pn(iv,:)   = real(sqrt(wn)) * 0.5 * (de(x,y,z,:) + dv(x,y,z,:))
                            pn(iv,:)   = pn(iv,:) / max(sqrt(sum(pn(iv,:)**2)), 1.e-30)
                        end do
                    end do
                end do
                cc = top_dirs(pn, 2)
                do it = 1, DFACT_ITS
                    asg = merge(1, 2, abs(matmul(pn, cc(1,:))) >= abs(matmul(pn, cc(2,:))))
                    do j = 1, 2
                        if( count(asg == j) > 2 ) cc(j:j,:) = top_dirs(pn(pack([(iv, iv=1,nv)], asg == j),:), 1)
                    end do
                end do
                cosv = abs(dot_product(cc(1,:), cc(2,:)))
                sh   = [(sum(vv, mask=asg == j) / sum(vv), j=1,2)]
                if( cosv < DFACT_COS .and. minval(sh) >= DFACT_SHARE .and. min(count(asg == 1), count(asg == 2)) >= DMINVOX )then
                    write(logfhandle,'(A,I0,A,F5.2,A,F5.2,A,F5.2)') '>>> SITE_DISCOVER component ', ic, &
                        &' splits into two factors: |cos| of their state profiles ', cosv, ', variance shares ', sh(1), ' / ', sh(2)
                    do j = 1, 2
                        allocate(b(n,n,n), source=.false.)
                        do iv = 1, nv
                            if( asg(iv) == j ) b(vidx(1,iv),vidx(2,iv),vidx(3,iv)) = .true.
                        end do
                        call label6(b, lb, nb)
                        do k = 1, nb
                            if( count(lb == k) < max(20, DMINVOX / 4) ) where( lb == k ) b = .false.
                        end do
                        call add(b, ic, j, cc(j,:) / real(sqrt(wn)))
                        deallocate(b)
                    end do
                    l_fac = .true.
                else
                    write(logfhandle,'(A,I0,A,F5.2,A,F5.2,A,F5.2,A)') '>>> SITE_DISCOVER component ', ic, ' is one factor (|cos| ', &
                        &cosv, ', shares ', sh(1), ' / ', sh(2), ')'
                endif
                deallocate(pn, vidx, vv)
            endif
            if( .not. l_fac ) call add(cm, ic, 0)
            if( all([(cands(k)%mol <= DMAXMOL, k=nc0+1,size(cands))]) ) cycle
            ! a component whose candidates are not local is re-cut inside itself at its cut plus MFRAC of its excess, while that splits it
            cut2 = regs(ir)%cut + MFRAC * (maxval(vs, mask=cm) - regs(ir)%cut)
            call label6(cm .and. vs > cut2, lb, nb)
            write(logfhandle,'(A,I0,A,F6.3,A,I0,A)') '>>> SITE_DISCOVER component ', ic, ' is not local: re-cut at ', cut2 / maxval(vs), &
                &' x max gives ', count([(count(lb == k) >= DMINVOX, k=1,nb)]), ' piece(s) of at least '//int2str(DMINVOX)//' voxels'
            if( count([(count(lb == k) >= DMINVOX, k=1,nb)]) < 2 ) cycle
            do k = 1, nb
                if( count(lb == k) < DMINVOX ) cycle
                np2 = np2 + 1
                call push(lb == k, cut2, .false., np2)
            end do
        end do
        nsite = 0
        ord   = pack([(i, i=1,size(cands))], [(cands(i)%mol <= DMAXMOL, i=1,size(cands))])
        write(logfhandle,'(A)') '>>> SITE_DISCOVER '//int2str(size(cands) - size(ord))//' of '//int2str(size(cands))//&
            &' candidates dropped as not local (mask > 0.5 over '//int2str(nint(100. * DMAXMOL))//' % of the molecule)'
        if( size(ord) == 0 ) return
        do i = 2, size(ord)
            j = ord(i)
            k = i - 1
            do while( k >= 1 )
                if( .not. ranks_above(j, ord(k)) ) exit
                ord(k+1) = ord(k)
                k = k - 1
            end do
            ord(k+1) = j
        end do
        allocate(rep(0))
        allocate(memb(size(cands)), source=0)
        rots = group_rotations(pgrp)
        do i = 1, size(ord)
            k = 0
            do j = 1, size(rep)
                dm = mate_dice(cands(ord(i))%m, cands(rep(j))%m, rots)
                if( dm >= DMATE )then
                    write(logfhandle,'(A,F6.3,A)') '>>> SITE_DISCOVER candidate '//ctag(ord(i))//' is a symmetry image of candidate '//&
                        &ctag(rep(j))//' (Dice ', dm, ' after a group rotation): folded into it, its mask not kept'
                    k = j
                    exit
                endif
                if( dice(cands(ord(i))%m, cands(rep(j))%m) < DDICE ) cycle
                l_sib = .false.
                do it = 1, size(cands)
                    if( memb(it) /= j ) cycle
                    l_sib = l_sib .or. (cands(ord(i))%fac > 0 .and. cands(it)%fac > 0 .and. cands(it)%comp == cands(ord(i))%comp)
                end do
                if( l_sib ) cycle
                k = j
                exit
            end do
            if( k == 0 )then
                rep = [rep, ord(i)]
                k   = size(rep)
            endif
            memb(ord(i)) = k
        end do
        nsite = size(rep)
        call img%new([n,n,n], smpd)
        call fopen(funit, string(outdir//'/sites.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('discover_sites: cannot write '//outdir//'/sites.txt', ios)
        write(funit,'(A)') '# site candidate repro voxels_above_0.5 molecule_share centre_radius_A seed_state members'
        do j = 1, nsite
            i = rep(j)
            call img%set_rmat(cands(i)%m, .false.)
            call img%write(string(outdir//'/site_'//pad2(j)//'.mrc'))
            mlist = ''
            do it = 1, size(cands)
                if( memb(it) == j ) mlist = mlist//' '//ctag(it)
            end do
            write(funit,'(I0,1X,A,1X,F6.3,1X,I0,1X,F6.3,1X,F6.1,1X,I0,A)') j, ctag(i), cands(i)%r, cands(i)%vox, cands(i)%mol, &
                &cands(i)%rad, cands(i)%seed, mlist
            write(logfhandle,'(A,I0,A,F6.3,A,F5.1,A,F5.1,A,I0)') '>>> SITE_DISCOVER SITE ', j, ': '//ctag(i)//' repro ', cands(i)%r, &
                &', mask > 0.5 ', 100. * cands(i)%mol, ' % of the molecule, centre ', cands(i)%rad, ' A, seed state ', cands(i)%seed
        end do
        call fclose(funit)
        call img%kill

    contains

        subroutine push( m, cut, top, comp )
            logical, intent(in) :: m(:,:,:)
            real,    intent(in) :: cut
            logical, intent(in) :: top
            integer, intent(in) :: comp
            integer :: q
            allocate(rtmp(size(regs)+1))
            rtmp(1:size(regs)) = regs
            call move_alloc(rtmp, regs)
            regs(size(regs))%v   = pack([(q, q=1,size(m))], reshape(m, [size(m)]))
            regs(size(regs))%cut = cut
            regs(size(regs))%top  = top
            regs(size(regs))%comp = comp
        end subroutine push

        subroutine set_voxels( m, v )
            logical, intent(inout) :: m(:,:,:)
            integer, intent(in)    :: v(:)
            integer :: q, nn
            nn = size(m,1)
            do q = 1, size(v)
                m(mod(v(q)-1, nn) + 1, mod((v(q)-1) / nn, nn) + 1, (v(q)-1) / (nn*nn) + 1) = .true.
            end do
        end subroutine set_voxels

        subroutine add( c, comp, fac, prof )
            logical,        intent(in) :: c(:,:,:)
            integer,        intent(in) :: comp, fac
            real, optional, intent(in) :: prof(:)
            real(dp) :: num, dsq, osq, ek(nk)
            real     :: cen(3)
            integer  :: kk, nc, xx, yy, zz
            allocate(tmp(size(cands)+1))
            tmp(1:size(cands)) = cands
            call move_alloc(tmp, cands)
            nc = size(cands)
            cands(nc)%m = dist_px(dilate6(c, nint(MDIL)), n, smpd)
            cands(nc)%m = exp(-cands(nc)%m**2 / (2. * MFALL**2))
            where( .not. sph ) cands(nc)%m = 0.
            num = 0._dp
            dsq = 0._dp
            osq = 0._dp
            do kk = 1, nk
                ek(kk) = sum(real(de(:,:,:,kk) * cands(nc)%m, dp) * real(dv(:,:,:,kk) * cands(nc)%m, dp))
                num    = num + wn(kk) * ek(kk)
                dsq    = dsq + wn(kk) * sum(real(de(:,:,:,kk) * cands(nc)%m, dp)**2)
                osq    = osq + wn(kk) * sum(real(dv(:,:,:,kk) * cands(nc)%m, dp)**2)
            end do
            cands(nc)%r     = real(num / max(sqrt(dsq * osq), 1.e-300_dp))
            cands(nc)%comp  = comp
            cands(nc)%fac   = fac
            cands(nc)%ncore = count(c)
            cands(nc)%vox   = count(cands(nc)%m > 0.5)
            cands(nc)%mol   = real(count(cands(nc)%m > 0.5 .and. sup)) / real(max(1, count(sup)))
            cands(nc)%share = sum(vs, mask=c) / tot
            cen = 0.
            do zz = 1, n
                do yy = 1, n
                    do xx = 1, n
                        if( c(xx,yy,zz) ) cen = cen + real([xx, yy, zz] - 1)
                    end do
                end do
            end do
            cen = cen / real(max(1, count(c)))
            cands(nc)%rad = sqrt(sum((cen - real(n / 2))**2)) * smpd
            if( present(prof) )then
                cands(nc)%seed = names(maxloc(merge(abs(prof), -huge(1.), okk), dim=1))
            else
                cands(nc)%seed = names(maxloc(merge(real(ek), -huge(1.), okk), dim=1))
            endif
            write(logfhandle,'(A,I0,A,F5.1,A,F5.1,A,F5.1,A,F6.3,A,I0)') '>>> SITE_DISCOVER candidate '//ctag(nc)//': core ', &
                &cands(nc)%ncore, ' voxels, ', 100. * cands(nc)%share, ' % of the node variance, mask > 0.5 ', 100. * cands(nc)%mol, &
                &' % of the molecule, centre ', cands(nc)%rad, ' A, repro ', cands(nc)%r, ', seed state ', cands(nc)%seed
        end subroutine add

        function ctag( i ) result( t )
            integer, intent(in) :: i
            character(len=:), allocatable :: t
            t = int2str(cands(i)%comp)
            if( cands(i)%fac > 0 ) t = t//'f'//int2str(cands(i)%fac)
        end function ctag

        logical function ranks_above( a, bb )
            integer, intent(in) :: a, bb
            ranks_above = cands(a)%share > cands(bb)%share .or. (cands(a)%share == cands(bb)%share .and. cands(a)%r > cands(bb)%r)
        end function ranks_above

        real function dice( m1, m2 )
            real, intent(in) :: m1(:,:,:), m2(:,:,:)
            dice = 2. * real(count(m1 > 0.5 .and. m2 > 0.5)) / real(max(1, count(m1 > 0.5) + count(m2 > 0.5)))
        end function dice

    end subroutine discover_sites

    function top_dirs( p, m ) result( c )
        real,    intent(in) :: p(:,:)
        integer, intent(in) :: m
        real, allocatable :: c(:,:)
        real(dp) :: g(size(p,2),size(p,2)), u(size(p,2),size(p,2)), lam(size(p,2))
        integer  :: nrot, j, ord(size(p,2))
        real     :: key(size(p,2))
        g   = matmul(transpose(real(p,dp)), real(p,dp))
        call jacobi(g, size(p,2), size(p,2), lam, u, nrot)
        key = real(-lam)
        ord = [(j, j=1,size(p,2))]
        call hpsort(key, ord)
        allocate(c(m,size(p,2)))
        do j = 1, m
            c(j,:) = real(u(:,ord(j)))
        end do
    end function top_dirs

    subroutine label6( mk, lab, nl )
        logical,              intent(in)  :: mk(:,:,:)
        integer, allocatable, intent(out) :: lab(:,:,:)
        integer,              intent(out) :: nl
        integer, allocatable :: stack(:,:)
        integer :: n1, n2, n3, x, y, z, top, c(3), d(3), q
        integer, parameter :: NB(3,6) = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], [3,6])
        n1 = size(mk,1)
        n2 = size(mk,2)
        n3 = size(mk,3)
        allocate(lab(n1,n2,n3), source=0)
        allocate(stack(3,max(1,count(mk))))
        nl = 0
        do z = 1, n3
            do y = 1, n2
                do x = 1, n1
                    if( .not. mk(x,y,z) .or. lab(x,y,z) > 0 ) cycle
                    nl  = nl + 1
                    top = 1
                    stack(:,1) = [x, y, z]
                    lab(x,y,z) = nl
                    do while( top > 0 )
                        c   = stack(:,top)
                        top = top - 1
                        do q = 1, 6
                            d = c + NB(:,q)
                            if( any(d < 1) .or. d(1) > n1 .or. d(2) > n2 .or. d(3) > n3 ) cycle
                            if( .not. mk(d(1),d(2),d(3)) .or. lab(d(1),d(2),d(3)) > 0 ) cycle
                            lab(d(1),d(2),d(3)) = nl
                            top = top + 1
                            stack(:,top) = d
                        end do
                    end do
                end do
            end do
        end do
    end subroutine label6

    function dilate6( mk, nit ) result( d )
        logical, intent(in) :: mk(:,:,:)
        integer, intent(in) :: nit
        logical, allocatable :: d(:,:,:), t(:,:,:)
        integer :: it, n1, n2, n3
        n1 = size(mk,1)
        n2 = size(mk,2)
        n3 = size(mk,3)
        d  = mk
        do it = 1, nit
            t = d
            t(2:n1,:,:) = t(2:n1,:,:) .or. d(1:n1-1,:,:)
            t(1:n1-1,:,:) = t(1:n1-1,:,:) .or. d(2:n1,:,:)
            t(:,2:n2,:) = t(:,2:n2,:) .or. d(:,1:n2-1,:)
            t(:,1:n2-1,:) = t(:,1:n2-1,:) .or. d(:,2:n2,:)
            t(:,:,2:n3) = t(:,:,2:n3) .or. d(:,:,1:n3-1)
            t(:,:,1:n3-1) = t(:,:,1:n3-1) .or. d(:,:,2:n3)
            d = t
        end do
    end function dilate6

    function filtered_vol( x, n, smpd, mode ) result( y )
        real,    intent(in) :: x(:,:,:), smpd
        integer, intent(in) :: n, mode
        real, allocatable :: y(:,:,:)
        type(image) :: vol
        call vol%new([n,n,n], smpd)
        call vol%set_rmat(x, .false.)
        call fourier_filter(vol, n, smpd, mode)
        y = vol%get_rmat()
        call vol%kill
    end function filtered_vol

    function dist_px( mk, n, smpd ) result( dpx )
        logical, intent(in) :: mk(:,:,:)
        integer, intent(in) :: n
        real,    intent(in) :: smpd
        real, allocatable :: dpx(:,:,:)
        type(image) :: vol
        call vol%new([n,n,n], smpd)
        call vol%set_rmat(merge(0., 1., mk), .false.)
        call vol%distance_transform
        dpx = vol%get_rmat() / smpd
        call vol%kill
    end function dist_px

    subroutine add_leaf( self, i, mem, l_junk, l_rest )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, mem(:)
        logical,            intent(in)    :: l_junk
        logical, optional,  intent(in)    :: l_rest
        character(len=11), parameter :: SFX(3) = [character(len=11) :: '', '_even_unfil', '_odd_unfil']
        type(image)          :: vol
        logical, allocatable :: inn(:)
        integer, allocatable :: rows(:), lab(:), sel(:), jk(:)
        real,    allocatable :: acc(:,:,:)
        character(len=:), allocatable :: d, name, b, fn, lrows
        real    :: wsum, smpd
        integer :: j, isfx, funit, ios, n, ldim(3), ifoo
        d    = self%ndir(i)
        rows = read_ints(d//'/rows.txt')
        lab  = project_states(d//'/stageB/b.simple')
        sel  = pack([(j, j=1,size(lab))], [(any(mem == lab(j)), j=1,size(lab))])
        if( present(l_rest) )then
            if( l_rest )then
                allocate(inn(size(lab)), source=local_rows(d))
                if( .not. local_rows(d) ) inn(rows) = .true.
                jk = junk_rows(d)
                if( size(jk) > 0 ) inn(jk) = .false.
                sel = pack([(j, j=1,size(lab))], [(any(mem == lab(j)) .or. (lab(j) == 0 .and. inn(j)), j=1,size(lab))])
            endif
        endif
        n    = size(sel)
        if( size(mem) == 1 )then
            name = self%nodes(i)%id//'_s'//int2str(mem(1))
            b    = d//'/stageB/recvol_state'//pad2(mem(1))
        else
            name = self%nodes(i)%id//'_g'//joined(mem, '_')
            call simple_mkdir(string(d//'/leafmaps'))
            b = d//'/leafmaps/'//name
            do isfx = 1, 3
                wsum = 0.
                do j = 1, size(mem)
                    fn = d//'/stageB/recvol_state'//pad2(mem(j))//trim(SFX(isfx))//'.mrc'
                    call find_ldim_nptcls(string(fn), ldim, ifoo)
                    smpd = find_img_smpd(string(fn))
                    call vol%new(ldim, smpd)
                    call vol%read(string(fn))
                    if( j == 1 )then
                        acc = real(count(lab == mem(j))) * vol%get_rmat()
                    else
                        acc = acc + real(count(lab == mem(j))) * vol%get_rmat()
                    endif
                    wsum = wsum + real(count(lab == mem(j)))
                    call vol%kill
                end do
                call vol%new(ldim, smpd)
                call vol%set_rmat(acc / max(wsum, 1.), .false.)
                call vol%write(string(b//trim(SFX(isfx))//'.mrc'))
                call vol%kill
            end do
        endif
        lrows = d//'/leaf_'//name//'_rows.txt'
        if( local_rows(d) ) sel = rows(sel)
        call write_ints(lrows, sel)
        call fopen(funit, string(self%root//'/'//self%ltab), status='unknown', position='append', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot append to '//self%ltab, ios)
        write(funit,'(A)') name//' '//merge('junk','leaf',l_junk)//' '//int2str(n)//' '//joined(mem, '+')//' '//b//'.mrc '//&
            &b//'_even_unfil.mrc '//b//'_odd_unfil.mrc '//lrows
        call fclose(funit)
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': states '//joined(mem)//' ('//int2str(n)//&
            &' particles) -> '//merge('junk leaf ','leaf      ',l_junk)//name
    end subroutine add_leaf

    subroutine common_keys( self, cl )
        class(hetero_tree), intent(in)    :: self
        class(cmdline),     intent(inout) :: cl
        integer :: np
        call cl%set('pgrp',    trim(self%p%pgrp))
        call cl%set('mskdiam', self%p%mskdiam)
        np = merge(self%npcur, self%p%nparts, self%npcur > 0)
        ! refine3D and reconstruct3D take any nparts as distributed: one part runs shared-memory without it
        if( np > 1 )then
            call cl%set('nparts', np)
        else
            call cl%delete('nparts')
        endif
        call cl%set('nthr',    self%p%nthr)
        if( self%p%rec_plane_reuse == 'yes' ) call cl%set('rec_plane_reuse', 'yes')
        if( self%p%state_assembly > 0 )then
            call cl%set('state_assembly', self%p%state_assembly)
            call cl%set('assembly_mem', self%p%assembly_mem)
        endif
    end subroutine common_keys

    subroutine crop_keys( self, cl, boxc )
        class(hetero_tree), intent(in)    :: self
        class(cmdline),     intent(inout) :: cl
        integer,            intent(in)    :: boxc
        if( self%box <= boxc ) return
        call cl%set('box_crop',  boxc)
        call cl%set('smpd_crop', real(self%box) * self%smpd / real(boxc))
    end subroutine crop_keys

    ! a reconstruction whose consumers read only its unfiltered halves needs neither ML regularization nor postprocessing
    subroutine halves_only_keys( cl )
        class(cmdline), intent(inout) :: cl
        call cl%set('ml_reg',      'no')
        call cl%set('postprocess', 'no')
    end subroutine halves_only_keys

    subroutine flex_keys( self, cl, mean, k )
        class(hetero_tree), intent(in)    :: self
        class(cmdline),     intent(inout) :: cl
        character(len=*),   intent(in)    :: mean
        integer,            intent(in)    :: k
        call cl%set('prg',            'flex_pca')
        call cl%set('projfile',       'f.simple')
        call cl%set('vol1',           mean)
        call cl%set('lp',             LPF)
        call cl%set('npreimages',     k)
        call cl%set('rec_states',     'yes')
        call cl%set('min_state_frac', 0.)
        call cl%set('umap',           'no')
        call cl%set('n_probe_iters',  NPROBE)
        call common_keys(self, cl)
    end subroutine flex_keys

    subroutine refine_keys( self, cl, lp, maxits, nspace, boxc, mask )
        class(hetero_tree), intent(in)    :: self
        class(cmdline),     intent(inout) :: cl
        real,               intent(in)    :: lp
        integer,            intent(in)    :: maxits, nspace, boxc
        character(len=*),   intent(in)    :: mask
        call cl%set('refine',          'prob_neigh')
        call cl%set('prob_neigh_mode', 'geom')
        call cl%set('prob_athres',     ATHR)
        call cl%set('lp',              lp)
        call cl%set('maxits',          maxits)
        call cl%set('center',          'no')
        call cl%set('rec_backend',     'pcg')
        call cl%set('pcg_mskfile',     mask)
        call cl%set('pcg_focus',       'yes')
        call cl%set('nspace',          nspace)
        call cl%set('box_crop',        boxc)
        call cl%set('smpd_crop',       real(self%box) * self%smpd / real(boxc))
        call common_keys(self, cl)
    end subroutine refine_keys

    subroutine fourier_filter( vol, n, smpd, mode )
        type(image), intent(inout) :: vol
        integer,     intent(in)    :: n, mode
        real,        intent(in)    :: smpd
        complex, allocatable :: cm(:,:,:)
        real    :: f2, fac
        integer :: h, k, l, nyq, ky, lz
        nyq = n / 2
        call vol%fft
        cm = vol%get_cmat()
        do l = 1, size(cm,3)
            lz = l - 1 - merge(n, 0, l > nyq)
            do k = 1, size(cm,2)
                ky = k - 1 - merge(n, 0, k > nyq)
                do h = 1, size(cm,1)
                    f2 = real((h-1)**2 + ky*ky + lz*lz) / real(n*n)
                    if( mode == 1 )then
                        fac = exp(-f2 / smpd**2 * MLP**2 * log(2.))
                    else
                        fac = exp(-2. * real(PI)**2 * MSIG**2 * f2)
                    endif
                    cm(h,k,l) = cm(h,k,l) * fac
                end do
            end do
        end do
        call vol%set_cmat(cm)
        call vol%ifft
    end subroutine fourier_filter

    function sphere( n, rad ) result( sph )
        integer, intent(in) :: n
        real,    intent(in) :: rad
        logical, allocatable :: sph(:,:,:)
        integer :: i, j, k, c
        c = n / 2 + 1
        allocate(sph(n,n,n))
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    sph(i,j,k) = real((i-c)**2 + (j-c)**2 + (k-c)**2) <= rad**2
                end do
            end do
        end do
    end function sphere

    function hv( s, half ) result( f )
        integer,          intent(in) :: s
        character(len=*), intent(in) :: half
        character(len=:), allocatable :: f
        type(string) :: t
        t = refine3D_state_halfvol_fname(s, half, .true.)
        f = t%to_char()
    end function hv

    function project_states( fname ) result( st )
        character(len=*), intent(in) :: fname
        integer, allocatable :: st(:)
        type(sp_project) :: spproj
        call spproj%read_segment('ptcl3D', string(fname))
        st = spproj%os_ptcl3D%get_all_asint('state')
        call spproj%kill
    end function project_states

    ! a copy of the project with only the particle state segments rewritten in place
    subroutine write_states( src, dst, st )
        character(len=*), intent(in) :: src, dst
        integer,          intent(in) :: st(:)
        call simple_copy_file(string(src), string(dst))
        call set_seg('ptcl3D')
        call set_seg('ptcl2D')

    contains

        subroutine set_seg( oritype )
            character(len=*), intent(in) :: oritype
            type(sp_project), target :: spproj
            class(oris),     pointer :: field
            integer :: j
            call spproj%read_segment(oritype, string(dst))
            call spproj%ptr2oritype(oritype, field)
            if( field%get_noris() /= size(st) ) THROW_HARD('state count does not match the project; write_states')
            do j = 1, size(st)
                call field%set_state(j, st(j))
            end do
            nullify(field)
            call spproj%write_segment_inside(oritype, string(dst))
            call spproj%kill
        end subroutine set_seg

    end subroutine write_states

    subroutine write_cut( src, dst, idx, l_direct )
        character(len=*), intent(in) :: src, dst
        integer,          intent(in) :: idx(:)
        logical,          intent(in) :: l_direct
        class(sp_project), allocatable :: selected
        if( any(idx(2:) <= idx(:size(idx)-1)) ) THROW_HARD('rows must ascend; write_cut')
        if( l_direct )then
            allocate(selected)
            call selected%write_selected_particles(string(src), string(dst), idx)
            call selected%kill
            deallocate(selected)
        else
            call simple_copy_file(string(src), string(dst))
        endif
        call cut_seg('ptcl3D')
        call cut_seg('ptcl2D')

    contains

        subroutine cut_seg( oritype )
            character(len=*), intent(in) :: oritype
            type(sp_project), target :: spproj
            class(oris),     pointer :: field
            logical, allocatable :: sel(:)
            integer :: j
            call spproj%read_segment(oritype, string(dst))
            call spproj%ptr2oritype(oritype, field)
            if( field%get_noris() == 0 )then
                nullify(field)
                call spproj%kill
                return
            endif
            if( .not. l_direct )then
                if( maxval(idx) > field%get_noris() ) THROW_HARD('rows beyond the project '//src//'; write_cut')
                allocate(sel(field%get_noris()), source=.false.)
                sel(idx) = .true.
                call field%compress(sel)
            endif
            do j = 1, field%get_noris()
                call field%set_state(j, 1)
            end do
            nullify(field)
            call spproj%write_segment_inside(oritype, string(dst))
            call spproj%kill
        end subroutine cut_seg

    end subroutine write_cut

    ! a node made with node_proj=yes: row j of its projects is row rows.txt(j) of the input project
    logical function local_rows( nd )
        character(len=*), intent(in) :: nd
        local_rows = file_exists(string(nd//'/'//NODE_LOCAL))
    end function local_rows

    function node_index( nrows, nd ) result( ix )
        integer,          intent(in) :: nrows
        character(len=*), intent(in) :: nd
        integer, allocatable :: ix(:), rows(:)
        integer :: j
        if( .not. local_rows(nd) )then
            ix = [(j, j=1,nrows)]
            return
        endif
        allocate(ix(nrows), source=0)
        rows     = read_ints(nd//'/rows.txt')
        ix(rows) = [(j, j=1,size(rows))]
    end function node_index

    ! node_rows=yes: the step's project keeps only its rows with state > 0 until expand_proj; NR_ROWS (the master size, then
    ! each kept row's master row) stays beside the step's per-row diagnostics
    subroutine compact_proj( sd, pf, l_direct )
        character(len=*), intent(in) :: sd, pf
        logical,          intent(in) :: l_direct
        type(sp_project)     :: spproj
        logical, allocatable :: sel(:)
        integer :: fu, io, j, nmics, nstks, nptcls
        call simple_copy_file(string(sd//'/'//pf), string(sd//'/'//NR_PROJ))
        call spproj%read_segment('ptcl3D', string(sd//'/'//pf))
        sel = spproj%os_ptcl3D%get_all_asint('state') > 0
        if( l_direct )then
            call spproj%read_data_info(string(sd//'/'//pf), nmics, nstks, nptcls)
            if( nptcls /= size(sel) ) THROW_HARD('ptcl2D and ptcl3D differ in size; compact_proj')
        else
            call spproj%read_segment('ptcl2D', string(sd//'/'//pf))
            if( spproj%os_ptcl2D%get_noris() /= size(sel) ) THROW_HARD('ptcl2D and ptcl3D differ in size; compact_proj')
        endif
        call write_ints(sd//'/'//NR_ROWS, [size(sel), pack([(j, j=1,size(sel))], sel)])
        if( l_direct )then
            call spproj%write_selected_particles(string(sd//'/'//NR_PROJ), string(sd//'/'//pf), &
                &pack([(j, j=1,size(sel))], sel))
        else
            call spproj%os_ptcl3D%compress(sel)
            call spproj%os_ptcl2D%compress(sel)
            call spproj%write_segment_inside('ptcl3D', string(sd//'/'//pf))
            call spproj%write_segment_inside('ptcl2D', string(sd//'/'//pf))
        endif
        call spproj%kill
        call fopen(fu, string(sd//'/'//NR_NAME), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write '//sd//'/'//NR_NAME, io)
        write(fu,'(A)') pf
        call fclose(fu)
    end subroutine compact_proj

    subroutine expand_proj( sd )
        character(len=*), intent(in) :: sd
        integer, allocatable :: rows(:)
        character(len=:), allocatable :: pf
        if( .not. file_exists(string(sd//'/'//NR_PROJ)) ) return
        pf = first_line(sd//'/'//NR_NAME, '')
        if( len(pf) == 0 ) THROW_HARD('no step project named in '//sd//'/'//NR_NAME//'; expand_proj')
        pf   = sd//'/'//pf
        rows = read_ints(sd//'/'//NR_ROWS)
        call expand_seg('ptcl3D')
        call expand_seg('ptcl2D')
        call del_file(sd//'/'//NR_PROJ)

    contains

        subroutine expand_seg( oritype )
            character(len=*), intent(in) :: oritype
            type(sp_project), target :: full, cpt
            class(oris),     pointer :: fo, co
            integer :: j
            call cpt%read_segment(oritype, string(pf))
            call cpt%ptr2oritype(oritype, co)
            if( co%get_noris() /= rows(1) )then
                if( co%get_noris() /= size(rows) - 1 ) THROW_HARD('rows of '//pf//' match neither its cut nor its master; expand_proj')
                call full%read_segment(oritype, string(sd//'/'//NR_PROJ))
                call full%ptr2oritype(oritype, fo)
                do j = 1, size(rows) - 1
                    call fo%transfer_ori(rows(j + 1), co, j)
                end do
                call full%write_segment_inside(oritype, string(pf))
                nullify(fo)
                call full%kill
            endif
            nullify(co)
            call cpt%kill
        end subroutine expand_seg

    end subroutine expand_proj

    subroutine write_contrast( proj, cfile )
        character(len=*), intent(in) :: proj, cfile
        type(sp_project), target :: spproj
        character(len=256) :: ln
        integer, allocatable :: r(:)
        real,    allocatable :: v(:)
        integer :: fu, io, n, j
        if( .not. file_exists(string(cfile)) ) return
        allocate(r(nlines(string(cfile))), v(nlines(string(cfile))))
        call fopen(fu, string(cfile), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        n = 0
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
            n = n + 1
            read(ln,*) r(n), v(n)
        end do
        call fclose(fu)
        call spproj%read_segment('ptcl3D', string(proj))
        do j = 1, n
            call spproj%os_ptcl3D%set(r(j), 'contrast', v(j))
        end do
        call spproj%write_segment_inside('ptcl3D', string(proj))
        call spproj%kill
    end subroutine write_contrast

    function present_labels( lab ) result( pres )
        integer, intent(in)  :: lab(:)
        integer, allocatable :: pres(:)
        integer :: s, smax
        smax = max(1, maxval([0,lab]))
        pres = pack([(s, s=1,smax)], [(count(lab == s) > 0, s=1,smax)])
    end function present_labels

    integer function mode_of( v )
        integer, intent(in) :: v(:)
        integer :: s, c, best
        mode_of = 0
        best    = 0
        do s = 1, maxval(v)
            c = count(v == s)
            if( c > best )then
                best    = c
                mode_of = s
            endif
        end do
    end function mode_of

    function read_ints( fname ) result( v )
        character(len=*), intent(in) :: fname
        integer, allocatable :: v(:)
        integer :: funit, ios, n, j
        if( .not. file_exists(string(fname)) )then
            allocate(v(0))
            return
        endif
        n = nlines(string(fname))
        allocate(v(n))
        if( n == 0 ) return
        call fopen(funit, string(fname), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//fname, ios)
        do j = 1, n
            read(funit,*,iostat=ios) v(j)
            if( ios /= 0 ) THROW_HARD('cannot read line '//int2str(j)//' of '//int2str(n)//' in '//fname//'; read_ints')
        end do
        call fclose(funit)
    end function read_ints

    subroutine stamp_local( fname )
        character(len=*), intent(in) :: fname
        type(sp_project) :: sp
        call sp%read_segment('compenv', string(fname))
        if( sp%compenv%get_noris() > 0 )then
            call sp%compenv%set(1, 'qsys_name', 'local')
            call sp%write_segment_inside('compenv', string(fname))
        endif
        call sp%kill
    end subroutine stamp_local

    integer function id_seed( id )
        character(len=*), intent(in) :: id
        integer :: k
        id_seed = 0
        do k = 1, len(id)
            id_seed = mod(31 * id_seed + ichar(id(k:k)), 1000003)
        end do
    end function id_seed

    subroutine write_ints( fname, v )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: v(:)
        integer :: funit, ios, j
        call fopen(funit, string(fname), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//fname, ios)
        do j = 1, size(v)
            write(funit,'(I0)') v(j)
        end do
        call fclose(funit)
    end subroutine write_ints

    integer function count_poles( fname )
        character(len=*), intent(in) :: fname
        character(len=4096) :: ln
        integer :: fu, io
        count_poles = 0
        if( .not. file_exists(string(fname)) ) return
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len(word(ln, 3)) == 0 ) cycle
            if( word(ln, 3) == 'rest' .and. len(word(ln, 4)) == 0 ) cycle
            count_poles = count_poles + 1
        end do
        call fclose(fu)
    end function count_poles

    function site_mask( self, j ) result( f )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: j
        character(len=:), allocatable :: f
        f = self%root//'/node_'//self%fs(j)%home//'/sites/site_'//pad2(self%fs(j)%k)//'.mrc'
    end function site_mask

    subroutine site_poles( sd, nrows, rows, pp, tp )
        character(len=*),     intent(in)  :: sd
        integer,              intent(in)  :: nrows
        integer, allocatable, intent(out) :: rows(:)
        real,    allocatable, intent(out) :: pp(:,:)
        logical, allocatable, intent(out) :: tp(:,:)
        type(string), allocatable :: ends(:)
        integer,      allocatable :: pole(:), it(:,:), top(:)
        real,         allocatable :: acc(:,:), v(:)
        logical,      allocatable :: inn(:)
        integer :: np, nr, e, a, fu, io
        if( file_exists(string(sd//'/poles.bin')) )then
            open(newunit=fu, file=sd//'/poles.bin', access='stream', form='unformatted', status='old', action='read', iostat=io)
            if( io == 0 )then
                read(fu, iostat=io) np, nr
                if( io == 0 )then
                    allocate(rows(nr), pp(nr,np), it(nr,np))
                    read(fu, iostat=io) rows, pp, it
                endif
                close(fu)
                if( io == 0 )then
                    tp = it > 0
                    return
                endif
                if( allocated(rows) ) deallocate(rows, pp, it)
            endif
        endif
        call pole_ends(sd, ends, pole, .false.)
        np = maxval([0, pole])
        allocate(acc(nrows,np), source=0.)
        allocate(it(nrows,np), source=0)
        allocate(inn(nrows), source=.false.)
        do e = 1, size(ends)
            v   = read_reals(sd//'/site_post_'//ends(e)%to_char()//'.txt')
            top = read_ints(sd//'/site_top_'//ends(e)%to_char()//'.txt')
            if( size(v) /= nrows .or. size(top) /= nrows ) THROW_HARD('site ranking files do not match the project: '//sd)
            inn = inn .or. v >= 0.
            acc(:,pole(e)) = acc(:,pole(e)) + v
            where( top > 0 ) it(:,pole(e)) = 1
        end do
        do a = 1, np
            acc(:,a) = acc(:,a) / real(count(pole == a))
        end do
        rows = pack([(e, e=1,nrows)], inn)
        nr   = size(rows)
        pp   = acc(rows,:)
        it   = it(rows,:)
        tp   = it > 0
        open(newunit=fu, file=sd//'/poles.bin', access='stream', form='unformatted', status='replace', action='write', iostat=io)
        if( io == 0 )then
            write(fu) np, nr
            write(fu) rows, pp, it
            close(fu)
        endif
    end subroutine site_poles

    subroutine pole_ends( sd, ends, pole, l_halves )
        character(len=*),          intent(in)  :: sd
        logical,                   intent(in)  :: l_halves
        type(string), allocatable, intent(out) :: ends(:)
        integer,      allocatable, intent(out) :: pole(:)
        type(string), allocatable :: nm(:)
        integer,      allocatable :: lab(:), pid(:)
        real,         allocatable :: z(:,:)
        character(len=4096) :: ln
        integer :: fu, io, nw, w, np, a, b, i
        allocate(ends(0), pole(0))
        if( file_exists(string(sd//'/starts.txt')) )then
            call fopen(fu, string(sd//'/starts.txt'), status='old', action='read', iostat=io)
            if( io /= 0 ) return
            np = 0
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                if( nwords(ln) < 2 ) cycle
                np   = np + 1
                ends = [ends, string(word(ln, 1))]
                pole = [pole, np]
            end do
            call fclose(fu)
            return
        endif
        if( l_halves )then
            if( .not. file_exists(string(sd//'/site_halves.txt')) ) return
            allocate(nm(0))
            call fopen(fu, string(sd//'/site_halves.txt'), status='old', action='read', iostat=io)
            if( io /= 0 ) return
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                if( ln(1:1) == '#' .or. index(ln, '|') == 0 ) cycle
                if( any(word(ln, 1) == [character(len=8) :: 'merge', 'accept', 'quarter', 'twin', 'rest']) ) cycle
                nm = [nm, string(word(ln, 1))]
            end do
            allocate(z(size(nm),2), source=-1.e6)
            allocate(lab(size(nm)), source=[(i, i=1,size(nm))])
            rewind(fu)
            do
                read(fu,'(A)',iostat=io) ln
                if( io /= 0 ) exit
                nw = nwords(ln)
                if( word(ln, 1) == 'accept' )then
                    do i = 1, size(nm)
                        if( nm(i)%to_char() /= word(ln, 2) ) cycle
                        z(i, merge(1, 2, word(ln, 3) == 'row')) = str2real(word(ln, nw - 1))
                    end do
                else if( word(ln, 1) == 'merge' .and. word(ln, nw) == 'yes' )then
                    a = 0
                    b = 0
                    do i = 1, size(nm)
                        if( nm(i)%to_char() == word(ln, 2) ) a = i
                        if( nm(i)%to_char() == word(ln, 3) ) b = i
                    end do
                    if( a > 0 .and. b > 0 ) where( lab == lab(b) ) lab = lab(a)
                endif
            end do
            call fclose(fu)
            np = 0
            allocate(pid(size(nm)), source=0)
            do i = 1, size(nm)
                if( lab(i) == 0 ) cycle
                if( any(lab(1:i-1) == lab(i)) )then
                    pid(i) = pid(findloc(lab(1:i-1), lab(i), dim=1))
                else
                    np     = np + 1
                    pid(i) = np
                endif
                pole = [pole, pid(i)]
                ends = [ends, nm(i)]
            end do
            return
        endif
        if( .not. file_exists(string(sd//'/site_poles.txt')) ) return
        call fopen(fu, string(sd//'/site_poles.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        np = 0
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            nw = nwords(ln)
            if( nw < 3 ) cycle
            if( word(ln, 3) == 'rest' .and. nw == 3 ) cycle
            np = np + 1
            do w = 3, nw
                ends = [ends, string(word(ln, w))]
                pole = [pole, np]
            end do
        end do
        call fclose(fu)
    end subroutine pole_ends

    integer function npoles_of( sd )
        character(len=*), intent(in) :: sd
        type(string), allocatable :: ends(:)
        integer,      allocatable :: pole(:)
        call pole_ends(sd, ends, pole, .false.)
        npoles_of = maxval([0, pole])
    end function npoles_of

    function pole_z( sd, np ) result( z )
        character(len=*), intent(in) :: sd
        integer,          intent(in) :: np
        real, allocatable :: z(:)
        type(string), allocatable :: ends(:)
        integer,      allocatable :: pole(:)
        character(len=4096) :: ln
        integer :: fu, io, e
        real    :: v
        allocate(z(np), source=huge(1.))
        call pole_ends(sd, ends, pole, .false.)
        if( file_exists(string(sd//'/site_halves.txt')) )then
            call fopen(fu, string(sd//'/site_halves.txt'), status='old', action='read', iostat=io)
            if( io == 0 )then
                do
                    read(fu,'(A)',iostat=io) ln
                    if( io /= 0 ) exit
                    if( word(ln, 1) /= 'accept' ) cycle
                    v = str2real(word(ln, 6))
                    do e = 1, size(ends)
                        if( ends(e)%to_char() == word(ln, 2) .and. pole(e) <= np ) z(pole(e)) = min(z(pole(e)), v)
                    end do
                end do
                call fclose(fu)
            endif
        endif
        where( z == huge(1.) ) z = 0.
    end function pole_z

    subroutine pole_map( fe, fo, fout, hstem, mskdiam )
        character(len=*), intent(in) :: fe, fo, fout, hstem
        real,             intent(in) :: mskdiam
        type(image) :: ve, vo
        complex, allocatable :: cm(:,:,:)
        integer :: ldim(3), ifoo, n, h, k, l, ky, lz
        real    :: smpd, f, k1, k2
        call simple_copy_file(string(fe), string(hstem//'_even.mrc'))
        call simple_copy_file(string(fo), string(hstem//'_odd.mrc'))
        call find_ldim_nptcls(string(fe), ldim, ifoo)
        n    = ldim(1)
        smpd = find_img_smpd(string(fe))
        call ve%new(ldim, smpd)
        call vo%new(ldim, smpd)
        call ve%read(string(fe))
        call vo%read(string(fo))
        call ve%set_rmat(0.5 * (ve%get_rmat() + vo%get_rmat()), .false.)
        call ve%fft
        cm = ve%get_cmat()
        k1 = 1. / (LPD + 1.)
        k2 = 1. / (LPD - 0.5)
        do l = 1, size(cm,3)
            lz = l - 1 - merge(n, 0, l > n / 2)
            do k = 1, size(cm,2)
                ky = k - 1 - merge(n, 0, k > n / 2)
                do h = 1, size(cm,1)
                    f = sqrt(real((h-1)**2 + ky*ky + lz*lz)) / (real(n) * smpd)
                    if( f >= k2 )then
                        cm(h,k,l) = cmplx(0.,0.)
                    else if( f > k1 )then
                        cm(h,k,l) = cm(h,k,l) * (0.5 + 0.5 * cos(real(PI) * (f - k1) / (k2 - k1)))
                    endif
                end do
            end do
        end do
        call ve%set_cmat(cm)
        call ve%ifft
        call ve%set_rmat(ve%get_rmat() * merge(1., 0., sphere(n, mskdiam / 2. / smpd)), .false.)
        call ve%write(string(fout))
        call ve%kill
        call vo%kill
    end subroutine pole_map

    real function mask_dice( f1, f2 )
        character(len=*), intent(in) :: f1, f2
        type(image) :: v
        real, allocatable :: m1(:,:,:), m2(:,:,:)
        integer :: ldim(3), ifoo
        mask_dice = 0.
        if( .not. (file_exists(string(f1)) .and. file_exists(string(f2))) ) return
        call find_ldim_nptcls(string(f1), ldim, ifoo)
        call v%new(ldim, find_img_smpd(string(f1)))
        call v%read(string(f1))
        m1 = v%get_rmat()
        call v%read(string(f2))
        m2 = v%get_rmat()
        call v%kill
        if( any(shape(m1) /= shape(m2)) ) return
        mask_dice = 2. * real(count(m1 > 0.5 .and. m2 > 0.5)) / real(max(1, count(m1 > 0.5) + count(m2 > 0.5)))
    end function mask_dice

    function junk_rows( nd ) result( jr )
        character(len=*), intent(in) :: nd
        integer, allocatable :: jr(:)
        jr = junk_of(seed_dir(nd))
    end function junk_rows

    function junk_of( sd ) result( jr )
        character(len=*), intent(in) :: sd
        integer, allocatable :: jr(:), buf(:)
        character(len=256) :: ln
        integer :: fu, io, r, fl, n
        allocate(jr(0))
        if( .not. file_exists(string(sd//'/site_junk.txt')) ) return
        allocate(buf(nlines(string(sd//'/site_junk.txt'))))
        call fopen(fu, string(sd//'/site_junk.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        n = 0
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
            read(ln,*) r, fl
            if( fl /= 1 ) cycle
            n      = n + 1
            buf(n) = r
        end do
        call fclose(fu)
        jr = buf(1:n)
    end function junk_of

    subroutine write_reals( fname, v )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: v(:)
        integer :: funit, ios, j
        call fopen(funit, string(fname), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//fname, ios)
        do j = 1, size(v)
            write(funit,'(F8.5)') v(j)
        end do
        call fclose(funit)
    end subroutine write_reals

    function read_reals( fname ) result( v )
        character(len=*), intent(in) :: fname
        real, allocatable :: v(:)
        integer :: funit, ios
        if( .not. file_exists(string(fname)) )then
            allocate(v(0))
            return
        endif
        allocate(v(nlines(string(fname))))
        call fopen(funit, string(fname), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//fname, ios)
        read(funit,*) v
        call fclose(funit)
    end function read_reals

    real function best_minz( sd )
        character(len=*), intent(in) :: sd
        type(string), allocatable :: ends(:)
        integer,      allocatable :: pole(:)
        real,         allocatable :: z(:,:)
        character(len=4096) :: ln
        integer :: fu, io, nw, i
        best_minz = -1.e6
        call pole_ends(sd, ends, pole, .true.)
        if( size(ends) == 0 ) return
        allocate(z(size(ends),2), source=-1.e6)
        call fopen(fu, string(sd//'/site_halves.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( word(ln, 1) /= 'accept' ) cycle
            nw = nwords(ln)
            do i = 1, size(ends)
                if( ends(i)%to_char() == word(ln, 2) ) z(i, merge(1, 2, word(ln, 3) == 'row')) = str2real(word(ln, nw - 1))
            end do
        end do
        call fclose(fu)
        best_minz = maxval(min(z(:,1), z(:,2)))
    end function best_minz

    real function curve_peak( fname )
        character(len=*), intent(in) :: fname
        character(len=256) :: ln
        integer :: fu, io, nn
        real    :: c
        curve_peak = -1.
        if( .not. file_exists(string(fname)) ) return
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
            read(ln,*) nn, c
            curve_peak = max(curve_peak, c)
        end do
        call fclose(fu)
    end function curve_peak

    logical function is_anc( a, u )
        character(len=*), intent(in) :: a, u
        is_anc = len(a) < len(u)
        if( is_anc ) is_anc = u(1:len(a)) == a
    end function is_anc

    logical function reentered_line( root, id )
        character(len=*), intent(in) :: root, id
        integer :: l
        reentered_line = .false.
        do l = len(id), 1, -1
            if( file_exists(string(root//'/node_'//id(1:l)//'/reentered.txt')) ) reentered_line = .true.
        end do
    end function reentered_line

    function parent_of( id ) result( p )
        character(len=*), intent(in) :: id
        character(len=:), allocatable :: p
        p = id(1:len(id)-1)
    end function parent_of

    function rstr( x ) result( s )
        real, intent(in) :: x
        character(len=:), allocatable :: s
        character(len=16) :: t
        write(t,'(F7.3)') x
        s = trim(adjustl(t))
    end function rstr

    function word( line, iw ) result( w )
        character(len=*), intent(in) :: line
        integer,          intent(in) :: iw
        character(len=:), allocatable :: w
        integer :: j, cnt, j0
        w   = ''
        cnt = 0
        j   = 1
        do while( j <= len_trim(line) )
            if( line(j:j) == ' ' )then
                j = j + 1
                cycle
            endif
            j0 = j
            do while( j <= len_trim(line) )
                if( line(j:j) == ' ' ) exit
                j = j + 1
            end do
            cnt = cnt + 1
            if( cnt == iw )then
                w = line(j0:j-1)
                return
            endif
        end do
    end function word

    integer function str2int( s )
        character(len=*), intent(in) :: s
        read(s,*) str2int
    end function str2int

    real function str2real( s )
        character(len=*), intent(in) :: s
        read(s,*) str2real
    end function str2real

    integer function nwords( line )
        character(len=*), intent(in) :: line
        integer :: j
        logical :: inw
        nwords = 0
        inw    = .false.
        do j = 1, len_trim(line)
            if( line(j:j) /= ' ' .and. .not. inw ) nwords = nwords + 1
            inw = line(j:j) /= ' '
        end do
    end function nwords

    function pad2( s ) result( str )
        integer, intent(in) :: s
        character(len=2) :: str
        write(str,'(I2.2)') s
    end function pad2

    function joined( vals, sep ) result( str )
        integer,                    intent(in) :: vals(:)
        character(len=*), optional, intent(in) :: sep
        character(len=:), allocatable :: str
        integer :: j
        str = ''
        do j = 1, size(vals)
            if( j > 1 )then
                if( present(sep) )then
                    str = str//sep
                else
                    str = str//' '
                endif
            endif
            str = str//int2str(vals(j))
        end do
    end function joined

    function failed_list( self ) result( str )
        class(hetero_tree), intent(in) :: self
        character(len=:), allocatable :: str
        integer :: j
        str = ''
        do j = 1, size(self%nodes)
            if( self%nodes(j)%step == S_FAILED ) str = str//' '//self%nodes(j)%id
        end do
    end function failed_list

    logical function live_step( s )
        integer, intent(in) :: s
        live_step = (s >= S_MEAN .and. s <= S_SPLIT) .or. s == S_MOTION .or. s == S_MBATCH .or. s == S_GATEW .or. s == S_SBATCH
    end function live_step

    function step_label( s ) result( lb )
        integer, intent(in) :: s
        character(len=:), allocatable :: lb
        select case(s)
            case(S_MEAN:S_SPLIT)
                lb = trim(STEP_DIRS(s))
            case(S_MOTION)
                lb = 'motion'
            case(S_MBATCH)
                lb = 'motion batch'
            case(S_GATEW)
                lb = 'gate floor wait'
            case(S_SBATCH)
                lb = 'site batch'
            case default
                lb = 'step '//int2str(s)
        end select
    end function step_label

    function step_dir( self, i, s ) result( sd )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i, s
        character(len=:), allocatable :: sd
        select case(s)
            case(S_SITE)
                sd = seed_dir(self%ndir(i))
            case(S_MOTION)
                sd = motion_dir(self%ndir(i), motion_cur(self%ndir(i)))
            case(S_MBATCH, S_GATEW)
                sd = self%ndir(i)//'/motion'
            case(S_SBATCH)
                sd = self%ndir(i)//'/sites'
            case(S_MEAN:S_RECB0, S_SEEDS:S_SPLIT)
                sd = self%ndir(i)//'/'//trim(STEP_DIRS(s))
            case default
                sd = self%ndir(i)
        end select
    end function step_dir

    function motion_dir( nd, c ) result( md )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        character(len=:), allocatable :: md
        md = nd//'/motion/c'//int2str(c)
    end function motion_dir

    function read_cands( nd ) result( cands )
        character(len=*), intent(in) :: nd
        integer, allocatable :: cands(:)
        allocate(cands(0))
        if( file_exists(string(nd//'/motion/cands.txt')) ) cands = read_ints(nd//'/motion/cands.txt')
    end function read_cands

    subroutine add_cand( nd, c )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        integer, allocatable :: cands(:)
        cands = read_cands(nd)
        if( any(cands == c) ) return
        call simple_mkdir(string(nd//'/motion'))
        call write_ints(nd//'/motion/cands.txt', [cands, c])
    end subroutine add_cand

    function motion_word( nd, c ) result( w )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        character(len=12) :: w
        w = word(first_line(motion_dir(nd, c)//'/motion_call.txt', '-'), 1)
    end function motion_word

    integer function motion_cur( nd )
        character(len=*), intent(in) :: nd
        integer, allocatable :: cands(:)
        integer :: ic
        motion_cur = -1
        cands = read_cands(nd)
        do ic = 1, size(cands)
            if( file_exists(string(motion_dir(nd, cands(ic))//'/motion_call.txt')) ) cycle
            motion_cur = cands(ic)
            return
        end do
    end function motion_cur

    function cand_name( c ) result( nm )
        integer, intent(in) :: c
        character(len=:), allocatable :: nm
        if( c > 0 )then
            nm = 'site '//int2str(c)
        else
            nm = 'the fflex fallback'
        endif
    end function cand_name

    function cand_label( nd, c ) result( nm )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        character(len=:), allocatable :: nm
        if( c > INH0 )then
            nm = 'inherited '//first_line(motion_dir(nd, c)//'/inherit.txt', '?')
        else
            nm = cand_name(c)
        endif
    end function cand_label

    function cand_origin( nd, c ) result( org )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        character(len=:), allocatable :: org
        if( c > INH0 )then
            org = 'inherited from '//first_line(motion_dir(nd, c)//'/inherit.txt', '?')
        else
            org = 'novel'
        endif
    end function cand_origin

    function group_rotations( pgrp ) result( rots )
        character(len=*), optional, intent(in) :: pgrp
        real, allocatable :: rots(:,:,:)
        real    :: th, c, s
        integer :: nsy, m, nr, r
        logical :: l_d
        allocate(rots(3,3,0))
        if( .not. present(pgrp) ) return
        if( len_trim(pgrp) < 2 ) return
        l_d = pgrp(1:1) == 'd'
        if( pgrp(1:1) /= 'c' .and. .not. l_d ) return
        nsy = nsym_of(trim(pgrp))
        nr  = nsy - 1 + merge(2 * nsy, 0, l_d)
        if( nr < 1 ) return
        deallocate(rots)
        allocate(rots(3,3,nr), source=0.)
        r = 0
        do m = 1, nsy - 1
            r  = r + 1
            th = 2. * PI * real(m) / real(nsy)
            c  = cos(th)
            s  = sin(th)
            rots(:,:,r) = reshape([c, s, 0., -s, c, 0., 0., 0., 1.], [3,3])
        end do
        if( .not. l_d ) return
        do m = 0, nsy - 1
            th = 2. * PI * real(m) / real(nsy)
            c  = cos(th)
            s  = sin(th)
            r  = r + 1
            rots(:,:,r) = matmul(reshape([1., 0., 0., 0., -1., 0., 0., 0., -1.], [3,3]), &
                &reshape([c, s, 0., -s, c, 0., 0., 0., 1.], [3,3]))
            r  = r + 1
            rots(:,:,r) = matmul(reshape([-1., 0., 0., 0., 1., 0., 0., 0., -1.], [3,3]), &
                &reshape([c, s, 0., -s, c, 0., 0., 0., 1.], [3,3]))
        end do
    end function group_rotations

    real function mate_dice( a, b, rots )
        real, intent(in) :: a(:,:,:), b(:,:,:), rots(:,:,:)
        integer, allocatable :: vox(:,:)
        real    :: cc(3), q(3), ctr
        integer :: n, h, k, l, r, t, na, nb, nhit, hh, kk, ll
        mate_dice = 0.
        if( size(rots,3) == 0 ) return
        n  = size(a,1)
        na = count(a > 0.5)
        nb = count(b > 0.5)
        if( na == 0 .or. nb == 0 ) return
        if( 2. * real(min(na, nb)) / real(na + nb) < DMATE ) return
        allocate(vox(3,na))
        t = 0
        do l = 1, n
            do k = 1, n
                do h = 1, n
                    if( a(h,k,l) <= 0.5 ) cycle
                    t = t + 1
                    vox(:,t) = [h, k, l]
                end do
            end do
        end do
        ctr = real(n / 2 + 1)
        do r = 1, size(rots,3)
            nhit = 0
            do t = 1, na
                cc = real(vox(:,t)) - ctr
                q  = matmul(rots(:,:,r), cc)
                hh = nint(q(1) + ctr)
                kk = nint(q(2) + ctr)
                ll = nint(q(3) + ctr)
                if( min(hh, kk, ll) < 1 .or. max(hh, kk, ll) > n ) cycle
                if( b(hh,kk,ll) > 0.5 ) nhit = nhit + 1
            end do
            mate_dice = max(mate_dice, 2. * real(min(nhit, nb)) / real(na + nb))
        end do
    end function mate_dice

    logical function off_open( nd )
        character(len=*), intent(in) :: nd
        off_open = file_exists(string(nd//'/motion/offpath.txt')) .and. .not. file_exists(string(nd//'/motion/offpath_done.txt'))
    end function off_open

    logical function floor_node( id )
        character(len=*), intent(in) :: id
        floor_node = id == 'r' .or. id == 'r1' .or. id == 'r2'
    end function floor_node

    integer(kind=8) function gate_seed( id )
        character(len=*), intent(in) :: id
        character(len=16) :: ev
        integer :: ios, es
        es = 0
        call get_environment_variable('SIMPLE_SEED', ev, status=ios)
        if( ios == 0 ) read(ev,*,iostat=ios) es
        if( ios /= 0 ) es = 0
        gate_seed = mod(1_8 + int(id_seed(id), 8) + 7919_8 * int(abs(es), 8), 2147483646_8) + 1_8
    end function gate_seed

    function lcg_perm( n, seed ) result( perm )
        integer,         intent(in)    :: n
        integer(kind=8), intent(inout) :: seed
        integer, allocatable :: perm(:)
        integer :: a, b, t
        perm = [(a, a=1,n)]
        do a = n, 2, -1
            seed    = mod(seed * 48271_8, 2147483647_8)
            b       = 1 + int(mod(seed, int(a, 8)))
            t       = perm(a)
            perm(a) = perm(b)
            perm(b) = t
        end do
    end function lcg_perm

    subroutine write_line( fname, text )
        character(len=*), intent(in) :: fname, text
        integer :: fu, io
        call fopen(fu, string(fname), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write '//fname, io)
        write(fu,'(A)') text
        call fclose(fu)
    end subroutine write_line

    subroutine write_coords( fname, hdr, lines )
        character(len=*), intent(in) :: fname
        type(string),     intent(in) :: hdr(:), lines(:)
        integer :: fu, io, j
        call fopen(fu, string(fname), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write '//fname, io)
        do j = 1, size(hdr)
            write(fu,'(A)') hdr(j)%to_char()
        end do
        do j = 1, size(lines)
            write(fu,'(A)') lines(j)%to_char()
        end do
        call fclose(fu)
    end subroutine write_coords

    subroutine write_lines( fname, lines )
        character(len=*), intent(in) :: fname
        type(string),     intent(in) :: lines(:)
        integer :: fu, io, j
        call fopen(fu, string(fname), status='replace', action='write', iostat=io)
        call fileiochk('hetero_tree: cannot write '//fname, io)
        do j = 1, size(lines)
            write(fu,'(A)') lines(j)%to_char()
        end do
        call fclose(fu)
    end subroutine write_lines

    function cand_mask( nd, kk ) result( mf )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: kk
        character(len=:), allocatable :: mf, sh
        character(len=8192) :: ln
        integer :: fu, io, j0, j1
        mf = nd//'/sites/site_'//pad2(kk)//'.mrc'
        if( file_exists(string(mf)) ) return
        mf = nd//'/mask.mrc'
        sh = site_dir_of(nd, kk)//'/hetero_step.sh'
        if( .not. file_exists(string(sh)) ) return
        call fopen(fu, string(sh), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            j0 = index(ln, 'pcg_mskfile=')
            if( j0 == 0 ) cycle
            j0 = j0 + len('pcg_mskfile=')
            j1 = index(ln(j0:), ' ')
            if( j1 == 0 ) j1 = len_trim(ln(j0:)) + 1
            mf = ln(j0:j0+j1-2)
            exit
        end do
        call fclose(fu)
    end function cand_mask

    function post_groups( p ) result( g )
        real, intent(in) :: p(:)
        integer, allocatable :: g(:)
        logical, allocatable :: ok(:)
        real    :: q1, q4
        integer :: nv
        allocate(g(size(p)), source=0)
        ok = p >= 0.
        nv = count(ok)
        if( nv < 2 * GATE_MINGRP ) return
        if( 2 * count(ok .and. abs(p - 0.5) > 0.45) > nv )then
            where( ok .and. p <  0.5 ) g = 1
            where( ok .and. p >= 0.5 ) g = 2
        else
            q1 = quant_of(pack(p, ok), 0.2)
            q4 = quant_of(pack(p, ok), 0.8)
            where( ok .and. p <= q1 ) g = 1
            where( ok .and. p >= q4 ) g = 2
        endif
    end function post_groups

    function groups4( ga, gb ) result( grp )
        integer, intent(in)  :: ga(:), gb(:)
        integer, allocatable :: grp(:)
        integer :: m
        allocate(grp(0))
        if( size(ga) == 0 .or. size(ga) /= size(gb) ) return
        grp = ga
        where( gb == 1 ) grp = 3
        where( gb == 2 ) grp = 4
        if( any([(count(grp == m) < GATE_MINGRP, m=1,4)]) )then
            deallocate(grp)
            allocate(grp(0))
        endif
    end function groups4

    function perm_groups( hr, seed0 ) result( grp )
        integer,         intent(in) :: hr(:)
        integer(kind=8), intent(in) :: seed0
        integer, allocatable :: grp(:), idx(:), perm(:)
        integer(kind=8) :: seed
        integer :: h, n5, m
        allocate(grp(size(hr)), source=0)
        seed = mod(seed0, 2147483646_8) + 1_8
        do h = 1, 2
            idx = pack([(m, m=1,size(hr))], hr == h)
            n5  = size(idx) / 5
            if( n5 < GATE_MINGRP )then
                deallocate(grp)
                allocate(grp(0))
                return
            endif
            perm = lcg_perm(size(idx), seed)
            grp(idx(perm(1:n5)))      = 2 * h - 1
            grp(idx(perm(n5+1:2*n5))) = 2 * h
        end do
    end function perm_groups

    function anc_index( ad, nd ) result( mx )
        character(len=*), intent(in) :: ad, nd
        integer, allocatable :: mx(:), st(:), rn(:), ra(:), pos(:)
        integer :: j
        allocate(mx(0))
        st = project_states(nd//'/'//NODE_PROJ)
        if( local_rows(nd) )then
            rn = read_ints(nd//'/rows.txt')
        else
            rn = [(j, j=1,size(st))]
        endif
        if( size(rn) /= size(st) .or. size(rn) == 0 ) return
        if( local_rows(ad) )then
            ra = read_ints(ad//'/rows.txt')
            allocate(pos(max(maxval([0, ra]), maxval(rn))), source=0)
            do j = 1, size(ra)
                pos(ra(j)) = j
            end do
            mx = [(pos(rn(j)), j=1,size(rn))]
        else
            mx = rn
        endif
        where( st <= 0 ) mx = 0
    end function anc_index

    real function gate_corr( rd, mskdiam )
        character(len=*), intent(in) :: rd
        real,             intent(in) :: mskdiam
        type(image) :: img
        real,    allocatable :: a(:,:,:), b(:,:,:), msk(:,:,:)
        logical, allocatable :: sel(:,:,:)
        integer :: ldim(3), ifoo, n, nh, h, k, l, j
        real    :: smpd, r, rr
        gate_corr = -2.
        do j = 1, 4
            if( .not. file_exists(string(rd//'/recvol_state'//pad2(j)//'.mrc')) ) return
        end do
        call find_ldim_nptcls(string(rd//'/recvol_state01.mrc'), ldim, ifoo)
        n    = ldim(1)
        nh   = n / 2
        smpd = find_img_smpd(string(rd//'/recvol_state01.mrc'))
        call img%new([n,n,n], smpd)
        rr = mskdiam / 2. / smpd
        allocate(msk(n,n,n), source=0.)
        do l = 1, n
            do k = 1, n
                do h = 1, n
                    r = sqrt(real((h-1-nh)**2 + (k-1-nh)**2 + (l-1-nh)**2))
                    if( r <= rr )then
                        msk(h,k,l) = 1.
                    else if( r <= rr + GATE_EDGE )then
                        msk(h,k,l) = 0.5 * (1. + cos(PI * (r - rr) / GATE_EDGE))
                    endif
                end do
            end do
        end do
        a = (vol(2) - vol(1)) * msk
        b = (vol(4) - vol(3)) * msk
        call lowpass(a)
        call lowpass(b)
        sel = msk > 0.5
        gate_corr = pearson(pack(a, sel), pack(b, sel))
        call img%kill

    contains

        function vol( jj ) result( v )
            integer, intent(in) :: jj
            real, allocatable :: v(:,:,:)
            call img%zero_and_unflag_ft
            call img%read(string(rd//'/recvol_state'//pad2(jj)//'.mrc'))
            v = img%get_rmat()
        end function vol

        subroutine lowpass( x )
            real, allocatable, intent(inout) :: x(:,:,:)
            complex, allocatable :: cm(:,:,:)
            real    :: f, f0, w
            integer :: h_, k_, l_, ky, lz
            f0 = smpd / GATE_LP
            w  = 0.1 * f0
            call img%zero_and_unflag_ft
            call img%set_rmat(x, .false.)
            call img%fft
            cm = img%get_cmat()
            do l_ = 1, size(cm,3)
                lz = l_ - 1 - merge(n, 0, l_ > nh)
                do k_ = 1, size(cm,2)
                    ky = k_ - 1 - merge(n, 0, k_ > nh)
                    do h_ = 1, size(cm,1)
                        f = sqrt(real((h_-1)**2 + ky*ky + lz*lz)) / real(n)
                        if( f > f0 + w )then
                            cm(h_,k_,l_) = cmplx(0., 0.)
                        else if( f > f0 )then
                            cm(h_,k_,l_) = cm(h_,k_,l_) * (0.5 * (1. + cos(PI * (f - f0) / w)))
                        endif
                    end do
                end do
            end do
            call img%set_cmat(cm)
            call img%ifft
            x = img%get_rmat()
        end subroutine lowpass

    end function gate_corr

    real function kept_margin( sd )
        character(len=*), intent(in) :: sd
        character(len=4096) :: ln
        integer :: fu, io, nw
        kept_margin = -huge(1.)
        if( .not. file_exists(string(sd//'/site_pole_rule.txt')) ) return
        call fopen(fu, string(sd//'/site_pole_rule.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            nw = nwords(ln)
            if( ln(1:1) == '#' .or. nw < 6 ) cycle
            if( word(ln, nw) /= 'kept' ) cycle
            kept_margin = max(kept_margin, min(str2real(word(ln, 4)), str2real(word(ln, 5))))
        end do
        call fclose(fu)
    end function kept_margin

    integer function strongest_site( nd, ks )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: ks(:)
        real    :: best, m
        integer :: k
        strongest_site = ks(1)
        best = -huge(best)
        do k = 1, size(ks)
            m = kept_margin(site_dir_of(nd, ks(k)))
            if( m > best )then
                best = m
                strongest_site = ks(k)
            endif
        end do
    end function strongest_site

    function kept_margins( nd, ks ) result( str )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: ks(:)
        character(len=:), allocatable :: str
        integer :: k
        str = ''
        do k = 1, size(ks)
            str = str//merge(', ', '  ', k > 1)//'site '//int2str(ks(k))//' '//rstr(kept_margin(site_dir_of(nd, ks(k))))
        end do
        str = trim(adjustl(str))
    end function kept_margins

    function kept_sites( nd ) result( ks )
        character(len=*), intent(in) :: nd
        integer, allocatable :: ks(:)
        integer :: k
        allocate(ks(0))
        do k = 1, max(1, nsites(nd))
            if( maxval([0, read_ints(site_dir_of(nd, k)//'/site_seeds.txt')]) >= 2 ) ks = [ks, k]
        end do
    end function kept_sites

    function kept_ends( sd ) result( ends )
        character(len=*), intent(in) :: sd
        type(string), allocatable :: ends(:)
        integer :: q
        call kept_pole(sd, q, ends)
    end function kept_ends

    subroutine kept_pole( sd, q, ends )
        character(len=*),          intent(in)  :: sd
        integer,                   intent(out) :: q
        type(string), allocatable, intent(out) :: ends(:)
        type(string), allocatable :: pe(:)
        integer,      allocatable :: sds(:), pid(:)
        integer :: w, nbest
        call site_pole_ends(sd, pe, pid)
        allocate(ends(0))
        q = 0
        if( size(pid) == 0 ) return
        sds   = read_ints(sd//'/site_seeds.txt')
        q     = pid(1)
        nbest = -1
        do w = 1, size(pid)
            if( count(sds == pid(w)) <= nbest ) cycle
            nbest = count(sds == pid(w))
            q     = pid(w)
        end do
        ends = pack(pe, pid == q)
    end subroutine kept_pole

    subroutine site_pole_ends( sd, pe, pid )
        character(len=*),          intent(in)  :: sd
        type(string), allocatable, intent(out) :: pe(:)
        integer,      allocatable, intent(out) :: pid(:)
        character(len=4096) :: ln
        integer :: fu, io, w, nw
        allocate(pe(0), pid(0))
        if( .not. file_exists(string(sd//'/site_poles.txt')) ) return
        call fopen(fu, string(sd//'/site_poles.txt'), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            nw = nwords(ln)
            if( ln(1:1) == '#' .or. nw < 3 ) cycle
            if( word(ln, 3) == 'rest' ) cycle
            do w = 3, nw
                pe  = [pe, string(word(ln, w))]
                pid = [pid, str2int(word(ln, 1))]
            end do
        end do
        call fclose(fu)
    end subroutine site_pole_ends

    function opp_end( e ) result( o )
        character(len=*), intent(in) :: e
        character(len=:), allocatable :: o
        o = ''
        if( len(e) < 3 ) return
        if( e(len(e)-1:) == 'lo' ) o = e(:len(e)-2)//'hi'
        if( e(len(e)-1:) == 'hi' ) o = e(:len(e)-2)//'lo'
    end function opp_end

    function mean_post( sd, ends ) result( p )
        character(len=*), intent(in) :: sd
        type(string),     intent(in) :: ends(:)
        real, allocatable :: p(:), v(:)
        integer :: e
        allocate(p(0))
        do e = 1, size(ends)
            v = read_reals(sd//'/site_post_'//ends(e)%to_char()//'.txt')
            if( size(v) == 0 )then
                deallocate(p)
                allocate(p(0))
                return
            endif
            if( e == 1 )then
                p = v / real(size(ends))
            else if( size(v) == size(p) )then
                p = p + v / real(size(ends))
            endif
        end do
    end function mean_post

    integer function pole_with_ends( sd, ends )
        character(len=*), intent(in) :: sd
        type(string),     intent(in) :: ends(:)
        type(string), allocatable :: pe(:)
        integer,      allocatable :: pid(:)
        integer :: w, u, v, q
        logical :: ok, found
        pole_with_ends = 0
        call site_pole_ends(sd, pe, pid)
        do w = 1, size(pid)
            q = pid(w)
            if( count(pid == q) /= size(ends) ) cycle
            ok = .true.
            do v = 1, size(ends)
                found = .false.
                do u = 1, size(pid)
                    if( pid(u) == q .and. pe(u)%to_char() == ends(v)%to_char() ) found = .true.
                end do
                ok = ok .and. found
            end do
            if( ok )then
                pole_with_ends = q
                return
            endif
        end do
    end function pole_with_ends

    logical function two_pole_seeds( nd, ka, minpop, ab, st0, msg, rule )
        character(len=*),                        intent(in)  :: nd
        integer,                                 intent(in)  :: ka, minpop
        integer,          allocatable,           intent(out) :: ab(:), st0(:)
        character(len=:), allocatable,           intent(out) :: msg
        character(len=:), allocatable, optional, intent(out) :: rule
        type(string), allocatable :: ea(:), eb(:)
        real,         allocatable :: pa(:), pb(:), sc(:)
        integer,      allocatable :: sa(:), sb(:), ks(:), ord(:)
        character(len=:), allocatable :: how, sda
        real    :: best, mk
        integer :: qa, qb, kb, k, q, e, na
        logical :: lopp
        two_pole_seeds = .false.
        msg = ''
        how = ''
        sda = site_dir_of(nd, ka)
        call kept_pole(sda, qa, ea)
        if( qa == 0 ) return
        sa = read_ints(sda//'/site_seeds.txt')
        pa = mean_post(sda, ea)
        if( size(pa) /= size(sa) ) return
        allocate(eb(size(ea)))
        lopp = .true.
        do e = 1, size(ea)
            eb(e) = string(opp_end(ea(e)%to_char()))
            if( len(eb(e)%to_char()) == 0 ) lopp = .false.
            if( any([(eb(e)%to_char() == ea(k)%to_char(), k=1,size(ea))]) ) lopp = .false.
        end do
        allocate(ab(size(sa)), source=0)
        where( sa == qa ) ab = 1
        allocate(pb(0))
        if( lopp )then
            kb = 0
            qb = 0
            best = -huge(best)
            ks = kept_sites(nd)
            do k = 1, size(ks)
                q = pole_with_ends(site_dir_of(nd, ks(k)), eb)
                if( q == 0 ) cycle
                mk = kept_margin(site_dir_of(nd, ks(k)))
                if( mk <= best ) cycle
                best = mk
                kb   = ks(k)
                qb   = q
            end do
            if( kb > 0 )then
                sb = read_ints(site_dir_of(nd, kb)//'/site_seeds.txt')
                pb = mean_post(site_dir_of(nd, kb), eb)
                if( size(sb) == size(sa) ) where( ab == 0 .and. sa > 0 .and. sb == qb ) ab = 2
                how = 'the kept pole of site '//int2str(kb)//' at the opposite endpoints'
            else
                pb = mean_post(sda, eb)
                if( size(pb) == size(sa) ) where( ab == 0 .and. sa > 0 .and. pb >= 0.5 ) ab = 2
                how = 'the opposite endpoints of the same site at posterior >= 0.5'
            endif
        endif
        if( size(pb) /= size(sa) ) pb = 1. - pa
        if( count(ab == 2) < max(minpop, PROB_MULTISTATE_MIN_POP + 1) )then
            ! no opposite end holds minpop: B mirrors A, the |A| rows lowest on the kept pole's score
            where( ab == 2 ) ab = 0
            na  = count(ab == 1)
            ord = pack([(k, k=1,size(sa))], ab == 0 .and. sa > 0)
            if( na == 0 .or. size(ord) < 2 * na ) return
            sc = pa(ord)
            call hpsort(sc, ord)
            ab(ord(1:na)) = 2
            how = 'mirrored, the '//int2str(na)//' rows lowest on the kept pole''s score (up to '//rstr(sc(na))//')'
        endif
        if( count(ab == 2) <= PROB_MULTISTATE_MIN_POP ) return
        st0 = ab
        where( st0 == 0 .and. sa > 0 .and. pa >= pb ) st0 = 1
        where( st0 == 0 .and. sa > 0 ) st0 = 2
        msg = 'pole A (state '//int2str(qa)//') '//int2str(count(ab == 1))//' seeds, the opposite end B ('//how//') '//&
            &int2str(count(ab == 2))//' seeds; the other '//int2str(count(ab == 0 .and. sa > 0))//' rows start at their nearer '//&
            &'end ('//int2str(count(st0 == 1 .and. ab == 0))//' / '//int2str(count(st0 == 2 .and. ab == 0))//') and stage A decides'
        if( present(rule) ) rule = how
        two_pole_seeds = .true.
    end function two_pole_seeds

    subroutine motion_factor( nd, c, l_sw, minpop, lv, cv, w )
        character(len=*),              intent(in)  :: nd
        integer,                       intent(in)  :: c, minpop
        logical,                       intent(in)  :: l_sw
        integer,          allocatable, intent(out) :: lv(:), cv(:)
        character(len=:), allocatable, intent(out) :: w
        integer, allocatable :: ab(:), st0(:), jk(:)
        character(len=:), allocatable :: msg, tag, brule
        integer :: k
        tag   = ''
        brule = 'as inherited'
        if( c > INH0 )then
            ab  = read_ints(motion_dir(nd, c)//'/inherit_ab.txt')
            st0 = read_ints(motion_dir(nd, c)//'/inherit_st0.txt')
            tag = 'inherited '
        else if( c == 0 )then
            lv = two_groups(nd)
            jk = junk_rows(nd)
            if( size(jk) > 0 .and. size(lv) > 0 ) lv(jk) = 0
            cv = lv
            w  = 'two fflex groups '//joined([count(lv == 1), count(lv == 2)], ' / ')
            call save_levels(lv, lv)
            return
        else if( .not. two_pole_seeds(nd, c, minpop, ab, st0, msg, brule) )then
            lv = cand_seeds(nd, c)
            cv = lv
            w  = 'no opposite end and no room to mirror A, its site seeds '//joined([(count(lv == k), k=1,maxval([0, lv]))], ' / ')
            call save_levels(lv, lv)
            return
        endif
        if( l_sw )then
            lv = st0
            where( ab == 1 ) lv = 1
            where( ab == 0 .and. st0 > 0 ) lv = 2
            where( ab == 2 ) lv = 3
            cv = lv
            w  = tag//'swing A / middle / B '//joined([count(lv == 1), count(lv == 2), count(lv == 3)], ' / ')//'; B: '//brule
        else
            lv = st0
            cv = ab
            w  = tag//'switch A / B seeds '//joined([count(ab == 1), count(ab == 2)], ' / ')//', start '//&
                &joined([count(st0 == 1), count(st0 == 2)], ' / ')//'; B: '//brule
        endif
        if( c <= INH0 ) call save_levels(ab, st0)

    contains

        subroutine save_levels( a, b )
            integer, intent(in) :: a(:), b(:)
            call write_ints(motion_dir(nd, c)//'/ab.txt', a)
            call write_ints(motion_dir(nd, c)//'/st0.txt', b)
        end subroutine save_levels

    end subroutine motion_factor

    function cand_score( nd, c ) result( sc )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        real,         allocatable :: sc(:), v(:)
        integer,      allocatable :: st(:), lab(:), jk(:)
        type(string), allocatable :: ends(:)
        character(len=:), allocatable :: sd
        integer :: e
        st = project_states(nd//'/'//NODE_PROJ)
        if( c > INH0 )then
            sc = read_reals(motion_dir(nd, c)//'/inherit_score.txt')
            if( size(sc) /= size(st) ) THROW_HARD('inherited score does not match the node project: '//nd)
            jk = junk_rows(nd)
        else if( c > 0 )then
            sd   = site_dir_of(nd, c)
            ends = kept_ends(sd)
            allocate(sc(size(st)), source=-1.)
            if( size(ends) > 0 ) sc = 0.
            do e = 1, size(ends)
                v = read_reals(sd//'/site_post_'//ends(e)%to_char()//'.txt')
                if( size(v) /= size(st) ) THROW_HARD('site posteriors do not match the node project: '//sd)
                sc = sc + v / real(size(ends))
            end do
            jk = junk_of(sd)
        else
            lab = two_groups(nd)
            jk  = junk_rows(nd)
            if( size(jk) > 0 .and. size(lab) > 0 ) lab(jk) = 0
            sc  = fisher_score(nd//'/fflex/flex_pca_coordinates.txt', lab, size(st))
        endif
        if( size(jk) > 0 ) sc(jk) = -1.
        where( st <= 0 ) sc = -1.
    end function cand_score

    function cand_seeds( nd, c ) result( sds )
        character(len=*), intent(in) :: nd
        integer,          intent(in) :: c
        integer, allocatable :: sds(:), jk(:)
        if( c > 0 )then
            sds = read_ints(site_dir_of(nd, c)//'/site_seeds.txt')
        else
            sds = two_groups(nd)
            jk  = junk_rows(nd)
            if( size(jk) > 0 .and. size(sds) > 0 ) sds(jk) = 0
        endif
    end function cand_seeds

    function fisher_score( fname, lab, nrows ) result( sc )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: lab(:), nrows
        real,     allocatable :: sc(:), key(:)
        real(dp), allocatable :: z(:,:), mu(:,:), sw(:,:), w(:)
        integer,  allocatable :: pinds(:), g(:), ord(:)
        integer :: n, q, a, k, r
        allocate(sc(nrows), source=-1.)
        call flex_coords(fname, pinds, z)
        n = size(pinds)
        q = size(z, 2)
        if( n < 2 .or. q < 1 ) return
        g = lab(pinds)
        if( count(g == 1) < 2 .or. count(g == 2) < 2 ) return
        allocate(mu(q,2), sw(q,q), source=0._dp)
        do a = 1, 2
            mu(:,a) = [(sum(z(:,k), mask=g == a) / real(count(g == a), dp), k=1,q)]
            do k = 1, q
                do r = 1, q
                    sw(k,r) = sw(k,r) + sum((z(:,k) - mu(k,a)) * (z(:,r) - mu(r,a)), mask=g == a) * &
                        &real(count(g == a), dp) / real(count(g == a) - 1, dp)
                end do
            end do
        end do
        sw = sw / real(count(g == 1) + count(g == 2), dp)
        do k = 1, q
            sw(k,k) = sw(k,k) + 1.e-9_dp
        end do
        w   = lin_solve(sw, mu(:,1) - mu(:,2))
        key = real(matmul(z, w))
        ord = [(r, r=1,n)]
        call hpsort(key, ord)
        do r = 1, n
            sc(pinds(ord(r))) = (real(r) - 0.5) / real(n)
        end do
    end function fisher_score

    function lin_solve( a, b ) result( x )
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), allocatable :: x(:), m(:,:), rowt(:)
        real(dp) :: f, xt
        integer  :: n, k, r, p
        n = size(b)
        m = a
        x = b
        do k = 1, n
            p = k - 1 + maxloc(abs(m(k:n,k)), dim=1)
            if( p /= k )then
                rowt   = m(k,:)
                m(k,:) = m(p,:)
                m(p,:) = rowt
                xt     = x(k)
                x(k)   = x(p)
                x(p)   = xt
            endif
            if( abs(m(k,k)) < 1.e-300_dp ) cycle
            do r = k + 1, n
                f        = m(r,k) / m(k,k)
                m(r,k:n) = m(r,k:n) - f * m(k,k:n)
                x(r)     = x(r) - f * x(k)
            end do
        end do
        do k = n, 1, -1
            if( abs(m(k,k)) < 1.e-300_dp )then
                x(k) = 0._dp
            else
                x(k) = (x(k) - sum(m(k,k+1:n) * x(k+1:n))) / m(k,k)
            endif
        end do
    end function lin_solve

    function quint_bins( sc ) result( bins )
        real, intent(in) :: sc(:)
        integer, allocatable :: bins(:), idx(:), ord(:)
        real,    allocatable :: key(:)
        integer :: n, r, q
        allocate(bins(size(sc)), source=0)
        idx = pack([(r, r=1,size(sc))], sc >= 0.)
        n   = size(idx)
        if( n < NMBIN ) return
        key = sc(idx)
        ord = [(r, r=1,n)]
        call hpsort(key, ord)
        do r = 1, n
            bins(idx(ord(r))) = 1 + count([(r > (n * q) / NMBIN, q=1,NMBIN-1)])
        end do
    end function quint_bins

    function product_cells( sv, minpop ) result( cell )
        integer, intent(in)  :: sv(:,:), minpop
        integer, allocatable :: cell(:), keys(:,:), cnt(:), rk(:)
        integer :: n, m, nk, r, k, j, small, best, dist, bd
        n  = size(sv,1)
        m  = size(sv,2)
        allocate(cell(n), source=0)
        allocate(keys(m,0))
        nk = 0
        do r = 1, n
            if( any(sv(r,:) <= 0) ) cycle
            k = 0
            do j = 1, nk
                if( all(keys(:,j) == sv(r,:)) )then
                    k = j
                    exit
                endif
            end do
            if( k == 0 )then
                nk   = nk + 1
                keys = reshape([keys, sv(r,:)], [m, nk])
                k    = nk
            endif
            cell(r) = k
        end do
        cnt = [(count(cell == k), k=1,nk)]
        do
            small = 0
            do k = 1, nk
                if( cnt(k) == 0 .or. cnt(k) >= minpop ) cycle
                if( small == 0 )then
                    small = k
                else if( cnt(k) < cnt(small) )then
                    small = k
                endif
            end do
            if( small == 0 .or. count(cnt > 0) < 2 ) exit
            best = 0
            bd   = huge(bd)
            do k = 1, nk
                if( k == small .or. cnt(k) == 0 ) cycle
                dist = count(keys(:,k) /= keys(:,small))
                if( best == 0 )then
                    best = k
                    bd   = dist
                else if( dist < bd .or. (dist == bd .and. cnt(k) > cnt(best)) )then
                    best = k
                    bd   = dist
                endif
            end do
            where( cell == small ) cell = best
            cnt(best)  = cnt(best) + cnt(small)
            cnt(small) = 0
        end do
        allocate(rk(0:nk), source=0)
        do k = 1, nk
            if( cnt(k) > 0 ) rk(k) = 1 + count(cnt > cnt(k)) + count(cnt(1:k-1) == cnt(k))
        end do
        cell = rk(cell)
    end function product_cells

    integer function nsym_of( pgrp )
        character(len=*), intent(in) :: pgrp
        integer :: ios
        nsym_of = 1
        if( len_trim(pgrp) < 2 ) return
        if( pgrp(1:1) /= 'c' .and. pgrp(1:1) /= 'd' ) return
        read(pgrp(2:), *, iostat=ios) nsym_of
        if( ios /= 0 .or. nsym_of < 1 ) nsym_of = 1
    end function nsym_of

    subroutine node_motion_summary( nd, cands, how )
        character(len=*), intent(in) :: nd, how
        integer,          intent(in) :: cands(:)
        integer :: fu, io, ic
        call fopen(fu, string(nd//'/motion_call.txt'), status='replace', action='write', iostat=io)
        if( io /= 0 ) return
        do ic = 1, size(cands)
            write(fu,'(A)') cand_label(nd, cands(ic))//': '//first_line(motion_dir(nd, cands(ic))//'/motion_call.txt', '-')//&
                &' | gate '//first_line(motion_dir(nd, cands(ic))//'/verdict.txt', '-')
        end do
        write(fu,'(A)') 'hand-back: '//how
        call fclose(fu)
    end subroutine node_motion_summary

    real function pearson( a, b )
        real, intent(in) :: a(:), b(:)
        real(dp) :: ma, mb, sab, saa, sbb
        pearson = 0.
        if( size(a) < 2 ) return
        ma  = sum(real(a,dp)) / real(size(a),dp)
        mb  = sum(real(b,dp)) / real(size(b),dp)
        sab = sum((real(a,dp) - ma) * (real(b,dp) - mb))
        saa = sum((real(a,dp) - ma)**2)
        sbb = sum((real(b,dp) - mb)**2)
        if( saa > 0._dp .and. sbb > 0._dp ) pearson = real(sab / sqrt(saa * sbb))
    end function pearson

    real function quant_of( v, q )
        real, intent(in) :: v(:), q
        real,    allocatable :: w(:)
        integer, allocatable :: ord(:)
        real    :: pos
        integer :: lo, r
        quant_of = 0.
        if( size(v) == 0 ) return
        w   = v
        ord = [(r, r=1,size(w))]
        call hpsort(w, ord)
        pos = real(size(w) - 1) * q
        lo  = int(pos)
        quant_of = w(lo+1) + (pos - real(lo)) * (w(min(lo+2, size(w))) - w(lo+1))
    end function quant_of

    function root_axes( self, i, c ) result( str )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i, c
        character(len=:), allocatable :: str, rd, d
        real,    allocatable :: s0(:), s1(:)
        integer, allocatable :: rc(:), mrow(:)
        logical, allocatable :: ok(:)
        integer :: ic, l
        str = ''
        if( self%nodes(i)%depth == 0 ) return
        rd = self%root//'/node_r'
        rc = read_cands(rd)
        if( size(rc) == 0 ) return
        d  = self%ndir(i)
        s1 = read_reals(motion_dir(d, c)//'/score.txt')
        if( local_rows(d) )then
            mrow = read_ints(d//'/rows.txt')
        else
            mrow = [(l, l=1,size(s1))]
        endif
        if( size(mrow) /= size(s1) ) return
        str = 'r with the root axes:'
        do ic = 1, size(rc)
            if( .not. file_exists(string(motion_dir(rd, rc(ic))//'/score.txt')) ) cycle
            s0 = read_reals(motion_dir(rd, rc(ic))//'/score.txt')
            if( maxval(mrow) > size(s0) ) cycle
            ok = s1 >= 0. .and. s0(mrow) >= 0.
            if( count(ok) < 3 ) cycle
            str = str//' '//cand_name(rc(ic))//' '//rstr(pearson(pack(s1, ok), pack(s0(mrow), ok)))
        end do
    end function root_axes

    subroutine motion_stats( md, maskf, nsym, cname, rootax, origin, line )
        character(len=*),              intent(in)  :: md, maskf, cname, rootax, origin
        integer,                       intent(in)  :: nsym
        character(len=:), allocatable, intent(out) :: line
        character(len=5), parameter :: HALF(0:2) = [character(len=5) :: '', '_even', '_odd']
        type(image)      :: img
        type(sp_project) :: sp
        complex,  allocatable :: cm(:,:,:), ce(:,:,:), qs(:,:,:), qf(:,:,:)
        real,     allocatable :: v(:,:,:), sm(:,:,:), smf(:,:,:), lpm(:,:,:), wlp(:,:,:), wg(:,:,:), dr(:), ndr(:), ad(:)
        real,     allocatable :: e1(:), e2(:), dfx(:)
        real(dp), allocatable :: fee(:,:), foo(:,:), feo(:,:), bw(:)
        integer,  allocatable :: ish(:,:,:), bidx(:,:), bsh(:), mvox(:,:), st(:)
        logical,  allocatable :: inm(:,:,:), foc(:), top(:)
        real(dp) :: rss(2:4), nbs(2:4), zrs(2:4), bts(2:4), rsf(2:4), nbf(2:4), zrf(2:4), btf(2:4), esns, zens, esnf, zenf
        real(dp) :: zloc, relloc, cosv, hmin, rel
        real     :: smpd, worst, rk(NMBIN), eo(NMBIN), rn(4), share, pmr, focf, cen(3), vexc, dexc, sdn, fc, f, wv, thr
        integer  :: ldim(3), mdim(3), ifoo, n, nh, nb, nmk, nfoc, smax, kk, j, h, k, l, ky, lz, m, t
        character(len=:), allocatable :: callw, flags, fname
        call find_ldim_nptcls(string(md//'/recvol_state01.mrc'), ldim, ifoo)
        call find_ldim_nptcls(string(maskf), mdim, ifoo)
        n    = ldim(1)
        nh   = n / 2
        smpd = find_img_smpd(string(md//'/recvol_state01.mrc'))
        if( any(mdim /= ldim) )then
            line = cname//': the node mask is not at the box of the bin maps, no call'
            call write_call('insufficient', line)
            return
        endif
        call img%new([n,n,n], smpd)
        call img%read(string(maskf))
        inm  = img%get_rmat() > 0.5
        nmk  = count(inm)
        allocate(mvox(3,nmk))
        t = 0
        do l = 1, n
            do k = 1, n
                do h = 1, n
                    if( .not. inm(h,k,l) ) cycle
                    t = t + 1
                    mvox(:,t) = [h, k, l]
                end do
            end do
        end do
        call img%set_rmat(merge(1., 0., inm), .false.)
        call img%fft
        cm = img%get_cmat()
        allocate(ish(size(cm,1),size(cm,2),size(cm,3)), wlp(size(cm,1),size(cm,2),size(cm,3)), &
            &wg(size(cm,1),size(cm,2),size(cm,3)))
        fc = smpd / MRES
        do l = 1, size(cm,3)
            lz = l - 1 - merge(n, 0, l > nh)
            do k = 1, size(cm,2)
                ky = k - 1 - merge(n, 0, k > nh)
                do h = 1, size(cm,1)
                    f  = sqrt(real((h-1)**2 + ky*ky + lz*lz))
                    ish(h,k,l) = nint(f)
                    f  = f / real(n)
                    wv = min(max((1.1 * fc - f) / (0.2 * fc), 0.), 1.)
                    wlp(h,k,l) = 0.5 - 0.5 * cos(real(PI) * wv)
                    wg(h,k,l)  = exp(-8. * real(PI)**2 * f * f)
                end do
            end do
        end do
        call img%set_cmat(cm * wg)
        call img%ifft
        sm   = img%get_rmat()
        smax = int(real(n) * smpd / MRES)
        nb   = count(ish >= 1 .and. ish <= smax)
        allocate(bidx(3,nb), bsh(nb), bw(nb))
        m = 0
        do l = 1, size(cm,3)
            do k = 1, size(cm,2)
                do h = 1, size(cm,1)
                    if( ish(h,k,l) < 1 .or. ish(h,k,l) > smax ) cycle
                    m = m + 1
                    bidx(:,m) = [h, k, l]
                    bsh(m)    = ish(h,k,l)
                    bw(m)     = merge(1._dp, 2._dp, h == 1 .or. h == size(cm,1))
                end do
            end do
        end do
        allocate(lpm(nmk,0:2,NMBIN), qs(nb,0:2,NMBIN), qf(nb,0:2,NMBIN))
        allocate(fee(nh,NMBIN), foo(nh,NMBIN), feo(nh,NMBIN), source=0._dp)
        do kk = 1, NMBIN
            do j = 0, 2
                call read_map(kk, j)
                call img%fft
                call img%set_cmat(img%get_cmat() * wlp)
                call img%ifft
                lpm(:,j,kk) = gather(img%get_rmat())
                call img%set_rmat(v * sm, .false.)
                call img%fft
                cm = img%get_cmat()
                qs(:,j,kk) = [(cm(bidx(1,m),bidx(2,m),bidx(3,m)), m=1,nb)]
                if( j == 1 ) ce = cm
                if( j == 2 ) call fsc_sums(kk)
            end do
        end do
        do kk = 1, NMBIN
            rk(kk) = 2. * smpd
            do k = 1, nh
                if( feo(k,kk) / sqrt(max(fee(k,kk) * foo(k,kk), 1.e-30_dp)) < 0.143_dp )then
                    rk(kk) = real(n) * smpd / real(k)
                    exit
                endif
            end do
        end do
        worst = maxval(rk)
        do kk = 1, NMBIN
            eo(kk) = 1. - pearson(lpm(:,1,kk), lpm(:,2,kk))
        end do
        do kk = 1, 4
            rn(kk) = (1. - pearson(lpm(:,0,kk), lpm(:,0,kk+1))) / max((eo(kk) + eo(kk+1)) / 4., 1.e-12)
        end do
        share = maxval(rn) / max(sum(rn), 1.e-30)
        dr    = lpm(:,0,NMBIN) - lpm(:,0,1)
        ndr   = ((lpm(:,1,NMBIN) - lpm(:,2,NMBIN)) - (lpm(:,1,1) - lpm(:,2,1))) / 2.
        sdn   = sqrt(sum((ndr - sum(ndr) / real(nmk))**2) / real(nmk))
        foc   = abs(dr) > 3. * sdn
        focf  = real(count(foc)) / real(nmk)
        if( real(count(foc)) < 0.005 * real(nmk) ) foc = abs(dr) >= quant_of(abs(dr), 0.95)
        nfoc  = count(foc)
        pmr   = real(sum(real(dr,dp)**2, mask=dr > 0.) / max(sum(real(dr,dp)**2, mask=dr < 0.), 1.e-30_dp))
        ad    = abs(dr)
        thr   = quant_of(ad, 0.99)
        top   = ad >= thr
        cen   = [(real(sum(real(mvox(k,:) - 1 - nh), mask=top)) / real(max(1, count(top))) * smpd, k=1,3)]
        v = 0.
        do t = 1, nmk
            if( foc(t) ) v(mvox(1,t),mvox(2,t),mvox(3,t)) = 1.
        end do
        call img%set_rmat(v, .false.)
        call img%fft
        call img%set_cmat(img%get_cmat() * wg)
        call img%ifft
        smf = img%get_rmat()
        do kk = 1, NMBIN
            do j = 0, 2
                call read_map(kk, j)
                call img%set_rmat(v * smf, .false.)
                call img%fft
                cm = img%get_cmat()
                qf(:,j,kk) = [(cm(bidx(1,m),bidx(2,m),bidx(3,m)), m=1,nb)]
            end do
        end do
        call img%kill
        call local_change
        call blend(qs, nmk, rss, nbs, zrs, bts, esns, zens)
        call blend(qf, nfoc, rsf, nbf, zrf, btf, esnf, zenf)
        call sp%read_segment('ptcl3D', string(md//'/m.simple'))
        st  = sp%os_ptcl3D%get_all_asint('state')
        e1  = sp%os_ptcl3D%get_all('e1')
        e2  = sp%os_ptcl3D%get_all('e2')
        dfx = sp%os_ptcl3D%get_all('dfx')
        call sp%kill
        call acquisition(vexc, dexc)
        flags = ''
        if( worst > MRES )then
            callw = 'insufficient'
        else if( zloc < ZLOC_MIN )then
            callw = 'none'
        else if( maxval(nbf) >= NB_SW .and. maxval(zrf) >= ZRES_SW )then
            callw = 'swing'
        else
            callw = 'switch'
            if( share < GRADED ) flags = flags//' graded'
            if( pmr > 1. / 3. .and. pmr < 3. ) flags = flags//' dipole'
        endif
        if( vexc >= VIEW_EXC .or. dexc >= DFX_EXC ) flags = flags//' acquisition'
        line = cname//' -> '//callw//flags//': worst-bin FSC 0.143 '//rstr(worst)//' A, z_loc '//rstr(real(zloc))//&
            &', nb_SF '//rstr(real(maxval(nbf)))//' at residual z '//rstr(real(maxval(zrf)))//', R_S '//&
            &rstr(real(maxval(rss)))//', step share '//rstr(share)//', +/- '//rstr(pmr)//', view / defocus excess '//&
            &rstr(vexc)//' / '//rstr(dexc)
        if( len(rootax) > 0 ) line = line//'; '//rootax
        line = line//'; '//origin
        call write_call(callw//flags, line)

    contains

        subroutine read_map( kk_, j_ )
            integer, intent(in) :: kk_, j_
            fname = md//'/recvol_state'//pad2(kk_)//trim(HALF(j_))//'.mrc'
            call img%zero_and_unflag_ft
            call img%read(string(fname))
            v = img%get_rmat()
        end subroutine read_map

        function gather( x ) result( y )
            real, intent(in) :: x(:,:,:)
            real, allocatable :: y(:)
            integer :: t_
            allocate(y(nmk))
            do t_ = 1, nmk
                y(t_) = x(mvox(1,t_),mvox(2,t_),mvox(3,t_))
            end do
        end function gather

        subroutine fsc_sums( kk_ )
            integer, intent(in) :: kk_
            real(dp) :: w_
            integer  :: h_, k_, l_, s_
            do l_ = 1, size(cm,3)
                do k_ = 1, size(cm,2)
                    do h_ = 1, size(cm,1)
                        s_ = ish(h_,k_,l_)
                        if( s_ < 1 .or. s_ > nh ) cycle
                        w_ = merge(1._dp, 2._dp, h_ == 1 .or. h_ == size(cm,1))
                        fee(s_,kk_) = fee(s_,kk_) + w_ * real(abs(ce(h_,k_,l_))**2, dp)
                        foo(s_,kk_) = foo(s_,kk_) + w_ * real(abs(cm(h_,k_,l_))**2, dp)
                        feo(s_,kk_) = feo(s_,kk_) + w_ * real(real(ce(h_,k_,l_) * conjg(cm(h_,k_,l_))), dp)
                    end do
                end do
            end do
        end subroutine fsc_sums

        real(dp) function cdot( x, y )
            complex, intent(in) :: x(:), y(:)
            cdot = sum(bw * (real(real(x),dp) * real(real(y),dp) + real(aimag(x),dp) * real(aimag(y),dp)))
        end function cdot

        subroutine local_change
            complex,  allocatable :: pm(:,:), dd(:,:), le(:), lo(:)
            real(dp), allocatable :: num(:), den(:), cs(:)
            real(dp) :: sig_m, sig_e, sig_l, rl, neff, hk
            integer  :: j_, k_, m_
            allocate(pm(nb,0:2), dd(nb,0:2))
            do j_ = 0, 2
                pm(:,j_) = sum(qs(:,j_,:), dim=2) / real(NMBIN)
                dd(:,j_) = qs(:,j_,NMBIN) - qs(:,j_,1)
            end do
            sig_m = cdot(pm(:,1), pm(:,2))
            sig_e = cdot(dd(:,1), dd(:,2))
            allocate(num(smax), den(smax), source=0._dp)
            do m_ = 1, nb
                num(bsh(m_)) = num(bsh(m_)) + bw(m_) * real(real(conjg(dd(m_,0)) * pm(m_,0)), dp)
                den(bsh(m_)) = den(bsh(m_)) + bw(m_) * real(abs(pm(m_,0))**2, dp)
            end do
            cs    = num / max(den, 1.e-30_dp)
            le    = dd(:,1) - real(cs(bsh)) * pm(:,1)
            lo    = dd(:,2) - real(cs(bsh)) * pm(:,2)
            sig_l = cdot(le, lo)
            rl    = sig_l / sqrt(max(cdot(le, le) * cdot(lo, lo), 1.e-300_dp))
            neff  = real(nmk,dp) * real(smpd,dp)**3 / (real(MRES,dp) / 2._dp)**3
            zloc  = rl * sqrt(neff)
            relloc = 0._dp
            rel    = 0._dp
            cosv   = 0._dp
            if( sig_m > 0._dp )then
                relloc = sig_l / sig_m
                rel    = sig_e / sig_m
            endif
            if( sig_e > 0._dp .and. sig_m > 0._dp ) cosv = cdot(dd(:,1), pm(:,2)) / sqrt(sig_e * sig_m)
            hmin = 1._dp
            do k_ = 1, NMBIN
                hk   = cdot(qs(:,1,k_), qs(:,2,k_)) / sqrt(max(cdot(qs(:,1,k_), qs(:,1,k_)) * cdot(qs(:,2,k_), qs(:,2,k_)), &
                    &1.e-300_dp))
                hmin = min(hmin, hk)
            end do
        end subroutine local_change

        subroutine blend( q, nvox, rr, nbv, zr, bt, esnr, zend )
            complex,  intent(in)  :: q(:,0:,:)
            integer,  intent(in)  :: nvox
            real(dp), intent(out) :: rr(2:4), nbv(2:4), zr(2:4), bt(2:4), esnr, zend
            complex,  allocatable :: de(:), dq(:), re(:), ro(:), r(:)
            real(dp), allocatable :: a11(:), a12(:), a22(:), y1(:), y2(:), a(:), b(:), pw(:), btv(:), nz1(:), nz5(:), nzk(:)
            real(dp) :: rep, ex, nf, det
            integer  :: k_, m_, s_
            nf   = sqrt(real(max(nvox, 1),dp) * real(smpd,dp)**3 / (real(MRES,dp) / 2._dp)**3)
            de   = q(:,1,NMBIN) - q(:,1,1)
            dq   = q(:,2,NMBIN) - q(:,2,1)
            rep  = cdot(de, dq)
            zend = rep / sqrt(max(cdot(de, de) * cdot(dq, dq), 1.e-300_dp)) * nf
            nz1  = real(abs(q(:,1,1) - q(:,2,1))**2, dp) / 4._dp
            nz5  = real(abs(q(:,1,NMBIN) - q(:,2,NMBIN))**2, dp) / 4._dp
            esnr = cdot(q(:,0,NMBIN) - q(:,0,1), q(:,0,NMBIN) - q(:,0,1)) / max(sum(bw * (nz1 + nz5)), 1.e-300_dp)
            do k_ = 2, 4
                allocate(a11(smax), a12(smax), a22(smax), y1(smax), y2(smax), pw(smax), source=0._dp)
                do m_ = 1, nb
                    s_      = bsh(m_)
                    a11(s_) = a11(s_) + bw(m_) * real(abs(q(m_,0,1))**2, dp)
                    a12(s_) = a12(s_) + bw(m_) * real(real(conjg(q(m_,0,1)) * q(m_,0,NMBIN)), dp)
                    a22(s_) = a22(s_) + bw(m_) * real(abs(q(m_,0,NMBIN))**2, dp)
                    y1(s_)  = y1(s_)  + bw(m_) * real(real(conjg(q(m_,0,1)) * q(m_,0,k_)), dp)
                    y2(s_)  = y2(s_)  + bw(m_) * real(real(conjg(q(m_,0,NMBIN)) * q(m_,0,k_)), dp)
                    pw(s_)  = pw(s_)  + bw(m_) * real(abs(q(m_,0,k_))**2, dp)
                end do
                allocate(a(smax), b(smax), btv(smax), source=0._dp)
                do s_ = 1, smax
                    det = a11(s_) * a22(s_) - a12(s_)**2
                    if( det <= 0._dp ) cycle
                    a(s_) = (a22(s_) * y1(s_) - a12(s_) * y2(s_)) / det
                    b(s_) = (a11(s_) * y2(s_) - a12(s_) * y1(s_)) / det
                end do
                r   = q(:,0,k_) - real(a(bsh)) * q(:,0,1) - real(b(bsh)) * q(:,0,NMBIN)
                nzk = real(abs(q(:,1,k_) - q(:,2,k_))**2, dp) / 4._dp
                ex  = sum(bw * (nzk + a(bsh)**2 * nz1 + b(bsh)**2 * nz5))
                rr(k_) = cdot(r, r) / max(ex, 1.e-300_dp)
                re  = q(:,1,k_) - (real(a(bsh)) * q(:,2,1) + real(b(bsh)) * q(:,2,NMBIN))
                ro  = q(:,2,k_) - (real(a(bsh)) * q(:,1,1) + real(b(bsh)) * q(:,1,NMBIN))
                zr(k_)  = cdot(re, ro) / sqrt(max(cdot(re, re) * cdot(ro, ro), 1.e-300_dp)) * nf
                nbv(k_) = 0._dp
                if( rep > 0._dp ) nbv(k_) = cdot(re, ro) / rep
                where( abs(a + b) > 1.e-9_dp ) btv = b / (a + b)
                bt(k_) = 0._dp
                if( sum(pw, mask=abs(a + b) > 1.e-9_dp) > 0._dp ) bt(k_) = sum(btv * pw, mask=abs(a + b) > 1.e-9_dp) / &
                    &sum(pw, mask=abs(a + b) > 1.e-9_dp)
                deallocate(a11, a12, a22, y1, y2, pw, a, b, btv)
            end do
        end subroutine blend

        subroutine acquisition( vx, dx )
            real, intent(out) :: vx, dx
            integer, allocatable :: sel(:), cellv(:), celld(:), lab(:)
            real,    allocatable :: dfs(:), qv(:)
            integer(kind=8) :: seed
            real    :: sec, t0(2), tsum(2)
            integer :: r_, ns, it, a_, b_, tmp
            vx = 0.
            dx = 0.
            sel = pack([(r_, r_=1,size(st))], st == 1 .or. st == NMBIN)
            ns  = size(sel)
            if( count(st(sel) == 1) < 2 .or. count(st(sel) == NMBIN) < 2 ) return
            sec = 360. / real(nsym)
            allocate(cellv(ns), celld(ns))
            do r_ = 1, ns
                cellv(r_) = min(int((cos(e2(sel(r_)) * real(PI) / 180.) + 1.) / 2. * 10.), 9) * 12 + &
                    &min(int(modulo(e1(sel(r_)), sec) / sec * 12.), 11)
            end do
            dfs = pack(dfx, st >= 1 .and. st <= NMBIN)
            qv  = [(quant_of(dfs, 0.1 * real(r_)), r_=1,9)]
            do r_ = 1, ns
                celld(r_) = count(qv < dfx(sel(r_)))
            end do
            lab  = st(sel)
            t0   = [tv(cellv, lab, 120), tv(celld, lab, 10)]
            tsum = 0.
            seed = 12345_8
            do it = 1, NVSHUF
                do a_ = ns, 2, -1
                    seed = mod(seed * 48271_8, 2147483647_8)
                    b_   = 1 + int(mod(seed, int(a_, 8)))
                    tmp     = lab(a_)
                    lab(a_) = lab(b_)
                    lab(b_) = tmp
                end do
                tsum = tsum + [tv(cellv, lab, 120), tv(celld, lab, 10)]
            end do
            vx = t0(1) - tsum(1) / real(NVSHUF)
            dx = t0(2) - tsum(2) / real(NVSHUF)
        end subroutine acquisition

        real function tv( cell, lab, nc )
            integer, intent(in) :: cell(:), lab(:), nc
            real(dp) :: p(0:nc-1), q(0:nc-1)
            integer  :: r_
            p = 0._dp
            q = 0._dp
            do r_ = 1, size(cell)
                if( lab(r_) == 1 )then
                    p(cell(r_)) = p(cell(r_)) + 1._dp
                else
                    q(cell(r_)) = q(cell(r_)) + 1._dp
                endif
            end do
            tv = real(0.5_dp * sum(abs(p / max(sum(p), 1._dp) - q / max(sum(q), 1._dp))))
        end function tv

        subroutine write_call( cw, ln )
            character(len=*), intent(in) :: cw, ln
            integer :: fu, io
            call fopen(fu, string(md//'/motion_call.txt'), status='replace', action='write', iostat=io)
            if( io /= 0 ) return
            write(fu,'(A)') trim(adjustl(cw))
            write(fu,'(A)') ln
            write(fu,'(A)') 'origin '//origin
            if( allocated(qs) )then
                write(fu,'(A)') 'worst_bin_fsc0143_A '//rstr(worst)
                write(fu,'(A)') 'z_loc '//rstr(real(zloc))//' rel_loc '//rstr(real(relloc))//' rel '//rstr(real(rel))//&
                    &' cos '//rstr(real(cosv))//' hmin '//rstr(real(hmin))
                write(fu,'(A)') 'nb_SF '//rstrs(real(nbf))//' residual_z_SF '//rstrs(real(zrf))//' beta_SF '//rstrs(real(btf))
                write(fu,'(A)') 'R_S '//rstrs(real(rss))//' nb_S '//rstrs(real(nbs))//' ends_snr_SF '//rstr(real(esnf))//&
                    &' z_ends_SF '//rstr(real(zenf))
                write(fu,'(A)') 'steps '//rstrs(rn)//' step_share '//rstr(share)//' sign_balance '//rstr(pmr)//&
                    &' focus_fraction '//rstr(focf)
                write(fu,'(A)') 'location_A '//rstrs(cen)//' view_excess '//rstr(vexc)//' defocus_excess '//rstr(dexc)
            endif
            call fclose(fu)
        end subroutine write_call

        function rstrs( x ) result( s_ )
            real, intent(in) :: x(:)
            character(len=:), allocatable :: s_
            integer :: i_
            s_ = ''
            do i_ = 1, size(x)
                s_ = s_//merge('/', ' ', i_ > 1)//rstr(x(i_))
            end do
            s_ = trim(adjustl(s_))
        end function rstrs

    end subroutine motion_stats

end module simple_hetero_tree
