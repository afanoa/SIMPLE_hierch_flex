!@descr: recursive heterogeneity tree: a resumable node queue whose mean, flex, focused refinement and split-test steps run as queue jobs, with cheap children, falloff-prior focus masks, node decisions and the leaves
module simple_hetero_tree
use simple_defs,            only: dp, logfhandle, PI, CWD_GLOB
use simple_defs_fname,      only: JOB_INFO_EXT
use simple_error,           only: simple_exception
use simple_string,          only: string
use simple_string_utils,    only: int2str
use simple_fileio,          only: fopen, fclose, fileiochk, nlines, read_exit_code, simple_copy_file
use simple_syslib,          only: file_exists, simple_mkdir, simple_chdir, simple_getcwd, simple_rename, simple_abspath
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
implicit none

public :: hetero_tree, discover_sites
private
#include "simple_local_flags.inc"

integer,  parameter :: KMAX = 8, NCHILD = 10, PERBIN = 1300, NBIN_MIN = 20, NKM_ITS = 40, POLL_SECS = 5, RND_SEED = 0
integer,  parameter :: DMINVOX = 200, DSEEDMIN = 1000, DFACT_ITS = 20
real,     parameter :: DMINSHARE = 0.05, DFACT_COS = 0.8, DFACT_SHARE = 0.1, DDICE = 0.9, DMAXMOL = 0.3
integer,  parameter :: ITA = 3, ITB = 3, BOXA = 64, BOXB = 88, NSPA = 2500, NSPB = 5000, NPROBE = 2
integer,  parameter :: NPARTS_MIN = 4, PTCLS_PER_PART = 20000, SITE_NTHR = 32
real,     parameter :: LPF = 15., LPA = 15., LPB = 10., MLP = 15., MFRAC = 0.2, MSIG = 1.5, MDIL = 3., MFALL = 6., MSUP = 0.1
real,     parameter :: ATHR = 0.01, LPD = 10.
integer,  parameter :: S_QUEUED = 0, S_MEAN = 1, S_FLEX = 2, S_CHEAP = 3, S_FFLEX = 4, S_STA0 = 5, S_RECB0 = 6, S_SITE = 7, &
                     & S_SEEDS = 8, S_STAGEA = 9, S_RECB = 10, S_STAGEB = 11, S_SPLIT = 12, S_DONE = 13, S_FAILED = -1
character(len=6), parameter :: STEP_DIRS(S_SPLIT) = [character(len=6) :: 'mean', 'flex', 'cheap', 'fflex', 'stA0', &
                                                     & 'recB0', 'site', 'seeds', 'stageA', 'recB', 'stageB', 'split']
character(len=*), parameter :: NODE_TABLE = 'hetero_tree_nodes.txt', LEAF_TABLE = 'hetero_tree_leaves.txt'
character(len=*), parameter :: STEP_EXIT = 'hetero_step.exit', NODE_PROJ = 'n.simple'

type :: tnode
    character(len=:), allocatable :: id, pproj, prows
    integer :: depth = 0, step = S_QUEUED
end type tnode

type :: hetero_tree
    private
    class(parameters), pointer :: p => null()
    type(tnode), allocatable :: nodes(:)
    type(string), allocatable :: subd(:)
    type(qsys_env) :: qenv, qsite
    character(len=:), allocatable :: root, sitemask
    real    :: smpd = 0.
    integer :: box = 0, nrows = 0, npcur = 0
    logical :: l_subtree = .false., l_lsdone = .false.
  contains
    procedure :: new
    procedure :: run
    procedure :: kill
    procedure, private :: start
    procedure, private :: advance
    procedure, private :: submit
    procedure, private :: submit_in
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
end type hetero_tree

contains

    subroutine new( self, params, l_rows )
        class(hetero_tree),        intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        logical,                   intent(in)    :: l_rows
        type(sp_project)     :: spproj
        type(string)         :: cwd, absproj
        integer, allocatable :: states(:), rows(:), sel(:)
        integer :: i
        call self%kill
        self%p         => params
        self%l_subtree = l_rows
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
        if( file_exists(string(NODE_TABLE)) )then
            call self%read_table
            write(logfhandle,'(A)') '>>> HETERO_TREE resumed from '//NODE_TABLE//': '//int2str(size(self%nodes))//' nodes'
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
        call self%qenv%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr)
        call self%qsite%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=max(params%nthr, SITE_NTHR))
    end subroutine new

    subroutine run( self )
        class(hetero_tree), intent(inout) :: self
        integer :: i, nact, code
        call seed_rnd_fixed(RND_SEED)
        do
            nact = 0
            i    = 0
            do while( i < size(self%nodes) )
                i = i + 1
                if( self%nodes(i)%step < S_MEAN .or. self%nodes(i)%step > S_SPLIT ) cycle
                code = self%step_status(i, self%nodes(i)%step)
                if( code == -1 )then
                    nact = nact + 1
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
                        &trim(STEP_DIRS(self%nodes(i)%step))//' FAILED (exit code '//int2str(code)//'), node stopped'
                    self%nodes(i)%step = S_FAILED
                    call self%write_table
                    cycle
                endif
                call self%advance(i)
                if( self%nodes(i)%step <= S_SPLIT ) nact = nact + 1
                call self%write_table
            end do
            do i = 1, size(self%nodes)
                if( nact >= max(1, self%p%nodepar) ) exit
                if( self%nodes(i)%step /= S_QUEUED ) cycle
                call self%start(i, S_MEAN)
                nact = nact + 1
                call self%write_table
            end do
            if( nact == 0 .and. .not. any(self%nodes(:)%step == S_QUEUED) )then
                if( self%final() ) exit
            endif
            call sleep(POLL_SECS)
        end do
    end subroutine run

    subroutine start( self, i, istep )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, istep
        type(cmdline)        :: cl
        integer, allocatable :: rows(:)
        character(len=:), allocatable :: d
        integer :: k
        if( istep > S_MEAN )then
            call self%advance(i, istep - 1)
            return
        endif
        d = self%ndir(i)
        if( .not. file_exists(string(d//'/'//NODE_PROJ)) ) call node_project(self, i)
        rows = read_ints(d//'/rows.txt')
        k    = min(KMAX, max(2, size(rows) / max(1, self%p%kfloor)))
        call self%set_parts(size(rows))
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' (depth '//int2str(self%nodes(i)%depth)//'): '//&
            &int2str(size(rows))//' particles, '//int2str(k)//' flex states, '//int2str(self%npcur)//' parts per step'
        call simple_mkdir(string(d//'/mean'))
        call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/mean/m.simple'))
        call cl%set('prg',      'reconstruct3D')
        call cl%set('projfile', 'm.simple')
        call cl%set('nstates',  1)
        call common_keys(self, cl)
        call self%submit(i, S_MEAN, cl)
    end subroutine start

    subroutine advance( self, i, done_step )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        integer, optional,  intent(in)    :: done_step
        type(cmdline)             :: cl
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: rows(:), st(:), pres(:), ch(:), seeds(:)
        real(dp),     allocatable :: w(:)
        character(len=:), allocatable :: d, mask, sitemask
        integer :: s, k, j, ks
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
                    call fflex_step
                else
                    call cheap_children(self, i, ch)
                    call simple_mkdir(string(d//'/cheap'))
                    call write_states(d//'/'//NODE_PROJ, d//'/cheap/c.simple', ch)
                    call cl%set('prg',      'reconstruct3D')
                    call cl%set('projfile', 'c.simple')
                    call cl%set('nstates',  NCHILD)
                    call common_keys(self, cl)
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
                call fflex_step
            case(S_FFLEX)
                if( self%p%site_disc .ne. 'yes' )then
                    call site_step(sitemask)
                    return
                endif
                st   = project_states(d//'/fflex/f.simple')
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
                call refine_keys(self, cl, LPA, ITA, NSPA, BOXA, mask)
                call self%submit(i, S_STA0, cl)
            case(S_STA0)
                ks = maxval(project_states(d//'/stA0/r.simple'))
                call simple_mkdir(string(d//'/recB0'))
                call simple_copy_file(string(d//'/stA0/r.simple'), string(d//'/recB0/b.simple'))
                call cl%set('prg',      'reconstruct3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                call common_keys(self, cl)
                call self%submit(i, S_RECB0, cl)
            case(S_RECB0)
                st   = project_states(d//'/recB0/b.simple')
                pres = present_labels(st)
                allocate(evens(size(pres)), odds(size(pres)), w(size(pres)))
                do j = 1, size(pres)
                    evens(j) = d//'/recB0/'//hv(pres(j), 'even')
                    odds(j)  = d//'/recB0/'//hv(pres(j), 'odd')
                    w(j)     = real(count(st == pres(j)), dp)
                end do
                call simple_mkdir(string(d//'/sites'))
                call discover_sites(evens, odds, w, pres, self%p%mskdiam, d//'/sites', ks)
                if( ks > 0 )then
                    sitemask = d//'/sites/site_01.mrc'
                else
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': no local site found, the node mask ranks the poles'
                    sitemask = mask
                endif
                call site_step(sitemask)
            case(S_SITE)
                seeds = read_ints(d//'/site/site_seeds.txt')
                ks    = maxval([0, seeds])
                call log_lines(d//'/site/site_poles.txt', 'pole (index, seeds, endpoints)')
                call log_lines(d//'/site/site_halves.txt', 'halves')
                call log_lines(d//'/site/site_spearman.txt', 'endpoint Spearman (merged at >= '//&
                    &trim(adjustl(rfmt2(self%p%pole_merge)))//')')
                call log_lines(d//'/site/site_junk_report.txt', 'junk screen')
                if( ks < 2 )then
                    write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': fewer than two seed states (accepted '//&
                        &'poles and the rest), the node is refined as one state and becomes a leaf'
                    if( ks < 1 ) seeds = merge(1, 0, project_states(d//'/'//NODE_PROJ) > 0)
                    ks    = 1
                endif
                call simple_mkdir(string(d//'/seeds'))
                call write_ints(d//'/seeds/seeds.txt', seeds)
                call write_states(d//'/'//NODE_PROJ, d//'/seeds/s.simple', seeds)
                call cl%set('prg',      'reconstruct3D')
                call cl%set('projfile', 's.simple')
                call cl%set('nstates',  ks)
                call common_keys(self, cl)
                call self%submit(i, S_SEEDS, cl)
            case(S_SEEDS)
                seeds = read_ints(d//'/seeds/seeds.txt')
                ks    = maxval(seeds)
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
                ks = maxval(project_states(d//'/stageA/r.simple'))
                call simple_mkdir(string(d//'/recB'))
                call simple_copy_file(string(d//'/stageA/r.simple'), string(d//'/recB/b.simple'))
                call cl%set('prg',      'reconstruct3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                call common_keys(self, cl)
                call self%submit(i, S_RECB, cl)
            case(S_RECB)
                ks = maxval(project_states(d//'/stageA/r.simple'))
                call simple_mkdir(string(d//'/stageB'))
                call simple_copy_file(string(d//'/stageA/r.simple'), string(d//'/stageB/b.simple'))
                call cl%set('prg',      'refine3D')
                call cl%set('projfile', 'b.simple')
                call cl%set('nstates',  ks)
                do j = 1, ks
                    call cl%set('vol'//int2str(j), d//'/recB/recvol_state'//pad2(j)//'.mrc')
                end do
                call refine_keys(self, cl, LPB, ITB, NSPB, BOXB, mask)
                call self%submit(i, S_STAGEB, cl)
            case(S_STAGEB)
                call simple_mkdir(string(d//'/split'))
                st = project_states(d//'/stageB/b.simple')
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
                call cl%set('nthr',        self%p%nthr)
                call self%submit(i, S_SPLIT, cl)
            case(S_SPLIT)
                call self%decide(i)
                self%nodes(i)%step = S_DONE
        end select

    contains

        subroutine site_step( smask )
            character(len=*), intent(in) :: smask
            call simple_mkdir(string(d//'/site'))
            call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/site/s.simple'))
            call cl%set('prg',         'site_rank')
            call cl%set('projfile',    's.simple')
            call cl%set('infile',      d//'/fflex/flex_pca_coordinates.txt')
            call cl%set('site_dirs',   trim(self%p%site_dirs))
            call cl%set('pole_merge',  self%p%pole_merge)
            call cl%set('junk_screen', trim(self%p%junk_screen))
            call cl%set('pcg_mskfile', smask)
            call cl%set('pgrp',        trim(self%p%pgrp))
            call cl%set('mskdiam',     self%p%mskdiam)
            call cl%set('nthr',        max(self%p%nthr, SITE_NTHR))
            call self%submit(i, S_SITE, cl)
        end subroutine site_step

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

        subroutine fflex_step
            call simple_mkdir(string(d//'/fflex'))
            call simple_copy_file(string(d//'/'//NODE_PROJ), string(d//'/fflex/f.simple'))
            call flex_keys(self, cl, d//'/mean/recvol_state01.mrc', k)
            call cl%set('pcg_mskfile', mask)
            call self%submit(i, S_FFLEX, cl)
        end subroutine fflex_step

    end subroutine advance

    subroutine submit( self, i, istep, cl )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, istep
        class(cmdline),     intent(inout) :: cl
        call self%submit_in(self%ndir(i)//'/'//trim(STEP_DIRS(istep)), cl)
        self%nodes(i)%step = istep
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': '//trim(STEP_DIRS(istep))//' submitted'
    end subroutine submit

    subroutine submit_in( self, sd, cl )
        class(hetero_tree), intent(inout) :: self
        character(len=*),   intent(in)    :: sd
        class(cmdline),     intent(inout) :: cl
        call cl%set('mkdir', 'no')
        call simple_chdir(string(sd))
        CWD_GLOB = sd
        if( cl%get_iarg('nthr') > self%p%nthr )then
            call self%qsite%exec_simple_prg_in_queue_async(cl, string('./hetero_step.sh'), string('hetero_step.log'), &
                &exit_code_fname=string(STEP_EXIT))
        else
            call self%qenv%exec_simple_prg_in_queue_async(cl, string('./hetero_step.sh'), string('hetero_step.log'), &
                &exit_code_fname=string(STEP_EXIT))
        endif
        call simple_chdir(string(self%root))
        CWD_GLOB = self%root
        if( .not. allocated(self%subd) ) allocate(self%subd(0))
        self%subd = [self%subd, string(sd)]
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
        final = .false.
        if( .not. self%plain() ) return
        if( .not. self%l_lsdone )then
            if( .not. self%landscape() ) return
            self%l_lsdone = .true.
        endif
        if( .not. file_exists(string(self%root//'/final_maps/choice.txt')) ) call final_maps(self)
        final = self%mergeup()
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
            np = count_poles(nd//'/site/site_poles.txt')
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
            call cl%set('infile',   nd//'/site/site_poles.txt')
            call cl%set('infile2',  nd//'/split/split_children.txt')
            call cl%set('mskdiam',  self%p%mskdiam)
            call cl%set('nthr',     self%p%nthr)
            call self%submit_in(ld, cl)
            write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': landscape submitted, '//int2str(np)//&
                &' accepted poles'
            nrun = nrun + 1
        end do
        if( .not. landscape ) return
        do i = 1, size(self%nodes)
            if( self%nodes(i)%step /= S_DONE ) cycle
            ld = self%ndir(i)//'/landscape'
            if( .not. file_exists(string(ld//'/'//STEP_EXIT)) ) cycle
            np = count_poles(self%ndir(i)//'/site/site_poles.txt')
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
        if( .not. file_exists(string(LEAF_TABLE)) ) return
        call leaf_lines(lines, nds)
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
        character(len=:), allocatable :: md, f
        logical :: err
        integer :: code, nn, j
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
        call cl%set('prg',      'state_mergeup')
        call cl%set('projfile', 'final.simple')
        call cl%set('infile',   self%root//'/leaf_tree/merge_tree.txt')
        call cl%set('infile2',  md//'/leaf_masks.txt')
        call cl%set('mskdiam',  self%p%mskdiam)
        call cl%set('nthr',     self%p%nthr)
        call self%submit_in(md, cl)
        write(logfhandle,'(A)') '>>> HETERO_TREE merge-up of the '//int2str(nn)//' leaves submitted'
    end function mergeup

    integer function step_status( self, i, istep )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i, istep
        character(len=:), allocatable :: f
        logical :: err
        f = self%ndir(i)//'/'//trim(STEP_DIRS(istep))//'/'//STEP_EXIT
        if( file_exists(string(f)) )then
            call read_exit_code(string(f), step_status, err)
            if( err ) step_status = -1
        else if( file_exists(string(f//JOB_INFO_EXT)) .or. self%pending(self%ndir(i)//'/'//trim(STEP_DIRS(istep))) )then
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
        character(len=:), allocatable :: d, name
        logical :: l_split
        integer :: s, g, ng, n, fg, big, j
        d      = self%ndir(i)
        rows   = read_ints(d//'/rows.txt')
        lab    = project_states(d//'/stageB/b.simple')
        ch     = read_ints(d//'/split/split_children.txt')
        junk   = read_ints(d//'/split/split_junk.txt')
        if( size(ch) /= size(lab) ) THROW_HARD('split children do not match the node; decide')
        states = present_labels(lab)
        allocate(to(size(states)), source=0)
        allocate(pop(max(1,maxval([0,ch]))), source=0)
        do g = 1, size(pop)
            pop(g) = count(ch == g)
        end do
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
        fg      = 0
        do g = 1, size(pop)
            if( .not. any(to == g) ) cycle
            mem = pack(states, to == g)
            n   = count([(any(mem == lab(j)), j=1,size(lab))])
            if( l_split .and. n >= self%p%minpop .and. self%nodes(i)%depth < self%p%maxdepth )then
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
                call write_ints(self%ndir(size(self%nodes))//'/rows.txt', pack([(j, j=1,size(lab))], &
                    &[(any(mem == lab(j)), j=1,size(lab))]))
                write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': states '//joined(mem)//' ('//&
                    &int2str(n)//' particles) -> node '//name
            else
                call add_leaf(self, i, mem, .false.)
            endif
        end do
        do j = 1, size(junk)
            if( count(lab == junk(j)) > 0 ) call add_leaf(self, i, [junk(j)], .true.)
        end do
    end subroutine decide

    subroutine finish( self, nn )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(out)   :: nn
        type(state_tree)          :: tree
        type(string), allocatable :: evens(:), odds(:)
        integer,      allocatable :: labels(:), rows(:), names(:)
        real(dp),     allocatable :: counts(:)
        character(len=4096) :: line
        integer :: funit, ios, nl, il
        if( any(self%nodes(:)%step == S_FAILED) ) write(logfhandle,'(A)') '>>> HETERO_TREE WARNING: failed nodes'//&
            &failed_list(self)
        allocate(labels(self%nrows), source=0)
        nl = 0
        if( file_exists(string(LEAF_TABLE)) ) nl = nlines(string(LEAF_TABLE))
        allocate(evens(nl), odds(nl), names(nl), counts(nl))
        nn = 0
        if( nl > 0 )then
            call fopen(funit, string(LEAF_TABLE), status='old', action='read', iostat=ios)
            call fileiochk('hetero_tree: cannot read '//LEAF_TABLE, ios)
            do il = 1, nl
                read(funit,'(A)') line
                rows = read_ints(word(line, 8))
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

    subroutine read_table( self )
        class(hetero_tree), intent(inout) :: self
        character(len=4096) :: line
        integer :: funit, ios, n, i
        n = nlines(string(NODE_TABLE))
        allocate(self%nodes(n))
        call fopen(funit, string(NODE_TABLE), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//NODE_TABLE, ios)
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
        call fopen(funit, string(NODE_TABLE//'.tmp'), status='replace', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot write '//NODE_TABLE, ios)
        do i = 1, size(self%nodes)
            write(funit,'(A)') self%nodes(i)%id//' '//int2str(self%nodes(i)%depth)//' '//int2str(self%nodes(i)%step)//&
                &' '//self%nodes(i)%pproj//' '//self%nodes(i)%prows
        end do
        call fclose(funit)
        call simple_rename(string(NODE_TABLE//'.tmp'), string(NODE_TABLE))
    end subroutine write_table

    function ndir( self, i ) result( d )
        class(hetero_tree), intent(in) :: self
        integer,            intent(in) :: i
        character(len=:), allocatable :: d
        d = self%root//'/node_'//self%nodes(i)%id
    end function ndir

    subroutine kill( self )
        class(hetero_tree), intent(inout) :: self
        if( allocated(self%nodes) ) deallocate(self%nodes)
        if( allocated(self%subd)  ) deallocate(self%subd)
        call self%qenv%kill
        call self%qsite%kill
        self%p => null()
    end subroutine kill

    subroutine final_project( self, fname )
        class(hetero_tree), intent(inout) :: self
        character(len=*),   intent(in)    :: fname
        type(sp_project)     :: master, nproj
        integer, allocatable :: rows(:)
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
                nd_read = nd
            endif
            rows = read_ints(lrows)
            do j = 1, size(rows)
                call master%os_ptcl3D%transfer_ori(rows(j), nproj%os_ptcl3D, rows(j))
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

    subroutine leaf_lines( lines, nds )
        type(string), allocatable, intent(out) :: lines(:), nds(:)
        character(len=4096) :: line
        character(len=:), allocatable :: lrows
        integer :: funit, ios, nl, il
        nl = nlines(string(LEAF_TABLE))
        allocate(lines(nl), nds(nl))
        call fopen(funit, string(LEAF_TABLE), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//LEAF_TABLE, ios)
        do il = 1, nl
            read(funit,'(A)') line
            lines(il) = trim(line)
            lrows     = word(line, 8)
            nds(il)   = lrows(1:index(lrows, '/leaf_') - 1)
        end do
        call fclose(funit)
    end subroutine leaf_lines

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

    subroutine node_project( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        integer, allocatable :: rows(:), st(:)
        rows = read_ints(self%ndir(i)//'/rows.txt')
        allocate(st(self%nrows), source=0)
        st(rows) = 1
        if( size(project_states(self%nodes(i)%pproj)) /= self%nrows )then
            THROW_HARD('the parent project does not hold every particle; node_project')
        endif
        call write_states(self%nodes(i)%pproj, self%ndir(i)//'/'//NODE_PROJ, st)
    end subroutine node_project

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

    subroutine discover_sites( evens, odds, counts, names, mskdiam, outdir, nsite )
        type(string),     intent(in)  :: evens(:), odds(:)
        real(dp),         intent(in)  :: counts(:)
        integer,          intent(in)  :: names(:)
        real,             intent(in)  :: mskdiam
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: nsite
        type :: dcand
            real, allocatable :: m(:,:,:)
            real    :: r = 0., share = 0., mol = 0., rad = 0.
            integer :: comp = 0, fac = 0, ncore = 0, seed = 0, vox = 0
        end type dcand
        type(dcand), allocatable :: cands(:), tmp(:)
        type(image)              :: img
        real,     allocatable :: vs(:,:,:), de(:,:,:,:), dv(:,:,:,:), pn(:,:), cc(:,:), vv(:)
        logical,  allocatable :: sph(:,:,:), sup(:,:,:), core(:,:,:), cm(:,:,:), b(:,:,:), okk(:)
        integer,  allocatable :: lab(:,:,:), lb(:,:,:), asg(:), vidx(:,:), ord(:), rep(:), memb(:)
        real(dp), allocatable :: wn(:)
        character(len=:), allocatable :: mlist
        real    :: smpd, tot, cosv, sh(2)
        integer :: n, nk, nl, ic, nv, iv, it, j, i, k, x, y, z, nb, funit, ios
        logical :: l_sib
        wn  = counts / sum(counts)
        nk  = size(evens)
        okk = counts >= real(DSEEDMIN,dp)
        call half_variance(evens, odds, counts, mskdiam, n, smpd, sph, sup, vs, de, dv)
        core = vs > MFRAC * maxval(vs)
        call label6(core, lab, nl)
        tot  = sum(vs, mask=core)
        allocate(cands(0))
        do ic = 1, nl
            cm = lab == ic
            if( count(cm) < DMINVOX ) cycle
            if( sum(vs, mask=cm) / tot < DMINSHARE ) cycle
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
                    deallocate(pn, vidx, vv)
                    cycle
                endif
                write(logfhandle,'(A,I0,A,F5.2,A,F5.2,A,F5.2,A)') '>>> SITE_DISCOVER component ', ic, ' is one factor (|cos| ', &
                    &cosv, ', shares ', sh(1), ' / ', sh(2), ')'
                deallocate(pn, vidx, vv)
            endif
            call add(cm, ic, 0)
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
        do i = 1, size(ord)
            k = 0
            do j = 1, size(rep)
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

    subroutine add_leaf( self, i, mem, l_junk )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i, mem(:)
        logical,            intent(in)    :: l_junk
        character(len=11), parameter :: SFX(3) = [character(len=11) :: '', '_even_unfil', '_odd_unfil']
        type(image)          :: vol
        integer, allocatable :: rows(:), lab(:), sel(:)
        real,    allocatable :: acc(:,:,:)
        character(len=:), allocatable :: d, name, b, fn, lrows
        real    :: wsum, smpd
        integer :: j, isfx, funit, ios, n, ldim(3), ifoo
        d    = self%ndir(i)
        rows = read_ints(d//'/rows.txt')
        lab  = project_states(d//'/stageB/b.simple')
        sel  = pack([(j, j=1,size(lab))], [(any(mem == lab(j)), j=1,size(lab))])
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
        call write_ints(lrows, sel)
        call fopen(funit, string(self%root//'/'//LEAF_TABLE), status='unknown', position='append', action='write', iostat=ios)
        call fileiochk('hetero_tree: cannot append to '//LEAF_TABLE, ios)
        write(funit,'(A)') name//' '//merge('junk','leaf',l_junk)//' '//int2str(n)//' '//joined(mem, '+')//' '//b//'.mrc '//&
            &b//'_even_unfil.mrc '//b//'_odd_unfil.mrc '//lrows
        call fclose(funit)
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//': states '//joined(mem)//' ('//int2str(n)//&
            &' particles) -> '//merge('junk leaf ','leaf      ',l_junk)//name
    end subroutine add_leaf

    subroutine common_keys( self, cl )
        class(hetero_tree), intent(in)    :: self
        class(cmdline),     intent(inout) :: cl
        call cl%set('pgrp',    trim(self%p%pgrp))
        call cl%set('mskdiam', self%p%mskdiam)
        call cl%set('nparts',  merge(self%npcur, self%p%nparts, self%npcur > 0))
        call cl%set('nthr',    self%p%nthr)
    end subroutine common_keys

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
            read(funit,*) v(j)
        end do
        call fclose(funit)
    end function read_ints

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

end module simple_hetero_tree
