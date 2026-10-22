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
use simple_linalg,          only: eigh
use simple_srch_sort_loc,   only: hpsort
use simple_parameters,      only: parameters
use simple_cmdline,         only: cmdline
use simple_sp_project,      only: sp_project
use simple_qsys_env,        only: qsys_env
use simple_image,           only: image
use simple_state_tree,      only: state_tree
use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
implicit none

public :: hetero_tree
private
#include "simple_local_flags.inc"

integer,  parameter :: KMAX = 8, NCHILD = 10, PERBIN = 1300, NBIN_MIN = 20, NKM_ITS = 40, POLL_SECS = 30, RND_SEED = 0
integer,  parameter :: ITA = 3, ITB = 3, BOXA = 64, BOXB = 88, NSPA = 2500, NSPB = 5000, NPROBE = 2
real,     parameter :: LPF = 15., LPA = 15., LPB = 10., MLP = 15., MFRAC = 0.2, MSIG = 1.5, MDIL = 3., MFALL = 6., MSUP = 0.1
integer,  parameter :: S_QUEUED = 0, S_MEAN = 1, S_FLEX = 2, S_CHEAP = 3, S_FFLEX = 4, S_STAGEA = 5, S_RECB = 6, &
                     & S_STAGEB = 7, S_SPLIT = 8, S_DONE = 9, S_FAILED = -1
character(len=6), parameter :: STEP_DIRS(S_SPLIT) = [character(len=6) :: 'mean', 'flex', 'cheap', 'fflex', 'stageA', &
                                                     & 'recB', 'stageB', 'split']
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
    type(qsys_env) :: qenv
    character(len=:), allocatable :: root
    real    :: smpd = 0.
    integer :: box = 0, nrows = 0
    logical :: l_subtree = .false.
  contains
    procedure :: new
    procedure :: run
    procedure :: kill
    procedure, private :: start
    procedure, private :: advance
    procedure, private :: submit
    procedure, private :: submit_in
    procedure, private :: mergeup
    procedure, private :: step_status
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
        self%root = cwd%to_char()
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
                if( self%mergeup() ) exit
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
        write(logfhandle,'(A)') '>>> HETERO_TREE node '//self%nodes(i)%id//' (depth '//int2str(self%nodes(i)%depth)//'): '//&
            &int2str(size(rows))//' particles, '//int2str(k)//' states'
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
        integer,      allocatable :: rows(:), st(:), pres(:), ch(:)
        real(dp),     allocatable :: w(:)
        character(len=:), allocatable :: d, mask
        integer :: s, k, j, ks
        s = self%nodes(i)%step
        if( present(done_step) ) s = done_step
        d    = self%ndir(i)
        mask = d//'/mask.mrc'
        rows = read_ints(d//'/rows.txt')
        k    = min(KMAX, max(2, size(rows) / max(1, self%p%kfloor)))
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
                st   = project_states(d//'/fflex/f.simple')
                pres = present_labels(st)
                ks   = size(pres)
                allocate(ch(size(st)), source=0)
                do j = 1, ks
                    where( st == pres(j) ) ch = j
                end do
                call simple_mkdir(string(d//'/stageA'))
                call write_states(d//'/'//NODE_PROJ, d//'/stageA/r.simple', ch)
                call cl%set('prg',      'refine3D')
                call cl%set('projfile', 'r.simple')
                call cl%set('nstates',  ks)
                do j = 1, ks
                    call cl%set('vol'//int2str(j), d//'/fflex/recvol_state'//pad2(pres(j))//'.mrc')
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
        call self%qenv%exec_simple_prg_in_queue_async(cl, string('./hetero_step.sh'), string('hetero_step.log'), &
            &exit_code_fname=string(STEP_EXIT))
        call simple_chdir(string(self%root))
        CWD_GLOB = self%root
    end subroutine submit_in

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
        if( file_exists(string(f//JOB_INFO_EXT)) ) return
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
        else if( file_exists(string(f//JOB_INFO_EXT)) )then
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
        call self%qenv%kill
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

    subroutine node_project( self, i )
        class(hetero_tree), intent(inout) :: self
        integer,            intent(in)    :: i
        type(sp_project)     :: spproj
        integer, allocatable :: rows(:)
        integer :: j, np
        rows = read_ints(self%ndir(i)//'/rows.txt')
        call spproj%read(string(self%nodes(i)%pproj))
        np = spproj%os_ptcl3D%get_noris()
        if( np /= self%nrows ) THROW_HARD('the parent project does not hold every particle; node_project')
        do j = 1, np
            call spproj%os_ptcl3D%set_state(j, 0)
            call spproj%os_ptcl2D%set_state(j, 0)
        end do
        do j = 1, size(rows)
            call spproj%os_ptcl3D%set_state(rows(j), 1)
            call spproj%os_ptcl2D%set_state(rows(j), 1)
        end do
        call spproj%write(string(self%ndir(i)//'/'//NODE_PROJ))
        call spproj%kill
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
        integer :: n, q, nb, it, j, b, funit, ios, nl, k, nsel, r, np, t
        fname = self%ndir(i)//'/flex/flex_pca_coordinates.txt'
        nl    = nlines(string(fname))
        call fopen(funit, string(fname), status='old', action='read', iostat=ios)
        call fileiochk('hetero_tree: cannot read '//fname, ios)
        n = 0
        q = 0
        do j = 1, nl
            read(funit,'(A)') line
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            n = n + 1
            if( q == 0 ) q = nwords(line) - 5
        end do
        rewind(funit)
        allocate(pinds(n), z(n,q))
        k = 0
        do j = 1, nl
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
        type(image)           :: img
        real,     allocatable :: e(:,:,:,:), o(:,:,:,:), em(:,:,:), om(:,:,:), v(:,:,:), vs(:,:,:), m(:,:,:)
        logical,  allocatable :: sph(:,:,:), sup(:,:,:), core(:,:,:)
        real(dp), allocatable :: wn(:)
        integer :: ldim(3), n, nm, j, ifoo
        real    :: smpd, rad, sc
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
        rad = self%p%mskdiam / 2. / smpd
        sph = sphere(n, rad)
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
        deallocate(e, o)
        vs = filtered(v, 2)
        where( .not. sph ) vs = 0.
        core = vs > MFRAC * maxval(vs)
        if( .not. any(core) )then
            m = merge(1., 0., sph)
        else
            m    = dist_to(core)
            core = m <= MDIL
            m    = dist_to(core)
            m    = exp(-m**2 / (2. * MFALL**2))
            where( .not. sph ) m = 0.
        endif
        call img%new([n,n,n], smpd)
        call img%set_rmat(m, .false.)
        call img%write(string(fname))
        call img%set_rmat(vs, .false.)
        call img%write(string(fname(1:len(fname)-4)//'_variance.mrc'))
        call img%kill
        write(logfhandle,'(A)') '>>> HETERO_TREE focus mask '//fname//' from '//int2str(nm)//' maps: '//&
            &int2str(count(m > 0.5))//' voxels above 0.5 ('//int2str(nint(100. * real(count(m > 0.5 .and. sup)) / &
            &real(max(1,count(sup)))))//' % of the molecule support)'

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

        function filtered( x, mode ) result( y )
            real,    intent(in) :: x(:,:,:)
            integer, intent(in) :: mode
            real, allocatable :: y(:,:,:)
            type(image) :: vol
            call vol%new([n,n,n], smpd)
            call vol%set_rmat(x, .false.)
            call fourier_filter(vol, n, smpd, mode)
            y = vol%get_rmat()
            call vol%kill
        end function filtered

        function dist_to( mk ) result( dpx )
            logical, intent(in) :: mk(:,:,:)
            real, allocatable :: dpx(:,:,:)
            type(image) :: vol
            call vol%new([n,n,n], smpd)
            call vol%set_rmat(merge(0., 1., mk), .false.)
            call vol%distance_transform
            dpx = vol%get_rmat() / smpd
            call vol%kill
        end function dist_to

    end subroutine falloff_mask

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
        call cl%set('nparts',  self%p%nparts)
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

    subroutine write_states( src, dst, st )
        character(len=*), intent(in) :: src, dst
        integer,          intent(in) :: st(:)
        type(sp_project) :: spproj
        integer :: j
        call spproj%read(string(src))
        if( spproj%os_ptcl3D%get_noris() /= size(st) ) THROW_HARD('state count does not match the project; write_states')
        do j = 1, size(st)
            call spproj%os_ptcl3D%set_state(j, st(j))
            call spproj%os_ptcl2D%set_state(j, st(j))
        end do
        call spproj%write(string(dst))
        call spproj%kill
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
