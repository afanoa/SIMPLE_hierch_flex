!@descr: the kernel-latent EM tree as one program (prg=klem_tree): the stages of the bash chain in Fortran, resumable
!! from existing outputs. stage=all: parent board, families, parent regions, one family job per family (async),
!! the landscape. stage=family: child board, k-means, refine3D_states (unless skip_refine), feedback board, regions,
!! one region job per region (async, parmax at a time). stage=region: focal board, node states, reconstructions
!! (unless norec), level 2 along the node's own axis; stage=zoom|pool|parents: zoom boards, merge, parent maps. Config file
!! (infile): "key value" lines; see klem_tree_cfg below. Boards run through prg=klem_board (distributed when nparts > 1).
module simple_commanders_klem_tree
use simple_core_module_api
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_qsys_async_job, only: qsys_async_job, ASYNC_JOB_RUNNING, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_klem_tree
implicit none
public :: commander_klem_tree
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_klem_tree
  contains
    procedure :: execute => exec_klem_tree
end type commander_klem_tree

type klem_tree_cfg
    character(len=PLEN) :: run = 'klem_tree', stage = 'all', meta = 'none', subset = 'none', docked = 'none', docked_all = 'none'
    character(len=PLEN) :: cons = 'none', closed = 'none', tag = '', mask = 'none', partition = 'none', refine_extra = 'ml_reg=no'
    character(len=PLEN) :: kind = 'region', rtag = 'klem'
    integer :: fam = 1, zstate = 0, nfam = 2, nparts = 20, nthr_part = 8, nthr = 8, np_refine = 20, maxits = 6, k_refine = 10
    integer :: nfree_parent = 2, nfree_sub = 1, maxreg = 6, nl2 = 3, parmax = 4, zparmax = 16, kgroups = 5, nmax = 30000
    integer :: ntop2 = 10000, maxk = 30, zmax = 10000, zmaxk = 10, zmin = 2000, l2min = 4000, minrec = 50, minhalf = 10
    integer :: nprec = 4, npz = 2, nthr_zoom = 12, nprec_cons = 6, chunk_cons = 60, np_parent = 10, minp = 1000
    integer :: minrec_leaf = 1000, uncap = 1, ctfabs = 1, pscale = 1, server = 1, kmax_load = 29, minvox = 40, ndilate = 2
    integer :: skip_refine = 0, norec = 0, mem_job = 64000, hours_job = 24
    real    :: sim = 0.7, frac = 0.8, minmass = 2000., minrep = 0.5, balance = 0.5, devw = 0.05, mskdiam = 232., mskdiam_refine = 220.
    real    :: lpstart = 20., lpstop = 8., region_thr = 0.5, cmerge = 0.8, floor = 12., fscmax = 12.
    character(len=256) :: lps_parent = '20 20 20 14 14 14 10 10 10', sig_parent = '1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7'
    character(len=256) :: lps_child  = '14 14 12 12 10 10 8 8 8',   sig_child  = '1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7'
    character(len=256) :: lps_region = '12 12 10 10 8 8 8 8 8',     sig_region = '1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7'
    character(len=256) :: lps_zoom   = '12 10 10 8 8 8',            sig_zoom   = '1.2 1.0 1.0 0.8 0.7 0.7'
end type klem_tree_cfg

contains

    subroutine exec_klem_tree( self, cline )
        class(commander_klem_tree), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)    :: params
        type(klem_tree_cfg) :: cfg
        type(sp_project)    :: spproj
        type(string)        :: cwd0
        character(len=PLEN) :: R, src
        integer :: nrows, box
        real    :: smpd
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        call params%new(cline)
        call read_cfg(params%infile%to_char(), cfg)
        call simple_getcwd(cwd0)
        R   = absolute(cfg%run, cwd0%to_char())
        src = absolute(params%projfile%to_char(), cwd0%to_char())
        box  = params%box
        smpd = params%smpd
        call simple_mkdir(trim(R))
        call spproj%read_segment('ptcl3D', params%projfile)
        nrows = spproj%os_ptcl3D%get_noris()
        call spproj%kill
        call simple_chdir(trim(R))
        call klog('tree', 'stage '//trim(cfg%stage)//' in '//trim(R)//' on '//trim(src))
        select case(trim(cfg%stage))
            case('all')
                call stage_all(params, cfg, R, src, nrows, box, smpd)
            case('family')
                call stage_family(params, cfg, R, src, nrows, box, smpd)
            case('region')
                call stage_region(cfg, R, src, nrows, smpd)
            case('zoom')
                call stage_zoom(cfg, R, src, nrows)
            case('pool')
                call stage_pool(cfg, R, src, nrows, smpd)
            case('parents')
                call stage_parents(cfg, R, src)
            case default
                THROW_HARD('unknown klem_tree stage '//trim(cfg%stage))
        end select
        call simple_chdir(cwd0)
        call simple_end('**** SIMPLE_KLEM_TREE NORMAL STOP ****')
    end subroutine exec_klem_tree

    function absolute( p, cwd ) result( a )
        character(len=*), intent(in) :: p, cwd
        character(len=PLEN) :: a
        if( p(1:1) == '/' )then
            a = trim(p)
        else
            a = trim(cwd)//'/'//trim(p)
        endif
    end function absolute

    subroutine klog( tag, msg )
        character(len=*), intent(in) :: tag, msg
        integer :: v(8)
        character(len=32) :: ts
        call date_and_time(values=v)
        write(ts,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') v(1), '-', v(2), '-', v(3), ' ', v(5), ':', v(6), ':', v(7)
        write(logfhandle,'(A)') '['//trim(ts)//' '//trim(tag)//'] '//trim(msg)
        call flush(logfhandle)
    end subroutine klog

    subroutine read_cfg( fname, c )
        character(len=*),    intent(in)    :: fname
        type(klem_tree_cfg), intent(inout) :: c
        character(len=1024) :: line, key, val
        integer :: u, io, p
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) THROW_HARD('cannot open the klem_tree config '//trim(fname))
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            line = adjustl(line)
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            p = index(trim(line), ' ')
            if( p <= 0 )then
                key = line
                val = ''
            else
                key = line(1:p-1)
                val = adjustl(line(p+1:))
            endif
            call set_key(c, trim(key), trim(val))
        enddo
        close(u)
    end subroutine read_cfg

    subroutine set_key( c, key, val )
        type(klem_tree_cfg), intent(inout) :: c
        character(len=*),    intent(in)    :: key, val
        select case(key)
            case('run');            c%run          = val
            case('rtag');           c%rtag         = val
            case('stage');          c%stage        = val
            case('kind');           c%kind         = val
            case('meta');           c%meta         = val
            case('subset');         c%subset       = val
            case('docked');         c%docked       = val
            case('docked_all');     c%docked_all   = val
            case('cons');           c%cons         = val
            case('closed');         c%closed       = val
            case('tag');            c%tag          = val
            case('mask');           c%mask         = val
            case('partition');      c%partition    = val
            case('refine_extra');   c%refine_extra = val
            case('lps_parent');     c%lps_parent   = val
            case('sig_parent');     c%sig_parent   = val
            case('lps_child');      c%lps_child    = val
            case('sig_child');      c%sig_child    = val
            case('lps_region');     c%lps_region   = val
            case('sig_region');     c%sig_region   = val
            case('lps_zoom');       c%lps_zoom     = val
            case('sig_zoom');       c%sig_zoom     = val
            case('fam');            read(val,*) c%fam
            case('zstate');         read(val,*) c%zstate
            case('nfam');           read(val,*) c%nfam
            case('nparts');         read(val,*) c%nparts
            case('nthr_part');      read(val,*) c%nthr_part
            case('nthr');           read(val,*) c%nthr
            case('np_refine');      read(val,*) c%np_refine
            case('maxits');         read(val,*) c%maxits
            case('k_refine');       read(val,*) c%k_refine
            case('nfree_parent');   read(val,*) c%nfree_parent
            case('nfree_sub');      read(val,*) c%nfree_sub
            case('maxreg');         read(val,*) c%maxreg
            case('nl2');            read(val,*) c%nl2
            case('parmax');         read(val,*) c%parmax
            case('zparmax');        read(val,*) c%zparmax
            case('kgroups');        read(val,*) c%kgroups
            case('nmax');           read(val,*) c%nmax
            case('ntop2');          read(val,*) c%ntop2
            case('maxk');           read(val,*) c%maxk
            case('zmax');           read(val,*) c%zmax
            case('zmaxk');          read(val,*) c%zmaxk
            case('zmin');           read(val,*) c%zmin
            case('l2min');          read(val,*) c%l2min
            case('minrec');         read(val,*) c%minrec
            case('minhalf');        read(val,*) c%minhalf
            case('nprec');          read(val,*) c%nprec
            case('npz');            read(val,*) c%npz
            case('nthr_zoom');      read(val,*) c%nthr_zoom
            case('nprec_cons');     read(val,*) c%nprec_cons
            case('chunk_cons');     read(val,*) c%chunk_cons
            case('np_parent');      read(val,*) c%np_parent
            case('minp');           read(val,*) c%minp
            case('minrec_leaf');    read(val,*) c%minrec_leaf
            case('uncap');          read(val,*) c%uncap
            case('ctfabs');         read(val,*) c%ctfabs
            case('pscale');         read(val,*) c%pscale
            case('server');         read(val,*) c%server
            case('kmax_load');      read(val,*) c%kmax_load
            case('minvox');         read(val,*) c%minvox
            case('ndilate');        read(val,*) c%ndilate
            case('skip_refine');    read(val,*) c%skip_refine
            case('norec');          read(val,*) c%norec
            case('mem_job');        read(val,*) c%mem_job
            case('hours_job');      read(val,*) c%hours_job
            case('sim');            read(val,*) c%sim
            case('frac');           read(val,*) c%frac
            case('minmass');        read(val,*) c%minmass
            case('minrep');         read(val,*) c%minrep
            case('balance');        read(val,*) c%balance
            case('devw');           read(val,*) c%devw
            case('mskdiam');        read(val,*) c%mskdiam
            case('mskdiam_refine'); read(val,*) c%mskdiam_refine
            case('lpstart');        read(val,*) c%lpstart
            case('lpstop');         read(val,*) c%lpstop
            case('region_thr');     read(val,*) c%region_thr
            case('cmerge');         read(val,*) c%cmerge
            case('floor');          read(val,*) c%floor
            case('fscmax');         read(val,*) c%fscmax
            case default
                THROW_HARD('unknown klem_tree config key '//trim(key))
        end select
    end subroutine set_key

    subroutine write_cfg( c, fname )
        type(klem_tree_cfg), intent(in) :: c
        character(len=*),    intent(in) :: fname
        integer :: u
        open(newunit=u, file=trim(fname), status='replace', action='write')
        call ws('run', c%run);               call ws('rtag', c%rtag);           call ws('stage', c%stage)
        call ws('kind', c%kind);             call ws('meta', c%meta);           call ws('subset', c%subset)
        call ws('docked', c%docked);         call ws('docked_all', c%docked_all); call ws('cons', c%cons)
        call ws('closed', c%closed);         call ws('tag', c%tag);             call ws('mask', c%mask)
        call ws('partition', c%partition);   call ws('refine_extra', c%refine_extra)
        call ws('lps_parent', c%lps_parent); call ws('sig_parent', c%sig_parent)
        call ws('lps_child', c%lps_child);   call ws('sig_child', c%sig_child)
        call ws('lps_region', c%lps_region); call ws('sig_region', c%sig_region)
        call ws('lps_zoom', c%lps_zoom);     call ws('sig_zoom', c%sig_zoom)
        call wi('fam', c%fam);               call wi('zstate', c%zstate);       call wi('nfam', c%nfam)
        call wi('nparts', c%nparts);         call wi('nthr_part', c%nthr_part); call wi('nthr', c%nthr)
        call wi('np_refine', c%np_refine);   call wi('maxits', c%maxits);       call wi('k_refine', c%k_refine)
        call wi('nfree_parent', c%nfree_parent); call wi('nfree_sub', c%nfree_sub); call wi('maxreg', c%maxreg)
        call wi('nl2', c%nl2);               call wi('parmax', c%parmax);       call wi('zparmax', c%zparmax)
        call wi('kgroups', c%kgroups);       call wi('nmax', c%nmax);           call wi('ntop2', c%ntop2)
        call wi('maxk', c%maxk);             call wi('zmax', c%zmax);           call wi('zmaxk', c%zmaxk)
        call wi('zmin', c%zmin);             call wi('l2min', c%l2min);         call wi('minrec', c%minrec)
        call wi('minhalf', c%minhalf);       call wi('nprec', c%nprec);         call wi('npz', c%npz)
        call wi('nthr_zoom', c%nthr_zoom);   call wi('nprec_cons', c%nprec_cons); call wi('chunk_cons', c%chunk_cons)
        call wi('np_parent', c%np_parent);   call wi('minp', c%minp);           call wi('minrec_leaf', c%minrec_leaf)
        call wi('uncap', c%uncap);           call wi('ctfabs', c%ctfabs);       call wi('pscale', c%pscale)
        call wi('server', c%server);         call wi('kmax_load', c%kmax_load); call wi('minvox', c%minvox)
        call wi('ndilate', c%ndilate);       call wi('skip_refine', c%skip_refine); call wi('norec', c%norec)
        call wi('mem_job', c%mem_job);       call wi('hours_job', c%hours_job)
        call wr('sim', c%sim);               call wr('frac', c%frac);           call wr('minmass', c%minmass)
        call wr('minrep', c%minrep);         call wr('balance', c%balance);     call wr('devw', c%devw)
        call wr('mskdiam', c%mskdiam);       call wr('mskdiam_refine', c%mskdiam_refine)
        call wr('lpstart', c%lpstart);       call wr('lpstop', c%lpstop);       call wr('region_thr', c%region_thr)
        call wr('cmerge', c%cmerge);         call wr('floor', c%floor);         call wr('fscmax', c%fscmax)
        close(u)

      contains

        subroutine ws( k, v )
            character(len=*), intent(in) :: k, v
            if( len_trim(v) > 0 ) write(u,'(A)') k//' '//trim(v)
        end subroutine ws

        subroutine wi( k, v )
            character(len=*), intent(in) :: k
            integer,          intent(in) :: v
            write(u,'(A)') k//' '//int2str(v)
        end subroutine wi

        subroutine wr( k, v )
            character(len=*), intent(in) :: k
            real,             intent(in) :: v
            write(u,'(A,1X,G0)') k, v
        end subroutine wr

    end subroutine write_cfg

    function simple_bin( name ) result( p )
        character(len=*), intent(in) :: name
        character(len=PLEN) :: p
        type(string) :: sp
        integer :: st
        sp = simple_getenv('SIMPLE_PATH', st)
        if( st /= 0 ) THROW_HARD('SIMPLE_PATH is not defined')
        p = sp%to_char()//'/bin/'//trim(name)
    end function simple_bin

    subroutine run_prog( dir, args, logname, nthr, ok )
        character(len=*), intent(in)  :: dir, args, logname
        integer,          intent(in)  :: nthr
        logical,          intent(out) :: ok
        integer :: st
        call exec_cmdline('( cd '//trim(dir)//' && OMP_NUM_THREADS='//int2str(nthr)//' '//trim(simple_bin('simple_exec'))//' '//&
            &trim(args)//' > '//trim(logname)//' 2>&1 )', suppress_errors=.true., exitstat=st)
        ok = st == 0
    end subroutine run_prog

    subroutine copyf( a, b )
        character(len=*), intent(in) :: a, b
        call simple_copy_file(string(trim(a)), string(trim(b)))
    end subroutine copyf

    subroutine touch( fname )
        character(len=*), intent(in) :: fname
        integer :: u
        open(newunit=u, file=trim(fname), status='replace', action='write')
        write(u,'(A)') 'ok'
        close(u)
    end subroutine touch

    function has_line( fname, pat ) result( yes )
        character(len=*), intent(in) :: fname, pat
        logical :: yes
        character(len=4096) :: line
        integer :: u, io
        yes = .false.
        if( .not. file_exists(trim(fname)) ) return
        open(newunit=u, file=trim(fname), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( index(line, pat) > 0 )then
                yes = .true.
                exit
            endif
        enddo
        close(u)
    end function has_line

    function first_error( fname ) result( msg )
        character(len=*), intent(in) :: fname
        character(len=512) :: msg
        character(len=4096) :: line
        integer :: u, io
        character(len=512) :: last
        msg  = 'see '//trim(fname)
        last = ''
        if( .not. file_exists(trim(fname)) ) return
        open(newunit=u, file=trim(fname), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) > 0 .and. index(line, 'NORMAL STOP') == 0 .and. index(line, 'Git Commit') == 0 .and. &
                &index(line, 'Execution time') == 0 ) last = adjustl(line(1:512))
            if( index(line, 'ERROR') > 0 .or. index(line, 'Error') > 0 .or. index(line, 'does not exist') > 0 .or. &
                &index(line, 'Segmentation') > 0 .or. index(line, 'FAILED') > 0 )then
                msg = line(1:512)
                close(u)
                return
            endif
        enddo
        close(u)
        if( len_trim(last) > 0 ) msg = last
    end function first_error

    function latest_dir( parent ) result( d )
        character(len=*), intent(in) :: parent
        character(len=PLEN) :: d
        type(string) :: last
        integer :: nxt
        d = ''
        if( .not. dir_exists(string(trim(parent))) ) return
        nxt = find_next_int_dir_prefix(string(trim(parent)), last)
        if( nxt > 1 ) d = trim(parent)//'/'//last%to_char()
    end function latest_dir

    subroutine set_partition( cfg, proj, dir )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: proj, dir
        logical :: ok
        if( trim(cfg%partition) == 'none' ) return
        call run_prog(dir, 'prg=update_project projfile='//trim(proj)//' qsys_partition='//trim(cfg%partition)//' mkdir=no', 'up_partition.log', 1, ok)
        if( .not. ok ) THROW_HARD('update_project qsys_partition failed on '//trim(proj))
    end subroutine set_partition

    subroutine select_states( proj, sfile, dir )
        character(len=*), intent(in) :: proj, sfile, dir
        logical :: ok
        call run_prog(dir, 'prg=selection projfile='//trim(proj)//' oritype=ptcl3D infile='//trim(sfile)//' mkdir=no', 'sel.log', 1, ok)
        if( .not. ok ) THROW_HARD('selection FAILED in '//trim(dir)//': '//trim(first_error(trim(dir)//'/sel.log')))
    end subroutine select_states

    subroutine readout( cfg, dir, smpd, save )
        type(klem_tree_cfg),        intent(in) :: cfg
        character(len=*),           intent(in) :: dir
        real,                       intent(in) :: smpd
        character(len=*), optional, intent(in) :: save
        character(len=PLEN), allocatable :: f1(:), f2(:), files(:)
        integer :: n1, n2
        if( trim(cfg%docked) == 'none' .or. trim(cfg%cons) == 'none' ) return
        call list_sorted(trim(dir)//'/cell_*.mrc', f1, n1)
        call list_sorted(trim(dir)//'/state_*.mrc', f2, n2)
        if( n1 + n2 == 0 ) return
        allocate(files(n1+n2))
        if( n1 > 0 ) files(1:n1) = f1
        if( n2 > 0 ) files(n1+1:n1+n2) = f2
        write(logfhandle,'(A)') '== '//trim(dir)//' vs the deposited conformers (readout)'
        if( present(save) )then
            call tree_deposited_match(files, n1+n2, cfg%docked, cfg%cons, smpd, 60., 20., save)
        else
            call tree_deposited_match(files, n1+n2, cfg%docked, cfg%cons, smpd, 60., 20.)
        endif
    end subroutine readout

    subroutine run_board( cfg, out, proj, subset, lps, sigmas, nfree, focus, nparts, nthr, mode, slist )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: out, proj, subset, lps, sigmas, focus, slist
        integer,             intent(in) :: nfree, nparts, nthr, mode
        character(len=PLEN) :: pcopy, srvexit
        logical :: ok
        integer :: u, waited, ec
        if( file_exists(trim(out)//'/klem_latent.txt') )then
            call klog('board', 'exists, reused: '//trim(out))
            return
        endif
        call simple_mkdir(trim(out))
        if( trim(cfg%meta) /= 'none' ) call tree_init_positions(cfg%meta, subset, 7, 7, trim(out)//'/init_positions.txt')
        open(newunit=u, file=trim(out)//'/board.txt', status='replace', action='write')
        write(u,'(A)') 'n1 7'
        write(u,'(A)') 'n2 7'
        write(u,'(A)') 'nfree '//int2str(nfree)
        write(u,'(A)') 'lps '//trim(lps)
        write(u,'(A)') 'sigmas '//trim(sigmas)
        write(u,'(A)') 'shrink 0.02'
        write(u,'(A)') 'eps0 0.05'
        write(u,'(A)') 'balance '//fr(real(cfg%balance,dp),1)
        write(u,'(A)') 'subset '//trim(subset)
        write(u,'(A)') 'nsub 0'
        write(u,'(A)') 'seed 1'
        write(u,'(A)') 'conv 0'
        write(u,'(A)') 'ss -1'
        write(u,'(A)') 'fs 1'
        write(u,'(A)') 'sign 1'
        write(u,'(A)') 'ctfabs '//int2str(cfg%ctfabs)
        write(u,'(A)') 'write_halves 0'
        write(u,'(A)') 'devw '//fr(real(cfg%devw,dp),2)
        write(u,'(A)') 'pnoise 1'
        if( trim(focus) /= 'none' ) write(u,'(A)') 'focus_mask '//trim(focus)
        if( trim(cfg%meta) /= 'none' ) write(u,'(A)') 'init_file '//trim(out)//'/init_positions.txt'
        if( nparts > 1 .and. mode /= 2 ) write(u,'(A)') 'nthr_part '//int2str(nthr)
        if( mode == 1 )then
            write(u,'(A)') 'server_list '//trim(slist)
            write(u,'(A)') 'kmax_load '//int2str(cfg%kmax_load)
        endif
        if( cfg%pscale /= 0 ) write(u,'(A)') 'pscale '//int2str(cfg%pscale)
        close(u)
        call klog('board', int2str(nlines(string(trim(subset))))//' particles -> '//trim(out))
        pcopy = trim(out)//'/p.simple'
        if( mode /= 2 )then
            call copyf(proj, pcopy)
            call set_partition(cfg, pcopy, out)
        endif
        call del_file(trim(out)//'/BOARD_DONE')
        select case(mode)
            case(0)
                call run_prog(out, 'prg=klem_board projfile='//trim(pcopy)//' infile='//trim(out)//'/board.txt dir='//trim(out)//&
                    &' nthr='//int2str(nthr)//' nparts='//int2str(max(nparts,1))//' mkdir=no', trim(out)//'/board.log', nthr, ok)
                if( .not. ok ) THROW_HARD('klem_board FAILED in '//trim(out)//': '//trim(first_error(trim(out)//'/board.log')))
            case(1)
                call del_file(trim(out)//'/board.exit')
                open(newunit=u, file=trim(slist), status='replace', action='write')
                close(u)
                call exec_cmdline('( cd '//trim(out)//' && OMP_NUM_THREADS='//int2str(nthr)//' '//trim(simple_bin('simple_exec'))//&
                    &' prg=klem_board projfile='//trim(pcopy)//' infile='//trim(out)//'/board.txt dir='//trim(out)//' nthr='//int2str(nthr)//&
                    &' nparts='//int2str(max(nparts,1))//' mkdir=no > '//trim(out)//'/board.log 2>&1; echo $? > '//trim(out)//&
                    &'/board.exit ) > /dev/null 2>&1 &', suppress_errors=.true.)
                call klog('board', 'board server started ('//int2str(nparts)//' workers, list '//trim(slist)//')')
                waited = 0
                do while( .not. file_exists(trim(out)//'/BOARD_DONE') )
                    if( file_exists(trim(out)//'/board.exit') )then
                        call sleep(5)
                        if( .not. file_exists(trim(out)//'/BOARD_DONE') ) THROW_HARD('board server died in '//trim(out)//': '//trim(first_error(trim(out)//'/board.log')))
                    endif
                    call sleep(10)
                    waited = waited + 10
                    if( waited > 6*3600 ) THROW_HARD('timed out waiting for the board server in '//trim(out))
                enddo
            case(2)
                srvexit = slist(1:index(slist, '/', back=.true.))//'klem_fb/board.exit'
                open(newunit=u, file=trim(slist), status='old', position='append', action='write')
                write(u,'(A)') trim(out)
                close(u)
                call klog('board', 'queued on the board server '//trim(slist))
                waited = 0
                do while( .not. file_exists(trim(out)//'/BOARD_DONE') )
                    if( file_exists(trim(srvexit)) ) THROW_HARD('the board server ended before serving '//trim(out))
                    call sleep(10)
                    waited = waited + 10
                    if( waited > 6*3600 ) THROW_HARD('timed out waiting for the board server for '//trim(out))
                enddo
        end select
        ec = 0
        if( .not. file_exists(trim(out)//'/klem_latent.txt') ) THROW_HARD('klem_board wrote no latent in '//trim(out))
        call klog('board', 'done: '//trim(out))
    end subroutine run_board

    subroutine eo_guard( cfg, proj, sfile, outfile, ns )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: proj, sfile, outfile
        integer,             intent(in) :: ns
        type(sp_project) :: spproj
        integer, allocatable :: st(:), c(:), ce(:), co(:), eo(:)
        integer :: i, n, s, noff
        call read_ints(sfile, st)
        n = size(st)
        call spproj%read_segment('ptcl3D', string(trim(proj)))
        if( spproj%os_ptcl3D%get_noris() /= n ) THROW_HARD('states file and project disagree on the particle count: '//trim(sfile))
        allocate(eo(n))
        do i = 1,n
            eo(i) = spproj%os_ptcl3D%get_eo(i)
        enddo
        call spproj%kill
        s = max(maxval(st), ns, 1)
        allocate(c(0:s), ce(0:s), co(0:s), source=0)
        do i = 1,n
            if( st(i) <= 0 ) cycle
            c(st(i)) = c(st(i)) + 1
            if( eo(i) == 0 ) ce(st(i)) = ce(st(i)) + 1
            if( eo(i) == 1 ) co(st(i)) = co(st(i)) + 1
        enddo
        noff = 0
        do s = 1,size(c)-1
            if( c(s) > 0 .and. (c(s) < cfg%minrec .or. ce(s) < cfg%minhalf .or. co(s) < cfg%minhalf) )then
                call klog('rec', 'state '//int2str(s)//' switched off ('//int2str(c(s))//' particles, '//int2str(ce(s))//' even, '//&
                    &int2str(co(s))//' odd)')
                noff = noff + 1
            endif
        enddo
        do i = 1,n
            s = st(i)
            if( s > 0 )then
                if( c(s) < cfg%minrec .or. ce(s) < cfg%minhalf .or. co(s) < cfg%minhalf ) st(i) = 0
            endif
        enddo
        call write_ints(outfile, st)
    end subroutine eo_guard

    subroutine rec_states( cfg, proj, sfile, ns, dir, np, nthr, local )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: proj, sfile, dir
        integer,             intent(in) :: ns, np, nthr
        logical,             intent(in) :: local
        character(len=PLEN) :: rdir, rr, f
        logical :: ok
        integer :: s, h
        character(len=8) :: suf(3)
        suf  = ['        ', '_even   ', '_odd    ']
        rdir = trim(dir)//'/rec'
        call simple_mkdir(trim(rdir))
        call simple_mkdir(trim(dir)//'/maps')
        call simple_mkdir(trim(dir)//'/states_dir')
        call copyf(proj, trim(rdir)//'/p.simple')
        call set_partition(cfg, trim(rdir)//'/p.simple', rdir)
        call eo_guard(cfg, proj, sfile, trim(rdir)//'/states_rec.txt', ns)
        if( local )then
            call run_prog(rdir, 'prg=update_project projfile=p.simple qsys_name=local nparts='//int2str(np)//' nthr='//int2str(nthr)//&
                &' mkdir=no', 'up.log', 1, ok)
            if( .not. ok ) THROW_HARD('update_project qsys_name=local FAILED in '//trim(rdir))
        endif
        call select_states('p.simple', trim(rdir)//'/states_rec.txt', rdir)
        call run_prog(rdir, 'prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam='//int2str(nint(cfg%mskdiam))//' nstates='//int2str(ns)//&
            &' nparts='//int2str(np)//' nthr='//int2str(nthr), 'rec.log', nthr, ok)
        if( .not. ok )then
            call klog('rec', 'reconstruct3D failed once ('//trim(dir)//'), retrying')
            call sleep(30)
            call run_prog(rdir, 'prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam='//int2str(nint(cfg%mskdiam))//' nstates='//int2str(ns)//&
                &' nparts='//int2str(np)//' nthr='//int2str(nthr), 'rec2.log', nthr, ok)
            if( .not. ok ) THROW_HARD('reconstruct3D FAILED twice in '//trim(dir)//': '//trim(first_error(trim(rdir)//'/rec2.log')))
        endif
        rr = latest_dir(rdir)
        if( len_trim(rr) == 0 ) THROW_HARD('no reconstruct3D output in '//trim(rdir))
        do s = 1,ns
            do h = 1,3
                f = trim(rr)//'/recvol_state'//int2str_pad(s,2)//trim(suf(h))//'.mrc'
                if( file_exists(trim(f)) ) call copyf(f, trim(dir)//'/maps/recvol_state'//int2str_pad(s,2)//trim(suf(h))//'.mrc')
            enddo
            f = trim(rr)//'/recvol_state'//int2str_pad(s,2)//'.mrc'
            if( file_exists(trim(f)) ) call copyf(f, trim(dir)//'/states_dir/state_'//int2str_pad(s,2)//'.mrc')
        enddo
    end subroutine rec_states

    function refined_proj( D, src ) result( ps )
        character(len=*), intent(in) :: D, src
        character(len=PLEN) :: ps, rd
        ps = src
        rd = latest_dir(trim(D)//'/refine')
        if( len_trim(rd) > 0 )then
            if( file_exists(trim(rd)//'/p.simple') ) ps = trim(rd)//'/p.simple'
        endif
    end function refined_proj

    function tag_mask( D, tag ) result( mk )
        character(len=*), intent(in) :: D, tag
        character(len=PLEN) :: mk
        if( index(tag, 'motion_') == 1 )then
            mk = trim(D)//'/regions_motion/region_'//trim(tag(8:))//'.mrc'
        else
            mk = trim(D)//'/regions/region'//trim(tag)//'.mrc'
        endif
    end function tag_mask

    subroutine stage_job( cfg, src, dir, label, stage, fam, tag, mask, zstate, ncpu, qenv, job )
        type(klem_tree_cfg),  intent(in)    :: cfg
        character(len=*),     intent(in)    :: src, dir, label, stage, tag, mask
        integer,              intent(in)    :: fam, zstate, ncpu
        type(qsys_env),       intent(inout) :: qenv
        type(qsys_async_job), intent(inout) :: job
        type(klem_tree_cfg) :: c
        type(cmdline)       :: cl
        c        = cfg
        c%stage  = stage
        c%fam    = fam
        c%tag    = tag
        c%mask   = mask
        c%zstate = zstate
        call simple_mkdir(trim(dir))
        call simple_mkdir(trim(dir)//'/stderrout')
        call write_cfg(c, trim(dir)//'/klem_tree_'//trim(label)//'.txt')
        call cl%set('prg', 'klem_tree')
        call cl%set('projfile', trim(src))
        call cl%set('infile', trim(dir)//'/klem_tree_'//trim(label)//'.txt')
        call cl%set('nthr', ncpu)
        call cl%set('mkdir', 'no')
        call qenv%qdescr%set('job_cpus_per_task', int2str(ncpu))
        call qenv%qdescr%set('job_memory_per_task', int2str(cfg%mem_job))
        call qenv%qdescr%set('job_time', '0-'//int2str(cfg%hours_job)//':0:0')
        call qenv%qdescr%set('job_name', 'klem_'//trim(cfg%rtag)//'_'//trim(label))
        call job%start(qenv, cl, string(trim(dir)), trim(label))
        call klog('jobs', trim(label)//' submitted (job dir '//trim(dir)//')')
        call cl%kill
    end subroutine stage_job

    subroutine run_stage_jobs( params, cfg, src, dirs, labels, stages, fams, tags, masks, zstates, nj, parmax, ncpu, what )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: src, what
        character(len=PLEN), intent(in)    :: dirs(:), labels(:), stages(:), tags(:), masks(:)
        integer,             intent(in)    :: fams(:), zstates(:), nj, parmax, ncpu
        type(qsys_env) :: qenv
        type(qsys_async_job), allocatable :: jobs(:)
        type(string) :: lg
        logical, allocatable :: started(:), finished(:), failed(:)
        integer :: i, st, nrun, nfail
        if( nj == 0 ) return
        if( trim(cfg%partition) /= 'none' )then
            call qenv%new(params, 1, qsys_nthr=ncpu, qsys_partition=string(trim(cfg%partition)))
        else
            call qenv%new(params, 1, qsys_nthr=ncpu)
        endif
        allocate(jobs(nj), started(nj), finished(nj), failed(nj))
        started  = .false.
        finished = .false.
        failed   = .false.
        do i = 1,nj
            if( file_exists(trim(dirs(i))//'/DONE') )then
                started(i)  = .true.
                finished(i) = .true.
            endif
        enddo
        call klog('jobs', what//': '//int2str(count(.not. finished))//' of '//int2str(nj)//' to run, '//int2str(parmax)//' at a time')
        do
            do i = 1,nj
                if( .not. started(i) .or. finished(i) ) cycle
                st = jobs(i)%status()
                if( st == ASYNC_JOB_DONE )then
                    finished(i) = .true.
                    if( .not. file_exists(trim(dirs(i))//'/DONE') )then
                        failed(i) = .true.
                        call klog('jobs', trim(labels(i))//' FAILED: exited without its DONE marker ('//trim(dirs(i))//')')
                    else
                        call klog('jobs', trim(labels(i))//' done')
                    endif
                else if( st == ASYNC_JOB_FAILED )then
                    finished(i) = .true.
                    failed(i)   = .true.
                    lg = jobs(i)%get_log()
                    call klog('jobs', trim(labels(i))//' FAILED (log '//lg%to_char()//')')
                endif
            enddo
            nrun = count(started .and. .not. finished)
            do i = 1,nj
                if( nrun >= parmax ) exit
                if( started(i) ) cycle
                call stage_job(cfg, src, dirs(i), labels(i), stages(i), fams(i), tags(i), masks(i), zstates(i), ncpu, qenv, jobs(i))
                started(i) = .true.
                nrun = nrun + 1
            enddo
            if( all(finished) ) exit
            call sleep(30)
        enddo
        nfail = count(failed)
        call qenv%kill
        if( nfail > 0 ) THROW_HARD(what//': '//int2str(nfail)//' job(s) FAILED')
        call klog('jobs', what//' finished')
    end subroutine run_stage_jobs

    subroutine stage_all( params, cfg, R, src, nrows, box, smpd )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src
        integer,             intent(in)    :: nrows, box
        real,                intent(in)    :: smpd
        type(board_t) :: B
        character(len=PLEN), allocatable :: dirs(:), labels(:), stages(:), tags(:), masks(:), leaves(:), rtags(:)
        integer,             allocatable :: fams(:), zst(:), nmem(:)
        character(len=PLEN) :: psub, D, RB
        character(len=4096) :: line
        integer :: j, n, nreg, nl, nt, i, k, u, io, nz, s
        real    :: sims(4)
        psub = trim(R)//'/subset_all.txt'
        if( trim(cfg%subset) /= 'none' )then
            call copyf(cfg%subset, psub)
            call klog('f692', 'parent subset from '//trim(cfg%subset)//': '//int2str(nlines(string(trim(psub))))//' particles')
        else
            call tree_parent_subset(cfg%meta, psub, n)
        endif
        call run_board(cfg, trim(R)//'/parent', src, psub, cfg%lps_parent, cfg%sig_parent, cfg%nfree_parent, 'none', cfg%nparts, cfg%nthr_part, 0, '')
        call readout(cfg, trim(R)//'/parent/cells', smpd)
        call simple_mkdir(trim(R)//'/expand')
        allocate(nmem(cfg%nfam))
        if( .not. file_exists(trim(R)//'/expand/family1_subset.txt') ) call tree_families(trim(R)//'/parent', cfg%nfam, trim(R)//'/expand', nmem)
        if( .not. file_exists(trim(R)//'/regions_parent/region_1.mrc') )then
            call read_board(trim(R)//'/parent', B)
            call tree_regions(B, 'parent', box, smpd, 1, trim(R)//'/regions_parent', 2, cfg%kgroups, cfg%region_thr, cfg%minvox, cfg%ndilate, cfg%cmerge, cfg%minrep, nreg)
            call kill_board(B)
        endif
        allocate(dirs(cfg%nfam), labels(cfg%nfam), stages(cfg%nfam), tags(cfg%nfam), masks(cfg%nfam), fams(cfg%nfam), zst(cfg%nfam))
        do j = 1,cfg%nfam
            dirs(j)   = trim(R)//'/fam'//int2str(j)
            labels(j) = 'fam'//int2str(j)
            stages(j) = 'family'
            tags(j)   = ''
            masks(j)  = 'none'
            fams(j)   = j
            zst(j)    = 0
        enddo
        call run_stage_jobs(params, cfg, src, dirs, labels, stages, fams, tags, masks, zst, cfg%nfam, cfg%nfam, max(cfg%nthr,8), 'families')
        deallocate(dirs, labels, stages, tags, masks, fams, zst)
        allocate(leaves(cfg%nfam * (cfg%maxreg + 1)))
        nl = 0
        do j = 1,cfg%nfam
            D = trim(R)//'/fam'//int2str(j)
            if( file_exists(trim(D)//'/klem_fb/klem_halves.bin') )then
                nl = nl + 1
                leaves(nl) = trim(D)//'/klem_fb'
            endif
        enddo
        do j = 1,cfg%nfam
            D = trim(R)//'/fam'//int2str(j)
            call read_tags(D, rtags, nt)
            do i = 1,nt
                RB = trim(D)//'/region'//trim(rtags(i))
                if( .not. file_exists(trim(RB)//'/board/klem_halves.bin') ) cycle
                nl = nl + 1
                leaves(nl) = trim(RB)//'/board'
                if( .not. file_exists(trim(RB)//'/board/cell_matches.txt') ) call readout(cfg, trim(RB)//'/board/cells', smpd, trim(RB)//'/board/cell_matches.txt')
            enddo
        enddo
        sims = [0.6, 0.7, 0.8, 0.9]
        if( nl > 0 )then
            write(logfhandle,'(A)') '== LANDSCAPE over the leaves:'
            call tree_count_states(leaves(1:nl), nl, cfg%minrep, sims, 4)
        endif
        nz = 0
        do k = 1,2
            do j = 1,cfg%nfam
                D = trim(R)//'/fam'//int2str(j)
                call read_tags(D, rtags, nt)
                do i = 1,nt
                    RB = trim(D)//'/region'//trim(rtags(i))
                    if( .not. file_exists(trim(RB)//'/level1/clusters.txt') ) cycle
                    open(newunit=u, file=trim(RB)//'/level1/clusters.txt', status='old', action='read')
                    do
                        read(u,'(A)',iostat=io) line
                        if( io /= 0 ) exit
                        if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                        read(line,*) s
                        nz = nz + 1
                        if( k == 2 )then
                            dirs(nz)   = trim(D)//'/zoomregion'//trim(rtags(i))//'_s'//int2str_pad(s,2)
                            labels(nz) = 'zoom_f'//int2str(j)//trim(rtags(i))//'_s'//int2str_pad(s,2)
                            stages(nz) = 'zoom'
                            tags(nz)   = rtags(i)
                            masks(nz)  = tag_mask(D, rtags(i))
                            fams(nz)   = j
                            zst(nz)    = s
                        endif
                    enddo
                    close(u)
                enddo
            enddo
            if( k == 1 )then
                allocate(dirs(max(nz,1)), labels(max(nz,1)), stages(max(nz,1)), tags(max(nz,1)), masks(max(nz,1)), fams(max(nz,1)), zst(max(nz,1)))
                nz = 0
            endif
        enddo
        call run_stage_jobs(params, cfg, src, dirs, labels, stages, fams, tags, masks, zst, nz, cfg%zparmax, cfg%nthr_zoom, 'zooms')
        deallocate(dirs, labels, stages, tags, masks, fams, zst)
        allocate(dirs(3), labels(3), stages(3), tags(3), masks(3), fams(3), zst(3))
        dirs(1) = trim(R)//'/pool_region'
        dirs(2) = trim(R)//'/pool_zoom'
        dirs(3) = trim(R)//'/parents'
        labels  = [character(len=PLEN) :: 'pool_region', 'pool_zoom', 'parents']
        stages  = [character(len=PLEN) :: 'pool', 'pool', 'parents']
        tags    = [character(len=PLEN) :: 'region', 'zoom', '']
        masks   = 'none'
        fams    = 0
        zst     = 0
        call run_stage_jobs(params, cfg, src, dirs, labels, stages, fams, tags, masks, zst, 3, 3, 16, 'consolidation and parents')
        do k = 1,2
            D = trim(R)//'/consolidated_'//fr(real(cfg%floor,dp),0)//'_'//pool_tag(cfg, trim(tags(k)))
            if( .not. file_exists(trim(D)//'/groups.txt') ) cycle
            write(logfhandle,'(A)') '== '//trim(D)//' composed with the per-family parents'
            call tree_compose(R, D, trim(R)//'/parents', smpd, trim(D)//'/composed')
        enddo
        call klog('f692', 'ALL DONE ('//trim(R)//')')
    end subroutine stage_all

    function pool_tag( cfg, kind ) result( t )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: kind
        character(len=PLEN) :: t
        if( cfg%uncap /= 0 )then
            t = 'fullw_'//trim(kind)
        else
            t = 'capped_'//trim(kind)
        endif
    end function pool_tag

    subroutine read_tags( D, tags, nt )
        character(len=*),                 intent(in)    :: D
        character(len=PLEN), allocatable, intent(inout) :: tags(:)
        integer,                          intent(out)   :: nt
        character(len=4096) :: line
        integer :: u, io
        nt = 0
        if( allocated(tags) ) deallocate(tags)
        allocate(tags(64))
        if( .not. file_exists(trim(D)//'/regions.txt') ) return
        open(newunit=u, file=trim(D)//'/regions.txt', status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            nt = nt + 1
            tags(nt) = line(1:index(line, ' ')-1)
        enddo
        close(u)
    end subroutine read_tags

    subroutine stage_family( params, cfg, R, src, nrows, box, smpd )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src
        integer,             intent(in)    :: nrows, box
        real,                intent(in)    :: smpd
        type(board_t) :: B
        character(len=PLEN), allocatable :: m1(:), m2(:), dirs(:), labels(:), stages(:), tags(:), masks(:)
        integer,             allocatable :: fams(:), zst(:)
        character(len=PLEN) :: D, sub, focus, ps, rd, slist, f, tg
        logical :: ok, served
        integer :: fam, nreg, n1, n2, nm, i, u, s, waited
        fam   = cfg%fam
        D     = trim(R)//'/fam'//int2str(fam)
        sub   = trim(R)//'/expand/family'//int2str(fam)//'_subset.txt'
        slist = trim(D)//'/server_boards.txt'
        call simple_mkdir(trim(D))
        call klog('fam'//int2str(fam), 'family '//int2str(fam)//' stage')
        focus = 'none'
        if( file_exists(trim(R)//'/regions_parent/region_1.mrc') ) focus = trim(R)//'/regions_parent/region_1.mrc'
        call run_board(cfg, trim(D)//'/child', src, sub, cfg%lps_child, cfg%sig_child, cfg%nfree_sub, focus, cfg%nparts, cfg%nthr_part, 0, '')
        call readout(cfg, trim(D)//'/child/cells', smpd)
        if( .not. file_exists(trim(D)//'/states.txt') )then
            call read_board(trim(D)//'/child', B, with_halves=.false.)
            call tree_kmeans_states(B, cfg%k_refine, nrows, trim(D)//'/states.txt', 'family '//int2str(fam))
            call kill_board(B)
        endif
        ps = refined_proj(D, src)
        if( trim(ps) == trim(src) .and. cfg%skip_refine == 0 )then
            call simple_mkdir(trim(D)//'/refine')
            call copyf(src, trim(D)//'/refine/p.simple')
            call set_partition(cfg, trim(D)//'/refine/p.simple', trim(D)//'/refine')
            call select_states('p.simple', trim(D)//'/states.txt', trim(D)//'/refine')
            call klog('fam'//int2str(fam), 'refine3D_states nstates='//int2str(cfg%k_refine)//' local '//fr(real(cfg%lpstart,dp),0)//'->'//&
                &fr(real(cfg%lpstop,dp),0)//' A maxits '//int2str(cfg%maxits)//' nparts '//int2str(cfg%np_refine)//' x '//int2str(cfg%nthr_part)//' threads')
            call run_prog(trim(D)//'/refine', 'prg=refine3D_states projfile=p.simple nstates='//int2str(cfg%k_refine)//' pose_policy=local lpstart='//&
                &fr(real(cfg%lpstart,dp),0)//' lpstop='//fr(real(cfg%lpstop,dp),0)//' maxits='//int2str(cfg%maxits)//' mskdiam='//&
                &int2str(nint(cfg%mskdiam_refine))//' pgrp=c1 nparts='//int2str(cfg%np_refine)//' nthr='//int2str(cfg%nthr_part)//' '//&
                &trim(cfg%refine_extra), 'refine.log', cfg%nthr_part, ok)
            if( .not. ok .or. .not. has_line(trim(D)//'/refine/refine.log', '>>> ITERATION') )then
                THROW_HARD('refine3D_states FAILED for family '//int2str(fam)//': '//trim(first_error(trim(D)//'/refine/refine.log')))
            endif
            ps = refined_proj(D, src)
            if( trim(ps) == trim(src) ) THROW_HARD('refine3D_states wrote no refined project for family '//int2str(fam))
        endif
        if( trim(ps) /= trim(src) )then
            rd = latest_dir(trim(D)//'/refine')
            call simple_mkdir(trim(D)//'/maps')
            call simple_mkdir(trim(D)//'/states_dir')
            do s = 1,cfg%k_refine
                f = trim(rd)//'/recvol_state'//int2str_pad(s,2)//'.mrc'
                if( .not. file_exists(trim(f)) ) cycle
                call copyf(f, trim(D)//'/maps/recvol_state'//int2str_pad(s,2)//'.mrc')
                call copyf(f, trim(D)//'/states_dir/state_'//int2str_pad(s,2)//'.mrc')
            enddo
            call klog('fam'//int2str(fam), 'refined project '//trim(ps))
            call readout(cfg, trim(D)//'/states_dir', smpd)
        else
            call klog('fam'//int2str(fam), 'skip_refine: feedback and region boards at the consensus poses')
        endif
        if( .not. file_exists(trim(D)//'/regions_child/region_1.mrc') )then
            call read_board(trim(D)//'/child', B)
            call tree_regions(B, 'child', box, smpd, 1, trim(D)//'/regions_child', 2, cfg%kgroups, cfg%region_thr, cfg%minvox, cfg%ndilate, cfg%cmerge, cfg%minrep, nreg)
            call kill_board(B)
        endif
        focus = 'none'
        if( file_exists(trim(D)//'/regions_child/region_1.mrc') ) focus = trim(D)//'/regions_child/region_1.mrc'
        served = .false.
        if( .not. file_exists(trim(D)//'/klem_fb/klem_latent.txt') .and. cfg%server /= 0 )then
            call run_board(cfg, trim(D)//'/klem_fb', ps, sub, cfg%lps_child, cfg%sig_child, cfg%nfree_sub, focus, cfg%nparts, cfg%nthr_part, 1, slist)
            served = .true.
        else
            call del_file(trim(slist))
            call run_board(cfg, trim(D)//'/klem_fb', ps, sub, cfg%lps_child, cfg%sig_child, cfg%nfree_sub, focus, cfg%nparts, cfg%nthr_part, 0, '')
        endif
        call readout(cfg, trim(D)//'/klem_fb/cells', smpd, trim(D)//'/klem_fb/cell_matches.txt')
        call read_board(trim(D)//'/klem_fb', B)
        if( .not. file_exists(trim(D)//'/regions/region_1.mrc') ) call tree_regions(B, 'klem_fb', box, smpd, 1, trim(D)//'/regions', 2, cfg%kgroups, &
            &cfg%region_thr, cfg%minvox, cfg%ndilate, cfg%cmerge, cfg%minrep, nreg)
        if( .not. file_exists(trim(D)//'/regions_motion/region_1.mrc') ) call tree_regions(B, 'klem_fb', box, smpd, 2, trim(D)//'/regions_motion', 6, &
            &cfg%kgroups, cfg%region_thr, cfg%minvox, cfg%ndilate, cfg%cmerge, cfg%minrep, nreg)
        call kill_board(B)
        call list_sorted(trim(D)//'/regions/region_?.mrc', m1, n1)
        call list_sorted(trim(D)//'/regions_motion/region_?.mrc', m2, n2)
        nm = min(n1 + n2, cfg%maxreg)
        allocate(dirs(max(nm,1)), labels(max(nm,1)), stages(max(nm,1)), tags(max(nm,1)), masks(max(nm,1)), fams(max(nm,1)), zst(max(nm,1)))
        open(newunit=u, file=trim(D)//'/regions.txt', status='replace', action='write')
        do i = 1,nm
            if( i <= n1 )then
                masks(i) = m1(i)
                f  = m1(i)(index(m1(i), '/region_', back=.true.)+7:)
                tg = f(1:index(f, '.mrc')-1)
            else
                masks(i) = m2(i-n1)
                f  = m2(i-n1)(index(m2(i-n1), '/region_', back=.true.)+8:)
                tg = 'motion_'//f(1:index(f, '.mrc')-1)
            endif
            tags(i)   = tg
            dirs(i)   = trim(D)//'/region'//trim(tg)
            labels(i) = 'reg'//int2str(fam)//trim(tg)
            stages(i) = 'region'
            fams(i)   = fam
            zst(i)    = 0
            write(u,'(A)') trim(tg)//' '//trim(masks(i))
        enddo
        close(u)
        call klog('regpar'//int2str(fam), int2str(nm)//' region masks (up to '//int2str(cfg%parmax)//' jobs at a time)')
        call run_stage_jobs(params, cfg, src, dirs, labels, stages, fams, tags, masks, zst, nm, cfg%parmax, max(cfg%nthr,8), 'regions of family '//int2str(fam))
        if( served )then
            open(newunit=u, file=trim(slist), status='old', position='append', action='write')
            write(u,'(A)') 'STOP'
            close(u)
            waited = 0
            do while( .not. file_exists(trim(D)//'/klem_fb/board.exit') )
                call sleep(10)
                waited = waited + 10
                if( waited > 3600 ) exit
            enddo
            call touch(trim(slist)//'.exit')
            call klog('fam'//int2str(fam), 'board server finished')
        endif
        call touch(trim(D)//'/DONE')
        call klog('fam'//int2str(fam), 'DONE')
    end subroutine stage_family

    subroutine stage_region( cfg, R, src, nrows, smpd )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: R, src
        integer,             intent(in) :: nrows
        real,                intent(in) :: smpd
        type(board_t) :: B
        integer, allocatable :: st(:)
        character(len=PLEN) :: D, RB, sub, ps, slist, L
        logical :: served
        integer :: ks, s, n, n2, i, u
        D     = trim(R)//'/fam'//int2str(cfg%fam)
        RB    = trim(D)//'/region'//trim(cfg%tag)
        sub   = trim(R)//'/expand/family'//int2str(cfg%fam)//'_subset.txt'
        ps    = refined_proj(D, src)
        slist = trim(D)//'/server_boards.txt'
        call simple_mkdir(trim(RB))
        served = file_exists(trim(slist)) .and. .not. file_exists(trim(slist)//'.exit') .and. .not. file_exists(trim(D)//'/klem_fb/board.exit')
        call klog('reg'//int2str(cfg%fam)//'/'//trim(cfg%tag), '(1) focal board at the refined poses, mask '//trim(cfg%mask))
        if( served )then
            call run_board(cfg, trim(RB)//'/board', ps, sub, cfg%lps_region, cfg%sig_region, cfg%nfree_sub, cfg%mask, cfg%nparts, cfg%nthr_part, 2, slist)
        else
            call run_board(cfg, trim(RB)//'/board', ps, sub, cfg%lps_region, cfg%sig_region, cfg%nfree_sub, cfg%mask, cfg%nparts, cfg%nthr_part, 0, '')
        endif
        call klog('reg'//int2str(cfg%fam)//'/'//trim(cfg%tag), '(2) states from the reproducible cell clusters')
        call read_board(trim(RB)//'/board', B)
        call tree_node_states(B, 'board', nrows, cfg%sim, cfg%nmax, cfg%frac, cfg%minmass, cfg%maxk, cfg%minrep, .true., trim(RB)//'/level1', ks)
        call kill_board(B)
        if( ks < 1 .or. cfg%norec /= 0 )then
            call touch(trim(RB)//'/DONE')
            return
        endif
        if( .not. file_exists(trim(RB)//'/level1/maps/recvol_state01.mrc') ) call rec_states(cfg, ps, trim(RB)//'/level1/states_clusters.txt', ks, &
            &trim(RB)//'/level1', cfg%nprec, 8, .false.)
        call readout(cfg, trim(RB)//'/level1/states_dir', smpd, trim(RB)//'/level1/state_matches.txt')
        call klog('reg'//int2str(cfg%fam)//'/'//trim(cfg%tag), '(3) level 2 for the largest '//int2str(cfg%nl2)//' states (inherited region, own main axis)')
        call read_ints(trim(RB)//'/level1/states_clusters.txt', st)
        do s = 1,min(cfg%nl2, ks)
            L = trim(RB)//'/level2_s'//int2str(s)
            call simple_mkdir(trim(L))
            n = count(st == s)
            if( n < cfg%l2min )then
                call klog('reg'//int2str(cfg%fam)//'/'//trim(cfg%tag), 'state '//int2str(s)//': '//int2str(n)//' particles, no level 2')
                cycle
            endif
            open(newunit=u, file=trim(L)//'/subset.txt', status='replace', action='write')
            do i = 1,size(st)
                if( st(i) == s ) write(u,'(I0)') i
            enddo
            close(u)
            call run_board(cfg, trim(L)//'/board', ps, trim(L)//'/subset.txt', cfg%lps_region, cfg%sig_region, 0, cfg%mask, 1, cfg%nthr, 0, '')
            n2 = min(n * 4 / 10, cfg%ntop2)
            call read_board(trim(L)//'/board', B)
            call tree_rank_pc1(B, n2, nrows, trim(L))
            call kill_board(B)
            if( .not. file_exists(trim(L)//'/maps/recvol_state01.mrc') ) call rec_states(cfg, ps, trim(L)//'/states_top.txt', 1, trim(L), cfg%nprec, 8, .false.)
            call readout(cfg, trim(L)//'/states_dir', smpd)
        enddo
        call touch(trim(RB)//'/DONE')
        call klog('reg'//int2str(cfg%fam)//'/'//trim(cfg%tag), 'DONE')
    end subroutine stage_region

    subroutine stage_zoom( cfg, R, src, nrows )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: R, src
        integer,             intent(in) :: nrows
        type(board_t) :: B
        integer, allocatable :: st(:)
        character(len=PLEN) :: D, RB, Z, L, ps, rr
        character(len=64)   :: zt
        type(string), allocatable :: parts(:)
        integer :: n, k, i, u
        D  = trim(R)//'/fam'//int2str(cfg%fam)
        RB = trim(D)//'/region'//trim(cfg%tag)
        Z  = trim(D)//'/zoomregion'//trim(cfg%tag)//'_s'//int2str_pad(cfg%zstate,2)
        L  = trim(Z)//'/level1'
        zt = 'zoom f'//int2str(cfg%fam)//'/'//trim(cfg%tag)//'/s'//int2str(cfg%zstate)
        if( file_exists(trim(Z)//'/DONE') ) return
        ps = refined_proj(D, src)
        call simple_mkdir(trim(L))
        call read_ints(trim(RB)//'/level1/states_clusters.txt', st)
        n = count(st == cfg%zstate)
        if( n < cfg%zmin )then
            call klog(zt, int2str(n)//' particles, no zoom')
            call touch(trim(Z)//'/DONE')
            return
        endif
        open(newunit=u, file=trim(Z)//'/subset.txt', status='replace', action='write')
        do i = 1,size(st)
            if( st(i) == cfg%zstate ) write(u,'(I0)') i
        enddo
        close(u)
        call klog(zt, 'zoom board on '//int2str(n)//' particles (focus '//trim(cfg%mask)//')')
        call run_board(cfg, trim(Z)//'/board', ps, trim(Z)//'/subset.txt', cfg%lps_zoom, cfg%sig_zoom, 1, cfg%mask, 1, cfg%nthr_zoom, 0, '')
        call read_board(trim(Z)//'/board', B)
        call tree_node_states(B, 'board', nrows, cfg%sim, cfg%zmax, cfg%frac, cfg%minmass, cfg%zmaxk, cfg%minrep, .true., trim(L), k)
        call kill_board(B)
        if( k < 1 )then
            call klog(zt, 'no supported sub-states')
            call touch(trim(Z)//'/DONE')
            return
        endif
        call klog(zt, int2str(k)//' sub-states')
        if( .not. file_exists(trim(L)//'/maps/recvol_state01_even.mrc') ) call rec_states(cfg, ps, trim(L)//'/states_clusters.txt', k, trim(L), &
            &cfg%npz, cfg%nthr_zoom, .true.)
        rr = latest_dir(trim(L)//'/rec')
        if( len_trim(rr) > 0 )then
            call simple_list_files(trim(rr)//'/*_part*.mrc', parts)
            if( allocated(parts) )then
                do i = 1,size(parts)
                    call del_file(parts(i))
                enddo
            endif
        endif
        call touch(trim(Z)//'/DONE')
        call klog(zt, 'DONE ('//int2str(k)//' sub-states)')
    end subroutine stage_zoom

    subroutine stage_pool( cfg, R, src, nrows, smpd )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: R, src
        integer,             intent(in) :: nrows
        real,                intent(in) :: smpd
        type(leaf_t), allocatable :: leaves(:)
        integer,      allocatable :: keepl(:), lab(:), st(:), sub(:), loc(:), nf(:,:)
        logical,      allocatable :: inf(:)
        character(len=PLEN), allocatable :: sfiles(:), names(:)
        character(len=PLEN) :: tag, C, Lc, out, ps, qd, rr, f, figs
        character(len=4096) :: line
        character(len=8)    :: suf(3), hdir(3)
        logical :: ok
        integer :: nl, nkeep, K, fam, q, nch, lo, hi, ns, i, s, u, io, best, h, nsf, sl, np_, cnt
        suf  = ['        ', '_even   ', '_odd    ']
        hdir = ['        ', 'halves/ ', 'halves/ ']
        tag = pool_tag(cfg, cfg%tag)
        out = trim(R)//'/merge_'//trim(tag)//'_'//fr(real(cfg%floor,dp),0)
        C   = trim(R)//'/consolidated_'//fr(real(cfg%floor,dp),0)//'_'//trim(tag)
        Lc  = trim(C)//'/fam0/regionall/level1'
        if( trim(cfg%cons) == 'none' ) THROW_HARD('the data-derived merge needs cons (the consensus map)')
        if( .not. file_exists(trim(C)//'/DONE') )then
            call build_leaves(R, cfg%tag, cfg%nfam, leaves, nl)
            write(logfhandle,'(A)') '== data-derived merge of the '//trim(cfg%tag)//' leaves of '//trim(R)//', floor '//fr(real(cfg%floor,dp),0)//' A'
            call tree_merge_leaves(R, cfg%cons, leaves, nl, cfg%floor, smpd, cfg%fscmax, cfg%minrec_leaf, out, keepl, lab, nkeep)
            if( nkeep == 0 .and. trim(cfg%tag) == 'zoom' )then
                call klog('pool', 'no zoom leaf states to merge')
                call touch(trim(R)//'/pool_'//trim(cfg%tag)//'/DONE')
                return
            endif
            if( nkeep == 0 ) THROW_HARD('no leaf states to merge ('//trim(cfg%tag)//')')
            call tree_consolidate(leaves, keepl, lab, nkeep, nrows, cfg%minp, cfg%uncap /= 0, cfg%floor, C, K)
            call simple_mkdir(trim(Lc)//'/maps/halves')
            call simple_mkdir(trim(Lc)//'/states_dir')
            call read_ints(trim(C)//'/states.txt', st)
            allocate(nf(cfg%nfam,K), source=0)
            nch = (K + cfg%chunk_cons - 1) / cfg%chunk_cons
            do fam = 1,cfg%nfam
                call read_ints(trim(R)//'/expand/family'//int2str(fam)//'_subset.txt', sub)
                allocate(inf(nrows), source=.false.)
                do i = 1,size(sub)
                    if( sub(i) >= 1 .and. sub(i) <= nrows ) inf(sub(i)) = .true.
                enddo
                do i = 1,nrows
                    if( inf(i) .and. st(i) > 0 ) nf(fam,st(i)) = nf(fam,st(i)) + 1
                enddo
                ps = refined_proj(trim(R)//'/fam'//int2str(fam), src)
                call klog('consolidate', 'family '//int2str(fam)//' poses from '//trim(ps)//'; '//int2str(sum(nf(fam,:)))//' particles in merged states')
                allocate(loc(nrows))
                do q = 0,nch-1
                    lo = q * cfg%chunk_cons
                    hi = (q + 1) * cfg%chunk_cons
                    ns = min(cfg%chunk_cons, K - lo)
                    qd = trim(C)//'/rec'//int2str(fam)//'/q'//int2str(q)
                    if( file_exists(trim(qd)//'/REC_DONE') ) cycle
                    do i = 1,nrows
                        loc(i) = 0
                        if( inf(i) .and. st(i) > lo .and. st(i) <= hi ) loc(i) = st(i) - lo
                    enddo
                    call simple_mkdir(trim(qd))
                    call write_ints(trim(C)//'/states_f'//int2str(fam)//'_q'//int2str(q)//'.txt', loc)
                    call copyf(ps, trim(qd)//'/p.simple')
                    call set_partition(cfg, trim(qd)//'/p.simple', qd)
                    call eo_guard(cfg, ps, trim(C)//'/states_f'//int2str(fam)//'_q'//int2str(q)//'.txt', trim(qd)//'/states_rec.txt', ns)
                    call select_states('p.simple', trim(qd)//'/states_rec.txt', qd)
                    call klog('consolidate', 'reconstruct3D family '//int2str(fam)//': states '//int2str(lo+1)//'-'//int2str(lo+ns)//' of '//&
                        &int2str(K)//' (cut '//fr(real(cfg%floor,dp),0)//')')
                    call run_prog(qd, 'prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam='//int2str(nint(cfg%mskdiam))//' nstates='//int2str(ns)//&
                        &' nparts='//int2str(cfg%nprec_cons)//' nthr=8', 'rec.log', 16, ok)
                    if( .not. ok ) THROW_HARD('reconstruct3D FAILED (family '//int2str(fam)//' chunk '//int2str(q)//'): '//trim(first_error(trim(qd)//'/rec.log')))
                    call touch(trim(qd)//'/REC_DONE')
                enddo
                deallocate(inf, loc)
            enddo
            do s = 1,K
                best = maxloc(nf(:,s), 1)
                q    = (s - 1) / cfg%chunk_cons
                sl   = s - q * cfg%chunk_cons
                rr   = latest_dir(trim(C)//'/rec'//int2str(best)//'/q'//int2str(q))
                do h = 1,3
                    f = trim(rr)//'/recvol_state'//int2str_pad(sl,2)//trim(suf(h))//'.mrc'
                    if( file_exists(trim(f)) ) call copyf(f, trim(Lc)//'/maps/'//trim(hdir(h))//'recvol_state'//int2str_pad(s,2)//trim(suf(h))//'.mrc')
                enddo
                f = trim(rr)//'/recvol_state'//int2str_pad(sl,2)//'.mrc'
                if( file_exists(trim(f)) ) call copyf(f, trim(Lc)//'/states_dir/state_'//int2str_pad(s,2)//'.mrc')
                line = 'merged state '//int2str(s)//': family '//int2str(best)//' map ('//int2str(nf(1,s))
                do fam = 2,cfg%nfam
                    line = trim(line)//' + '//int2str(nf(fam,s))
                enddo
                call klog('consolidate', trim(line)//' particles)')
            enddo
            open(newunit=u, file=trim(C)//'/groups.txt', status='old', action='read')
            open(newunit=io, file=trim(Lc)//'/clusters.txt', status='replace', action='write')
            do
                read(u,'(A)',iostat=i) line
                if( i /= 0 ) exit
                if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                read(line,*) s, np_
                write(io,'(A)') int2str(s)//' 0 '//int2str(np_)//' '//int2str(np_)
            enddo
            close(u)
            close(io)
            call readout(cfg, trim(Lc)//'/states_dir', smpd, trim(Lc)//'/state_matches.txt')
            call touch(trim(C)//'/DONE')
            call klog('consolidate', 'DONE: '//int2str(K)//' merged states, maps in '//trim(Lc)//'/states_dir, groups '//trim(C)//'/groups.txt')
        endif
        call list_sorted(trim(Lc)//'/states_dir/state_*.mrc', sfiles, nsf)
        figs = trim(R)//'/figs'
        call simple_mkdir(trim(figs))
        write(logfhandle,'(A)') '== gallery'
        call tree_gallery(C, cfg%cons, trim(figs)//'/gallery_'//trim(cfg%rtag)//'_consolidated_'//fr(real(cfg%floor,dp),0)//'_'//trim(tag), 8, 7., smpd)
        if( trim(cfg%closed) /= 'none' )then
            write(logfhandle,'(A)') '== arm readout (closed/open occupancy per state)'
            call tree_arm_readout(cfg%cons, cfg%closed, sfiles, nsf, trim(C)//'/arm_states.txt', smpd)
        endif
        if( trim(cfg%docked_all) /= 'none' .and. nsf > 0 )then
            write(logfhandle,'(A)') '== strict recovery of the deposited conformers'
            call simple_mkdir(trim(C)//'/readout55')
            call tree_deposited_match(sfiles, nsf, cfg%docked_all, cfg%cons, smpd, 60., 20., trim(C)//'/readout55/matches.txt')
            allocate(names(nsf))
            do i = 1,nsf
                names(i) = sfiles(i)(index(sfiles(i), '/', back=.true.)+1:)
            enddo
            call tree_recover(trim(C)//'/readout55/matches.txt', names, nsf, 0.4, 0.05)
        endif
        cnt = 0
        call touch(trim(R)//'/pool_'//trim(cfg%tag)//'/DONE')
        write(logfhandle,'(A)') 'POOLED DONE '//trim(C)
    end subroutine stage_pool

    subroutine stage_parents( cfg, R, src )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: R, src
        type(sp_project) :: spproj
        integer, allocatable :: st(:)
        character(len=PLEN) :: O, D, ps, rr, f, part
        character(len=8)    :: suf(3), hdir(3)
        logical :: ok
        integer :: fam, i, n, h
        type(string), allocatable :: parts(:)
        suf  = ['_even   ', '_odd    ', '        ']
        hdir = ['halves/ ', 'halves/ ', '        ']
        O    = trim(R)//'/parents'
        call simple_mkdir(trim(O)//'/halves')
        do fam = 1,cfg%nfam
            if( file_exists(trim(O)//'/fam'//int2str(fam)//'.mrc') ) cycle
            D  = trim(O)//'/rec'//int2str(fam)
            ps = refined_proj(trim(R)//'/fam'//int2str(fam), src)
            call simple_mkdir(trim(D))
            call copyf(ps, trim(D)//'/p.simple')
            call spproj%read_segment('ptcl3D', string(trim(D)//'/p.simple'))
            n = spproj%os_ptcl3D%get_noris()
            allocate(st(n))
            do i = 1,n
                st(i) = merge(1, 0, spproj%os_ptcl3D%get_state(i) > 0)
            enddo
            call spproj%kill
            call write_ints(trim(D)//'/states.txt', st)
            call klog('parentfam', 'family '//int2str(fam)//': '//int2str(count(st == 1))//' particles with a state')
            deallocate(st)
            call select_states('p.simple', trim(D)//'/states.txt', D)
            if( trim(cfg%partition) /= 'none' )then
                part = ' qsys_partition='//trim(cfg%partition)
            else
                part = ''
            endif
            call run_prog(D, 'prg=update_project projfile=p.simple qsys_name=slurm'//trim(part)//' job_memory_per_task=32000 walltime=20000 nparts='//&
                &int2str(cfg%np_parent)//' nthr=8 mkdir=no', 'up.log', 1, ok)
            if( .not. ok ) THROW_HARD('update_project FAILED in '//trim(D))
            call klog('parentfam', 'reconstruct3D family '//int2str(fam)//' at its refined poses')
            call run_prog(D, 'prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam='//int2str(nint(cfg%mskdiam))//' nstates=1 nparts='//&
                &int2str(cfg%np_parent)//' nthr=8', 'rec.log', 8, ok)
            if( .not. ok ) THROW_HARD('reconstruct3D FAILED for the family '//int2str(fam)//' parent: '//trim(first_error(trim(D)//'/rec.log')))
            rr = latest_dir(D)
            do h = 1,3
                f = trim(rr)//'/recvol_state01'//trim(suf(h))//'.mrc'
                if( .not. file_exists(trim(f)) ) THROW_HARD('parent map missing: '//trim(f))
                call copyf(f, trim(O)//'/'//trim(hdir(h))//'fam'//int2str(fam)//trim(suf(h))//'.mrc')
            enddo
            call simple_list_files(trim(rr)//'/*_part*.mrc', parts)
            if( allocated(parts) )then
                do i = 1,size(parts)
                    call del_file(parts(i))
                enddo
            endif
        enddo
        call touch(trim(O)//'/DONE')
        call klog('parentfam', 'DONE')
    end subroutine stage_parents

end module simple_commanders_klem_tree
