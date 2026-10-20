!@descr: the kernel-latent EM tree in Fortran: what the bash/Python chain did around the board. A finished board
!! (klem_latent.txt, klem_states.txt, klem_resp.bin, klem_halves.bin) is read back; families by spectral bisection of
!! the cell-deviation graph; data-derived regions from the cross-half reproducible variability (spatial components and
!! motion-pattern groups) written as full-box soft masks; states = average-linkage clusters of surviving cells with the
!! top particles by responsibility mass; ranking along the node's own main deviation axis; identification of cells
!! against deposited maps (readout only); the landscape count over the leaves.
module simple_klem_tree
use simple_core_module_api
use iso_fortran_env, only: int8, int16, int32, int64
use simple_rnd,      only: seed_rnd_fixed
use simple_jpg,      only: write_rgb_jpeg
use simple_fftw3,    only: fftw_plan_dft_3d, fftw_execute_dft, fftw_destroy_plan, FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
!$ use omp_lib
implicit none
public :: board_t, leaf_t, read_board, kill_board, read_mrc, write_mrc, read_ints, write_ints, list_sorted, sort_names
public :: tree_init_positions, tree_parent_subset, tree_families, tree_regions, tree_node_states, tree_rank_pc1
public :: tree_kmeans_states, tree_deposited_match, tree_count_states, tree_merge_leaves, tree_consolidate
public :: tree_compose, tree_arm_readout, tree_recover, tree_gallery, leaf_ok, leaf_mask_path, build_leaves, fr, PLEN
private
#include "simple_local_flags.inc"

integer, parameter :: PLEN = 512

type board_t
    integer :: n = 0, ncl = 0, ngrid = 0, kmax = 0, n3 = 0, n1 = 0, n2 = 0, nb = 0
    integer,  allocatable :: pinds(:), am(:), i1(:), i2(:), bh(:), bk(:), bl(:)
    real,     allocatable :: r(:,:), crep(:), w(:)
    real(dp), allocatable :: z1(:), z2(:), gm(:), maxr(:)
    complex,  allocatable :: D0(:,:), D1(:,:)
end type board_t

type leaf_t
    integer :: fam = 0, state = 0, nrec = 1
    character(len=PLEN) :: name = '', dir = ''
end type leaf_t

type fft3_t
    integer     :: n = 0
    type(c_ptr) :: pf = c_null_ptr, pb = c_null_ptr
    complex(dp), allocatable :: a(:,:,:)
end type fft3_t

interface
    subroutine dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info )
        character        :: jobz, uplo
        integer          :: n, lda, lwork, liwork, info, iwork(*)
        double precision :: a(lda,*), w(*), work(*)
    end subroutine dsyevd
    subroutine cgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info )
        character :: jobz
        integer   :: m, n, lda, ldu, ldvt, lwork, info, iwork(*)
        complex   :: a(lda,*), u(ldu,*), vt(ldvt,*), work(*)
        real      :: s(*), rwork(*)
    end subroutine cgesdd
    subroutine cgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
        character :: transa, transb
        integer   :: m, n, k, lda, ldb, ldc
        complex   :: alpha, beta, a(lda,*), b(ldb,*), c(ldc,*)
    end subroutine cgemm
end interface

contains

    subroutine read_mrc( fname, v, n )
        character(len=*),  intent(in)    :: fname
        real, allocatable, intent(inout) :: v(:,:,:)
        integer,           intent(out)   :: n(3)
        integer :: h(256), u, io
        if( allocated(v) ) deallocate(v)
        open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='old', action='read', iostat=io)
        if( io /= 0 ) THROW_HARD('cannot open '//trim(fname))
        read(u) h
        n = h(1:3)
        allocate(v(n(1),n(2),n(3)))
        read(u, pos=1025+h(24)) v
        close(u)
    end subroutine read_mrc

    subroutine write_mrc( fname, v, smpd )
        character(len=*), intent(in) :: fname
        real(dp),         intent(in) :: v(:,:,:)
        real,             intent(in) :: smpd
        real, allocatable :: f(:,:,:)
        integer  :: h(256), u, n(3)
        real(dp) :: mu, sd
        n = shape(v)
        allocate(f(n(1),n(2),n(3)))
        f  = real(v)
        mu = sum(real(f,dp)) / real(size(f),dp)
        sd = sqrt(sum((real(f,dp) - mu)**2) / real(size(f),dp))
        h        = 0
        h(1:3)   = n
        h(4)     = 2
        h(8:10)  = n
        h(11:13) = transfer(real(n) * smpd, 0, 3)
        h(14:16) = transfer([90.,90.,90.], 0, 3)
        h(17:19) = [1,2,3]
        h(20)    = transfer(minval(f), 0)
        h(21)    = transfer(maxval(f), 0)
        h(22)    = transfer(real(mu), 0)
        h(53)    = transfer('MAP ', 0)
        h(54)    = 16708
        h(55)    = transfer(real(sd), 0)
        open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='replace', action='write')
        write(u) h, f
        close(u)
    end subroutine write_mrc

    subroutine read_ints( fname, v )
        character(len=*),     intent(in)    :: fname
        integer, allocatable, intent(inout) :: v(:)
        character(len=4096) :: line
        integer :: u, io, n, i
        if( allocated(v) ) deallocate(v)
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) THROW_HARD('cannot open '//trim(fname))
        n = 0
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            n = n + 1
        enddo
        allocate(v(n))
        rewind(u)
        i = 0
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1
            read(line,*) v(i)
        enddo
        close(u)
    end subroutine read_ints

    subroutine write_ints( fname, v )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: v(:)
        integer :: u, i
        open(newunit=u, file=trim(fname), status='replace', action='write')
        do i = 1,size(v)
            write(u,'(I0)') v(i)
        enddo
        close(u)
    end subroutine write_ints

    subroutine list_sorted( pattern, names, n )
        character(len=*),                 intent(in)    :: pattern
        character(len=PLEN), allocatable, intent(inout) :: names(:)
        integer,                          intent(out)   :: n
        type(string), allocatable :: l(:)
        integer :: i
        call simple_list_files(pattern, l)
        n = 0
        if( allocated(l) ) n = size(l)
        if( allocated(names) ) deallocate(names)
        allocate(names(n))
        do i = 1,n
            names(i) = l(i)%to_char()
        enddo
        call sort_names(names, n)
    end subroutine list_sorted

    subroutine sort_names( names, n )
        integer,             intent(in)    :: n
        character(len=PLEN), intent(inout) :: names(n)
        character(len=PLEN) :: t
        integer :: i, j
        do i = 2,n
            t = names(i)
            j = i - 1
            do while( j >= 1 )
                if( names(j) <= t ) exit
                names(j+1) = names(j)
                j = j - 1
            enddo
            names(j+1) = t
        enddo
    end subroutine sort_names

    function fr( x, d ) result( s )
        real(dp), intent(in) :: x
        integer,  intent(in) :: d
        character(len=:), allocatable :: s
        character(len=64) :: b
        character(len=16) :: fm
        integer :: l
        if( x /= x )then
            s = 'nan'
            return
        endif
        write(fm,'(A,I0,A)') '(F0.', d, ')'
        write(b, fm) x
        b = adjustl(b)
        if( b(1:1) == '.' ) b = '0'//trim(b)
        if( b(1:2) == '-.' ) b = '-0'//trim(b(2:))
        l = len_trim(b)
        if( d == 0 .and. b(l:l) == '.' ) b(l:l) = ' '
        if( trim(b) == '-0' ) b = '0'
        s = trim(b)
    end function fr

    function frs( x, d ) result( s )
        real(dp), intent(in) :: x
        integer,  intent(in) :: d
        character(len=:), allocatable :: s
        s = fr(x, d)
        if( s(1:1) /= '-' ) s = '+'//s
    end function frs

    function leaf_mask_path( run, name ) result( p )
        character(len=*), intent(in) :: run, name
        character(len=PLEN) :: p, rest
        integer :: q
        logical :: mot
        rest = name(4:)
        mot  = index(rest, 'motion_') == 1
        if( mot ) rest = rest(8:)
        if( rest(1:1) == '_' ) rest = rest(2:)
        q = 1
        do while( q <= len_trim(rest) )
            if( rest(q:q) < '0' .or. rest(q:q) > '9' ) exit
            q = q + 1
        enddo
        if( mot )then
            p = trim(run)//'/fam'//name(2:2)//'/regions_motion/region_'//rest(1:q-1)//'.mrc'
        else
            p = trim(run)//'/fam'//name(2:2)//'/regions/region_'//rest(1:q-1)//'.mrc'
        endif
    end function leaf_mask_path

    subroutine npz_locate( fname, key, prefix, found, dpos )
        character(len=*),    intent(in)  :: fname, key
        logical,             intent(in)  :: prefix
        character(len=PLEN), intent(out) :: found
        integer(int64),      intent(out) :: dpos
        character(len=:), allocatable :: nm
        integer(int64) :: pos, csize, usize, q, qend
        integer(int32) :: sig, cs32, us32
        integer(int16) :: method, nlen16, xlen16, hid, hsz
        integer        :: u, io, nlen, xlen
        logical        :: hit
        found = ''
        dpos  = -1
        open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='old', action='read', iostat=io)
        if( io /= 0 ) THROW_HARD('cannot open '//trim(fname))
        pos = 1
        do
            read(u, pos=pos, iostat=io) sig
            if( io /= 0 ) exit
            if( sig /= 67324752_int32 ) exit
            read(u, pos=pos+8)  method
            read(u, pos=pos+18) cs32, us32
            read(u, pos=pos+26) nlen16, xlen16
            nlen = iand(int(nlen16), 65535)
            xlen = iand(int(xlen16), 65535)
            if( allocated(nm) ) deallocate(nm)
            allocate(character(len=nlen) :: nm)
            read(u, pos=pos+30) nm
            csize = iand(int(cs32, int64), 4294967295_int64)
            usize = iand(int(us32, int64), 4294967295_int64)
            if( cs32 == -1_int32 .or. us32 == -1_int32 )then
                q    = pos + 30 + nlen
                qend = q + xlen
                do while( q + 4 <= qend )
                    read(u, pos=q) hid, hsz
                    if( hid == 1_int16 )then
                        q = q + 4
                        if( us32 == -1_int32 )then
                            read(u, pos=q) usize
                            q = q + 8
                        endif
                        if( cs32 == -1_int32 ) read(u, pos=q) csize
                        exit
                    endif
                    q = q + 4 + iand(int(hsz, int64), 65535_int64)
                enddo
            endif
            if( method /= 0_int16 ) THROW_HARD('compressed npz entries are not supported: '//trim(fname))
            if( prefix )then
                hit = index(nm, trim(key)) == 1
            else
                hit = nm == trim(key)//'.npy'
            endif
            if( hit )then
                found = nm
                dpos  = pos + 30 + nlen + xlen
                exit
            endif
            pos = pos + 30 + nlen + xlen + csize
        enddo
        close(u)
    end subroutine npz_locate

    subroutine npz_read( fname, key, a, nrow, ncol, prefix, keyname )
        character(len=*),              intent(in)    :: fname, key
        real(dp), allocatable,         intent(inout) :: a(:,:)
        integer,                       intent(out)   :: nrow, ncol
        logical,             optional, intent(in)    :: prefix
        character(len=PLEN), optional, intent(out)   :: keyname
        character(len=PLEN) :: found
        character(len=:), allocatable :: hdr, inner
        character(len=6)  :: magic
        character(len=16) :: descr
        integer(int64)    :: dpos, hpos
        integer(int8)     :: major, minor
        integer(int16)    :: hl16
        integer(int32)    :: hl32
        real(dp),       allocatable :: bd(:)
        real,           allocatable :: bf(:)
        integer(int64), allocatable :: bl(:)
        integer(int32), allocatable :: bi(:)
        integer(int8),  allocatable :: bb(:)
        integer :: u, hlen, p1, p2, nel, sh(2), nsh, i, io
        logical :: lpref, intok
        lpref = .false.
        if( present(prefix) ) lpref = prefix
        call npz_locate(fname, key, lpref, found, dpos)
        if( dpos < 0 ) THROW_HARD('no entry '//trim(key)//' in '//trim(fname))
        if( present(keyname) ) keyname = found(1:len_trim(found)-4)
        open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='old', action='read')
        read(u, pos=dpos) magic, major, minor
        if( major == 1_int8 )then
            read(u, pos=dpos+8) hl16
            hlen = iand(int(hl16), 65535)
            hpos = dpos + 10
        else
            read(u, pos=dpos+8) hl32
            hlen = int(hl32)
            hpos = dpos + 12
        endif
        allocate(character(len=hlen) :: hdr)
        read(u, pos=hpos) hdr
        if( index(hdr, "'fortran_order': True") > 0 ) THROW_HARD('Fortran-ordered npy arrays are not supported: '//trim(found))
        p1    = index(hdr, "'descr': '") + 10
        p2    = p1 + index(hdr(p1:), "'") - 2
        descr = hdr(p1:p2)
        p1    = index(hdr, "'shape': (") + 10
        p2    = p1 + index(hdr(p1:), ")") - 2
        inner = hdr(p1:p2)//' '
        do i = 1,len(inner)
            if( inner(i:i) == ',' ) inner(i:i) = ' '
        enddo
        nsh   = 0
        sh    = 1
        intok = .false.
        do i = 1,len(inner)
            if( inner(i:i) /= ' ' .and. .not. intok )then
                nsh   = nsh + 1
                intok = .true.
                read(inner(i:),*,iostat=io) sh(min(nsh,2))
            else if( inner(i:i) == ' ' )then
                intok = .false.
            endif
        enddo
        nrow = sh(1)
        ncol = 1
        if( nsh >= 2 ) ncol = sh(2)
        nel = nrow * ncol
        if( allocated(a) ) deallocate(a)
        allocate(a(ncol,nrow))
        select case(trim(descr))
            case('<f8')
                allocate(bd(nel))
                read(u, pos=hpos+hlen) bd
                a = reshape(bd, [ncol,nrow])
            case('<f4')
                allocate(bf(nel))
                read(u, pos=hpos+hlen) bf
                a = reshape(real(bf,dp), [ncol,nrow])
            case('<i8')
                allocate(bl(nel))
                read(u, pos=hpos+hlen) bl
                a = reshape(real(bl,dp), [ncol,nrow])
            case('<i4')
                allocate(bi(nel))
                read(u, pos=hpos+hlen) bi
                a = reshape(real(bi,dp), [ncol,nrow])
            case('|b1','|u1','|i1')
                allocate(bb(nel))
                read(u, pos=hpos+hlen) bb
                a = reshape(real(bb,dp), [ncol,nrow])
            case default
                THROW_HARD('unsupported npy dtype '//trim(descr)//' in '//trim(found))
        end select
        close(u)
    end subroutine npz_read

    subroutine kt_sort( v )
        real(dp), intent(inout) :: v(:)
        integer  :: n, i, ir, j, l
        real(dp) :: rra
        n = size(v)
        if( n < 2 ) return
        l  = n/2 + 1
        ir = n
        do
            if( l > 1 )then
                l   = l - 1
                rra = v(l)
            else
                rra   = v(ir)
                v(ir) = v(1)
                ir    = ir - 1
                if( ir == 1 )then
                    v(1) = rra
                    return
                endif
            endif
            i = l
            j = l + l
            do while( j <= ir )
                if( j < ir )then
                    if( v(j) < v(j+1) ) j = j + 1
                endif
                if( rra < v(j) )then
                    v(i) = v(j)
                    i    = j
                    j    = j + j
                else
                    j = ir + 1
                endif
            enddo
            v(i) = rra
        enddo
    end subroutine kt_sort

    function kt_percentile( v, q ) result( p )
        real(dp), intent(in) :: v(:)
        real(dp), intent(in) :: q
        real(dp) :: p, idx, t, a, b, d
        real(dp), allocatable :: s(:)
        integer :: n, lo
        n = size(v)
        if( n == 0 )then
            p = 0.d0
            return
        endif
        allocate(s(n), source=v)
        call kt_sort(s)
        idx = real(n - 1, dp) * (q / 100.d0)
        lo  = min(max(int(floor(idx)), 0), n - 1)
        t   = idx - real(lo, dp)
        a   = s(lo + 1)
        b   = s(min(lo + 2, n))
        d   = b - a
        if( t >= 0.5d0 )then
            p = b - d * (1.d0 - t)
        else
            p = a + d * t
        endif
    end function kt_percentile

    function kt_median( v ) result( m )
        real(dp), intent(in) :: v(:)
        real(dp) :: m
        real(dp), allocatable :: s(:)
        integer :: n
        n = size(v)
        if( n == 0 )then
            m = 0.d0
            return
        endif
        allocate(s(n), source=v)
        call kt_sort(s)
        if( mod(n,2) == 1 )then
            m = s((n+1)/2)
        else
            m = (s(n/2) + s(n/2+1)) / 2.d0
        endif
    end function kt_median

    subroutine kt_argsort_desc( v, idx )
        real(dp), intent(in)  :: v(:)
        integer,  intent(out) :: idx(:)
        integer, allocatable :: tmp(:)
        integer :: n, w, i, lo, mid, hi, a, b, k
        n = size(v)
        do i = 1,n
            idx(i) = i
        enddo
        allocate(tmp(n))
        w = 1
        do while( w < n )
            lo = 1
            do while( lo <= n )
                mid = min(lo + w - 1, n)
                hi  = min(lo + 2*w - 1, n)
                if( mid < hi )then
                    a = lo
                    b = mid + 1
                    k = lo
                    do while( a <= mid .and. b <= hi )
                        if( v(idx(b)) > v(idx(a)) )then
                            tmp(k) = idx(b)
                            b = b + 1
                        else
                            tmp(k) = idx(a)
                            a = a + 1
                        endif
                        k = k + 1
                    enddo
                    do while( a <= mid )
                        tmp(k) = idx(a)
                        a = a + 1
                        k = k + 1
                    enddo
                    do while( b <= hi )
                        tmp(k) = idx(b)
                        b = b + 1
                        k = k + 1
                    enddo
                    idx(lo:hi) = tmp(lo:hi)
                endif
                lo = lo + 2*w
            enddo
            w = 2*w
        enddo
    end subroutine kt_argsort_desc

    function kt_pairwise_sum( a ) result( s )
        real(dp), intent(in) :: a(:)
        real(dp) :: s, r(8)
        integer  :: n, i, j
        n = size(a)
        if( n < 8 )then
            s = 0.d0
            do i = 1,n
                s = s + a(i)
            enddo
            return
        endif
        r = a(1:8)
        i = 9
        do while( i <= n - mod(n,8) )
            do j = 1,8
                r(j) = r(j) + a(i+j-1)
            enddo
            i = i + 8
        enddo
        s = ((r(1) + r(2)) + (r(3) + r(4))) + ((r(5) + r(6)) + (r(7) + r(8)))
        do while( i <= n )
            s = s + a(i)
            i = i + 1
        enddo
    end function kt_pairwise_sum

    subroutine fft3_new( F, n )
        type(fft3_t), intent(inout) :: F
        integer,      intent(in)    :: n
        call fft3_kill(F)
        F%n = n
        allocate(F%a(n,n,n))
        !$omp critical (klem_fftw_planner)
        F%pf = fftw_plan_dft_3d(int(n,c_int), int(n,c_int), int(n,c_int), F%a, F%a, FFTW_FORWARD,  FFTW_ESTIMATE)
        F%pb = fftw_plan_dft_3d(int(n,c_int), int(n,c_int), int(n,c_int), F%a, F%a, FFTW_BACKWARD, FFTW_ESTIMATE)
        !$omp end critical (klem_fftw_planner)
    end subroutine fft3_new

    subroutine fft3_fwd( F )
        type(fft3_t), intent(inout) :: F
        call fftw_execute_dft(F%pf, F%a, F%a)
    end subroutine fft3_fwd

    subroutine fft3_bwd( F )
        type(fft3_t), intent(inout) :: F
        call fftw_execute_dft(F%pb, F%a, F%a)
    end subroutine fft3_bwd

    subroutine fft3_kill( F )
        type(fft3_t), intent(inout) :: F
        !$omp critical (klem_fftw_planner)
        if( c_associated(F%pf) ) call fftw_destroy_plan(F%pf)
        if( c_associated(F%pb) ) call fftw_destroy_plan(F%pb)
        !$omp end critical (klem_fftw_planner)
        F%pf = c_null_ptr
        F%pb = c_null_ptr
        if( allocated(F%a) ) deallocate(F%a)
        F%n = 0
    end subroutine fft3_kill

    pure function kt_freq( i, n ) result( f )
        integer, intent(in) :: i, n
        integer :: f
        f = i - 1
        if( f >= n - n/2 ) f = f - n
    end function kt_freq

    subroutine kt_shells( n, sh )
        integer,              intent(in)    :: n
        integer, allocatable, intent(inout) :: sh(:,:,:)
        integer :: i, j, k, fi, fj, fk
        if( allocated(sh) ) deallocate(sh)
        allocate(sh(n,n,n))
        do k = 1,n
            fk = kt_freq(k, n)
            do j = 1,n
                fj = kt_freq(j, n)
                do i = 1,n
                    fi = kt_freq(i, n)
                    sh(i,j,k) = min(nint(sqrt(real(fi*fi + fj*fj + fk*fk, dp))), n/2)
                enddo
            enddo
        enddo
    end subroutine kt_shells

    subroutine kt_lowpass( v, F, smpd, mode, res )
        real(dp),     intent(inout) :: v(:,:,:)
        type(fft3_t), intent(inout) :: F
        real,         intent(in)    :: smpd, res
        integer,      intent(in)    :: mode
        integer  :: n, i, j, k
        real(dp) :: k2, g, fct, sm, rs
        n   = size(v,1)
        fct = 1.d0 / real(n,dp)**3
        sm  = real(smpd,dp)
        rs  = real(res,dp)
        F%a = cmplx(v, 0.d0, dp)
        call fft3_fwd(F)
        do k = 1,n
            do j = 1,n
                do i = 1,n
                    k2 = (real(kt_freq(i,n),dp)**2 + real(kt_freq(j,n),dp)**2 + real(kt_freq(k,n),dp)**2) / real(n,dp)**2
                    select case(mode)
                        case(1)
                            g = exp(-3.d0 * k2 / (sm / rs)**2)
                        case(2)
                            g = exp(-0.5d0 * (k2 / sm**2) * (rs / 1.2d0)**2)
                        case default
                            g = exp(-0.5d0 * k2 * (sm * rs / 1.2d0)**2 / sm**2)
                    end select
                    F%a(i,j,k) = F%a(i,j,k) * g
                enddo
            enddo
        enddo
        call fft3_bwd(F)
        v = real(F%a, dp) * fct
    end subroutine kt_lowpass

    function kt_fsc_res( e, o, m, sh, smpd ) result( res )
        real(dp), intent(in) :: e(:,:,:), o(:,:,:), m(:,:,:)
        integer,  intent(in) :: sh(:,:,:)
        real,     intent(in) :: smpd
        real(dp) :: res
        type(fft3_t) :: FE, FO
        real(dp), allocatable :: num(:), d1(:), d2(:), fsc(:)
        integer :: n, i, j, k, s, kr
        n = size(e,1)
        call fft3_new(FE, n)
        call fft3_new(FO, n)
        FE%a = cmplx(e * m, 0.d0, dp)
        FO%a = cmplx(o * m, 0.d0, dp)
        call fft3_fwd(FE)
        call fft3_fwd(FO)
        allocate(num(0:n/2), d1(0:n/2), d2(0:n/2), fsc(0:n/2), source=0.d0)
        do k = 1,n
            do j = 1,n
                do i = 1,n
                    s = sh(i,j,k)
                    num(s) = num(s) + real(conjg(FE%a(i,j,k)) * FO%a(i,j,k), dp)
                    d1(s)  = d1(s)  + abs(FE%a(i,j,k))**2
                    d2(s)  = d2(s)  + abs(FO%a(i,j,k))**2
                enddo
            enddo
        enddo
        fsc = num / max(sqrt(d1 * d2), 1.d-30)
        kr  = n/2 - 1
        do s = 1,n/2-1
            if( fsc(s) < 0.143d0 )then
                kr = s
                exit
            endif
        enddo
        res = real(n,dp) * real(smpd,dp) / real(max(kr,1),dp)
        call fft3_kill(FE)
        call fft3_kill(FO)
    end function kt_fsc_res

    function kt_sphere( n ) result( m )
        integer, intent(in) :: n
        real(dp), allocatable :: m(:,:,:)
        integer  :: i, j, k
        real(dp) :: rr
        allocate(m(n,n,n))
        do k = 1,n
            do j = 1,n
                do i = 1,n
                    rr = sqrt(real((i-1-n/2)**2 + (j-1-n/2)**2 + (k-1-n/2)**2, dp))
                    m(i,j,k) = min(max((0.45d0 * real(n,dp) - rr) / (0.05d0 * real(n,dp)), 0.d0), 1.d0)
                enddo
            enddo
        enddo
    end function kt_sphere

    subroutine read_board( dir, B, with_halves )
        character(len=*),  intent(in)    :: dir
        type(board_t),     intent(inout) :: B
        logical, optional, intent(in)    :: with_halves
        character(len=1024)  :: line
        complex, allocatable :: cube(:,:,:)
        integer  :: u, io, i, ic, ih, isg, n1, n2, m, ib, h, k, l, kx, nc2
        real(dp) :: a1, a2, a5, a6, a7
        real     :: pop, rep
        logical  :: lh
        lh = .true.
        if( present(with_halves) ) lh = with_halves
        call kill_board(B)
        open(newunit=u, file=trim(dir)//'/klem_resp.bin', access='stream', form='unformatted', status='old', action='read', iostat=io)
        if( io /= 0 ) THROW_HARD('no board in '//trim(dir))
        read(u) B%n, B%ncl, B%ngrid, n1, n2
        B%n1 = n1
        B%n2 = n2
        allocate(B%pinds(B%n), B%i1(B%ngrid), B%i2(B%ngrid), B%r(B%ncl,B%n))
        read(u) B%pinds
        read(u) B%i1
        read(u) B%i2
        read(u) B%r
        close(u)
        allocate(B%am(B%n), B%z1(B%n), B%z2(B%n), B%gm(B%n), B%maxr(B%n))
        i = 0
        open(newunit=u, file=trim(dir)//'/klem_latent.txt', status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1
            if( i > B%n ) exit
            read(line,*) m, a1, a2, ic, a5, a6, a7
            if( m /= B%pinds(i) ) THROW_HARD('klem_latent.txt and klem_resp.bin disagree on the particle order in '//trim(dir))
            B%z1(i)   = a1
            B%z2(i)   = a2
            B%am(i)   = ic
            B%maxr(i) = a5
            B%gm(i)   = a6
        enddo
        close(u)
        if( i /= B%n ) THROW_HARD('klem_latent.txt row count differs from klem_resp.bin in '//trim(dir))
        allocate(B%crep(B%ncl), source=0.)
        open(newunit=u, file=trim(dir)//'/klem_states.txt', status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            read(line,*) ic, pop, rep, isg
            if( ic >= 1 .and. ic <= B%ncl ) B%crep(ic) = rep
        enddo
        close(u)
        if( .not. lh ) return
        open(newunit=u, file=trim(dir)//'/klem_halves.bin', access='stream', form='unformatted', status='old', action='read')
        read(u) ic, m, B%kmax, B%n3
        if( ic /= B%ncl .or. m /= B%ngrid ) THROW_HARD('klem_halves.bin header disagrees with klem_resp.bin in '//trim(dir))
        kx  = B%kmax
        nc2 = (B%n3/2)**2
        allocate(B%w(B%ngrid))
        read(u) B%w
        B%nb = 0
        do l = -kx,kx
            do k = -kx,kx
                do h = -kx,kx
                    if( h*h + k*k + l*l >= 1 .and. h*h + k*k + l*l <= nc2 ) B%nb = B%nb + 1
                enddo
            enddo
        enddo
        allocate(B%bh(B%nb), B%bk(B%nb), B%bl(B%nb))
        ib = 0
        do l = -kx,kx
            do k = -kx,kx
                do h = -kx,kx
                    if( h*h + k*k + l*l >= 1 .and. h*h + k*k + l*l <= nc2 )then
                        ib = ib + 1
                        B%bh(ib) = h
                        B%bk(ib) = k
                        B%bl(ib) = l
                    endif
                enddo
            enddo
        enddo
        allocate(B%D0(B%nb,B%ngrid), B%D1(B%nb,B%ngrid), cube(-kx:kx,-kx:kx,-kx:kx))
        do ih = 1,2
            do ic = 1,B%ngrid
                read(u) cube
                do ib = 1,B%nb
                    if( ih == 1 )then
                        B%D0(ib,ic) = cube(B%bh(ib),B%bk(ib),B%bl(ib))
                    else
                        B%D1(ib,ic) = cube(B%bh(ib),B%bk(ib),B%bl(ib))
                    endif
                enddo
            enddo
        enddo
        close(u)
    end subroutine read_board

    subroutine kill_board( B )
        type(board_t), intent(inout) :: B
        if( allocated(B%pinds) ) deallocate(B%pinds, B%am, B%i1, B%i2, B%r, B%z1, B%z2, B%gm, B%maxr, B%crep)
        if( allocated(B%w)     ) deallocate(B%w)
        if( allocated(B%D0)    ) deallocate(B%D0, B%D1, B%bh, B%bk, B%bl)
        B%n     = 0
        B%ncl   = 0
        B%ngrid = 0
        B%kmax  = 0
        B%n3    = 0
        B%nb    = 0
    end subroutine kill_board

    subroutine band_to_real( B, d, F, v )
        type(board_t), intent(in)    :: B
        complex,       intent(in)    :: d(:)
        type(fft3_t),  intent(inout) :: F
        real(dp),      intent(out)   :: v(:,:,:)
        integer  :: ib, n, i, j, k
        real(dp) :: fct
        n   = B%n3
        fct = 1.d0 / real(n,dp)**3
        F%a = cmplx(0.d0, 0.d0, dp)
        do ib = 1,B%nb
            F%a(modulo(B%bh(ib),n)+1, modulo(B%bk(ib),n)+1, modulo(B%bl(ib),n)+1) = cmplx(d(ib), kind=dp)
        enddo
        call fft3_bwd(F)
        do k = 1,n
            do j = 1,n
                do i = 1,n
                    v(i,j,k) = real(F%a(modulo(i-1-B%kmax,n)+1, modulo(j-1-B%kmax,n)+1, modulo(k-1-B%kmax,n)+1), dp) * fct
                enddo
            enddo
        enddo
    end subroutine band_to_real

    subroutine kt_gauss( v, sigma )
        real(dp), intent(inout) :: v(:,:,:)
        real,     intent(in)    :: sigma
        real(dp), allocatable :: w(:), ext(:), res(:)
        real(dp) :: s2
        integer  :: r, n(3), ax, i, j, q, m
        r  = int(4.d0 * real(sigma,dp) + 0.5d0)
        s2 = real(sigma,dp) * real(sigma,dp)
        allocate(w(-r:r))
        do q = -r,r
            w(q) = exp(-0.5d0 / s2 * real(q*q,dp))
        enddo
        w = w / kt_pairwise_sum(w)
        n = shape(v)
        do ax = 3,1,-1
            m = n(ax)
            allocate(ext(1-r:m+r), res(m))
            select case(ax)
                case(3)
                    do j = 1,n(2)
                        do i = 1,n(1)
                            ext(1:m) = v(i,j,:)
                            call filt
                            v(i,j,:) = res
                        enddo
                    enddo
                case(2)
                    do j = 1,n(3)
                        do i = 1,n(1)
                            ext(1:m) = v(i,:,j)
                            call filt
                            v(i,:,j) = res
                        enddo
                    enddo
                case(1)
                    do j = 1,n(3)
                        do i = 1,n(2)
                            ext(1:m) = v(:,i,j)
                            call filt
                            v(:,i,j) = res
                        enddo
                    enddo
            end select
            deallocate(ext, res)
        enddo

      contains

        subroutine filt
            integer :: p, jj
            do p = 1,r
                ext(1-p) = ext(p)
                ext(m+p) = ext(m+1-p)
            enddo
            do p = 1,m
                res(p) = ext(p) * w(0)
                do jj = -r,-1
                    res(p) = res(p) + (ext(p+jj) + ext(p-jj)) * w(jj)
                enddo
            enddo
        end subroutine filt

    end subroutine kt_gauss

    subroutine kt_erode( m )
        logical, intent(inout) :: m(:,:,:)
        logical, allocatable :: t(:,:,:)
        integer :: n(3), i, j, k
        n = shape(m)
        allocate(t(n(1),n(2),n(3)))
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    t(i,j,k) = m(i,j,k)
                    if( .not. t(i,j,k) ) cycle
                    if( i == 1 .or. i == n(1) .or. j == 1 .or. j == n(2) .or. k == 1 .or. k == n(3) )then
                        t(i,j,k) = .false.
                    else
                        t(i,j,k) = m(i-1,j,k) .and. m(i+1,j,k) .and. m(i,j-1,k) .and. m(i,j+1,k) .and. m(i,j,k-1) .and. m(i,j,k+1)
                    endif
                enddo
            enddo
        enddo
        m = t
    end subroutine kt_erode

    subroutine kt_dilate( m, iters )
        logical, intent(inout) :: m(:,:,:)
        integer, intent(in)    :: iters
        logical, allocatable :: t(:,:,:)
        integer :: n(3), i, j, k, it
        n = shape(m)
        allocate(t(n(1),n(2),n(3)))
        do it = 1,iters
            do k = 1,n(3)
                do j = 1,n(2)
                    do i = 1,n(1)
                        t(i,j,k) = m(i,j,k)
                        if( t(i,j,k) ) cycle
                        if( i > 1    ) t(i,j,k) = t(i,j,k) .or. m(i-1,j,k)
                        if( i < n(1) ) t(i,j,k) = t(i,j,k) .or. m(i+1,j,k)
                        if( j > 1    ) t(i,j,k) = t(i,j,k) .or. m(i,j-1,k)
                        if( j < n(2) ) t(i,j,k) = t(i,j,k) .or. m(i,j+1,k)
                        if( k > 1    ) t(i,j,k) = t(i,j,k) .or. m(i,j,k-1)
                        if( k < n(3) ) t(i,j,k) = t(i,j,k) .or. m(i,j,k+1)
                    enddo
                enddo
            enddo
            m = t
        enddo
    end subroutine kt_dilate

    subroutine kt_label6( m, lab, nl )
        logical, intent(in)    :: m(:,:,:)
        integer, intent(inout) :: lab(:,:,:)
        integer, intent(out)   :: nl
        integer, allocatable :: st(:,:)
        integer :: n(3), i, j, k, top, x, y, z, q, nx, ny, nz, dx(6), dy(6), dz(6)
        dx  = [1,-1,0,0,0,0]
        dy  = [0,0,1,-1,0,0]
        dz  = [0,0,0,0,1,-1]
        n   = shape(m)
        lab = 0
        nl  = 0
        allocate(st(3,product(n)))
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    if( .not. m(i,j,k) .or. lab(i,j,k) /= 0 ) cycle
                    nl  = nl + 1
                    top = 1
                    st(:,1)    = [i,j,k]
                    lab(i,j,k) = nl
                    do while( top > 0 )
                        x   = st(1,top)
                        y   = st(2,top)
                        z   = st(3,top)
                        top = top - 1
                        do q = 1,6
                            nx = x + dx(q)
                            ny = y + dy(q)
                            nz = z + dz(q)
                            if( nx < 1 .or. nx > n(1) .or. ny < 1 .or. ny > n(2) .or. nz < 1 .or. nz > n(3) ) cycle
                            if( m(nx,ny,nz) .and. lab(nx,ny,nz) == 0 )then
                                lab(nx,ny,nz) = nl
                                top = top + 1
                                st(:,top) = [nx,ny,nz]
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine kt_label6

    subroutine kt_zoom( src, m, dst )
        real(dp), intent(in)  :: src(:,:,:)
        integer,  intent(in)  :: m
        real(dp), intent(out) :: dst(:,:,:)
        real(dp), allocatable :: t(:)
        integer,  allocatable :: s0(:), s1(:)
        integer  :: n, i, j, k, a, b, c
        real(dp) :: f, cc, wz, wy, wx, acc
        n = size(src,1)
        allocate(t(m), s0(m), s1(m))
        f = real(n - 1, dp) / real(m - 1, dp)
        do i = 1,m
            cc    = real(i - 1, dp) * f
            s0(i) = int(floor(cc))
            t(i)  = cc - real(s0(i), dp)
            s1(i) = min(s0(i) + 1, n - 1)
        enddo
        do k = 1,m
            do j = 1,m
                do i = 1,m
                    acc = 0.d0
                    do c = 0,1
                        wz = merge(1.d0 - t(k), t(k), c == 0)
                        do b = 0,1
                            wy = merge(1.d0 - t(j), t(j), b == 0)
                            do a = 0,1
                                wx  = merge(1.d0 - t(i), t(i), a == 0)
                                acc = acc + src(merge(s0(i), s1(i), a == 0) + 1, merge(s0(j), s1(j), b == 0) + 1, &
                                    &merge(s0(k), s1(k), c == 0) + 1) * (wz * wy * wx)
                            enddo
                        enddo
                    enddo
                    dst(i,j,k) = acc
                enddo
            enddo
        enddo
    end subroutine kt_zoom

    subroutine kt_avg_linkage( S, t, lab, K )
        real(dp), intent(in)  :: S(:,:)
        real(dp), intent(in)  :: t
        integer,  intent(out) :: lab(:), K
        real(dp), allocatable :: D(:,:)
        integer,  allocatable :: sz(:), newl(:)
        logical,  allocatable :: act(:)
        integer  :: n, i, j, a, b, c
        real(dp) :: dmin, dnew
        n = size(S,1)
        allocate(D(n,n), sz(n), act(n), newl(n))
        do j = 1,n
            do i = 1,n
                if( i < j )then
                    D(i,j) = min(max(1.d0 - S(i,j), 0.d0), 2.d0)
                else if( i > j )then
                    D(i,j) = min(max(1.d0 - S(j,i), 0.d0), 2.d0)
                else
                    D(i,j) = 0.d0
                endif
            enddo
        enddo
        sz  = 1
        act = .true.
        do i = 1,n
            lab(i) = i
        enddo
        do
            dmin = huge(dmin)
            a    = 0
            b    = 0
            do i = 1,n
                if( .not. act(i) ) cycle
                do j = i+1,n
                    if( .not. act(j) ) cycle
                    if( D(i,j) < dmin )then
                        dmin = D(i,j)
                        a    = i
                        b    = j
                    endif
                enddo
            enddo
            if( a == 0 ) exit
            if( dmin > t ) exit
            do c = 1,n
                if( .not. act(c) .or. c == a .or. c == b ) cycle
                dnew   = (real(sz(a),dp) * D(a,c) + real(sz(b),dp) * D(b,c)) / real(sz(a) + sz(b), dp)
                D(a,c) = dnew
                D(c,a) = dnew
            enddo
            act(b) = .false.
            sz(a)  = sz(a) + sz(b)
            do i = 1,n
                if( lab(i) == b ) lab(i) = a
            enddo
        enddo
        newl = 0
        K    = 0
        do i = 1,n
            if( newl(lab(i)) == 0 )then
                K = K + 1
                newl(lab(i)) = K
            endif
        enddo
        do i = 1,n
            lab(i) = newl(lab(i))
        enddo
    end subroutine kt_avg_linkage

    subroutine kt_kmeans_pp( X, k, nrestart, niter, lab, cent )
        real(dp), intent(in)  :: X(:,:)
        integer,  intent(in)  :: k, nrestart, niter
        integer,  intent(out) :: lab(:)
        real(dp), intent(out) :: cent(:,:)
        real(dp), allocatable :: C(:,:), Cn(:,:), d2(:), cum(:)
        integer,  allocatable :: l(:), cnt(:)
        real(dp) :: best, inert, dd, dmin, rr, tot
        integer  :: n, s, i, j, it, jb
        n = size(X,2)
        allocate(C(size(X,1),k), Cn(size(X,1),k), d2(n), cum(n), l(n), cnt(k))
        best = huge(best)
        do s = 0,nrestart-1
            call seed_rnd_fixed(s)
            i = min(n, 1 + int(real(ran3(),dp) * real(n,dp)))
            C(:,1) = X(:,i)
            do j = 2,k
                do i = 1,n
                    dmin = huge(dmin)
                    do jb = 1,j-1
                        dd   = sum((X(:,i) - C(:,jb))**2)
                        dmin = min(dmin, dd)
                    enddo
                    d2(i) = dmin
                enddo
                tot    = max(sum(d2), 1.d-300)
                rr     = real(ran3(), dp)
                cum(1) = d2(1) / tot
                do i = 2,n
                    cum(i) = cum(i-1) + d2(i) / tot
                enddo
                jb = n
                do i = 1,n
                    if( cum(i) >= rr )then
                        jb = i
                        exit
                    endif
                enddo
                C(:,j) = X(:,jb)
            enddo
            do it = 1,niter
                do i = 1,n
                    dmin = huge(dmin)
                    jb   = 1
                    do j = 1,k
                        dd = sum((X(:,i) - C(:,j))**2)
                        if( dd < dmin )then
                            dmin = dd
                            jb   = j
                        endif
                    enddo
                    l(i) = jb
                enddo
                Cn  = 0.d0
                cnt = 0
                do i = 1,n
                    Cn(:,l(i)) = Cn(:,l(i)) + X(:,i)
                    cnt(l(i))  = cnt(l(i)) + 1
                enddo
                do j = 1,k
                    if( cnt(j) > 0 )then
                        Cn(:,j) = Cn(:,j) / real(cnt(j),dp)
                    else
                        Cn(:,j) = C(:,j)
                    endif
                enddo
                C = Cn
            enddo
            inert = 0.d0
            do i = 1,n
                inert = inert + sum((X(:,i) - C(:,l(i)))**2)
            enddo
            if( inert < best )then
                best = inert
                lab  = l
                cent = C
            endif
        enddo
    end subroutine kt_kmeans_pp

    subroutine tree_init_positions( meta, subset, n1, n2, out )
        character(len=*), intent(in) :: meta, subset, out
        integer,          intent(in) :: n1, n2
        real(dp), allocatable :: Z(:,:), fl(:,:), zz1(:), zz2(:), g1(:), g2(:), mu(:), sd(:)
        integer,  allocatable :: sub(:), pos(:), okrow(:)
        character(len=PLEN) :: zkey
        character(len=256)  :: sdl
        integer :: nz, dz, nf, df, nref, i, j, u, nok, row
        call npz_read(meta, 'z', Z, nz, dz, prefix=.true., keyname=zkey)
        call npz_read(meta, 'flags', fl, nf, df)
        nref = count(fl(1,:) == 1.d0)
        allocate(pos(nf), source=-1)
        j = 0
        do i = 1,nf
            if( fl(1,i) == 1.d0 )then
                if( nz == nref )then
                    pos(i) = j
                else
                    pos(i) = i - 1
                endif
                j = j + 1
            endif
        enddo
        allocate(mu(dz), sd(dz), source=0.d0)
        do i = 1,nz
            mu = mu + Z(:,i)
        enddo
        mu = mu / real(nz,dp)
        do i = 1,nz
            sd = sd + (Z(:,i) - mu)**2
        enddo
        sd = sqrt(sd / real(nz,dp))
        call read_ints(subset, sub)
        allocate(okrow(size(sub)))
        nok = 0
        do i = 1,size(sub)
            if( sub(i) < 1 .or. sub(i) > nf ) cycle
            if( nz == nref .and. pos(sub(i)) < 0 ) cycle
            nok = nok + 1
            okrow(nok) = i
        enddo
        allocate(zz1(nok), zz2(nok), g1(nok), g2(nok))
        do i = 1,nok
            row = sub(okrow(i))
            if( nz == nref )then
                j = pos(row) + 1
            else
                j = row
            endif
            zz1(i) = (Z(1,j) - mu(1)) / max(sd(1), 1.d-12)
            zz2(i) = (Z(2,j) - mu(2)) / max(sd(2), 1.d-12)
        enddo
        call scale_to_grid(zz1, n1, g1)
        call scale_to_grid(zz2, n2, g2)
        open(newunit=u, file=trim(out), status='replace', action='write')
        do i = 1,nok
            write(u,'(I0,1X,F0.4,1X,F0.4)') sub(okrow(i)), g1(i), g2(i)
        enddo
        close(u)
        sdl = ''
        do i = 1,min(4,dz)
            sdl = trim(sdl)//' '//fr(sd(i),3)
        enddo
        write(logfhandle,'(A)') trim(out)//': '//int2str(nok)//' of '//int2str(size(sub))//' particles placed from '//trim(zkey)//&
            &' PCs 1-2 (std ['//trim(adjustl(sdl))//'])'
        call flush(logfhandle)

      contains

        subroutine scale_to_grid( v, nn, g )
            real(dp), intent(in)  :: v(:)
            integer,  intent(in)  :: nn
            real(dp), intent(out) :: g(:)
            real(dp) :: p5, p95
            p5  = kt_percentile(v, 5.d0)
            p95 = kt_percentile(v, 95.d0)
            g   = min(max((v - p5) / max(p95 - p5, 1.d-9), 0.d0), 1.d0) * real(nn - 1, dp)
        end subroutine scale_to_grid

    end subroutine tree_init_positions

    subroutine tree_parent_subset( meta, out, n )
        character(len=*), intent(in)  :: meta, out
        integer,          intent(out) :: n
        real(dp), allocatable :: fl(:,:)
        integer :: nf, df, i, u
        call npz_read(meta, 'flags', fl, nf, df)
        n = 0
        open(newunit=u, file=trim(out), status='replace', action='write')
        do i = 1,nf
            if( nint(fl(1,i)) == 1 )then
                write(u,'(I0)') i
                n = n + 1
            endif
        enddo
        close(u)
        write(logfhandle,'(A)') 'parent subset: '//int2str(n)//' particles ('//int2str(n)//' flagged)'
    end subroutine tree_parent_subset

    subroutine tree_families( dir, nfam, outdir, nmem )
        character(len=*), intent(in)  :: dir, outdir
        integer,          intent(in)  :: nfam
        integer,          intent(out) :: nmem(nfam)
        type(board_t) :: B
        character(len=PLEN), allocatable :: fc(:), fs(:), files(:)
        real,     allocatable :: v(:,:,:)
        real(dp), allocatable :: mean(:), Dc(:,:), S(:,:), W(:,:), L(:,:), ev(:), dg(:), work(:), emb(:,:), cent(:,:), col(:)
        integer,  allocatable :: fam(:), iwork(:), lab(:)
        integer  :: nc, ns, K, i, j, npix, n(3), info, lwork, liwork, u, kk
        real(dp) :: med, m
        call list_sorted(trim(dir)//'/cells/cell_*.mrc', fc, nc)
        call list_sorted(trim(dir)//'/cells/state_*.mrc', fs, ns)
        K = nc + ns
        if( K < 2 ) THROW_HARD('tree_families: fewer than two cell maps in '//trim(dir))
        allocate(files(K))
        files(1:nc) = fc
        if( ns > 0 ) files(nc+1:K) = fs
        call read_mrc(files(1), v, n)
        npix = product(n)
        allocate(mean(npix), source=0.d0)
        do i = 1,K
            call read_mrc(files(i), v, n)
            mean = mean + real(reshape(v, [npix]), dp)
        enddo
        mean = mean / real(K,dp)
        allocate(Dc(npix,K))
        do i = 1,K
            call read_mrc(files(i), v, n)
            Dc(:,i) = real(reshape(v, [npix]), dp) - mean
            m = sum(Dc(:,i)) / real(npix,dp)
            Dc(:,i) = Dc(:,i) - m
        enddo
        allocate(S(K,K), col(K))
        do i = 1,K
            col(i) = sqrt(sum(Dc(:,i)**2))
        enddo
        !$omp parallel do default(shared) private(i,j) schedule(dynamic)
        do j = 1,K
            do i = 1,j
                S(i,j) = dot_product(Dc(:,i), Dc(:,j)) / max(col(i) * col(j), 1.d-30)
                S(j,i) = S(i,j)
            enddo
        enddo
        !$omp end parallel do
        deallocate(Dc, mean)
        allocate(W(K,K), L(K,K), ev(K), dg(K))
        W = max(S, 0.d0)
        do i = 1,K
            W(i,i) = 0.d0
        enddo
        do i = 1,K
            dg(i) = sum(W(i,:)) + 1.d-9
        enddo
        do j = 1,K
            do i = 1,K
                L(i,j) = merge(1.d0, 0.d0, i == j) - W(i,j) / sqrt(dg(i) * dg(j))
            enddo
        enddo
        lwork  = 1 + 6*K + 2*K*K
        liwork = 3 + 5*K
        allocate(work(lwork), iwork(liwork))
        call dsyevd('V', 'L', K, L, K, ev, work, lwork, iwork, liwork, info)
        if( info /= 0 ) THROW_HARD('tree_families: dsyevd failed')
        allocate(fam(K))
        if( nfam == 2 )then
            med = kt_median(L(:,2))
            do i = 1,K
                fam(i) = merge(2, 1, L(i,2) > med)
            enddo
        else
            allocate(emb(nfam-1,K), lab(K), cent(nfam-1,nfam))
            emb = transpose(L(:,2:nfam))
            call kt_kmeans_pp(emb, nfam, 1, 100, lab, cent)
            fam = lab
        endif
        call read_board(dir, B, with_halves=.false.)
        nmem = 0
        do kk = 1,nfam
            open(newunit=u, file=trim(outdir)//'/family'//int2str(kk)//'_subset.txt', status='replace', action='write')
            do i = 1,B%n
                if( B%gm(i) > 0.5d0 .and. B%am(i) >= 1 .and. B%am(i) <= K )then
                    if( fam(B%am(i)) == kk )then
                        write(u,'(I0)') B%pinds(i)
                        nmem(kk) = nmem(kk) + 1
                    endif
                endif
            enddo
            close(u)
            write(logfhandle,'(A)') 'family '//int2str(kk)//': '//int2str(count(fam == kk))//' states, '//int2str(nmem(kk))//&
                &' particles; within-family mean state similarity '//fr(mean_sim(kk, kk), 2)
        enddo
        if( nfam == 2 ) write(logfhandle,'(A)') 'between-family mean state similarity '//fr(mean_sim(1, 2), 2)
        open(newunit=u, file=trim(outdir)//'/families.txt', status='replace', action='write')
        write(u,'(A)') '# state family'
        do i = 1,K
            write(u,'(I0,1X,I0)') i, fam(i)
        enddo
        close(u)
        call kill_board(B)
        call flush(logfhandle)

      contains

        function mean_sim( a, b ) result( ms )
            integer, intent(in) :: a, b
            real(dp) :: ms
            integer  :: i_, j_, cnt
            ms  = 0.d0
            cnt = 0
            do j_ = 1,K
                do i_ = 1,K
                    if( fam(i_) == a .and. fam(j_) == b )then
                        ms  = ms + S(i_,j_)
                        cnt = cnt + 1
                    endif
                enddo
            enddo
            ms = ms / real(max(cnt,1),dp)
        end function mean_sim

    end subroutine tree_families

    subroutine tree_regions( B, bname, box, smpd, mode, outdir, maxreg, kgroups, thr, minvox, ndilate, cmerge, minrep, nreg )
        type(board_t),    intent(in)  :: B
        character(len=*), intent(in)  :: bname, outdir
        integer,          intent(in)  :: box, mode, maxreg, kgroups, minvox, ndilate
        real,             intent(in)  :: smpd, thr, cmerge, minrep
        integer,          intent(out) :: nreg
        type(fft3_t) :: F
        real(dp), allocatable :: v(:,:,:), r0(:,:,:), r1(:,:,:), fm(:,:,:), devm(:,:), prof(:,:), cent(:,:), Sc(:,:)
        real(dp), allocatable :: soft(:,:,:), fullm(:,:,:), uni(:,:,:), mass(:), pmass(:)
        logical,  allocatable :: m(:,:,:), comp(:,:,:)
        integer,  allocatable :: lab(:,:,:), sizes(:), order(:), keep(:), idx(:), labv(:), parent(:), pgrp(:), pptr(:), pvox(:)
        integer,  allocatable :: porder(:), gl(:,:,:)
        integer  :: n, nk, k, kk, i, j, l, nl, nmov, c, g, np, nsel, nv, ngl
        real(dp) :: p99, smpd_c, cx, cy, cz, rad, radmax, mx, vmov
        n      = B%n3
        smpd_c = real(box,dp) * real(smpd,dp) / real(n,dp)
        nk     = count(B%crep(1:B%ngrid) >= minrep)
        allocate(keep(nk))
        kk = 0
        do k = 1,B%ngrid
            if( B%crep(k) >= minrep )then
                kk = kk + 1
                keep(kk) = k
            endif
        enddo
        allocate(v(n,n,n), source=0.d0)
        allocate(r0(n,n,n), r1(n,n,n), fm(n,n,n), soft(n,n,n), fullm(box,box,box), uni(box,box,box))
        if( mode == 2 ) allocate(devm(n*n*n,nk))
        call fft3_new(F, n)
        do kk = 1,nk
            call band_to_real(B, B%D0(:,keep(kk)), F, r0)
            call band_to_real(B, B%D1(:,keep(kk)), F, r1)
            v = v + r0 * r1
            if( mode == 2 ) devm(:,kk) = reshape(0.5d0 * (r0 + r1), [n*n*n])
        enddo
        call fft3_kill(F)
        v = max(v, 0.d0)
        call kt_gauss(v, 1.5)
        p99 = kt_percentile(reshape(v, [n*n*n]), 99.d0)
        fm  = sqrt(min(max(v / max(p99, 1.d-30), 0.d0), 1.d0))
        uni = 0.d0
        call simple_mkdir(trim(outdir))
        allocate(m(n,n,n), comp(n,n,n), lab(n,n,n))
        nreg = 0
        if( mode == 1 )then
            m = fm > real(thr,dp)
            call kt_erode(m)
            call kt_dilate(m, 1)
            call kt_label6(m, lab, nl)
            allocate(sizes(max(nl,1)), mass(max(nl,1)), order(max(nl,1)))
            sizes = 0
            mass  = 0.d0
            do l = 1,n
                do j = 1,n
                    do i = 1,n
                        if( lab(i,j,l) > 0 )then
                            sizes(lab(i,j,l)) = sizes(lab(i,j,l)) + 1
                            mass(lab(i,j,l))  = mass(lab(i,j,l)) + v(i,j,l)
                        endif
                    enddo
                enddo
            enddo
            if( nl > 0 ) call kt_argsort_desc(mass(1:nl), order(1:nl))
            nsel = 0
            do c = 1,nl
                if( sizes(order(c)) >= minvox .and. nsel < maxreg ) nsel = nsel + 1
            enddo
            write(logfhandle,'(A)') 'board '//trim(bname)//': '//int2str(nk)//' of '//int2str(B%ngrid)//&
                &' cells survive; coarse grid '//int2str(n)//'^3 ('//fr(smpd_c,1)//' A/px); variability above '//fr(real(thr,dp),2)//&
                &': '//int2str(count(m))//' voxels in '//int2str(nl)//' components, '//int2str(nsel)//' kept (>= '//int2str(minvox)//' voxels)'
            do c = 1,nl
                if( nreg >= maxreg ) exit
                k = order(c)
                if( sizes(k) < minvox ) cycle
                nreg = nreg + 1
                comp = lab == k
                cx = 0.d0
                cy = 0.d0
                cz = 0.d0
                radmax = 0.d0
                do l = 1,n
                    do j = 1,n
                        do i = 1,n
                            if( comp(i,j,l) )then
                                cx = cx + real(i-1,dp)
                                cy = cy + real(j-1,dp)
                                cz = cz + real(l-1,dp)
                                radmax = max(radmax, sqrt((real(i-1,dp) - real(n,dp)/2.d0)**2 + (real(j-1,dp) - real(n,dp)/2.d0)**2 + &
                                    &(real(l-1,dp) - real(n,dp)/2.d0)**2) * smpd_c)
                            endif
                        enddo
                    enddo
                enddo
                cx  = cx / real(sizes(k),dp)
                cy  = cy / real(sizes(k),dp)
                cz  = cz / real(sizes(k),dp)
                rad = sqrt((cx - real(n,dp)/2.d0)**2 + (cy - real(n,dp)/2.d0)**2 + (cz - real(n,dp)/2.d0)**2) * smpd_c
                call write_region(nreg)
                write(logfhandle,'(A)') '  region '//int2str(nreg)//': '//int2str(sizes(k))//' coarse voxels ('//&
                    &fr(100.d0 * mass(k) / sum(mass(1:nl)), 1)//' % of the variability mass), centroid at '//fr(rad,0)//&
                    &' A radius, max radius '//fr(radmax,0)//' A'
            enddo
            call write_mrc(trim(outdir)//'/all_moving.mrc', uni, smpd)
            call write_mrc(trim(outdir)//'/variability_coarse.mrc', fm, real(smpd_c))
            write(logfhandle,'(A)') 'masks -> '//trim(outdir)//'/region_*.mrc ('//int2str(nreg)//'), all_moving.mrc, variability_coarse.mrc'
        else
            m    = fm > real(thr,dp)
            nmov = count(m)
            allocate(idx(nmov), prof(nk,nmov), labv(nmov), cent(nk,kgroups), parent(kgroups), Sc(kgroups,kgroups))
            c = 0
            do l = 1,n
                do j = 1,n
                    do i = 1,n
                        if( m(i,j,l) )then
                            c = c + 1
                            idx(c) = i + (j-1)*n + (l-1)*n*n
                        endif
                    enddo
                enddo
            enddo
            do c = 1,nmov
                prof(:,c) = devm(idx(c),:)
                prof(:,c) = prof(:,c) - sum(prof(:,c)) / real(nk,dp)
                prof(:,c) = prof(:,c) / max(sqrt(sum(prof(:,c)**2)), 1.d-30)
            enddo
            call kt_kmeans_pp(prof, kgroups, 12, 10, labv, cent)
            do g = 1,kgroups
                cent(:,g) = cent(:,g) / max(sqrt(sum(cent(:,g)**2)), 1.d-30)
            enddo
            Sc = matmul(transpose(cent), cent)
            do g = 1,kgroups
                parent(g) = g
            enddo
            do i = 1,kgroups
                do j = i+1,kgroups
                    if( Sc(i,j) > real(cmerge,dp) ) parent(j) = parent(i)
                enddo
            enddo
            do c = 1,nmov
                labv(c) = parent(labv(c))
            enddo
            mx = -huge(mx)
            do i = 1,kgroups
                do j = 1,kgroups
                    if( i /= j ) mx = max(mx, Sc(i,j))
                enddo
            enddo
            ngl = 0
            do g = 1,kgroups
                if( any(labv == g) ) ngl = ngl + 1
            enddo
            write(logfhandle,'(A)') 'board '//trim(bname)//': '//int2str(nk)//' cells; '//int2str(nmov)//&
                &' moving coarse voxels (fm > '//fr(real(thr,dp),2)//'); k-means k='//int2str(kgroups)//' on motion profiles -> '//&
                &int2str(ngl)//' groups after merging centroids above '//fr(real(cmerge,dp),2)//' (centroid corr matrix off-diagonal max '//fr(mx,2)//')'
            allocate(pgrp(max(nmov,1)), pptr(nmov+1), pvox(max(nmov,1)), pmass(max(nmov,1)), gl(n,n,n))
            np = 0
            pptr(1) = 1
            do g = 1,kgroups
                if( .not. any(labv == g) ) cycle
                m = .false.
                do c = 1,nmov
                    if( labv(c) == g ) m(1 + mod(idx(c)-1, n), 1 + mod((idx(c)-1)/n, n), 1 + (idx(c)-1)/(n*n)) = .true.
                enddo
                call kt_label6(m, gl, nl)
                do k = 1,nl
                    nv = count(gl == k)
                    if( nv < minvox ) cycle
                    np = np + 1
                    pgrp(np)   = g
                    pmass(np)  = 0.d0
                    pptr(np+1) = pptr(np)
                    do l = 1,n
                        do j = 1,n
                            do i = 1,n
                                if( gl(i,j,l) == k )then
                                    pvox(pptr(np+1)) = i + (j-1)*n + (l-1)*n*n
                                    pptr(np+1) = pptr(np+1) + 1
                                    pmass(np)  = pmass(np) + v(i,j,l)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            allocate(porder(max(np,1)))
            if( np > 0 ) call kt_argsort_desc(pmass(1:np), porder(1:np))
            vmov = 0.d0
            do l = 1,n
                do j = 1,n
                    do i = 1,n
                        if( fm(i,j,l) > real(thr,dp) ) vmov = vmov + v(i,j,l)
                    enddo
                enddo
            enddo
            do c = 1,min(np, maxreg)
                k = porder(c)
                comp = .false.
                cx = 0.d0
                cy = 0.d0
                cz = 0.d0
                radmax = 0.d0
                do i = pptr(k),pptr(k+1)-1
                    l = pvox(i)
                    comp(1 + mod(l-1, n), 1 + mod((l-1)/n, n), 1 + (l-1)/(n*n)) = .true.
                    cx = cx + real(mod(l-1, n),dp)
                    cy = cy + real(mod((l-1)/n, n),dp)
                    cz = cz + real((l-1)/(n*n),dp)
                    radmax = max(radmax, sqrt((real(mod(l-1, n),dp) - real(n,dp)/2.d0)**2 + (real(mod((l-1)/n, n),dp) - real(n,dp)/2.d0)**2 + &
                        &(real((l-1)/(n*n),dp) - real(n,dp)/2.d0)**2) * smpd_c)
                enddo
                nv  = pptr(k+1) - pptr(k)
                cx  = cx / real(nv,dp)
                cy  = cy / real(nv,dp)
                cz  = cz / real(nv,dp)
                rad = sqrt((cx - real(n,dp)/2.d0)**2 + (cy - real(n,dp)/2.d0)**2 + (cz - real(n,dp)/2.d0)**2) * smpd_c
                nreg = nreg + 1
                call write_region(nreg)
                write(logfhandle,'(A)') '  region '//int2str(nreg)//' (motion group '//int2str(pgrp(k)-1)//'): '//int2str(nv)//&
                    &' coarse voxels, '//fr(100.d0 * pmass(k) / max(vmov, 1.d-30), 1)//' % of the variability mass, centroid radius '//&
                    &fr(rad,0)//' A, max radius '//fr(radmax,0)//' A'
            enddo
            call write_mrc(trim(outdir)//'/all_moving.mrc', uni, smpd)
            write(logfhandle,'(A)') 'masks -> '//trim(outdir)//'/region_*.mrc ('//int2str(nreg)//')'
        endif
        call flush(logfhandle)

      contains

        subroutine write_region( ireg )
            integer, intent(in) :: ireg
            call kt_dilate(comp, ndilate)
            soft = merge(1.d0, 0.d0, comp)
            call kt_gauss(soft, 1.0)
            mx   = maxval(soft)
            soft = min(max(soft / max(mx, 1.d-30), 0.d0), 1.d0)
            call kt_zoom(soft, box, fullm)
            fullm = min(max(fullm, 0.d0), 1.d0)
            call write_mrc(trim(outdir)//'/region_'//int2str(ireg)//'.mrc', fullm, smpd)
            uni = max(uni, fullm)
        end subroutine write_region

    end subroutine tree_regions

    subroutine tree_node_states( B, bname, nrows, sim, nmax, frac, minmass, maxk, minrep, support, outdir, nstates )
        type(board_t),    intent(in)  :: B
        character(len=*), intent(in)  :: bname, outdir
        integer,          intent(in)  :: nrows, nmax, maxk
        real,             intent(in)  :: sim, frac, minmass, minrep
        logical,          intent(in)  :: support
        integer,          intent(out) :: nstates
        integer,  allocatable :: keep(:), lab(:), lab2(:), order(:), states(:), am(:), mine(:), top(:), cells(:)
        real(dp), allocatable :: S(:,:), massc(:), M(:,:), col(:), okrep(:), oksep(:), med(:)
        logical,  allocatable :: lok(:)
        character(len=8192) :: line
        integer  :: nk, i, j, ii, jj, K, c, d, e, rank, nsel, nmine, u, p
        real(dp) :: rep, sep, mrep, msep, best, x
        nk = count(B%crep(1:B%ngrid) >= minrep)
        allocate(keep(nk))
        j = 0
        do i = 1,B%ngrid
            if( B%crep(i) >= minrep )then
                j = j + 1
                keep(j) = i
            endif
        enddo
        allocate(S(nk,nk), source=0.d0)
        !$omp parallel do default(shared) private(ii,jj) schedule(dynamic)
        do jj = 1,nk
            do ii = 1,jj
                S(ii,jj) = 0.5d0 * (cc_c(B%D0(:,keep(ii)), B%D1(:,keep(jj))) + cc_c(B%D1(:,keep(ii)), B%D0(:,keep(jj))))
                S(jj,ii) = S(ii,jj)
            enddo
        enddo
        !$omp end parallel do
        allocate(lab(nk))
        call kt_avg_linkage(S, 1.d0 - real(sim,dp), lab, K)
        allocate(massc(K))
        call calc_mass
        do while( support .and. K > 1 )
            allocate(lok(K), okrep(K), oksep(K))
            do c = 1,K
                call supported(c, lok(c), okrep(c), oksep(c))
            enddo
            if( all(lok) )then
                deallocate(lok, okrep, oksep)
                exit
            endif
            c    = 0
            best = huge(best)
            do e = 1,K
                if( lok(e) ) cycle
                if( massc(e) < best )then
                    best = massc(e)
                    c    = e
                endif
            enddo
            rep  = okrep(c)
            sep  = oksep(c)
            d    = 0
            best = -huge(best)
            do e = 1,K
                if( e == c ) cycle
                x = xh(cells_of(lab, c), cells_of(lab, e))
                if( x > best )then
                    best = x
                    d    = e
                endif
            enddo
            allocate(lab2(nk), source=lab)
            where( lab2 == c ) lab2 = d
            mrep = xh(cells_of(lab2, d), cells_of(lab2, d))
            msep = -1.d0
            do e = 1,K
                if( e == c .or. e == d ) cycle
                msep = max(msep, xh(cells_of(lab2, d), cells_of(lab2, e)))
            enddo
            if( mrep >= real(sim,dp) .and. msep < real(sim,dp) )then
                write(logfhandle,'(A)') '  cluster of '//int2str(count(lab == c))//' cells (mass '//fr(massc(c),0)//') unsupported: reproducibility '//&
                    &fr(rep,2)//' separability '//fr(sep,2)//' -> merged into its neighbour ('//int2str(count(lab == d))//' cells, mass '//&
                    &fr(massc(d),0)//'); merged: '//fr(mrep,2)//' / '//fr(msep,2)
                lab = lab2
            else
                write(logfhandle,'(A)') '  cluster of '//int2str(count(lab == c))//' cells (mass '//fr(massc(c),0)//') unsupported: reproducibility '//&
                    &fr(rep,2)//' separability '//fr(sep,2)//' -> no state (merge would give '//fr(mrep,2)//' / '//fr(msep,2)//')'
                keep = pack(keep, lab /= c)
                lab  = pack(lab,  lab /= c)
                nk   = size(keep)
            endif
            where( lab > c ) lab = lab - 1
            K = K - 1
            deallocate(lab2, lok, okrep, oksep, massc)
            allocate(massc(K))
            call calc_mass
        enddo
        allocate(order(K))
        call kt_argsort_desc(massc, order)
        line = ''
        do c = 1,min(K,12)
            line = trim(line)//' '//fr(massc(order(c)),0)
        enddo
        write(logfhandle,'(A)') 'board '//trim(bname)//': '//int2str(nk)//' of '//int2str(B%ngrid)//' cells survive, '//int2str(K)//&
            &' clusters at similarity '//fr(real(sim,dp),1)//'; cluster masses (particles):'//trim(line)
        allocate(M(B%n,K), source=0.d0)
        do c = 1,K
            do ii = 1,nk
                if( lab(ii) /= c ) cycle
                do p = 1,B%n
                    M(p,c) = M(p,c) + real(B%r(keep(ii),p), dp)
                enddo
            enddo
        enddo
        allocate(am(B%n))
        do p = 1,B%n
            am(p) = maxloc(M(p,:), 1)
        enddo
        allocate(states(nrows), source=0)
        call simple_mkdir(trim(outdir))
        open(newunit=u, file=trim(outdir)//'/clusters.txt', status='replace', action='write')
        write(u,'(A)') '# state ncells mass nselected cells'
        nstates = 0
        do rank = 1,min(K, maxk)
            c = order(rank)
            if( massc(c) < real(minmass,dp) ) cycle
            nmine = count(am == c)
            mine  = pack([(p, p=1,B%n)], am == c)
            nsel  = int(min(real(nmax,dp), real(frac,dp) * massc(c), real(nmine,dp)))
            allocate(col(nmine), top(nmine))
            do i = 1,nmine
                col(i) = M(mine(i),c)
            enddo
            call kt_argsort_desc(col, top)
            do i = 1,nsel
                states(B%pinds(mine(top(i)))) = rank
            enddo
            cells = pack(keep, lab == c)
            line  = ''
            do i = 1,size(cells)
                line = trim(line)//' '//int2str(cells(i))
            enddo
            write(u,'(A)') int2str(rank)//' '//int2str(size(cells))//' '//fr(massc(c),0)//' '//int2str(nsel)//trim(line)
            allocate(med(nsel))
            do i = 1,nsel
                med(i) = col(top(i))
            enddo
            write(logfhandle,'(A)') '  state '//int2str(rank)//': '//int2str(size(cells))//' cells, mass '//fr(massc(c),0)//', selected '//&
                &int2str(nsel)//' particles (mass in the cluster of the selected: median '//fr(kt_median(med),2)//')'
            deallocate(col, top, med)
            nstates = nstates + 1
        enddo
        close(u)
        call write_ints(trim(outdir)//'/states_clusters.txt', states)
        write(logfhandle,'(A)') int2str(nstates)//' states -> '//trim(outdir)
        call flush(logfhandle)

      contains

        subroutine calc_mass
            integer :: c_, ii_, p_
            massc = 0.d0
            do ii_ = 1,nk
                c_ = lab(ii_)
                do p_ = 1,B%n
                    massc(c_) = massc(c_) + real(B%r(keep(ii_),p_), dp)
                enddo
            enddo
        end subroutine calc_mass

        function cells_of( lb, c_ ) result( cl )
            integer, intent(in) :: lb(:), c_
            integer, allocatable :: cl(:)
            cl = pack(keep, lb == c_)
        end function cells_of

        function xh( ca, cb ) result( val )
            integer, intent(in) :: ca(:), cb(:)
            real(dp) :: val
            complex(dp), allocatable :: a0(:), a1(:), b0(:), b1(:)
            integer :: i_
            allocate(a0(B%nb), a1(B%nb), b0(B%nb), b1(B%nb), source=cmplx(0.d0,0.d0,dp))
            do i_ = 1,size(ca)
                a0 = a0 + B%D0(:,ca(i_))
                a1 = a1 + B%D1(:,ca(i_))
            enddo
            do i_ = 1,size(cb)
                b0 = b0 + B%D0(:,cb(i_))
                b1 = b1 + B%D1(:,cb(i_))
            enddo
            a0  = a0 / real(size(ca),dp)
            a1  = a1 / real(size(ca),dp)
            b0  = b0 / real(size(cb),dp)
            b1  = b1 / real(size(cb),dp)
            val = 0.5d0 * (cc_z(a0, b1) + cc_z(a1, b0))
        end function xh

        subroutine supported( c_, ok, rp, sp )
            integer,  intent(in)  :: c_
            logical,  intent(out) :: ok
            real(dp), intent(out) :: rp, sp
            integer :: e_
            rp = xh(cells_of(lab, c_), cells_of(lab, c_))
            sp = -1.d0
            do e_ = 1,K
                if( e_ == c_ ) cycle
                sp = max(sp, xh(cells_of(lab, c_), cells_of(lab, e_)))
            enddo
            ok = rp >= real(sim,dp) .and. sp < real(sim,dp)
        end subroutine supported

    end subroutine tree_node_states

    function cc_c( x, y ) result( c )
        complex, intent(in) :: x(:), y(:)
        real(dp) :: c, num, nx, ny
        integer  :: i
        num = 0.d0
        nx  = 0.d0
        ny  = 0.d0
        do i = 1,size(x)
            num = num + real(conjg(cmplx(x(i),kind=dp)) * cmplx(y(i),kind=dp), dp)
            nx  = nx + abs(cmplx(x(i),kind=dp))**2
            ny  = ny + abs(cmplx(y(i),kind=dp))**2
        enddo
        c = num / max(sqrt(nx * ny), 1.d-30)
    end function cc_c

    function cc_z( x, y ) result( c )
        complex(dp), intent(in) :: x(:), y(:)
        real(dp) :: c
        c = real(sum(conjg(x) * y), dp) / max(sqrt(sum(abs(x)**2) * sum(abs(y)**2)), 1.d-30)
    end function cc_z

    subroutine tree_rank_pc1( B, ntop, nrows, outdir )
        type(board_t),    intent(in) :: B
        integer,          intent(in) :: ntop, nrows
        character(len=*), intent(in) :: outdir
        complex,  allocatable :: A(:,:), Dc(:,:), U(:,:), VT(:,:), work(:), mu(:)
        real,     allocatable :: sv(:), rwork(:), cs(:)
        real(dp), allocatable :: s(:)
        integer,  allocatable :: iwork(:), order(:), st(:)
        character(len=256) :: svl
        integer :: K, nb, i, kk, info, lwork, mn, mx, uu
        complex :: wq(1)
        K  = B%ngrid
        nb = B%nb
        allocate(A(K,nb), mu(nb))
        do kk = 1,K
            A(kk,:) = 0.5 * (B%D0(:,kk) + B%D1(:,kk))
        enddo
        mu = cmplx(0.,0.)
        do kk = 1,K
            mu = mu + A(kk,:)
        enddo
        mu = mu / real(K)
        do kk = 1,K
            A(kk,:) = A(kk,:) - mu
        enddo
        allocate(Dc(K,nb), source=A)
        mn = min(K,nb)
        mx = max(K,nb)
        allocate(sv(mn), U(K,mn), VT(mn,nb), iwork(8*mn), rwork(max(1, mn * max(5*mn + 7, 2*mx + 2*mn + 1))))
        call cgesdd('S', K, nb, A, K, sv, U, K, VT, mn, wq, -1, rwork, iwork, info)
        lwork = max(1, int(real(wq(1))))
        allocate(work(lwork))
        call cgesdd('S', K, nb, A, K, sv, U, K, VT, mn, work, lwork, rwork, iwork, info)
        if( info /= 0 ) THROW_HARD('tree_rank_pc1: cgesdd failed')
        allocate(cs(K))
        do kk = 1,K
            cs(kk) = real(sum(cmplx(Dc(kk,:),kind=dp) * conjg(cmplx(VT(1,:),kind=dp))))
        enddo
        cs = cs - sum(cs) / real(K)
        svl = ''
        do i = 1,min(4,mn)
            svl = trim(svl)//' '//fr(real(sv(i) / sv(1), dp), 2)
        enddo
        write(logfhandle,'(A)') 'own main axis: singular values (first 4, relative) ['//trim(adjustl(svl))//']'
        allocate(s(B%n))
        do i = 1,B%n
            s(i) = sum(real(B%r(1:K,i), dp) * real(cs, dp))
        enddo
        allocate(order(B%n))
        call kt_argsort_desc(s, order)
        allocate(st(nrows), source=0)
        do i = 1,min(ntop, B%n)
            st(B%pinds(order(i))) = 1
        enddo
        do i = max(1, B%n - ntop + 1),B%n
            st(B%pinds(order(i))) = 2
        enddo
        call simple_mkdir(trim(outdir))
        call write_ints(trim(outdir)//'/states_top.txt', st)
        open(newunit=uu, file=trim(outdir)//'/score_per_particle.txt', status='replace', action='write')
        do i = 1,B%n
            write(uu,'(A)') int2str(B%pinds(i))//' '//fr(s(i), 4)
        enddo
        close(uu)
        write(logfhandle,'(A)') 'cells scored against the own main axis (pc1, both ends): cell score range '//frs(real(minval(cs),dp),2)//' .. '//&
            &frs(real(maxval(cs),dp),2)//'; particle score quantiles 10/50/90 = '//frs(kt_percentile(s,10.d0),3)//' '//frs(kt_percentile(s,50.d0),3)//&
            &' '//frs(kt_percentile(s,90.d0),3)//'; top '//int2str(ntop)//' >= '//frs(s(order(min(ntop,B%n))),3)
        call flush(logfhandle)
    end subroutine tree_rank_pc1

    subroutine tree_kmeans_states( B, K, nrows, outfile, label )
        type(board_t),    intent(in) :: B
        integer,          intent(in) :: K, nrows
        character(len=*), intent(in) :: outfile, label
        real(dp), allocatable :: X(:,:), C(:,:), Cn(:,:), d(:,:)
        integer,  allocatable :: a(:), bestlab(:), perm(:), cnt(:), st(:), sel(:)
        character(len=1024) :: line
        real(dp) :: best, cost
        integer  :: nx, i, j, s, it
        nx = count(B%gm > 0.5d0)
        allocate(X(2,nx), sel(nx))
        j = 0
        do i = 1,B%n
            if( B%gm(i) > 0.5d0 )then
                j = j + 1
                X(1,j) = B%z1(i)
                X(2,j) = B%z2(i)
                sel(j) = i
            endif
        enddo
        allocate(C(2,K), Cn(2,K), d(nx,K), a(nx), bestlab(nx), perm(nx), cnt(K))
        call seed_rnd_fixed(0)
        best = huge(best)
        do s = 1,10
            do i = 1,nx
                perm(i) = i
            enddo
            call partial_shuffle(perm, K)
            do j = 1,K
                C(:,j) = X(:,perm(j))
            enddo
            do it = 1,100
                do j = 1,K
                    do i = 1,nx
                        d(i,j) = (X(1,i) - C(1,j))**2 + (X(2,i) - C(2,j))**2
                    enddo
                enddo
                do i = 1,nx
                    a(i) = minloc(d(i,:), 1)
                enddo
                Cn  = 0.d0
                cnt = 0
                do i = 1,nx
                    Cn(:,a(i)) = Cn(:,a(i)) + X(:,i)
                    cnt(a(i))  = cnt(a(i)) + 1
                enddo
                do j = 1,K
                    if( cnt(j) > 0 )then
                        Cn(:,j) = Cn(:,j) / real(cnt(j),dp)
                    else
                        Cn(:,j) = C(:,j)
                    endif
                enddo
                if( all(abs(Cn - C) <= 1.d-8 + 1.d-5 * abs(C)) ) exit
                C = Cn
            enddo
            cost = 0.d0
            do i = 1,nx
                cost = cost + d(i,a(i))
            enddo
            if( cost < best )then
                best    = cost
                bestlab = a
            endif
        enddo
        allocate(st(nrows), source=0)
        do i = 1,nx
            st(B%pinds(sel(i))) = bestlab(i)
        enddo
        call write_ints(outfile, st)
        line = ''
        do j = 1,K
            line = trim(line)//' '//int2str(count(bestlab == j))
        enddo
        write(logfhandle,'(A)') trim(label)//': k-means K='//int2str(K)//' on the child board posterior means: sizes ['//trim(adjustl(line))//']'
        call flush(logfhandle)
    end subroutine tree_kmeans_states

    subroutine tree_deposited_match( files, nf, docked, consfile, smpd, lo, hi, save )
        character(len=PLEN),        intent(in) :: files(:)
        integer,                    intent(in) :: nf
        character(len=*),           intent(in) :: docked, consfile
        real,                       intent(in) :: smpd, lo, hi
        character(len=*), optional, intent(in) :: save
        type(fft3_t) :: F
        character(len=PLEN), allocatable :: dep(:)
        real,     allocatable :: v(:,:,:), cons(:,:,:)
        real(dp), allocatable :: mean(:,:,:), msk(:,:,:), X(:,:), Y(:,:), S(:,:)
        integer,  allocatable :: ids(:), bh(:), bk(:), bl(:), best(:), nbest(:)
        character(len=8192) :: line
        integer  :: nd, n(3), nn, i, j, k, ib, nb, u, nid, p
        real(dp) :: rr, rlo, rhi, bc
        if( nf < 1 ) return
        call list_sorted(trim(docked)//'/conf*_m160.mrc', dep, nd)
        if( nd < 1 ) return
        allocate(ids(nd))
        do j = 1,nd
            p = index(dep(j), '/conf', back=.true.) + 5
            k = p
            do while( k <= len_trim(dep(j)) )
                if( dep(j)(k:k) < '0' .or. dep(j)(k:k) > '9' ) exit
                k = k + 1
            enddo
            read(dep(j)(p:k-1),*) ids(j)
        enddo
        call read_mrc(consfile, cons, n)
        nn  = n(1)
        rlo = real(nn,dp) * real(smpd,dp) / real(lo,dp)
        rhi = real(nn,dp) * real(smpd,dp) / real(hi,dp)
        nb  = 0
        do k = 1,nn
            do j = 1,nn
                do i = 1,nn
                    rr = sqrt(real(kt_freq(i,nn)**2 + kt_freq(j,nn)**2 + kt_freq(k,nn)**2, dp))
                    if( rr >= rlo .and. rr <= rhi ) nb = nb + 1
                enddo
            enddo
        enddo
        allocate(bh(nb), bk(nb), bl(nb))
        ib = 0
        do k = 1,nn
            do j = 1,nn
                do i = 1,nn
                    rr = sqrt(real(kt_freq(i,nn)**2 + kt_freq(j,nn)**2 + kt_freq(k,nn)**2, dp))
                    if( rr >= rlo .and. rr <= rhi )then
                        ib = ib + 1
                        bh(ib) = i
                        bk(ib) = j
                        bl(ib) = k
                    endif
                enddo
            enddo
        enddo
        allocate(mean(nn,nn,nn), source=0.d0)
        if( nf > 1 )then
            do i = 1,nf
                call read_mrc(files(i), v, n)
                mean = mean + real(v,dp)
            enddo
            mean = mean / real(nf,dp)
        else
            mean = real(cons,dp)
        endif
        allocate(msk(nn,nn,nn))
        do k = 1,nn
            do j = 1,nn
                do i = 1,nn
                    rr = sqrt(real((i-1-nn/2)**2 + (j-1-nn/2)**2 + (k-1-nn/2)**2, dp))
                    msk(i,j,k) = 1.d0 / (1.d0 + exp((rr - 100.d0 / real(smpd,dp)) / 2.d0))
                enddo
            enddo
        enddo
        call fft3_new(F, nn)
        allocate(X(2*nb,nf), Y(2*nb,nd))
        do i = 1,nf
            call read_mrc(files(i), v, n)
            F%a = cmplx(real(v,dp) - mean, 0.d0, dp)
            call fft3_fwd(F)
            call pack_band(X(:,i))
        enddo
        do j = 1,nd
            call read_mrc(dep(j), v, n)
            F%a = cmplx((real(v,dp) - real(cons,dp)) * msk, 0.d0, dp)
            call fft3_fwd(F)
            call pack_band(Y(:,j))
        enddo
        call fft3_kill(F)
        allocate(S(nf,nd), best(nf), nbest(nd))
        nbest = 0
        do i = 1,nf
            do j = 1,nd
                S(i,j) = dot_product(X(:,i), Y(:,j))
            enddo
            best(i) = maxloc(S(i,:), 1)
            nbest(best(i)) = nbest(best(i)) + 1
        enddo
        write(logfhandle,'(A)') 'identification vs the '//int2str(nd)//' deposited conformers (deviation corr at '//fr(real(lo,dp),0)//'-'//&
            &fr(real(hi,dp),0)//' A on the coarse grid), '//int2str(nf)//' states:'
        write(logfhandle,'(A)') '  conformer  best state corr  states whose best match it is  (coverage: a conformer with >= 1 state and corr >= 0.3 counts as identified)'
        nid = 0
        do j = 1,nd
            bc = maxval(S(:,j))
            if( nbest(j) >= 1 .and. bc >= 0.3d0 ) nid = nid + 1
            write(line,'(A,I5,A,A,A,I2,A,A)') '  ', ids(j), '      ', frs(bc,2), '            ', nbest(j), '   ', &
                &merge('identified', '-         ', nbest(j) >= 1 .and. bc >= 0.3d0)
            write(logfhandle,'(A)') trim(line)
        enddo
        write(logfhandle,'(A)') '  identified '//int2str(nid)//' / '//int2str(nd)//'; mean best-match corr over states '//fr(sum(maxval(S,2)) / real(nf,dp), 2)
        if( present(save) )then
            open(newunit=u, file=trim(save), status='replace', action='write')
            line = '# state best_conformer best_corr | corr per conformer'
            do j = 1,nd
                line = trim(line)//' '//int2str(ids(j))
            enddo
            write(u,'(A)') trim(line)
            do i = 1,nf
                line = int2str(i)//' '//int2str(ids(best(i)))//' '//fr(S(i,best(i)),4)
                do j = 1,nd
                    line = trim(line)//' '//fr(S(i,j),3)
                enddo
                write(u,'(A)') trim(line)
            enddo
            close(u)
            write(logfhandle,'(A)') '  per-state matches -> '//trim(save)
        endif
        call flush(logfhandle)

      contains

        subroutine pack_band( z )
            real(dp), intent(out) :: z(:)
            integer  :: q
            real(dp) :: nrm
            do q = 1,nb
                z(2*q-1) = real(F%a(bh(q),bk(q),bl(q)), dp)
                z(2*q)   = aimag(F%a(bh(q),bk(q),bl(q)))
            enddo
            nrm = sqrt(sum(z**2))
            z   = z / max(nrm, 1.d-30)
        end subroutine pack_band

    end subroutine tree_deposited_match

    subroutine tree_count_states( dirs, nd, minrep, sims, nsims )
        character(len=PLEN), intent(in) :: dirs(:)
        integer,             intent(in) :: nd, nsims
        real,                intent(in) :: minrep, sims(:)
        type(board_t) :: Bd
        complex,  allocatable :: C0(:,:), C1(:,:), Cm(:,:)
        real(dp), allocatable :: S(:,:), W(:), wsz(:), nrm0(:), nrm1(:), offd(:), selfr(:), depc(:)
        integer,  allocatable :: lab(:), dep(:), ddom(:), cnt(:)
        character(len=8192) :: line
        character(len=PLEN) :: mf
        integer  :: j, k, ntot, c, nbb, i, is, nK, ic, idd, u, io, q, nid, maxid
        real(dp) :: corr
        ntot = 0
        nbb  = -1
        do j = 1,nd
            call read_board(dirs(j), Bd)
            ntot = ntot + count(Bd%crep(1:Bd%ngrid) >= minrep)
            if( nbb < 0 ) nbb = Bd%nb
            if( Bd%nb /= nbb ) THROW_HARD('count_states: the boards have different coarse grids')
            call kill_board(Bd)
        enddo
        write(logfhandle,'(A)') int2str(ntot)//' surviving cells (reproducibility >= '//fr(real(minrep,dp),2)//') from '//int2str(nd)//' boards'
        if( ntot == 0 ) return
        allocate(C0(ntot,nbb), C1(ntot,nbb), W(ntot), dep(ntot), depc(ntot))
        dep  = 0
        depc = 0.d0
        c    = 0
        do j = 1,nd
            call read_board(dirs(j), Bd)
            mf = ''
            do q = 1,4
                select case(q)
                    case(1)
                        mf = trim(dirs(j))//'/cell_matches.txt'
                    case(2)
                        mf = trim(dirs(j))//'/matched_rec/cell_matches.txt'
                    case(3)
                        mf = trim(dirs(j))//'/../matched_rec/cell_matches.txt'
                    case(4)
                        mf = trim(dirs(j))//'/../cell_matches.txt'
                end select
                if( file_exists(trim(mf)) ) exit
                mf = ''
            enddo
            do k = 1,Bd%ngrid
                if( Bd%crep(k) < minrep ) cycle
                c = c + 1
                C0(c,:) = Bd%D0(:,k)
                C1(c,:) = Bd%D1(:,k)
                W(c)    = real(Bd%w(k), dp)
                if( len_trim(mf) > 0 )then
                    open(newunit=u, file=trim(mf), status='old', action='read')
                    do
                        read(u,'(A)',iostat=io) line
                        if( io /= 0 ) exit
                        if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                        read(line,*) ic, idd, corr
                        if( ic == k )then
                            dep(c)  = idd
                            depc(c) = corr
                            exit
                        endif
                    enddo
                    close(u)
                endif
            enddo
            call kill_board(Bd)
        enddo
        allocate(nrm0(ntot), nrm1(ntot), Cm(ntot,ntot), S(ntot,ntot))
        do i = 1,ntot
            nrm0(i) = sqrt(sum(abs(cmplx(C0(i,:),kind=dp))**2))
            nrm1(i) = sqrt(sum(abs(cmplx(C1(i,:),kind=dp))**2))
        enddo
        call cgemm('N', 'C', ntot, ntot, nbb, cmplx(1.,0.), C0, ntot, C1, ntot, cmplx(0.,0.), Cm, ntot)
        do j = 1,ntot
            do i = 1,ntot
                S(i,j) = real(real(Cm(i,j)),dp) / max(nrm0(i) * nrm1(j), 1.d-30)
            enddo
        enddo
        allocate(selfr(ntot))
        do i = 1,ntot
            selfr(i) = S(i,i)
        enddo
        S = 0.5d0 * (S + transpose(S))
        allocate(offd(max(1,ntot*(ntot-1)/2)))
        k = 0
        do j = 2,ntot
            do i = 1,j-1
                k = k + 1
                offd(k) = S(i,j)
            enddo
        enddo
        if( k > 0 ) write(logfhandle,'(A)') 'cross-half similarity between surviving cells: median '//fr(kt_median(offd(1:k)),2)//&
            &', 90th pct '//fr(kt_percentile(offd(1:k),90.d0),2)//'; self-reproducibility median '//fr(kt_median(selfr),2)
        allocate(lab(ntot))
        maxid = max(1, maxval(dep))
        do is = 1,nsims
            call kt_avg_linkage(S, 1.d0 - real(sims(is),dp), lab, nK)
            allocate(wsz(nK), source=0.d0)
            do i = 1,ntot
                wsz(lab(i)) = wsz(lab(i)) + W(i)
            enddo
            line = '  similarity '//fr(real(sims(is),dp),1)//': '//int2str(nK)//' states; particles per state median '//&
                &int2str(int(kt_median(wsz)))//' (min '//int2str(int(minval(wsz)))//', max '//int2str(int(maxval(wsz)))//')'
            if( any(dep > 0) )then
                allocate(ddom(nK), cnt(0:maxid))
                do c = 1,nK
                    cnt = 0
                    do i = 1,ntot
                        if( lab(i) == c .and. depc(i) >= 0.3d0 .and. dep(i) > 0 ) cnt(dep(i)) = cnt(dep(i)) + 1
                    enddo
                    ddom(c) = 0
                    if( maxval(cnt) > 0 ) ddom(c) = maxloc(cnt, 1) - 1
                enddo
                nid = 0
                mf  = ''
                do q = 1,maxid
                    if( .not. any(ddom == q) ) cycle
                    nid = nid + 1
                    mf  = trim(mf)//' '//int2str(q)//':'//int2str(count(ddom == q))
                enddo
                line = trim(line)//' | deposited conformers identified '//int2str(nid)//' (states per conformer: '//trim(adjustl(mf))//&
                    &'; unmatched states '//int2str(count(ddom == 0))//')'
                deallocate(ddom, cnt)
            endif
            write(logfhandle,'(A)') trim(line)
            deallocate(wsz)
        enddo
        call flush(logfhandle)
    end subroutine tree_count_states

    function leaf_ok( L1, s, smpd, fscmax, minrec, sph, sh ) result( ok )
        character(len=*), intent(in) :: L1
        integer,          intent(in) :: s, minrec
        real,             intent(in) :: smpd, fscmax
        real(dp),         intent(in) :: sph(:,:,:)
        integer,          intent(in) :: sh(:,:,:)
        logical :: ok
        character(len=PLEN) :: fe, fo
        character(len=4096) :: line
        real, allocatable :: e(:,:,:), o(:,:,:)
        integer :: nrec, u, io, st, nc, nm, n(3)
        ok = .true.
        fe = trim(L1)//'/maps/recvol_state'//int2str_pad(s,2)//'_even.mrc'
        fo = trim(L1)//'/maps/recvol_state'//int2str_pad(s,2)//'_odd.mrc'
        if( .not. file_exists(trim(fe)) .or. .not. file_exists(trim(fo)) ) return
        nrec = 0
        if( file_exists(trim(L1)//'/clusters.txt') )then
            open(newunit=u, file=trim(L1)//'/clusters.txt', status='old', action='read')
            do
                read(u,'(A)',iostat=io) line
                if( io /= 0 ) exit
                if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                read(line,*) st, nc, nm, nrec
                if( st == s ) exit
                nrec = 0
            enddo
            close(u)
        endif
        if( nrec > 0 .and. nrec < minrec )then
            ok = .false.
            return
        endif
        call read_mrc(fe, e, n)
        call read_mrc(fo, o, n)
        ok = kt_fsc_res(real(e,dp), real(o,dp), sph, sh, smpd) <= real(fscmax,dp)
    end function leaf_ok

    subroutine build_leaves( R, kind, nfam, leaves, nl )
        character(len=*),          intent(in)    :: R, kind
        integer,                   intent(in)    :: nfam
        type(leaf_t), allocatable, intent(inout) :: leaves(:)
        integer,                   intent(out)   :: nl
        type(string), allocatable :: sub(:)
        character(len=PLEN), allocatable :: names(:), evens(:)
        character(len=PLEN) :: L1, nm, sfx
        character(len=4096) :: line
        integer :: f, i, j, ne, nd, s, u, io, st, nc, nm_, nrec, cap, p
        character(len=16) :: pre
        cap = 4096
        if( allocated(leaves) ) deallocate(leaves)
        allocate(leaves(cap))
        nl = 0
        pre = 'region'
        if( trim(kind) == 'zoom' ) pre = 'zoomregion'
        do f = 1,nfam
            sub = simple_list_dirs(trim(R)//'/fam'//int2str(f))
            if( .not. allocated(sub) ) cycle
            nd = size(sub)
            allocate(names(nd))
            j = 0
            do i = 1,nd
                nm = sub(i)%to_char()
                if( index(nm, trim(pre)) /= 1 ) cycle
                if( trim(kind) == 'region' .and. index(nm, 'zoom') == 1 ) cycle
                if( .not. dir_exists(string(trim(R)//'/fam'//int2str(f)//'/'//trim(nm)//'/level1/maps')) ) cycle
                j = j + 1
                names(j) = nm
            enddo
            call sort_names(names, j)
            do i = 1,j
                L1  = trim(R)//'/fam'//int2str(f)//'/'//trim(names(i))//'/level1'
                p   = index(names(i), 'region') + 6
                sfx = names(i)(p:)
                call list_sorted(trim(L1)//'/maps/recvol_state??_even.mrc', evens, ne)
                do s = 1,ne
                    p = index(evens(s), 'recvol_state', back=.true.) + 12
                    read(evens(s)(p:p+1),*) st
                    if( .not. file_exists(trim(L1)//'/maps/recvol_state'//int2str_pad(st,2)//'_odd.mrc') ) cycle
                    nrec = 1
                    if( file_exists(trim(L1)//'/clusters.txt') )then
                        open(newunit=u, file=trim(L1)//'/clusters.txt', status='old', action='read')
                        do
                            read(u,'(A)',iostat=io) line
                            if( io /= 0 ) exit
                            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                            read(line,*) nc, nm_, nm_, nm_
                            if( nc == st )then
                                read(line,*) nc, nm_, nm_, nrec
                                exit
                            endif
                        enddo
                        close(u)
                    endif
                    nl = nl + 1
                    if( nl > cap ) THROW_HARD('build_leaves: more than 4096 leaves')
                    leaves(nl)%fam   = f
                    leaves(nl)%state = st
                    leaves(nl)%nrec  = nrec
                    leaves(nl)%dir   = L1
                    leaves(nl)%name  = 'f'//int2str(f)//'_'//trim(sfx)//'_s'//int2str_pad(st,2)
                enddo
            enddo
            deallocate(names)
        enddo
    end subroutine build_leaves

    subroutine tree_merge_leaves( run, consfile, leaves, nl, floor, smpd, fscmax, minrec, out, keepl, lab, nkeep )
        character(len=*),     intent(in)    :: run, consfile, out
        type(leaf_t),         intent(in)    :: leaves(:)
        integer,              intent(in)    :: nl, minrec
        real,                 intent(in)    :: floor, smpd, fscmax
        integer, allocatable, intent(inout) :: keepl(:), lab(:)
        integer,              intent(out)   :: nkeep
        type(fft3_t) :: F, Fc
        real,     allocatable :: cons(:,:,:), vin(:,:,:)
        real(dp), allocatable :: lpc(:,:,:), sph(:,:,:), GE(:,:,:,:), GO(:,:,:,:), GM(:,:,:,:), GMd(:,:,:,:), envc(:,:,:)
        real(dp), allocatable :: mtmp(:,:,:), mc(:,:,:), dmap(:,:,:), R(:,:), self_res(:), gn(:), cnt(:), vals(:), cov(:)
        integer,  allocatable :: sh(:,:,:), shc(:,:,:), sl(:), grp(:), gorder(:), dmk(:), pi_(:), pj_(:), ord(:)
        logical,  allocatable :: supp(:,:,:), act(:)
        character(len=PLEN), allocatable :: mpath(:)
        character(len=8192) :: line
        integer  :: n, c, i, k, ii, jj, ib, jb, nmerge, klo, khi, ng, u, nmk, q, n3c(3), npairs
        real(dp) :: fct, sxx, sx, sn, worst, wi, wj, lastw(5)
        call read_mrc(consfile, cons, n3c)
        n = n3c(1)
        c = 80
        call kt_shells(n, sh)
        sph = kt_sphere(n)
        allocate(keepl(nl), cnt(nl))
        nkeep = 0
        do i = 1,nl
            if( .not. leaf_ok(leaves(i)%dir, leaves(i)%state, smpd, fscmax, minrec, sph, sh) ) cycle
            nkeep = nkeep + 1
            keepl(nkeep) = i
            cnt(nkeep)   = real(leaves(i)%nrec, dp)
        enddo
        write(logfhandle,'(A)') int2str(nkeep)//' leaf states (floor '//fr(real(floor,dp),0)//' A)'
        call flush(logfhandle)
        allocate(lab(max(nkeep,1)))
        if( nkeep == 0 ) return
        call fft3_new(F, n)
        call fft3_new(Fc, c)
        fct = 1.d0 / real(c,dp)**3
        allocate(lpc(n,n,n), supp(n,n,n), dmap(n,n,n), mtmp(n,n,n), mc(c,c,c))
        lpc = real(cons, dp)
        call kt_lowpass(lpc, F, smpd, 3, 10.)
        supp = lpc > 0.2d0 * maxval(lpc)
        sn   = real(count(supp), dp)
        sx   = sum(lpc, mask=supp)
        sxx  = sum(lpc**2, mask=supp)
        allocate(sl(c))
        do i = 1,c/2
            sl(i)     = i
            sl(c/2+i) = n - c/2 + i
        enddo
        allocate(shc(c,c,c))
        do k = 1,c
            do jj = 1,c
                do ii = 1,c
                    shc(ii,jj,k) = min(nint(sqrt(real(kt_freq(ii,c)**2 + kt_freq(jj,c)**2 + kt_freq(k,c)**2, dp))), n/2)
                enddo
            enddo
        enddo
        allocate(GE(c,c,c,nkeep), GO(c,c,c,nkeep))
        do ii = 1,nkeep
            i = keepl(ii)
            call dev_crop(trim(leaves(i)%dir)//'/maps/recvol_state'//int2str_pad(leaves(i)%state,2)//'_even.mrc', GE(:,:,:,ii))
            call dev_crop(trim(leaves(i)%dir)//'/maps/recvol_state'//int2str_pad(leaves(i)%state,2)//'_odd.mrc',  GO(:,:,:,ii))
        enddo
        allocate(envc(c,c,c))
        mtmp = merge(1.d0, 0.d0, lpc > 0.1d0 * maxval(lpc))
        call kt_lowpass(mtmp, F, smpd, 1, 15.)
        call crop_real(mtmp, envc)
        envc = min(max(envc, 0.d0), 1.d0)
        allocate(mpath(nkeep), dmk(nkeep))
        nmk = 0
        do ii = 1,nkeep
            dmk(ii) = 0
            do q = 1,nmk
                if( mpath(q) == leaf_mask_path(run, leaves(keepl(ii))%name) ) dmk(ii) = q
            enddo
            if( dmk(ii) == 0 )then
                nmk = nmk + 1
                mpath(nmk) = leaf_mask_path(run, leaves(keepl(ii))%name)
                dmk(ii) = nmk
            endif
        enddo
        allocate(GMd(c,c,c,nmk), GM(c,c,c,nkeep), cov(nkeep))
        do q = 1,nmk
            if( file_exists(trim(mpath(q))) )then
                call read_mrc(trim(mpath(q)), vin, n3c)
                mtmp = real(vin, dp)
                call kt_lowpass(mtmp, F, smpd, 1, 15.)
                call crop_real(mtmp, mc)
                GMd(:,:,:,q) = real(real(min(max(mc, 0.d0), 1.d0)) * real(envc), dp)
            else
                write(logfhandle,'(A)') '  region mask '//trim(mpath(q))//' MISSING: whole envelope'
                GMd(:,:,:,q) = envc
            endif
        enddo
        do ii = 1,nkeep
            GM(:,:,:,ii) = GMd(:,:,:,dmk(ii))
            cov(ii) = sum(GM(:,:,:,ii)) / max(sum(envc), 1.d-30)
        enddo
        deallocate(GMd)
        call fft3_kill(F)
        write(logfhandle,'(A)') 'region masks: '//int2str(nmk)//' distinct, mean coverage of the envelope '//fr(sum(cov) / real(nkeep,dp), 2)
        allocate(self_res(nkeep))
        do ii = 1,nkeep
            self_res(ii) = band_res(GE(:,:,:,ii) * GM(:,:,:,ii), GO(:,:,:,ii) * GM(:,:,:,ii))
        enddo
        write(logfhandle,'(A)') 'own resolution of the leaves: median '//fr(kt_median(self_res),1)//' A, worst '//fr(maxval(self_res),1)//' A'
        call flush(logfhandle)
        khi = nint(real(n,dp) * real(smpd,dp) / real(floor,dp))
        klo = nint(real(n,dp) * real(smpd,dp) / (1.25d0 * real(floor,dp)))
        allocate(R(nkeep,nkeep), source=huge(1.d0))
        npairs = nkeep * (nkeep - 1) / 2
        !$omp parallel do default(shared) private(ii,jj) schedule(dynamic)
        do jj = 2,nkeep
            do ii = 1,jj-1
                R(ii,jj) = dres(ii, jj)
            enddo
        enddo
        !$omp end parallel do
        if( npairs > 0 )then
            allocate(vals(npairs), pi_(npairs), pj_(npairs), ord(npairs))
            k = 0
            do jj = 2,nkeep
                do ii = 1,jj-1
                    k = k + 1
                    vals(k) = R(ii,jj)
                    pi_(k)  = ii
                    pj_(k)  = jj
                enddo
            enddo
            call kt_argsort_desc(-vals, ord)
            do q = 1,min(6, npairs)
                write(logfhandle,'(A)') '  least-resolved difference at the start: '//trim(leaves(keepl(pi_(ord(q))))%name)//' vs '//&
                    &trim(leaves(keepl(pj_(ord(q))))%name)//': band FSC '//fr(vals(ord(q)),2)
            enddo
            write(logfhandle,'(A)') 'pairs at the start: difference band FSC median '//fr(kt_median(vals),2)//', '//int2str(count(vals < 0.143d0))//&
                &' of '//int2str(npairs)//' pairs unresolved (< 0.143) at '//fr(real(floor,dp),0)//' A'
            call flush(logfhandle)
        endif
        allocate(act(nkeep), gn(nkeep), grp(nkeep))
        act = .true.
        gn  = cnt(1:nkeep)
        do ii = 1,nkeep
            grp(ii) = ii
        enddo
        nmerge = 0
        lastw  = 0.d0
        do
            worst = huge(worst)
            ib    = 0
            jb    = 0
            do jj = 2,nkeep
                if( .not. act(jj) ) cycle
                do ii = 1,jj-1
                    if( .not. act(ii) ) cycle
                    if( R(ii,jj) < worst )then
                        worst = R(ii,jj)
                        ib    = ii
                        jb    = jj
                    endif
                enddo
            enddo
            if( ib == 0 ) exit
            if( worst >= 0.143d0 ) exit
            wi = gn(ib)
            wj = gn(jb)
            GE(:,:,:,ib) = (wi * GE(:,:,:,ib) + wj * GE(:,:,:,jb)) / (wi + wj)
            GO(:,:,:,ib) = (wi * GO(:,:,:,ib) + wj * GO(:,:,:,jb)) / (wi + wj)
            GM(:,:,:,ib) = max(GM(:,:,:,ib), GM(:,:,:,jb))
            gn(ib)  = wi + wj
            act(jb) = .false.
            where( grp == jb ) grp = ib
            nmerge = nmerge + 1
            lastw  = [lastw(2:5), worst]
            !$omp parallel do default(shared) private(jj) schedule(dynamic)
            do jj = 1,nkeep
                if( .not. act(jj) .or. jj == ib ) cycle
                if( jj < ib )then
                    R(jj,ib) = dres(jj, ib)
                else
                    R(ib,jj) = dres(ib, jj)
                endif
            enddo
            !$omp end parallel do
        enddo
        ng = count(act)
        allocate(gorder(ng))
        call order_groups
        do ii = 1,nkeep
            do q = 1,ng
                if( gorder(q) == grp(ii) ) lab(ii) = q
            enddo
        enddo
        line = ''
        do q = 6 - min(nmerge,5),5
            line = trim(line)//' '//fr(lastw(q),2)
        enddo
        write(logfhandle,'(A)') int2str(nmerge)//' merges -> '//int2str(ng)//' states; the last merges happened at difference band FSC'//trim(line)
        if( ng > 1 )then
            if( allocated(vals) ) deallocate(vals)
            allocate(vals(ng*(ng-1)/2))
            k = 0
            do jj = 2,nkeep
                if( .not. act(jj) ) cycle
                do ii = 1,jj-1
                    if( .not. act(ii) ) cycle
                    k = k + 1
                    vals(k) = R(ii,jj)
                enddo
            enddo
            write(logfhandle,'(A)') 'remaining pairs: difference band FSC median '//fr(kt_median(vals(1:k)),2)//', least-resolved pair '//&
                &fr(minval(vals(1:k)),2)//' (threshold 0.143 at '//fr(real(floor,dp),0)//' A)'
        else
            write(logfhandle,'(A)') 'one state'
        endif
        open(newunit=u, file=trim(out)//'.txt', status='replace', action='write')
        write(u,'(A)') '# data-derived merge: states whose difference half-maps have FSC < 0.143 in the '//fr(real(floor,dp),0)//&
            &' A band (region mask, whiten 0) are one state; '//int2str(nkeep)//' leaves -> '//int2str(ng)//' states'
        do q = 1,ng
            line = int2str(q)//' '//int2str(nint(sum(cnt(1:nkeep), mask=lab(1:nkeep) == q)))
            do ii = 1,nkeep
                if( lab(ii) == q ) line = trim(line)//' '//trim(leaves(keepl(ii))%name)
            enddo
            write(u,'(A)') trim(line)
        enddo
        close(u)
        open(newunit=u, file=trim(out)//'_labels.txt', status='replace', action='write')
        write(u,'(A)') '# leaf label nrec own_resolution_A'
        do ii = 1,nkeep
            write(u,'(A)') trim(leaves(keepl(ii))%name)//' '//int2str(lab(ii))//' '//int2str(leaves(keepl(ii))%nrec)//' '//fr(self_res(ii),2)
        enddo
        close(u)
        write(logfhandle,'(A)') '-> '//trim(out)//'.txt'
        call fft3_kill(Fc)
        call flush(logfhandle)

      contains

        subroutine dev_crop( fname, outv )
            character(len=*), intent(in)  :: fname
            real(dp),         intent(out) :: outv(:,:,:)
            real(dp) :: sy_, sxy_, det_, a_, b_
            call read_mrc(fname, vin, n3c)
            sy_  = sum(real(vin,dp), mask=supp)
            sxy_ = sum(real(vin,dp) * lpc, mask=supp)
            det_ = sn * sxx - sx * sx
            a_   = (sn * sxy_ - sx * sy_) / det_
            b_   = (sxx * sy_ - sx * sxy_) / det_
            dmap = real(real((real(vin,dp) - b_) / max(a_, 1.d-12) - lpc), dp)
            F%a  = cmplx(dmap * sph, 0.d0, dp)
            call fft3_fwd(F)
            Fc%a = F%a(sl,sl,sl)
            call fft3_bwd(Fc)
            outv = real(real(real(Fc%a, dp) * fct), dp)
        end subroutine dev_crop

        subroutine crop_real( vv, outc )
            real(dp), intent(in)  :: vv(:,:,:)
            real(dp), intent(out) :: outc(:,:,:)
            F%a = cmplx(vv, 0.d0, dp)
            call fft3_fwd(F)
            Fc%a = F%a(sl,sl,sl)
            call fft3_bwd(Fc)
            outc = real(real(real(Fc%a, dp) * fct), dp) * (real(c,dp) / real(n,dp))**3
        end subroutine crop_real

        function band_res( e, o ) result( res )
            real(dp), intent(in) :: e(:,:,:), o(:,:,:)
            real(dp) :: res
            real(dp), allocatable :: num(:), d1(:), d2(:), fsc(:)
            complex(dp), allocatable :: E_(:,:,:)
            integer :: s_, kr, i_, j_, k_
            allocate(num(0:n/2), d1(0:n/2), d2(0:n/2), fsc(0:n/2), source=0.d0)
            Fc%a = cmplx(e, 0.d0, dp)
            call fft3_fwd(Fc)
            allocate(E_(c,c,c), source=Fc%a)
            Fc%a = cmplx(o, 0.d0, dp)
            call fft3_fwd(Fc)
            do k_ = 1,c
                do j_ = 1,c
                    do i_ = 1,c
                        s_ = shc(i_,j_,k_)
                        num(s_) = num(s_) + real(conjg(E_(i_,j_,k_)) * Fc%a(i_,j_,k_), dp)
                        d1(s_)  = d1(s_)  + abs(E_(i_,j_,k_))**2
                        d2(s_)  = d2(s_)  + abs(Fc%a(i_,j_,k_))**2
                    enddo
                enddo
            enddo
            fsc = num / max(sqrt(d1 * d2), 1.d-30)
            kr  = n/2 - 1
            do s_ = 1,n/2-1
                if( fsc(s_) < 0.143d0 )then
                    kr = s_
                    exit
                endif
            enddo
            res = real(n,dp) * real(smpd,dp) / real(kr,dp)
        end function band_res

        function dres( i_, j_ ) result( val )
            integer, intent(in) :: i_, j_
            real(dp) :: val
            type(fft3_t) :: W
            real(dp),    allocatable :: mm(:,:,:), num(:), d1(:), d2(:)
            complex(dp), allocatable :: E_(:,:,:)
            integer :: s_, a_, b_, c_
            allocate(mm(c,c,c))
            allocate(num(klo:khi), d1(klo:khi), d2(klo:khi), source=0.d0)
            call fft3_new(W, c)
            mm  = max(GM(:,:,:,i_), GM(:,:,:,j_))
            W%a = cmplx((GE(:,:,:,i_) - GE(:,:,:,j_)) * mm, 0.d0, dp)
            call fft3_fwd(W)
            allocate(E_(c,c,c), source=W%a)
            W%a = cmplx((GO(:,:,:,i_) - GO(:,:,:,j_)) * mm, 0.d0, dp)
            call fft3_fwd(W)
            do c_ = 1,c
                do b_ = 1,c
                    do a_ = 1,c
                        s_ = shc(a_,b_,c_)
                        if( s_ < klo .or. s_ > khi ) cycle
                        num(s_) = num(s_) + real(conjg(E_(a_,b_,c_)) * W%a(a_,b_,c_), dp)
                        d1(s_)  = d1(s_)  + abs(E_(a_,b_,c_))**2
                        d2(s_)  = d2(s_)  + abs(W%a(a_,b_,c_))**2
                    enddo
                enddo
            enddo
            call fft3_kill(W)
            val = sum(num / max(sqrt(d1 * d2), 1.d-30)) / real(khi - klo + 1, dp)
        end function dres

        subroutine order_groups
            integer :: q_, g_
            real(dp), allocatable :: key(:)
            integer,  allocatable :: ids_(:), ord_(:)
            allocate(key(ng), ids_(ng), ord_(ng))
            q_ = 0
            do g_ = 1,nkeep
                if( .not. act(g_) ) cycle
                q_ = q_ + 1
                ids_(q_) = g_
                key(q_)  = gn(g_)
            enddo
            call kt_argsort_desc(key, ord_)
            do q_ = 1,ng
                gorder(q_) = ids_(ord_(q_))
            enddo
        end subroutine order_groups

    end subroutine tree_merge_leaves

    subroutine tree_consolidate( leaves, keepl, lab, nkeep, nrows, minp, uncap, cut, outdir, K )
        type(leaf_t),     intent(in)  :: leaves(:)
        integer,          intent(in)  :: keepl(:), lab(:), nkeep, nrows, minp
        logical,          intent(in)  :: uncap
        real,             intent(in)  :: cut
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: K
        type(board_t) :: B
        real,     allocatable :: votes(:,:), wrow(:)
        integer,  allocatable :: labrow(:), states(:), sizes(:), remap0(:), lab2(:), order(:), remap(:), c2s(:), st(:)
        real(dp), allocatable :: key(:)
        character(len=PLEN)  :: cur
        character(len=8192)  :: line, gline
        integer :: Kg, ii, i, row, g, u, io, c, nsmall
        Kg = maxval(lab(1:nkeep))
        allocate(votes(nrows,Kg), source=0.)
        allocate(labrow(nrows), wrow(nrows))
        cur = ''
        do ii = 1,nkeep
            i = keepl(ii)
            if( leaves(i)%dir /= cur )then
                cur    = leaves(i)%dir
                labrow = 0
                wrow   = 0.
                if( uncap )then
                    call read_board(trim(cur)//'/../board', B, with_halves=.false.)
                    allocate(c2s(B%ncl), source=0)
                    open(newunit=u, file=trim(cur)//'/clusters.txt', status='old', action='read')
                    do
                        read(u,'(A)',iostat=io) line
                        if( io /= 0 ) exit
                        if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                        call parse_cluster_line(line, c2s)
                    enddo
                    close(u)
                    do row = 1,B%n
                        if( B%am(row) >= 1 .and. B%am(row) <= B%ncl ) labrow(B%pinds(row)) = c2s(B%am(row))
                        wrow(B%pinds(row)) = real(B%maxr(row))
                    enddo
                    deallocate(c2s)
                    call kill_board(B)
                else
                    call read_ints(trim(cur)//'/states_clusters.txt', st)
                    labrow = st
                    wrow   = 1.
                endif
            endif
            do row = 1,nrows
                if( labrow(row) == leaves(i)%state ) votes(row,lab(ii)) = votes(row,lab(ii)) + wrow(row)
            enddo
        enddo
        allocate(states(nrows), source=0)
        do row = 1,nrows
            if( sum(votes(row,:)) > 0. ) states(row) = maxloc(votes(row,:), 1)
        enddo
        allocate(sizes(Kg), source=0)
        do row = 1,nrows
            if( states(row) > 0 ) sizes(states(row)) = sizes(states(row)) + 1
        enddo
        nsmall = count(sizes < minp)
        if( nsmall > 0 )then
            line = ''
            do g = 1,Kg
                if( sizes(g) < minp ) line = trim(line)//' '//int2str(sizes(g))
            enddo
            write(logfhandle,'(A)') 'dropped '//int2str(nsmall)//' merged groups below '//int2str(minp)//' particles:'//trim(line)
        endif
        allocate(remap0(0:Kg), source=0)
        K = 0
        do g = 1,Kg
            if( sizes(g) >= minp )then
                K = K + 1
                remap0(g) = K
            endif
        enddo
        do row = 1,nrows
            states(row) = remap0(states(row))
        enddo
        allocate(lab2(nkeep))
        do ii = 1,nkeep
            lab2(ii) = remap0(lab(ii))
        enddo
        deallocate(sizes)
        allocate(sizes(max(K,1)), source=0)
        do row = 1,nrows
            if( states(row) > 0 ) sizes(states(row)) = sizes(states(row)) + 1
        enddo
        allocate(key(max(K,1)), order(max(K,1)), remap(0:K))
        key   = real(sizes, dp)
        if( K > 0 ) call kt_argsort_desc(key(1:K), order(1:K))
        remap = 0
        do c = 1,K
            remap(order(c)) = c
        enddo
        do row = 1,nrows
            states(row) = remap(states(row))
        enddo
        call simple_mkdir(trim(outdir))
        call write_ints(trim(outdir)//'/states.txt', states)
        open(newunit=u, file=trim(outdir)//'/groups.txt', status='replace', action='write')
        write(u,'(A)') '# merged state (cut '//fr(real(cut,dp),2)//'): particles, members (leaf states)'
        line = ''
        do c = 1,K
            g     = order(c)
            gline = int2str(c)//' '//int2str(sizes(g))
            do ii = 1,nkeep
                if( lab2(ii) == g ) gline = trim(gline)//' '//trim(leaves(keepl(ii))%name)
            enddo
            write(u,'(A)') trim(gline)
            line = trim(line)//' '//int2str(sizes(g))
        enddo
        close(u)
        if( uncap )then
            write(logfhandle,'(A)') 'cut '//fr(real(cut,dp),2)//': '//int2str(K)//' merged states from '//int2str(nkeep)//&
                &' leaf states (ALL particles of each state); particles per merged state:'//trim(line)
        else
            write(logfhandle,'(A)') 'cut '//fr(real(cut,dp),2)//': '//int2str(K)//' merged states from '//int2str(nkeep)//&
                &' leaf states (capped leaf selections); particles per merged state:'//trim(line)
        endif
        call flush(logfhandle)

      contains

        subroutine parse_cluster_line( l, cs )
            character(len=*), intent(in)    :: l
            integer,          intent(inout) :: cs(:)
            integer :: p, q, val, ntok, io2, s
            s    = 0
            ntok = 0
            p    = 1
            do
                do while( p <= len_trim(l) )
                    if( l(p:p) /= ' ' ) exit
                    p = p + 1
                enddo
                if( p > len_trim(l) ) exit
                q = p
                do while( q <= len_trim(l) )
                    if( l(q:q) == ' ' ) exit
                    q = q + 1
                enddo
                ntok = ntok + 1
                read(l(p:q-1),*,iostat=io2) val
                if( ntok == 1 ) s = val
                if( ntok >= 5 )then
                    if( val >= 1 .and. val <= size(cs) ) cs(val) = s
                endif
                p = q
            enddo
        end subroutine parse_cluster_line

    end subroutine tree_consolidate

    subroutine tree_compose( run, cdir, parents_dir, smpd, outdir )
        character(len=*), intent(in) :: run, cdir, parents_dir, outdir
        real,             intent(in) :: smpd
        type(fft3_t) :: F
        real,     allocatable :: vin(:,:,:)
        real(dp), allocatable :: P(:,:,:,:), Pl(:,:,:), Se(:,:,:,:), M(:,:,:,:), mk(:,:,:), sph(:,:,:), wm(:,:,:), lps(:,:,:)
        integer,  allocatable :: sh(:,:,:)
        logical,  allocatable :: supp(:,:,:)
        character(len=PLEN) :: leaf, mf, fe, fo, ff, pfam
        character(len=8192) :: line
        character(len=8)    :: suf(3)
        integer  :: n, n3(3), u, io, s, np_, h, curfam, fam
        real(dp) :: sn, sx, sxx, sy, sxy, det, a, b
        suf = ['_even   ', '_odd    ', '        ']
        call simple_mkdir(trim(outdir))
        curfam = -1
        n      = 0
        pfam   = ''
        write(logfhandle,'(A)') 'state  particles  leaf  region_mask_fraction  FSC_global(state->composed)  FSC_in_region(state->composed)  '//&
            &'FSC_core_outside(state->composed)'
        open(newunit=u, file=trim(cdir)//'/groups.txt', status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            read(line,*) s, np_, leaf
            read(leaf(2:2),*) fam
            mf = leaf_mask_path(run, leaf)
            if( .not. file_exists(trim(mf)) )then
                write(logfhandle,'(A)') int2str(s)//' '//trim(leaf)//': region mask missing ('//trim(mf)//')'
                cycle
            endif
            fe = state_map(s, '_even')
            fo = state_map(s, '_odd')
            ff = state_map(s, '')
            if( .not. file_exists(trim(fe)) .or. .not. file_exists(trim(fo)) .or. .not. file_exists(trim(ff)) ) cycle
            if( fam /= curfam )then
                curfam = fam
                pfam   = trim(parents_dir)//'/fam'//int2str(fam)
                do h = 1,3
                    call read_mrc(trim(pfam)//trim(suf(h))//'.mrc', vin, n3)
                    if( .not. allocated(P) )then
                        n = n3(1)
                        allocate(P(n,n,n,3), Pl(n,n,n), Se(n,n,n,3), M(n,n,n,3), mk(n,n,n), wm(n,n,n), lps(n,n,n), supp(n,n,n))
                        call fft3_new(F, n)
                        call kt_shells(n, sh)
                        sph = kt_sphere(n)
                    endif
                    P(:,:,:,h) = real(vin, dp)
                enddo
                Pl = P(:,:,:,3)
                call kt_lowpass(Pl, F, smpd, 1, 10.)
                Pl   = real(real(Pl), dp)
                supp = Pl > 0.2d0 * maxval(Pl)
                sn   = real(count(supp), dp)
                sx   = sum(Pl, mask=supp)
                sxx  = sum(Pl**2, mask=supp)
            endif
            call read_mrc(trim(mf), vin, n3)
            mk = real(vin, dp)
            call kt_lowpass(mk, F, smpd, 1, 15.)
            mk = min(max(real(real(mk), dp), 0.d0), 1.d0)
            do h = 1,3
                select case(h)
                    case(1)
                        call read_mrc(trim(fe), vin, n3)
                    case(2)
                        call read_mrc(trim(fo), vin, n3)
                    case(3)
                        call read_mrc(trim(ff), vin, n3)
                end select
                Se(:,:,:,h) = real(vin, dp)
                lps = Se(:,:,:,h)
                call kt_lowpass(lps, F, smpd, 1, 10.)
                lps = real(real(lps), dp)
                sy  = sum(lps, mask=supp)
                sxy = sum(lps * Pl, mask=supp)
                det = sn * sxx - sx * sx
                a   = (sn * sxy - sx * sy) / det
                b   = (sxx * sy - sx * sxy) / det
                Se(:,:,:,h) = (Se(:,:,:,h) - b) / max(a, 1.d-12)
                M(:,:,:,h)  = P(:,:,:,h) + mk * (Se(:,:,:,h) - P(:,:,:,h))
            enddo
            call write_mrc(trim(outdir)//'/state_'//int2str_pad(s,2)//'.mrc',      M(:,:,:,3), smpd)
            call write_mrc(trim(outdir)//'/state_'//int2str_pad(s,2)//'_even.mrc', M(:,:,:,1), smpd)
            call write_mrc(trim(outdir)//'/state_'//int2str_pad(s,2)//'_odd.mrc',  M(:,:,:,2), smpd)
            wm   = mk**2 * sph
            line = int2str(s)//'  '//int2str(np_)//'  '//trim(leaf)//' '//fr(sum(mk)/real(size(mk),dp),2)//'   '//&
                &fr(kt_fsc_res(Se(:,:,:,1), Se(:,:,:,2), sph, sh, smpd),1)//' -> '//fr(kt_fsc_res(M(:,:,:,1), M(:,:,:,2), sph, sh, smpd),1)//'   '//&
                &fr(kt_fsc_res(Se(:,:,:,1), Se(:,:,:,2), wm, sh, smpd),1)//' -> '//fr(kt_fsc_res(M(:,:,:,1), M(:,:,:,2), wm, sh, smpd),1)
            wm   = (1.d0 - mk)**2 * sph
            line = trim(line)//'   '//fr(kt_fsc_res(Se(:,:,:,1), Se(:,:,:,2), wm, sh, smpd),1)//' -> '//&
                &fr(kt_fsc_res(M(:,:,:,1), M(:,:,:,2), wm, sh, smpd),1)
            write(logfhandle,'(A)') trim(line)
            call flush(logfhandle)
        enddo
        close(u)
        if( allocated(P) )then
            write(logfhandle,'(A)') 'parent '//trim(pfam)//': global FSC '//fr(kt_fsc_res(P(:,:,:,1), P(:,:,:,2), sph, sh, smpd),1)//' A'
            call fft3_kill(F)
        endif

      contains

        function state_map( s_, sf ) result( p )
            integer,          intent(in) :: s_
            character(len=*), intent(in) :: sf
            character(len=PLEN) :: p, q
            p = trim(cdir)//'/fam0/regionall/level1/maps/recvol_state'//int2str_pad(s_,2)//trim(sf)//'.mrc'
            q = trim(cdir)//'/fam0/regionall/level1/maps/halves/recvol_state'//int2str_pad(s_,2)//trim(sf)//'.mrc'
            if( .not. file_exists(trim(p)) .and. file_exists(trim(q)) ) p = q
        end function state_map

    end subroutine tree_compose

    subroutine tree_arm_readout( consfile, closedfile, files, nf, out, smpd )
        character(len=*),    intent(in) :: consfile, closedfile, out
        character(len=PLEN), intent(in) :: files(:)
        integer,             intent(in) :: nf
        real,                intent(in) :: smpd
        type(fft3_t) :: F
        real,     allocatable :: vin(:,:,:)
        real(dp), allocatable :: cons(:,:,:), closed(:,:,:), d(:,:,:), v(:,:,:), occ(:,:), vals(:)
        logical,  allocatable :: supp(:,:,:), outside(:,:,:), mC(:,:,:), mO(:,:,:)
        character(len=16), allocatable :: calls(:)
        integer,  allocatable :: order(:)
        integer  :: n3(3), n, i, u
        real(dp) :: thr, refC, refO, cmax
        call read_mrc(consfile, vin, n3)
        n = n3(1)
        call fft3_new(F, n)
        allocate(cons(n,n,n), closed(n,n,n), d(n,n,n), v(n,n,n), supp(n,n,n), outside(n,n,n), mC(n,n,n), mO(n,n,n))
        cons = real(vin, dp)
        call kt_lowpass(cons, F, smpd, 2, 10.)
        supp    = cons > 0.2d0  * maxval(cons)
        outside = cons < 0.02d0 * maxval(cons)
        call standardise(cons)
        call read_mrc(closedfile, vin, n3)
        closed = real(vin, dp)
        call kt_lowpass(closed, F, smpd, 2, 10.)
        call standardise(closed)
        d    = closed - cons
        cmax = maxval(closed)
        vals = pack(abs(d), supp .or. closed > 0.2d0 * cmax)
        thr  = 0.5d0 * kt_percentile(vals, 99.d0)
        mC   = d >  thr
        mO   = d < -thr
        refC = sum(closed, mask=mC) / real(max(count(mC),1), dp)
        refO = sum(cons,   mask=mO) / real(max(count(mO),1), dp)
        write(logfhandle,'(A)') 'closed-arm mask '//int2str(count(mC))//' voxels ('//fr(100.d0 * real(count(mC),dp) / real(max(count(supp),1),dp),1)//&
            &'% of the support), open-arm mask '//int2str(count(mO))//' voxels; reference occupancies '//fr(refC,2)//' / '//fr(refO,2)
        allocate(occ(2,max(nf,1)), calls(max(nf,1)))
        do i = 1,nf
            call read_mrc(files(i), vin, n3)
            v = real(vin, dp)
            call kt_lowpass(v, F, smpd, 2, 10.)
            call standardise(v)
            occ(1,i) = sum(v, mask=mC) / real(max(count(mC),1), dp) / refC
            occ(2,i) = sum(v, mask=mO) / real(max(count(mO),1), dp) / refO
            if( occ(1,i) > 0.6d0 .and. occ(2,i) < 0.45d0 )then
                calls(i) = 'closed'
            else if( occ(2,i) > 0.6d0 .and. occ(1,i) < 0.45d0 )then
                calls(i) = 'open'
            else
                calls(i) = 'intermediate'
            endif
        enddo
        call fft3_kill(F)
        open(newunit=u, file=trim(out), status='replace', action='write')
        write(u,'(A)') '# map closed_arm_occupancy open_arm_occupancy call (occupancies relative to the deposited closed map / the consensus)'
        do i = 1,nf
            write(u,'(A)') trim(files(i))//' '//fr(occ(1,i),3)//' '//fr(occ(2,i),3)//' '//trim(calls(i))
        enddo
        close(u)
        if( nf > 0 )then
            allocate(order(nf))
            call kt_argsort_desc(occ(1,1:nf), order)
            write(logfhandle,'(A)') 'arm readout, maps by closed-arm occupancy (first 12):'
            do i = 1,min(12,nf)
                write(logfhandle,'(A)') trim(files(order(i)))//' closed '//fr(occ(1,order(i)),2)//' open '//fr(occ(2,order(i)),2)//' -> '//trim(calls(order(i)))
            enddo
        endif
        write(logfhandle,'(A)') 'summary: '//int2str(count(calls(1:nf) == 'closed'))//' closed, '//int2str(count(calls(1:nf) == 'open'))//' open, '//&
            &int2str(count(calls(1:nf) == 'intermediate'))//' intermediate of '//int2str(nf)//' maps'
        call flush(logfhandle)

      contains

        subroutine standardise( x )
            real(dp), intent(inout) :: x(:,:,:)
            real(dp) :: mo, rms
            mo  = sum(x, mask=outside) / real(max(count(outside),1), dp)
            x   = x - mo
            rms = sqrt(sum(x**2, mask=supp) / real(max(count(supp),1), dp))
            x   = x / max(rms, 1.d-12)
        end subroutine standardise

    end subroutine tree_arm_readout

    subroutine tree_recover( matches, names, nn, cmin, margin )
        character(len=*),    intent(in) :: matches
        character(len=PLEN), intent(in) :: names(:)
        integer,             intent(in) :: nn
        real,                intent(in) :: cmin, margin
        character(len=16384) :: line
        real(dp), allocatable :: C(:,:), bestc(:), bestm(:)
        integer,  allocatable :: ids(:), bestn(:), order(:)
        integer  :: u, io, nd, ns, i, j, p, k, i1, i2
        real(dp) :: c1, c2
        open(newunit=u, file=trim(matches), status='old', action='read')
        read(u,'(A)') line
        p    = index(line, 'conformer') + 9
        line = line(p:)
        nd   = count_tokens(line)
        allocate(ids(nd))
        read(line,*) ids
        ns = 0
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            ns = ns + 1
        enddo
        allocate(C(ns,nd))
        rewind(u)
        read(u,'(A)') line
        i = 0
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1
            read(line,*) k, k, c1, C(i,:)
        enddo
        close(u)
        allocate(bestc(nd), bestm(nd), bestn(nd))
        bestc = -huge(1.d0)
        bestm = 0.d0
        bestn = 0
        do i = 1,ns
            i1 = maxloc(C(i,:), 1)
            c1 = C(i,i1)
            c2 = -huge(1.d0)
            do j = 1,nd
                if( j /= i1 ) c2 = max(c2, C(i,j))
            enddo
            if( c1 >= real(cmin,dp) .and. c1 - c2 >= real(margin,dp) )then
                if( bestn(i1) == 0 .or. c1 > bestc(i1) )then
                    bestc(i1) = c1
                    bestm(i1) = c1 - c2
                    bestn(i1) = i
                endif
            endif
        enddo
        allocate(order(nd))
        do j = 1,nd
            order(j) = j
        enddo
        do i = 2,nd
            k = order(i)
            j = i - 1
            do while( j >= 1 )
                if( ids(order(j)) <= ids(k) ) exit
                order(j+1) = order(j)
                j = j - 1
            enddo
            order(j+1) = k
        enddo
        line = ''
        do j = 1,nd
            if( bestn(order(j)) > 0 ) line = trim(line)//' '//int2str(ids(order(j)))
        enddo
        write(logfhandle,'(A)') 'recovered '//int2str(count(bestn > 0))//' of '//int2str(nd)//' deposited conformers (corr >= '//&
            &fr(real(cmin,dp),2)//', margin >= '//fr(real(margin,dp),2)//'):'//trim(line)
        do j = 1,nd
            i2 = order(j)
            if( bestn(i2) == 0 ) cycle
            if( bestn(i2) <= nn )then
                write(logfhandle,'(A)') '  conformer '//int2str(ids(i2))//': corr '//fr(bestc(i2),2)//' margin '//fr(bestm(i2),2)//'  by '//trim(names(bestn(i2)))
            else
                write(logfhandle,'(A)') '  conformer '//int2str(ids(i2))//': corr '//fr(bestc(i2),2)//' margin '//fr(bestm(i2),2)//'  by state '//int2str(bestn(i2))
            endif
        enddo
        call flush(logfhandle)

      contains

        function count_tokens( l ) result( nt )
            character(len=*), intent(in) :: l
            integer :: nt, q
            logical :: intok
            nt    = 0
            intok = .false.
            do q = 1,len_trim(l)
                if( l(q:q) /= ' ' .and. .not. intok )then
                    nt    = nt + 1
                    intok = .true.
                else if( l(q:q) == ' ' )then
                    intok = .false.
                endif
            enddo
        end function count_tokens

    end subroutine tree_recover

    subroutine tree_gallery( cdir, consfile, outprefix, nrow, lp, smpd )
        character(len=*), intent(in) :: cdir, consfile, outprefix
        integer,          intent(in) :: nrow
        real,             intent(in) :: lp, smpd
        type(fft3_t) :: F
        character(len=PLEN), allocatable :: files(:)
        real,     allocatable :: vin(:,:,:), rgb(:,:,:)
        real(dp), allocatable :: cons(:,:,:), v(:,:,:), proj(:,:,:,:), dev(:,:,:,:), mass(:), nsel(:), corr(:), res(:), allp(:), alld(:)
        real(dp), allocatable :: sph(:,:,:), e(:,:,:), o(:,:,:)
        integer,  allocatable :: order(:), sh(:,:,:), conf(:), snum(:)
        logical,  allocatable :: supp(:,:,:), outside(:,:,:)
        character(len=PLEN) :: L1, cname
        character(len=8192) :: line
        integer  :: nf, n3(3), n, i, j, k, u, io, npage, pg, r0, nr, ip, x, y, st, nc, nsl, cf
        real(dp) :: vmin, vmax, dmax, val, ms, cr
        L1 = trim(cdir)//'/fam0/regionall/level1'
        call list_sorted(trim(L1)//'/states_dir/state_*.mrc', files, nf)
        if( nf < 1 ) return
        allocate(mass(nf), nsel(nf), corr(nf), res(nf), conf(nf), snum(nf))
        mass = 0.d0
        nsel = 0.d0
        corr = 0.d0
        conf = 0
        do i = 1,nf
            j = index(files(i), 'state_', back=.true.) + 6
            read(files(i)(j:j+1),*) snum(i)
        enddo
        if( file_exists(trim(L1)//'/clusters.txt') )then
            open(newunit=u, file=trim(L1)//'/clusters.txt', status='old', action='read')
            do
                read(u,'(A)',iostat=io) line
                if( io /= 0 ) exit
                if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                read(line,*) st, nc, ms, nsl
                do i = 1,nf
                    if( snum(i) == st )then
                        mass(i) = ms
                        nsel(i) = real(nsl, dp)
                    endif
                enddo
            enddo
            close(u)
        endif
        if( file_exists(trim(L1)//'/state_matches.txt') )then
            open(newunit=u, file=trim(L1)//'/state_matches.txt', status='old', action='read')
            do
                read(u,'(A)',iostat=io) line
                if( io /= 0 ) exit
                if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                read(line,*) st, cf, cr
                do i = 1,nf
                    if( snum(i) == st )then
                        conf(i) = cf
                        corr(i) = cr
                    endif
                enddo
            enddo
            close(u)
        endif
        call read_mrc(consfile, vin, n3)
        n = n3(1)
        call fft3_new(F, n)
        call kt_shells(n, sh)
        sph = kt_sphere(n)
        res = -1.d0
        do i = 1,nf
            if( file_exists(trim(L1)//'/maps/recvol_state'//int2str_pad(snum(i),2)//'_even.mrc') .and. &
               &file_exists(trim(L1)//'/maps/recvol_state'//int2str_pad(snum(i),2)//'_odd.mrc') )then
                call read_mrc(trim(L1)//'/maps/recvol_state'//int2str_pad(snum(i),2)//'_even.mrc', vin, n3)
                e = real(vin, dp)
                call read_mrc(trim(L1)//'/maps/recvol_state'//int2str_pad(snum(i),2)//'_odd.mrc', vin, n3)
                o = real(vin, dp)
                res(i) = kt_fsc_res(e, o, sph, sh, smpd)
            endif
        enddo
        allocate(order(nf))
        call kt_argsort_desc(mass, order)
        call read_mrc(consfile, vin, n3)
        cons = real(vin, dp)
        call kt_lowpass(cons, F, smpd, 2, lp)
        allocate(supp(n,n,n), outside(n,n,n))
        supp    = cons > 0.2d0  * maxval(cons)
        outside = cons < 0.02d0 * maxval(cons)
        call standardise(cons)
        allocate(proj(n,n,3,nf), dev(n,n,3,nf), v(n,n,n))
        do i = 1,nf
            call read_mrc(files(order(i)), vin, n3)
            v = real(vin, dp)
            call kt_lowpass(v, F, smpd, 2, lp)
            call standardise(v)
            proj(:,:,1,i) = sum(v, 3)
            proj(:,:,2,i) = sum(v, 2)
            proj(:,:,3,i) = sum(v, 1)
            dev(:,:,1,i)  = sum(v - cons, 3)
            dev(:,:,2,i)  = sum(v - cons, 2)
            dev(:,:,3,i)  = sum(v - cons, 1)
        enddo
        call fft3_kill(F)
        allp = reshape(proj, [size(proj)])
        alld = abs(reshape(dev, [size(dev)]))
        vmin = kt_percentile(allp, 1.d0)
        vmax = kt_percentile(allp, 99.8d0)
        dmax = kt_percentile(alld, 99.5d0)
        npage = (nf + nrow - 1) / nrow
        do pg = 1,npage
            r0 = (pg - 1) * nrow
            nr = min(nrow, nf - r0)
            allocate(rgb(3, 6*n, nr*n))
            do i = 1,nr
                do ip = 1,6
                    do y = 1,n
                        do x = 1,n
                            if( ip <= 3 )then
                                val = (proj(x,y,ip,r0+i) - vmin) / max(vmax - vmin, 1.d-30)
                                rgb(:,(ip-1)*n + x, (i-1)*n + n - y + 1) = real(min(max(val, 0.d0), 1.d0))
                            else
                                val = min(max(dev(x,y,ip-3,r0+i) / max(dmax, 1.d-30), -1.d0), 1.d0)
                                if( val >= 0.d0 )then
                                    rgb(:,(ip-1)*n + x, (i-1)*n + n - y + 1) = real([1.d0, 1.d0 - val, 1.d0 - val])
                                else
                                    rgb(:,(ip-1)*n + x, (i-1)*n + n - y + 1) = real([1.d0 + val, 1.d0 + val, 1.d0])
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
            call write_rgb_jpeg(trim(outprefix)//'_p'//int2str_pad(pg,2)//'.jpg', rgb, 90)
            deallocate(rgb)
        enddo
        open(newunit=u, file=trim(outprefix)//'_table.txt', status='replace', action='write')
        write(u,'(A)') '# rank label mass nrec fsc_A rep own_rep sep best_conformer corr file   (pages: 6 panels per row = projections '//&
            &'of the low-passed map and of its deviation from the consensus)'
        do k = 1,nf
            i = order(k)
            if( conf(i) > 0 )then
                cname = int2str(conf(i))//' '//fr(corr(i),3)
            else
                cname = '? nan'
            endif
            line = int2str(k)//' f0_all_L1_s'//int2str_pad(snum(i),2)//' '//fr(mass(i),0)//' '//int2str(int(nsel(i)))//' '//fr(res(i),1)//&
                &' nan nan nan '//trim(cname)//' '//trim(files(i))
            write(u,'(A)') trim(line)
        enddo
        close(u)
        write(logfhandle,'(A)') 'pages: '//int2str(npage)//' -> '//trim(outprefix)//'_pNN.jpg, table '//trim(outprefix)//'_table.txt'
        call flush(logfhandle)

      contains

        subroutine standardise( xx )
            real(dp), intent(inout) :: xx(:,:,:)
            real(dp) :: mo, rms
            mo  = sum(xx, mask=outside) / real(max(count(outside),1), dp)
            xx  = xx - mo
            rms = sqrt(sum(xx**2, mask=supp) / real(max(count(supp),1), dp))
            xx  = xx / max(rms, 1.d-12)
        end subroutine standardise

    end subroutine tree_gallery

end module simple_klem_tree
