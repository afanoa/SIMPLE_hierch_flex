!@descr: the add-only per-site split of a would-be leaf: at every site of the node (the tree's sites and the islands) a cross-validated two-state classification at fixed poses, then the GAPw discreteness call on the two groups' half maps; labels and statistics only, the project is never written
!! The call (doc/implementation_notes/hetero_tree_site_split_2026_11.md): G1 the first cross-half factor of the two
!! groups' site deviations is signal-dominated (z1 >= 3, r1 >= 0.5 after the nuisance terms), G2 its lobes move
!! density (centroids >= 10 A apart), G3 the even and odd lobe axes agree within 30 deg, and the consensus density
!! between the lobes (the hull minus the lobes dilated by 10 A, >= 10 voxels) is nearly empty: GAPw < 0.28.
module simple_site_split
use simple_defs,          only: dp, logfhandle, PI
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_string_utils,  only: int2str
use simple_fileio,        only: fopen, fclose, fileiochk, nlines
use simple_syslib,        only: file_exists, simple_mkdir
use simple_parameters,    only: parameters
use simple_builder,       only: builder
use simple_cmdline,       only: cmdline
use simple_sp_project,    only: sp_project
use simple_image,         only: image
use simple_eo_partition,  only: ensure_field_halves
use simple_state_mergeup, only: state_mergeup, EM_LIM
implicit none

public :: run_site_split, gap_call, gap_stats, hull_mask
private
#include "simple_local_flags.inc"

integer, parameter :: KSPLIT = 2, NITER = 10, SEED0 = 1, NB_MIN = 10, LOBE_MIN = 4
real,    parameter :: LPG = 15., TAU = 0.5, DIL_A = 10., LSEP_MIN = 10., ANG_MAX = 30., Z1_MIN = 3., R1_MIN = 0.5, GAP_THR = 0.28

type :: gap_stats
    real    :: z1 = 0., r1 = 0., lsep = 0., ang = 0., gapw = -1., gape = -1., gapo = -1.
    integer :: nb = 0, nplus = 0, nminus = 0
    character(len=12) :: call = 'no_call'
    character(len=24) :: fails = ''
end type gap_stats

contains

    !> the worker: cline keys projfile (the node's stage-B project), dir (the node directory with sites/), mskdiam, pgrp,
    !! nthr; writes summary.txt and per site <tag>/labels.txt and stats.txt in the current directory
    subroutine run_site_split( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(sp_project)     :: spproj
        type(state_mergeup)  :: split
        type(string), allocatable :: masks(:), tags(:)
        integer, allocatable :: labels(:), rows(:), out(:)
        complex, allocatable :: ma(:,:,:,:), mb(:,:,:,:)
        real,    allocatable :: e1(:,:,:), e2(:,:,:), o1(:,:,:), o2(:,:,:), sm(:,:,:)
        type(gap_stats) :: st
        type(image)     :: vol
        character(len=:), allocatable :: d, fl
        character(len=8) :: envv
        real    :: smpd
        integer :: box, j, i, funit, ios, n1, n2
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp')    ) call cline%set('hp',    40.)
        if( .not. cline%defined('lp')    ) call cline%set('lp',    12.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('box_crop') )then
            if( .not. cline%defined('projfile') ) THROW_HARD('site_split needs projfile (or box_crop)')
            call spproj%read_segment('stk', cline%get_carg('projfile'))
            box  = spproj%get_box()
            smpd = spproj%get_smpd()
            call spproj%kill
            call cline%set('box_crop', min(box, max(26, 2*ceiling(real(box)*smpd/cline%get_rarg('lp')))))
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        d = params%dir%to_char()
        call site_masks(d, masks, tags)
        call fopen(funit, string('summary.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('site_split: cannot write summary.txt', ios)
        write(funit,'(A)') '# site n1 n2 z1 r1 lsep_A axis_deg nB GAPw GAPw_even GAPw_odd call fails'
        if( size(masks) == 0 )then
            write(logfhandle,'(A)') '>>> SITE_SPLIT no site mask under '//d//'/sites; nothing to test'
            call fclose(funit)
            call build%kill_general_tbox()
            return
        endif
        call get_environment_variable('SIMPLE_SPLIT_CONTRAST', envv, status=ios)
        if( ios /= 0 ) envv = ''
        call split%new(params, build, EM_LIM, contrast=any(trim(envv) == ['yes  ', 'maps ']), &
            &amplitude=any(trim(envv) == ['yes  ', 'estep']))
        rows = split%get_rows()
        do j = 1, size(masks)
            write(logfhandle,'(A)') '>>> SITE_SPLIT site '//tags(j)%to_char()//': '//masks(j)%to_char()
            call split%site_focus(masks(j))
            call split%site_split_em(KSPLIT, NITER, SEED0, labels, ma, mb)
            n1 = count(labels == 1)
            n2 = count(labels == 2)
            call to_real(ma(:,:,:,1), params%box_crop, params%smpd_crop, e1)
            call to_real(ma(:,:,:,2), params%box_crop, params%smpd_crop, e2)
            call to_real(mb(:,:,:,1), params%box_crop, params%smpd_crop, o1)
            call to_real(mb(:,:,:,2), params%box_crop, params%smpd_crop, o2)
            call vol%read_and_crop(masks(j), params%smpd, params%box_crop, params%smpd_crop)
            sm = vol%get_rmat()
            call vol%kill
            call gap_call(e1, e2, o1, o2, real(n1), real(n2), sm > 0.5, params%smpd_crop, st)
            call simple_mkdir(string(tags(j)%to_char()))
            allocate(out(build%spproj_field%get_noris()), source=0)
            out(rows) = labels
            call write_ints(tags(j)%to_char()//'/labels.txt', out)
            deallocate(out)
            call write_stats(tags(j)%to_char()//'/stats.txt', st, n1, n2)
            fl = '-'
            if( len_trim(st%fails) > 0 ) fl = trim(st%fails)
            write(funit,'(A,1X,I0,1X,I0,1X,F7.2,1X,F6.3,1X,F6.1,1X,F6.1,1X,I0,1X,F6.3,1X,F6.3,1X,F6.3,1X,A,1X,A)') &
                &tags(j)%to_char(), n1, n2, st%z1, st%r1, st%lsep, st%ang, st%nb, st%gapw, st%gape, st%gapo, &
                &trim(st%call), fl
            flush(funit)
            write(logfhandle,'(A,F6.2,A,F5.2,A,F5.1,A,F5.1,A,I0,A,F6.3,A)') '>>> SITE_SPLIT site '//tags(j)%to_char()//': groups '//&
                &int2str(n1)//' / '//int2str(n2)//', z1 ', st%z1, ', r1 ', st%r1, ', lobes ', st%lsep, ' A apart, axes ', st%ang, &
                &' deg, B ', st%nb, ' voxels, GAPw ', st%gapw, ' -> '//trim(st%call)
            deallocate(labels, ma, mb, e1, e2, o1, o2, sm)
        end do
        call fclose(funit)
        call split%kill
        call build%kill_general_tbox()
    end subroutine run_site_split

    !> the node's site masks: sites/sites.txt (site_XX.mrc) and sites/islands.txt (island_XX.mrc), tags sN / iN
    subroutine site_masks( d, masks, tags )
        character(len=*), intent(in) :: d
        type(string), allocatable, intent(out) :: masks(:), tags(:)
        character(len=1024) :: line
        integer :: funit, ios, k
        allocate(masks(0), tags(0))
        call read_table(d//'/sites/sites.txt',   'site_',   's')
        call read_table(d//'/sites/islands.txt', 'island_', 'i')

        contains

            subroutine read_table( fname, stem, tag )
                character(len=*), intent(in) :: fname, stem, tag
                if( .not. file_exists(string(fname)) ) return
                call fopen(funit, string(fname), status='old', action='read', iostat=ios)
                call fileiochk('site_split: cannot read '//fname, ios)
                do
                    read(funit,'(A)',iostat=ios) line
                    if( ios /= 0 ) exit
                    if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
                    read(line,*,iostat=ios) k
                    if( ios /= 0 ) cycle
                    if( .not. file_exists(string(d//'/sites/'//stem//pad2(k)//'.mrc')) ) cycle
                    masks = [masks, string(d//'/sites/'//stem//pad2(k)//'.mrc')]
                    tags  = [tags,  string(tag//int2str(k))]
                end do
                call fclose(funit)
            end subroutine read_table

    end subroutine site_masks

    function pad2( i ) result( str )
        integer, intent(in) :: i
        character(len=:), allocatable :: str
        str = int2str(i)
        if( i < 10 ) str = '0'//str
    end function pad2

    !> a reconstructed map (the EM's Fourier layout) to real space, 15 A Gaussian low-passed
    subroutine to_real( cm, box, smpd, rmat )
        complex, intent(in)  :: cm(:,:,:)
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd
        real, allocatable, intent(out) :: rmat(:,:,:)
        type(image) :: img
        complex, allocatable :: c(:,:,:)
        integer :: h, k, l, kk, ll
        real    :: f2
        call img%new([box,box,box], smpd)
        call img%zero_and_flag_ft
        c = img%get_cmat()
        do l = 1, size(c,3)
            ll = l - 1 - merge(box, 0, l > box/2 + 1)
            do k = 1, size(c,2)
                kk = k - 1 - merge(box, 0, k > box/2 + 1)
                do h = 1, size(c,1)
                    f2 = real((h-1)**2 + kk*kk + ll*ll) / real(box*box)
                    c(h,k,l) = cm(h,k,l) * exp(-f2 / smpd**2 * LPG**2 * log(2.))
                end do
            end do
        end do
        call img%set_cmat(c)
        call img%ifft
        rmat = img%get_rmat()
        call img%kill
    end subroutine to_real

    subroutine write_ints( fname, v )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: v(:)
        integer :: funit, ios, i
        call fopen(funit, string(fname), status='replace', action='write', iostat=ios)
        call fileiochk('site_split: cannot write '//fname, ios)
        do i = 1, size(v)
            write(funit,'(I0)') v(i)
        end do
        call fclose(funit)
    end subroutine write_ints

    subroutine write_stats( fname, st, n1, n2 )
        character(len=*), intent(in) :: fname
        type(gap_stats),  intent(in) :: st
        integer,          intent(in) :: n1, n2
        integer :: funit, ios
        call fopen(funit, string(fname), status='replace', action='write', iostat=ios)
        call fileiochk('site_split: cannot write '//fname, ios)
        write(funit,'(A,I0)')   'n1=', n1
        write(funit,'(A,I0)')   'n2=', n2
        write(funit,'(A,F8.3)') 'z1=', st%z1
        write(funit,'(A,F8.3)') 'r1=', st%r1
        write(funit,'(A,F8.2)') 'lobe_half_separation_A=', st%lsep
        write(funit,'(A,F8.2)') 'lobe_axis_angle_deg=', st%ang
        write(funit,'(A,I0)')   'lobe_plus_voxels=', st%nplus
        write(funit,'(A,I0)')   'lobe_minus_voxels=', st%nminus
        write(funit,'(A,I0)')   'between_voxels=', st%nb
        write(funit,'(A,F8.3)') 'GAPw=', st%gapw
        write(funit,'(A,F8.3)') 'GAPw_even=', st%gape
        write(funit,'(A,F8.3)') 'GAPw_odd=', st%gapo
        write(funit,'(A)')      'call='//trim(st%call)
        write(funit,'(A)')      'fails='//trim(st%fails)
        call fclose(funit)
    end subroutine write_stats

    ! ------------------------------------------------------------------ the GAPw call

    !> the two groups' half maps (15 A low-passed, the crop box), their sizes, the site mask and the sampling
    subroutine gap_call( e1, e2, o1, o2, n1, n2, m, smpd, st )
        real,            intent(in)  :: e1(:,:,:), e2(:,:,:), o1(:,:,:), o2(:,:,:), n1, n2, smpd
        logical,         intent(in)  :: m(:,:,:)
        type(gap_stats), intent(out) :: st
        real,    allocatable :: ve(:,:,:), vo(:,:,:), f1e(:,:,:), f1o(:,:,:), w(:,:,:)
        real(dp),allocatable :: fa(:), fb(:)
        logical, allocatable :: lpe(:,:,:), lne(:,:,:), lpo(:,:,:), lno(:,:,:), lpw(:,:,:), lnw(:,:,:), lob(:,:,:), hull(:,:,:), bb(:,:,:)
        real(dp) :: cep(3), cen(3), cop(3), con(3), cwp(3), cwn(3), ae(3), ao(3), neff, r
        integer  :: n, nm, nrd
        real     :: de, dn
        logical  :: g1, g2, g3, gb
        st = gap_stats()
        n  = size(e1,1)
        nm = count(m)
        if( nm < LOBE_MIN )then
            st%fails = 'mask'
            return
        endif
        ve = (n1 * e1 + n2 * e2) / max(n1 + n2, 1.)
        vo = (n1 * o1 + n2 * o2) / max(n1 + n2, 1.)
        ! G1: the groups' site deviation after the nuisance terms (mean, consensus, its Laplacian and rotational blur)
        call nuisance_projected(e1 - e2, ve, m, fa)
        call nuisance_projected(o1 - o2, vo, m, fb)
        neff  = real(nm,dp) * 4._dp / 3._dp * real(PI,dp) * (real(smpd,dp) / real(LPG,dp))**3 / 2._dp
        r     = pearson(fa, fb)
        st%z1 = real(atanh(min(max(r, -0.999999_dp), 0.999999_dp)) * sqrt(max(neff - 3._dp, 1._dp)))
        st%r1 = real(tanh(real(st%z1,dp) / sqrt(max(neff - 3._dp, 1._dp))))
        ! the first factor on either half is the group difference; the cross-half lobes from their geometric mean
        f1e = (e1 - e2) / sqrt(2.)
        f1o = (o1 - o2) / sqrt(2.)
        w   = sign(1., f1e + f1o) * sqrt(max(f1e * f1o, 0.))
        call lobes(m, f1e, lpe, lne)
        call lobes(m, f1o, lpo, lno)
        call lobes(m, w,   lpw, lnw)
        st%nplus  = count(lpw)
        st%nminus = count(lnw)
        if( min(count(lpe), count(lne), count(lpo), count(lno), count(lpw), count(lnw)) < LOBE_MIN )then
            st%fails = merge('lobes   ', 'G1+lobes', st%z1 >= Z1_MIN .and. st%r1 >= R1_MIN)
            return
        endif
        cep = centroid(lpe, f1e); cen = centroid(lne, f1e)
        cop = centroid(lpo, f1o); con = centroid(lno, f1o)
        cwp = centroid(lpw, w);   cwn = centroid(lnw, w)
        ae  = cep - cen
        ao  = cop - con
        st%ang  = real(acos(min(max(dot_product(ae, ao) / max(norm2(ae) * norm2(ao), 1.e-300_dp), -1._dp), 1._dp)) * 180._dp / real(PI,dp))
        st%lsep = real(norm2(cwp - cwn) / 2._dp) * smpd
        ! B: the convex hull of the lobes minus the lobes dilated by DIL_A
        lob = lpw .or. lnw
        call hull_mask(lob, hull)
        nrd = max(1, nint(DIL_A / smpd))
        bb  = hull .and. .not. dilate6(lob, nrd)
        st%nb = count(bb)
        ! GAPw stays -1 (not available) when the weaker lobe lies on solvent, where the ratio means nothing
        if( st%nb >= 3 )then
            de = min(real(sum(ve, mask=lpw)) / real(st%nplus), real(sum(ve, mask=lnw)) / real(st%nminus))
            dn = min(real(sum(vo, mask=lpw)) / real(st%nplus), real(sum(vo, mask=lnw)) / real(st%nminus))
            if( de > 0. .and. dn > 0. )then
                st%gape = min(real(sum(ve, mask=bb) / real(st%nb)) / de, 99.999)
                st%gapo = min(real(sum(vo, mask=bb) / real(st%nb)) / dn, 99.999)
                st%gapw = 0.5 * (st%gape + st%gapo)
            endif
        endif
        g1 = st%z1 >= Z1_MIN .and. st%r1 >= R1_MIN
        g2 = st%lsep >= LSEP_MIN
        g3 = st%ang <= ANG_MAX
        gb = st%nb >= NB_MIN
        st%fails = ''
        if( .not. g1 ) st%fails = trim(st%fails)//'+G1'
        if( .not. g2 ) st%fails = trim(st%fails)//'+G2'
        if( .not. g3 ) st%fails = trim(st%fails)//'+G3'
        if( .not. gb ) st%fails = trim(st%fails)//'+B'
        if( len_trim(st%fails) > 0 ) st%fails = st%fails(2:)
        if( g1 .and. g2 .and. g3 .and. gb )then
            st%call = merge('DISCRETE    ', 'not_discrete', st%gapw >= 0. .and. st%gapw < GAP_THR)
        else
            st%call = 'no_call'
        endif

        contains

            !> the per-half group difference inside the site, the nuisance columns (1, consensus, its Laplacian and rotational blur) projected out
            subroutine nuisance_projected( dmap, v, mk, x )
                real,    intent(in) :: dmap(:,:,:), v(:,:,:)
                logical, intent(in) :: mk(:,:,:)
                real(dp), allocatable, intent(out) :: x(:)
                real(dp), allocatable :: q(:,:)
                real,     allocatable :: lap(:,:,:), rot(:,:,:)
                real(dp) :: nrm
                integer  :: jc, jp
                call blur_terms(v, lap, rot)
                allocate(q(nm,4))
                q(:,1) = 1._dp
                q(:,2) = real(pack(v,   mk), dp)
                q(:,3) = real(pack(lap, mk), dp)
                q(:,4) = real(pack(rot, mk), dp)
                x = real(pack(dmap, mk), dp)
                ! Gram-Schmidt on the four columns, then the projection
                do jc = 1, 4
                    do jp = 1, jc - 1
                        q(:,jc) = q(:,jc) - dot_product(q(:,jp), q(:,jc)) * q(:,jp)
                    end do
                    nrm = norm2(q(:,jc))
                    if( nrm > 1.e-12_dp )then
                        q(:,jc) = q(:,jc) / nrm
                    else
                        q(:,jc) = 0._dp
                    endif
                end do
                do jc = 1, 4
                    x = x - dot_product(q(:,jc), x) * q(:,jc)
                end do
            end subroutine nuisance_projected

            function centroid( mk, f ) result( c )
                logical, intent(in) :: mk(:,:,:)
                real,    intent(in) :: f(:,:,:)
                real(dp) :: c(3), ws, aw
                integer  :: i, j, k
                c  = 0._dp
                ws = 0._dp
                do k = 1, n
                    do j = 1, n
                        do i = 1, n
                            if( .not. mk(i,j,k) ) cycle
                            aw = abs(real(f(i,j,k),dp))
                            c  = c + aw * real([i, j, k] - 1 - n/2, dp)
                            ws = ws + aw
                        end do
                    end do
                end do
                c = c / max(ws, 1.e-300_dp)
            end function centroid

    end subroutine gap_call

    !> the lobes of a factor map inside the site: above TAU of its maximum and below TAU of its minimum
    subroutine lobes( m, f, lp, ln )
        logical, intent(in) :: m(:,:,:)
        real,    intent(in) :: f(:,:,:)
        logical, allocatable, intent(out) :: lp(:,:,:), ln(:,:,:)
        real :: fmax, fmin
        fmax = maxval(f, mask=m)
        fmin = minval(f, mask=m)
        lp = m .and. f > TAU * fmax .and. fmax > 0.
        ln = m .and. f < TAU * fmin .and. fmin < 0.
    end subroutine lobes

    !> the Laplacian of a map and its rotational blur (each angular-momentum generator applied twice, summed)
    subroutine blur_terms( v, lap, rot )
        real, intent(in) :: v(:,:,:)
        real, allocatable, intent(out) :: lap(:,:,:), rot(:,:,:)
        real, allocatable :: t(:,:,:)
        integer :: n, a
        n = size(v,1)
        allocate(lap(n,n,n), rot(n,n,n), source=0.)
        lap(2:n-1,:,:) = v(3:n,:,:) + v(1:n-2,:,:) - 2. * v(2:n-1,:,:)
        lap(:,2:n-1,:) = lap(:,2:n-1,:) + v(:,3:n,:) + v(:,1:n-2,:) - 2. * v(:,2:n-1,:)
        lap(:,:,2:n-1) = lap(:,:,2:n-1) + v(:,:,3:n) + v(:,:,1:n-2) - 2. * v(:,:,2:n-1)
        do a = 1, 3
            t   = generator(a, generator(a, v))
            rot = rot + t
        end do
    end subroutine blur_terms

    !> the a-th angular-momentum generator of f (x, y, z = the voxel indices relative to the box centre)
    function generator( a, f ) result( out )
        integer, intent(in) :: a
        real,    intent(in) :: f(:,:,:)
        real, allocatable :: out(:,:,:), g1(:,:,:), g2(:,:,:), g3(:,:,:)
        integer :: n, i, j, k
        n = size(f,1)
        call gradient(f, g1, g2, g3)
        allocate(out(n,n,n))
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    select case(a)
                        case(1)
                            out(i,j,k) = real(j-1-n/2) * g3(i,j,k) - real(k-1-n/2) * g2(i,j,k)
                        case(2)
                            out(i,j,k) = real(k-1-n/2) * g1(i,j,k) - real(i-1-n/2) * g3(i,j,k)
                        case default
                            out(i,j,k) = real(i-1-n/2) * g2(i,j,k) - real(j-1-n/2) * g1(i,j,k)
                    end select
                end do
            end do
        end do
    end function generator

    !> central differences inside, one-sided at the faces
    subroutine gradient( f, g1, g2, g3 )
        real, intent(in) :: f(:,:,:)
        real, allocatable, intent(out) :: g1(:,:,:), g2(:,:,:), g3(:,:,:)
        integer :: n
        n = size(f,1)
        allocate(g1(n,n,n), g2(n,n,n), g3(n,n,n), source=0.)
        g1(2:n-1,:,:) = 0.5 * (f(3:n,:,:) - f(1:n-2,:,:)); g1(1,:,:) = f(2,:,:) - f(1,:,:); g1(n,:,:) = f(n,:,:) - f(n-1,:,:)
        g2(:,2:n-1,:) = 0.5 * (f(:,3:n,:) - f(:,1:n-2,:)); g2(:,1,:) = f(:,2,:) - f(:,1,:); g2(:,n,:) = f(:,n,:) - f(:,n-1,:)
        g3(:,:,2:n-1) = 0.5 * (f(:,:,3:n) - f(:,:,1:n-2)); g3(:,:,1) = f(:,:,2) - f(:,:,1); g3(:,:,n) = f(:,:,n) - f(:,:,n-1)
    end subroutine gradient

    !> 6-connected dilation by nit voxels
    function dilate6( mk, nit ) result( d )
        logical, intent(in) :: mk(:,:,:)
        integer, intent(in) :: nit
        logical, allocatable :: d(:,:,:), t(:,:,:)
        integer :: it, n
        n = size(mk,1)
        d = mk
        do it = 1, nit
            t = d
            t(2:n,:,:)   = t(2:n,:,:)   .or. d(1:n-1,:,:)
            t(1:n-1,:,:) = t(1:n-1,:,:) .or. d(2:n,:,:)
            t(:,2:n,:)   = t(:,2:n,:)   .or. d(:,1:n-1,:)
            t(:,1:n-1,:) = t(:,1:n-1,:) .or. d(:,2:n,:)
            t(:,:,2:n)   = t(:,:,2:n)   .or. d(:,:,1:n-1)
            t(:,:,1:n-1) = t(:,:,1:n-1) .or. d(:,:,2:n)
            d = t
        end do
    end function dilate6

    real(dp) function pearson( a, b )
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: ma, mb
        ma = sum(a) / real(max(size(a), 1),dp)
        mb = sum(b) / real(max(size(b), 1),dp)
        pearson = sum((a - ma) * (b - mb)) / max(sqrt(sum((a - ma)**2) * sum((b - mb)**2)), 1.e-300_dp)
    end function pearson

    ! ------------------------------------------------------------------ the convex hull of a voxel set

    !> inside = the voxels of the convex hull of the set mk (incremental hull of its boundary voxels); a flat set gives none
    subroutine hull_mask( mk, inside )
        logical, intent(in) :: mk(:,:,:)
        logical, allocatable, intent(out) :: inside(:,:,:)
        real(dp), allocatable :: p(:,:), fn(:,:), fo(:)
        integer,  allocatable :: fv(:,:)
        logical,  allocatable :: fdel(:), vis(:), bnd(:,:,:)
        real(dp) :: cen(3), nrm(3), d, eps, dmax, scale
        integer  :: n, np, nf, i, j, k, q, i0, i1, i2, i3, f, lo(3), hi(3), nvis, ne, a, b, c, cnt
        integer,  allocatable :: edges(:,:)
        n = size(mk,1)
        allocate(inside(n,n,n), source=.false.)
        if( count(mk) < 4 ) return
        ! the boundary voxels carry the hull
        bnd = mk .and. .not. (eoshift(mk, 1, .false., 1) .and. eoshift(mk, -1, .false., 1) .and. eoshift(mk, 1, .false., 2) .and. &
            &eoshift(mk, -1, .false., 2) .and. eoshift(mk, 1, .false., 3) .and. eoshift(mk, -1, .false., 3))
        np = count(bnd)
        allocate(p(3,np))
        q  = 0
        lo = n
        hi = 1
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if( mk(i,j,k) )then
                        lo = min(lo, [i, j, k])
                        hi = max(hi, [i, j, k])
                    endif
                    if( .not. bnd(i,j,k) ) cycle
                    q = q + 1
                    p(:,q) = real([i, j, k], dp)
                end do
            end do
        end do
        scale = real(maxval(hi - lo) + 1, dp)
        eps   = 1.e-9_dp * scale
        ! the initial simplex: extreme in x, farthest from it, farthest from that line, farthest from that plane
        i0 = minloc(p(1,:), dim=1)
        i1 = maxloc([(norm2(p(:,q) - p(:,i0)), q=1,np)], dim=1)
        if( norm2(p(:,i1) - p(:,i0)) < eps ) return
        dmax = -1._dp
        i2   = 0
        do q = 1, np
            d = norm2(cross(p(:,i1) - p(:,i0), p(:,q) - p(:,i0))) / norm2(p(:,i1) - p(:,i0))
            if( d > dmax )then
                dmax = d
                i2   = q
            endif
        end do
        if( dmax < eps ) return
        nrm  = cross(p(:,i1) - p(:,i0), p(:,i2) - p(:,i0))
        nrm  = nrm / norm2(nrm)
        dmax = -1._dp
        i3   = 0
        do q = 1, np
            d = abs(dot_product(nrm, p(:,q) - p(:,i0)))
            if( d > dmax )then
                dmax = d
                i3   = q
            endif
        end do
        if( dmax < eps ) return
        cen = 0.25_dp * (p(:,i0) + p(:,i1) + p(:,i2) + p(:,i3))
        allocate(fv(3,8*np+16), fn(3,8*np+16), fo(8*np+16), fdel(8*np+16))
        nf = 0
        call add_face(i0, i1, i2)
        call add_face(i0, i1, i3)
        call add_face(i0, i2, i3)
        call add_face(i1, i2, i3)
        allocate(vis(size(fdel)), edges(2,6*np+12))
        do q = 1, np
            if( any(q == [i0, i1, i2, i3]) ) cycle
            nvis = 0
            do f = 1, nf
                vis(f) = .false.
                if( fdel(f) ) cycle
                if( dot_product(fn(:,f), p(:,q)) - fo(f) > eps )then
                    vis(f) = .true.
                    nvis   = nvis + 1
                endif
            end do
            if( nvis == 0 ) cycle
            ! the horizon: edges of visible faces not shared by another visible face
            ne = 0
            do f = 1, nf
                if( .not. vis(f) ) cycle
                do c = 1, 3
                    a = fv(c,f)
                    b = fv(mod(c,3)+1,f)
                    cnt = 0
                    do k = 1, nf
                        if( .not. vis(k) .or. k == f ) cycle
                        if( any(fv(:,k) == a) .and. any(fv(:,k) == b) ) cnt = cnt + 1
                    end do
                    if( cnt == 0 )then
                        ne = ne + 1
                        edges(:,ne) = [a, b]
                    endif
                end do
            end do
            do f = 1, nf
                if( vis(f) ) fdel(f) = .true.
            end do
            do c = 1, ne
                call add_face(edges(1,c), edges(2,c), q)
            end do
        end do
        ! the voxels inside every face of the hull
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    inside(i,j,k) = .true.
                    do f = 1, nf
                        if( fdel(f) ) cycle
                        if( dot_product(fn(:,f), real([i, j, k], dp)) - fo(f) > eps )then
                            inside(i,j,k) = .false.
                            exit
                        endif
                    end do
                end do
            end do
        end do

        contains

            function cross( u, v ) result( w )
                real(dp), intent(in) :: u(3), v(3)
                real(dp) :: w(3)
                w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
            end function cross

            !> a face with its normal pointing away from the simplex centre
            subroutine add_face( a1, a2, a3 )
                integer, intent(in) :: a1, a2, a3
                real(dp) :: nn(3), ln
                nn = cross(p(:,a2) - p(:,a1), p(:,a3) - p(:,a1))
                ln = norm2(nn)
                if( ln < 1.e-12_dp ) return
                nn = nn / ln
                if( dot_product(nn, cen - p(:,a1)) > 0._dp ) nn = -nn
                if( nf >= size(fdel) ) return
                nf = nf + 1
                fv(:,nf) = [a1, a2, a3]
                fn(:,nf) = nn
                fo(nf)   = dot_product(nn, p(:,a1))
                fdel(nf) = .false.
            end subroutine add_face

    end subroutine hull_mask

end module simple_site_split
