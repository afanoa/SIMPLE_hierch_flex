!@descr: tree-weighted final maps of a tree's units: per Fourier shell the units' cross-half covariance fitted with the tree (NNLS) or taken as measured (PSD), inside and outside a focus mask where the units differ between halves, and every unit's posterior-mean map
module simple_tree_wiener
use simple_defs,         only: dp, logfhandle, PI
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: fopen, fclose, fileiochk
use simple_syslib,       only: simple_mkdir
use simple_imghead,      only: find_ldim_nptcls, find_img_smpd
use simple_image,        only: image
use simple_linalg,       only: jacobi, matinv
use simple_rnd,          only: seed_rnd_fixed, ran3
implicit none

public :: tree_wiener_maps
private
#include "simple_local_flags.inc"

real,    parameter :: HPFIT = 40., LPFIT = 8., R2MIN = 0.7, LPMASK = 20., MTHR = 0.15, MSIG = 2., MDIL = 3., EDGE = 6.
integer, parameter :: NEST = 4, RND_SEED = 1
character(len=5), parameter :: ESTS(NEST) = [character(len=5) :: 'raw', 'own', 'tree2', 'emp2']

contains

    subroutine tree_wiener_maps( evens, odds, counts, names, ljunk, steps, mskdiam, lp, outdir )
        type(string),     intent(in) :: evens(:), odds(:)
        real(dp),         intent(in) :: counts(:)
        integer,          intent(in) :: names(:), steps(:,:)
        logical,          intent(in) :: ljunk(:)
        real,             intent(in) :: mskdiam, lp
        character(len=*), intent(in) :: outdir
        type(image)           :: vol
        complex,  allocatable :: ef(:,:), of(:,:), ea(:,:), oa(:,:), est(:,:,:,:)
        real,     allocatable :: msk(:,:,:), mv(:,:,:), eb(:,:,:), ob(:,:,:), vv(:,:,:), el(:,:,:), ol(:,:,:), x(:,:,:)
        real,     allocatable :: wg(:), filt(:), fres(:,:)
        real(dp), allocatable :: wl(:), atr(:,:), ash(:,:), ast(:,:), gtr(:,:), gsh(:,:), gst(:,:), r2(:,:)
        integer,  allocatable :: phys(:,:), shell(:), lo(:), hi(:), perm(:), cnt(:), ord(:)
        logical,  allocatable :: inm(:,:,:), core(:,:,:)
        character(len=:), allocatable :: choice
        real(dp) :: r2m(3)
        real     :: smpd, r0, rr, fq, rsh
        integer  :: ldim(3), ifoo, n, nu, nsh, nvox, i, j, k, h, l, s, u, lims(3,2), ie, nfit, funit, ios, c
        nu = size(evens)
        if( nu < 2 ) THROW_HARD('tree-weighted maps need at least two units')
        call find_ldim_nptcls(evens(1), ldim, ifoo)
        n    = ldim(1)
        smpd = find_img_smpd(evens(1))
        nsh  = n / 2
        wl   = counts / sum(counts)
        allocate(msk(n,n,n), eb(n,n,n), ob(n,n,n), vv(n,n,n), source=0.)
        r0 = mskdiam / 2. / smpd
        c  = n / 2 + 1
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    rr = sqrt(real((i-c)**2 + (j-c)**2 + (k-c)**2))
                    if( rr <= r0 )then
                        msk(i,j,k) = 1.
                    else if( rr <= r0 + EDGE )then
                        msk(i,j,k) = 0.5 + 0.5 * cos(PI * (rr - r0) / EDGE)
                    endif
                end do
            end do
        end do
        inm = msk > 0.5
        call vol%new([n,n,n], smpd)
        do u = 1, nu
            el = filtered(read_masked(evens(u)), 1, LPMASK)
            ol = filtered(read_masked(odds(u)),  1, LPMASK)
            eb = eb + real(wl(u)) * el
            ob = ob + real(wl(u)) * ol
            vv = vv + real(wl(u)) * el * ol
        end do
        vv   = filtered(vv - eb * ob, 2, MSIG)
        core = (vv > MTHR * maxval(vv, mask=inm)) .and. inm
        if( any(core) )then
            mv = merge(1., 0., dist_to(core) <= MDIL)
            mv = min(max(filtered(mv, 2, MSIG), 0.), 1.) * msk
        else
            allocate(mv(n,n,n), source=0.)
        endif
        deallocate(eb, ob, vv, el, ol, core)
        lims = vol%loop_lims(2)
        nvox = 0
        allocate(cnt(nsh), source=0)
        do l = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    s = nint(sqrt(real(h*h + k*k + l*l)))
                    if( s < 1 .or. s >= nsh ) cycle
                    nvox   = nvox + 1
                    cnt(s) = cnt(s) + 1
                end do
            end do
        end do
        allocate(lo(nsh), hi(nsh))
        lo(1) = 1
        do s = 1, nsh - 1
            hi(s) = lo(s) + cnt(s) - 1
            if( s < nsh ) lo(s+1) = hi(s) + 1
        end do
        allocate(phys(3,nvox), shell(nvox), source=0)
        allocate(wg(nvox), filt(nvox), source=0.)
        ord = lo
        do l = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    s = nint(sqrt(real(h*h + k*k + l*l)))
                    if( s < 1 .or. s >= nsh ) cycle
                    j          = ord(s)
                    ord(s)     = ord(s) + 1
                    phys(:,j)  = vol%comp_addr_phys(h,k,l)
                    shell(j)   = s
                    wg(j)      = merge(1., 2., h == 0 .or. h == n/2)
                    fq         = sqrt(real(h*h + k*k + l*l)) / (real(n) * smpd)
                    filt(j)    = raised_cosine(fq)
                end do
            end do
        end do
        allocate(ef(nvox,nu), of(nvox,nu), ea(nvox,nu), oa(nvox,nu))
        do u = 1, nu
            x = read_masked(evens(u))
            ef(:,u) = coefs(x)
            ea(:,u) = coefs(x * mv)
            x = read_masked(odds(u))
            of(:,u) = coefs(x)
            oa(:,u) = coefs(x * mv)
        end do
        call designs
        gtr = matmul(transpose(atr), atr)
        gsh = matmul(transpose(ash), ash)
        gst = matmul(transpose(ast), ast)
        allocate(est(nvox,nu,3,3), source=cmplx(0.,0.))
        allocate(r2(3,nsh), source=0._dp)
        !$omp parallel do default(shared) private(s) schedule(dynamic) proc_bind(close)
        do s = 1, nsh - 1
            if( hi(s) < lo(s) ) cycle
            call shell_pass(s)
        end do
        !$omp end parallel do
        nfit = 0
        r2m  = 0._dp
        do s = 1, nsh - 1
            rsh = real(n) * smpd / real(s)
            if( rsh < LPFIT .or. rsh > HPFIT ) cycle
            nfit = nfit + 1
            r2m  = r2m + r2(:,s)
        end do
        r2m    = r2m / real(max(1, nfit), dp)
        choice = merge('tree2', 'emp2 ', r2m(1) >= real(R2MIN,dp))
        choice = trim(choice)
        call simple_mkdir(string(outdir))
        allocate(fres(nu,NEST))
        do ie = 1, NEST
            call simple_mkdir(string(outdir//'/'//trim(ESTS(ie))))
            call simple_mkdir(string(outdir//'/'//trim(ESTS(ie))//'/halves'))
            do u = 1, nu
                call write_unit(ie, u)
            end do
        end do
        call summary
        ! choice.txt goes last: the tree's restart guard reads its presence as the final maps being complete
        call fopen(funit, string(outdir//'/choice.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('tree_wiener_maps: cannot write choice.txt', ios)
        write(funit,'(A)') choice//' '//outdir//'/'//choice
        call fclose(funit)
        call vol%kill

    contains

        function read_masked( fname ) result( y )
            type(string), intent(in) :: fname
            real, allocatable :: y(:,:,:)
            call vol%set_ft(.false.)
            call vol%read(fname)
            y = vol%get_rmat() * msk
        end function read_masked

        function filtered( y, mode, par ) result( z )
            real,    intent(in) :: y(:,:,:), par
            integer, intent(in) :: mode
            real, allocatable :: z(:,:,:)
            complex, allocatable :: cm(:,:,:)
            real    :: f2, fac
            integer :: hh, kk, ll, ky, lz
            call vol%set_ft(.false.)
            call vol%set_rmat(y, .false.)
            call vol%fft
            cm = vol%get_cmat()
            do ll = 1, size(cm,3)
                lz = ll - 1 - merge(n, 0, ll > n/2)
                do kk = 1, size(cm,2)
                    ky = kk - 1 - merge(n, 0, kk > n/2)
                    do hh = 1, size(cm,1)
                        f2 = real((hh-1)**2 + ky*ky + lz*lz) / real(n*n)
                        if( mode == 1 )then
                            fac = exp(-f2 / smpd**2 * par**2 * log(2.))
                        else
                            fac = exp(-2. * PI**2 * par**2 * f2)
                        endif
                        cm(hh,kk,ll) = cm(hh,kk,ll) * fac
                    end do
                end do
            end do
            call vol%set_cmat(cm)
            call vol%ifft
            z = vol%get_rmat()
        end function filtered

        function dist_to( mk ) result( dpx )
            logical, intent(in) :: mk(:,:,:)
            real, allocatable :: dpx(:,:,:)
            call vol%set_ft(.false.)
            call vol%set_rmat(merge(0., 1., mk), .false.)
            call vol%distance_transform
            dpx = vol%get_rmat() / smpd
        end function dist_to

        function coefs( y ) result( cv )
            real, intent(in) :: y(:,:,:)
            complex, allocatable :: cv(:), cm(:,:,:)
            integer :: jj
            call vol%set_ft(.false.)
            call vol%set_rmat(y, .false.)
            call vol%fft
            cm = vol%get_cmat()
            allocate(cv(nvox))
            do jj = 1, nvox
                cv(jj) = cm(phys(1,jj),phys(2,jj),phys(3,jj))
            end do
        end function coefs

        real function raised_cosine( f )
            real, intent(in) :: f
            real :: k1, k2
            raised_cosine = 1.
            if( lp <= 0. ) return
            k1 = 1. / (lp + 1.)
            k2 = 1. / (lp - 0.5)
            if( f >= k2 )then
                raised_cosine = 0.
            else if( f > k1 )then
                raised_cosine = 0.5 + 0.5 * cos(PI * (f - k1) / (k2 - k1))
            endif
        end function raised_cosine

        subroutine designs
            integer, allocatable :: par(:), paths(:,:), npath(:), pshuf(:,:), nshuf(:)
            integer :: nn, ist, root, p, q, r, ii, col, tmp
            logical :: l_shared
            nn   = 2 * nu - 1
            root = steps(3, nu - 1)
            allocate(par(0:nn-1), source=-1)
            do ist = 1, nu - 1
                par(steps(1,ist)) = steps(3,ist)
                par(steps(2,ist)) = steps(3,ist)
            end do
            allocate(paths(nn,nu), npath(nu), source=0)
            do ii = 1, nu
                p = ii - 1
                do while( p /= root .and. p >= 0 )
                    npath(ii) = npath(ii) + 1
                    paths(npath(ii),ii) = p
                    p = par(p)
                end do
            end do
            call seed_rnd_fixed(RND_SEED)
            perm = [(ii, ii=1,nu)]
            do ii = nu, 2, -1
                q = 1 + int(ran3() * real(ii))
                q = min(q, ii)
                tmp = perm(ii)
                perm(ii) = perm(q)
                perm(q)  = tmp
            end do
            pshuf = paths(:,perm)
            nshuf = npath(perm)
            allocate(atr(nu*(nu+1)/2, nn), ash(nu*(nu+1)/2, nn), ast(nu*(nu+1)/2, nn), source=0._dp)
            r = 0
            do p = 1, nu
                do q = p, nu
                    r = r + 1
                    atr(r,1) = 1._dp
                    ash(r,1) = 1._dp
                    ast(r,1) = 1._dp
                    do ii = 1, npath(p)
                        l_shared = any(paths(1:npath(q),q) == paths(ii,p))
                        col = 1 + paths(ii,p) + 1
                        if( paths(ii,p) > root ) col = col - 1
                        if( l_shared ) atr(r,col) = 1._dp
                    end do
                    do ii = 1, nshuf(p)
                        l_shared = any(pshuf(1:nshuf(q),q) == pshuf(ii,p))
                        col = 1 + pshuf(ii,p) + 1
                        if( pshuf(ii,p) > root ) col = col - 1
                        if( l_shared ) ash(r,col) = 1._dp
                    end do
                    if( p == q ) ast(r, 1 + p) = 1._dp
                end do
            end do
        end subroutine designs

        subroutine shell_pass( sh )
            integer, intent(in) :: sh
            complex, allocatable :: e(:,:), o(:,:)
            real(dp) :: ce(nu,nu), nv(nu), cf(nu,nu), d(nu), r2s(3)
            integer  :: a1, a2, aa
            a1 = lo(sh)
            a2 = hi(sh)
            e = ef(a1:a2,:)
            o = of(a1:a2,:)
            call cov(e, o, wg(a1:a2), ce, nv)
            do aa = 1, nu
                d(aa) = max(ce(aa,aa), 0._dp)
                est(a1:a2,aa,1,1) = real(d(aa) / (d(aa) + nv(aa))) * 0.5 * (e(:,aa) + o(:,aa))
                est(a1:a2,aa,2,1) = real(d(aa) / (d(aa) + 2._dp * nv(aa))) * e(:,aa)
                est(a1:a2,aa,3,1) = real(d(aa) / (d(aa) + 2._dp * nv(aa))) * o(:,aa)
            end do
            call treefit(atr, gtr, ce, cf, r2s(1))
            call treefit(ash, gsh, ce, cf, r2s(2))
            call treefit(ast, gst, ce, cf, r2s(3))
            r2(:,sh) = r2s
            e = ea(a1:a2,:)
            o = oa(a1:a2,:)
            call part_pass(e, o, a1, a2, .true.)
            e = ef(a1:a2,:) - ea(a1:a2,:)
            o = of(a1:a2,:) - oa(a1:a2,:)
            call part_pass(e, o, a1, a2, .false.)
        end subroutine shell_pass

        subroutine part_pass( ep, op, a1, a2, l_first )
            complex, intent(in) :: ep(:,:), op(:,:)
            integer, intent(in) :: a1, a2
            logical, intent(in) :: l_first
            complex, allocatable :: xo(:,:), eo(:,:), oo(:,:)
            real(dp) :: ce(nu,nu), nv(nu), cf(nu,nu), r2d
            integer  :: ke
            call cov(ep, op, wg(a1:a2), ce, nv)
            do ke = 2, 3
                if( ke == 2 )then
                    call treefit(atr, gtr, ce, cf, r2d)
                else
                    call psd(ce, cf)
                endif
                call wiener(cf, nv, ep, op, xo, eo, oo)
                if( l_first )then
                    est(a1:a2,:,1,ke) = xo
                    est(a1:a2,:,2,ke) = eo
                    est(a1:a2,:,3,ke) = oo
                else
                    est(a1:a2,:,1,ke) = est(a1:a2,:,1,ke) + xo
                    est(a1:a2,:,2,ke) = est(a1:a2,:,2,ke) + eo
                    est(a1:a2,:,3,ke) = est(a1:a2,:,3,ke) + oo
                endif
            end do
        end subroutine part_pass

        subroutine cov( e, o, w, ce, nv )
            complex,  intent(in)  :: e(:,:), o(:,:)
            real,     intent(in)  :: w(:)
            real(dp), intent(out) :: ce(nu,nu), nv(nu)
            real(dp) :: g(nu,nu), sw
            integer  :: aa, bb
            sw = sum(real(w,dp))
            do bb = 1, nu
                do aa = 1, nu
                    g(aa,bb) = sum(real(w,dp) * (real(real(e(:,aa)),dp) * real(real(o(:,bb)),dp) + &
                        &real(aimag(e(:,aa)),dp) * real(aimag(o(:,bb)),dp)))
                end do
                nv(bb) = max(sum(real(w,dp) * real(abs(e(:,bb) - o(:,bb)),dp)**2) / (4._dp * sw), 1.e-300_dp)
            end do
            ce = (g + transpose(g)) / (2._dp * sw)
        end subroutine cov

        subroutine psd( ce, cf )
            real(dp), intent(in)  :: ce(nu,nu)
            real(dp), intent(out) :: cf(nu,nu)
            real(dp) :: a(nu,nu), v(nu,nu), dv(nu)
            integer  :: nrot, aa
            a = ce
            call jacobi(a, nu, nu, dv, v, nrot)
            cf = 0._dp
            do aa = 1, nu
                cf = cf + max(dv(aa), 0._dp) * spread(v(:,aa), 2, nu) * spread(v(:,aa), 1, nu)
            end do
        end subroutine psd

        subroutine treefit( a, g, ce, cf, r2o )
            real(dp), intent(in)  :: a(:,:), g(:,:), ce(nu,nu)
            real(dp), intent(out) :: cf(nu,nu), r2o
            real(dp), allocatable :: y(:), yf(:), th(:)
            real(dp) :: ym, sres, stot
            integer  :: p, q, r, noff
            allocate(y(size(a,1)))
            r = 0
            do p = 1, nu
                do q = p, nu
                    r    = r + 1
                    y(r) = ce(p,q)
                end do
            end do
            th = nnls(a, g, y)
            yf = matmul(a, th)
            r  = 0
            ym = 0._dp
            noff = 0
            do p = 1, nu
                do q = p, nu
                    r = r + 1
                    cf(p,q) = yf(r)
                    cf(q,p) = yf(r)
                    if( p /= q )then
                        ym   = ym + y(r)
                        noff = noff + 1
                    endif
                end do
            end do
            ym   = ym / real(max(1,noff), dp)
            sres = 0._dp
            stot = 0._dp
            r    = 0
            do p = 1, nu
                do q = p, nu
                    r = r + 1
                    if( p == q ) cycle
                    sres = sres + (y(r) - yf(r))**2
                    stot = stot + (y(r) - ym)**2
                end do
            end do
            r2o = 1._dp - sres / max(stot, 1.e-300_dp)
        end subroutine treefit

        subroutine wiener( cf, nv, e, o, xo, eo, oo )
            real(dp), intent(in) :: cf(nu,nu), nv(nu)
            complex,  intent(in) :: e(:,:), o(:,:)
            complex, allocatable, intent(out) :: xo(:,:), eo(:,:), oo(:,:)
            real(dp) :: m1(nu,nu), inv(nu,nu)
            complex  :: w1(nu,nu), w2(nu,nu)
            integer  :: aa, err
            m1 = cf
            do aa = 1, nu
                m1(aa,aa) = m1(aa,aa) + nv(aa)
            end do
            call matinv(m1, inv, nu, err)
            w1 = cmplx(real(matmul(cf, inv)))
            m1 = cf
            do aa = 1, nu
                m1(aa,aa) = m1(aa,aa) + 2._dp * nv(aa)
            end do
            call matinv(m1, inv, nu, err)
            w2 = cmplx(real(matmul(cf, inv)))
            xo = matmul(0.5 * (e + o), transpose(w1))
            eo = matmul(e, transpose(w2))
            oo = matmul(o, transpose(w2))
        end subroutine wiener

        subroutine write_unit( ie, uu )
            integer, intent(in) :: ie, uu
            complex, allocatable :: cx(:), ce_(:), co(:)
            character(len=:), allocatable :: dir, stem
            if( ie == 1 )then
                ce_ = ef(:,uu)
                co  = of(:,uu)
                cx  = 0.5 * (ce_ + co)
            else
                cx  = est(:,uu,1,ie-1)
                ce_ = est(:,uu,2,ie-1)
                co  = est(:,uu,3,ie-1)
            endif
            dir  = outdir//'/'//trim(ESTS(ie))
            stem = trim(merge('junk_', '     ', ljunk(uu)))//'unit_'//pad2(names(uu))
            call put(cx * filt, dir//'/'//stem//'.mrc')
            call put(ce_, dir//'/halves/'//stem//'_even.mrc')
            call put(co,  dir//'/halves/'//stem//'_odd.mrc')
            fres(uu,ie) = fsc_res(ce_, co)
        end subroutine write_unit

        subroutine put( cv, fn )
            complex,          intent(in) :: cv(:)
            character(len=*), intent(in) :: fn
            complex, allocatable :: cm(:,:,:)
            integer :: jj
            call vol%zero_and_flag_ft
            cm = vol%get_cmat()
            do jj = 1, nvox
                cm(phys(1,jj),phys(2,jj),phys(3,jj)) = cv(jj)
            end do
            call vol%set_cmat(cm)
            call vol%ifft
            call vol%write(string(fn))
        end subroutine put

        real function fsc_res( e1, o1 )
            complex, intent(in) :: e1(:), o1(:)
            real(dp) :: num, de, do_
            integer  :: ss, jj
            fsc_res = real(n) * smpd / real(nsh - 1)
            do ss = 1, nsh - 1
                num = 0._dp
                de  = 0._dp
                do_ = 0._dp
                do jj = lo(ss), hi(ss)
                    num = num + wg(jj) * real(real(e1(jj) * conjg(o1(jj))),dp)
                    de  = de  + wg(jj) * real(abs(e1(jj)),dp)**2
                    do_ = do_ + wg(jj) * real(abs(o1(jj)),dp)**2
                end do
                if( num / max(sqrt(de * do_), 1.e-300_dp) < 0.143_dp )then
                    if( ss > 1 ) fsc_res = real(n) * smpd / real(ss - 1)
                    return
                endif
            end do
        end function fsc_res

        subroutine summary
            character(len=1024) :: line
            character(len=64)   :: piece
            integer, allocatable :: real_units(:)
            integer :: iee
            real_units = pack([(i, i=1,nu)], .not. ljunk)
            call fopen(funit, string(outdir//'/summary.txt'), status='replace', action='write', iostat=ios)
            call fileiochk('tree_wiener_maps: cannot write summary.txt', ios)
            write(line,'(I0,A,I0,A,I0,A,F6.3,A,F5.1,A,F5.1,A)') nu, ' units (', count(ljunk), ' junk), box ', n, ' at ', smpd, &
                &' A/px, low-pass ', lp, ' A, unsharpened; focus mask ', 100. * sum(mv) / sum(msk), ' % of the sphere'
            call emit(trim(line))
            write(line,'(A,F4.0,A,F4.0,A,F6.3,A,F6.3,A,F6.3,A)') 'covariance fit, off-diagonal R^2 over ', HPFIT, '-', LPFIT, &
                &' A: tree ', r2m(1), ', shuffled tree ', r2m(2), ', star ', r2m(3), ' -> deliver '//choice//&
                &' (rule: tree2 when the tree R^2 >= 0.70, else the unstructured emp2)'
            call emit(trim(line))
            line = 'half-map FSC 0.143 resolution, median / worst over the '//int2str(size(real_units))//' non-junk units:'
            do iee = 1, NEST
                if( size(real_units) == 0 ) exit
                write(piece,'(A,A,F5.1,A,F5.1,A)') ' ', trim(ESTS(iee)), median(fres(real_units,iee)), ' / ', &
                    &maxval(fres(real_units,iee)), ' A'
                line = trim(line)//trim(piece)
            end do
            call emit(trim(line))
            call fclose(funit)
        end subroutine summary

        subroutine emit( str )
            character(len=*), intent(in) :: str
            write(funit,'(A)') str
            write(logfhandle,'(A)') '>>> TREE_WIENER '//str
        end subroutine emit

    end subroutine tree_wiener_maps

    function nnls( a, g, y ) result( th )
        real(dp), intent(in) :: a(:,:), g(:,:), y(:)
        real(dp), allocatable :: th(:), b(:), w(:), sp(:), s(:), gp(:,:), gi(:,:)
        integer,  allocatable :: idx(:)
        logical,  allocatable :: pas(:)
        real(dp) :: tol, alpha, rg
        integer  :: nc, it, inner, t, i, np, err
        nc = size(a,2)
        b  = matmul(transpose(a), y)
        allocate(th(nc), s(nc), source=0._dp)
        allocate(pas(nc), source=.false.)
        tol = 1.e-10_dp * max(maxval(abs(b)), 1.e-300_dp)
        do it = 1, 3 * nc
            w = b - matmul(g, th)
            if( all(pas) ) exit
            t = maxloc(w, dim=1, mask=.not. pas)
            if( w(t) <= tol ) exit
            pas(t) = .true.
            do inner = 1, 3 * nc
                idx = pack([(i, i=1,nc)], pas)
                np  = size(idx)
                gp  = g(idx,idx)
                rg  = 1.e-12_dp * max(sum([(gp(i,i), i=1,np)]) / real(np,dp), 1.e-300_dp)
                do i = 1, np
                    gp(i,i) = gp(i,i) + rg
                end do
                allocate(gi(np,np))
                call matinv(gp, gi, np, err)
                sp = matmul(gi, b(idx))
                deallocate(gi)
                s      = 0._dp
                s(idx) = sp
                if( all(sp > 0._dp) )then
                    th = s
                    exit
                endif
                alpha = huge(1._dp)
                do i = 1, np
                    if( sp(i) <= 0._dp ) alpha = min(alpha, th(idx(i)) / max(th(idx(i)) - sp(i), 1.e-300_dp))
                end do
                th = th + alpha * (s - th)
                where( pas .and. th <= 1.e-14_dp * max(maxval(th), 1.e-300_dp) )
                    pas = .false.
                    th  = 0._dp
                end where
                if( .not. any(pas) ) exit
            end do
        end do
    end function nnls

    real function median( v )
        real, intent(in) :: v(:)
        real, allocatable :: srt(:)
        integer :: i, j, nv
        real    :: tmp
        srt = v
        nv  = size(srt)
        do i = 2, nv
            tmp = srt(i)
            j   = i - 1
            do while( j >= 1 )
                if( srt(j) <= tmp ) exit
                srt(j+1) = srt(j)
                j = j - 1
            end do
            srt(j+1) = tmp
        end do
        median = merge(srt((nv+1)/2), 0.5 * (srt(nv/2) + srt(nv/2+1)), mod(nv,2) == 1)
    end function median

    function pad2( i ) result( str )
        integer, intent(in) :: i
        character(len=:), allocatable :: str
        str = int2str(i)
        if( i < 10 ) str = '0'//str
    end function pad2

end module simple_tree_wiener
