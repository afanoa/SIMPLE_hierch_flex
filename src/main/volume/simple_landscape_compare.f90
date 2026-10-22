!@descr: reproducibility of a node's conformational landscape across independent runs: per run the state maps' masked, low-passed deviations and their half-map reliability; across runs each map explained by the other run's maps against its noise ceiling, the matched pairs' cross-run deviation FSC, the principal angles of the top deviation subspaces, the occupancy along the pooled principal components, the extreme maps at both ends of each pooled component matched against half-map noise, and optionally pole rankings grouped across runs
module simple_landscape_compare
use simple_defs,          only: dp, logfhandle
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_string_utils,  only: int2str
use simple_fileio,        only: fopen, fclose, fileiochk, nlines
use simple_imghead,       only: find_ldim_nptcls, find_img_smpd
use simple_image,         only: image
use simple_linalg,        only: jacobi
use simple_srch_sort_loc, only: hpsort
use simple_stat,          only: avg_ranks
use simple_estimate_ssnr, only: get_resolution
implicit none

public :: landscape_compare
private
#include "simple_local_flags.inc"

integer,  parameter :: NSUB = 3, NPC = 2, NQ = 99, LNAME = 7, NTOP(2) = [1000, 2500]
real,     parameter :: SPEAR_GROUP = 0.99
real(dp), parameter :: EIG_TOL = 1.e-9_dp

type :: lrun
    character(len=:), allocatable :: name
    integer                       :: k = 0
    integer,          allocatable :: ks(:), lab(:)
    real(dp),         allocatable :: d(:,:), p(:), rel(:), he(:,:), ho(:,:), rh(:)
    real,             allocatable :: res_own(:)
    type(image),      allocatable :: fe(:), fo(:)
end type lrun

contains

    subroutine landscape_compare( runfile, maskfile, lp, outfile, polefile )
        class(string),           intent(in) :: runfile, maskfile
        real,                    intent(in) :: lp
        character(len=*),        intent(in) :: outfile
        class(string), optional, intent(in) :: polefile
        character(len=4096), allocatable :: specs(:)
        type(lrun),  allocatable :: runs(:)
        type(image)              :: vol
        real(dp),    allocatable :: w(:), va(:,:), vb(:,:), ea(:), eb(:), ball(:,:), pcs(:,:), shares(:), coords(:,:)
        real(dp),    allocatable :: qa(:), qb(:), c1(:), xc(:)
        real,        allocatable :: msk(:,:,:), res(:)
        logical,     allocatable :: sup(:,:,:)
        integer,     allocatable :: ext(:)
        character(len=:), allocatable :: msg
        real(dp) :: wd, spread, noise, tot, e2
        real     :: smpd
        integer  :: fout, ios, nrun, ir, ia, ib, k, j, n, nsup, ntot, off, iend, kb, nrank
        call read_specs(runfile, specs)
        nrun = size(specs)
        if( nrun < 2 ) THROW_HARD('landscape_compare needs at least two runs')
        call fopen(fout, string(outfile), status='replace', action='write', iostat=ios)
        call fileiochk('landscape_compare: cannot write '//outfile, ios)
        allocate(runs(nrun))
        do ir = 1, nrun
            call load_run(specs(ir), runs(ir))
        end do
        do ir = 1, nrun
            call emit(padname(runs(ir)%name)//' '//int2str(runs(ir)%k)//' leaves; leaf fractions'//fmt_list(runs(ir)%p, 2)//&
                &'; half-map reliability of the deviation maps'//fmt_list(runs(ir)%rel, 2))
        end do
        call emit('== (1) each map of run A explained by run B''s maps: fit R^2 / own reliability (1 = within noise), and '//&
            &'best single match / noise ceiling')
        do ia = 1, nrun
            do ib = 1, nrun
                if( ia == ib ) cycle
                msg = '  '//runs(ia)%name//' by '//runs(ib)%name//':'
                do k = 1, runs(ia)%k
                    if( k > 1 ) msg = msg//' '
                    msg = msg//' '//explained(runs(ia), k, runs(ib))
                end do
                call emit(msg)
            end do
        end do
        call emit('== (1b) matched-pair deviation FSC: each map of run A with its best single match in run B, the FSC 0.143 '//&
            &'resolution (A) of the pair''s cross-run even x odd deviations / of the map''s own half-maps')
        do ia = 1, nrun
            do ib = 1, nrun
                if( ia == ib ) cycle
                msg = '  '//runs(ia)%name//' by '//runs(ib)%name//':'
                do k = 1, runs(ia)%k
                    kb = best_match(runs(ia), k, runs(ib))
                    if( k > 1 ) msg = msg//' '
                    msg = msg//' '//resfmt(fsc_res(runs(ia)%fe(k), runs(ia)%fo(k), runs(ib)%fe(kb), runs(ib)%fo(kb)))//&
                        &'/'//resfmt(runs(ia)%res_own(k))
                end do
                call emit(msg)
            end do
        end do
        call emit('== (1c) cross-half match: each map of run A with its best match in run B by the correlation of A''s even with '//&
            &'B''s odd deviation (and odd with even) over both maps'' own even x odd correlation; 1 = within noise, unbiased by '//&
            &'particles the two maps share')
        do ia = 1, nrun
            do ib = 1, nrun
                if( ia == ib ) cycle
                msg = '  '//runs(ia)%name//' by '//runs(ib)%name//':'
                do k = 1, runs(ia)%k
                    if( k > 1 ) msg = msg//' '
                    msg = msg//' '//rfmt(maxval([(xmatch(runs(ia), k, runs(ib), kb), kb=1,runs(ib)%k)]), 2)
                end do
                call emit(msg)
            end do
        end do
        call emit('== (2) principal angles (cosines) between the runs'' top-3 weighted deviation subspaces')
        do ia = 1, nrun - 1
            do ib = ia + 1, nrun
                call subspace(runs(ia), va, ea)
                call subspace(runs(ib), vb, eb)
                call emit('  '//runs(ia)%name//' vs '//runs(ib)%name//':'//fmt_list(cosines(va, vb), 2)//'   (variance shares'//&
                    &fmt_list(ea, 2)//' |'//fmt_list(eb, 2)//')')
            end do
        end do
        call emit('== (3) occupancy on the pooled top-2 PCs (leaf-size weighted): mean / sd / 5-95 % range, and the PC1 '//&
            &'Wasserstein distance')
        ntot = sum(runs(:)%k)
        allocate(ball(nsup,ntot))
        off = 0
        do ir = 1, nrun
            do k = 1, runs(ir)%k
                ball(:,off+k) = runs(ir)%d(:,k) * sqrt(runs(ir)%p(k))
            end do
            off = off + runs(ir)%k
        end do
        call top_vectors(ball, NPC, pcs, shares)
        deallocate(ball)
        do ir = 1, nrun
            coords = matmul(transpose(runs(ir)%d), pcs)
            msg = '  '//padname(runs(ir)%name)//' '
            do j = 1, size(pcs,2)
                if( j > 1 ) msg = msg//' | '
                msg = msg//'PC'//int2str(j)//' '//occ(coords(:,j), runs(ir)%p)
            end do
            call emit(msg)
        end do
        do ia = 1, nrun - 1
            do ib = ia + 1, nrun
                c1     = matmul(transpose(runs(ia)%d), pcs(:,1))
                qa     = quantiles(c1, runs(ia)%p)
                qb     = quantiles(matmul(transpose(runs(ib)%d), pcs(:,1)), runs(ib)%p)
                wd     = sum(abs(qa - qb)) / real(NQ,dp)
                spread = sqrt(sum(runs(ia)%p * c1**2))
                call emit('  '//runs(ia)%name//' vs '//runs(ib)%name//': PC1 Wasserstein '//rfmt(wd, 3)//' (= '//&
                    &rfmt(wd / max(spread, 1.e-12_dp), 2)//' of run '//runs(ia)%name//'''s PC1 spread)')
            end do
        end do
        call emit('== (4) extrema: along each pooled PC that carries more than the pooled noise share of the variance, each '//&
            &'run''s most extreme map at both ends (state, coordinate, fraction); across runs that map''s best cross-half match '//&
            &'in the other run with its rank toward the same end (1 = the other run''s extreme too), deviation match / noise '//&
            &'ceiling, cross-half match, membership Jaccard and cross-run deviation FSC 0.143 (A) against the maps'' own')
        allocate(ext(nrun), xc(nrun))
        noise = 0._dp
        tot   = 0._dp
        do ir = 1, nrun
            do k = 1, runs(ir)%k
                e2    = runs(ir)%p(k) * sum(runs(ir)%d(:,k)**2)
                tot   = tot + e2
                noise = noise + (1._dp - runs(ir)%rel(k)) * e2
            end do
        end do
        noise = noise / max(tot, 1.e-300_dp)
        do j = 1, size(pcs,2)
            call emit('  PC'//int2str(j)//': pooled variance share '//rfmt(shares(j), 2)//', noise share '//rfmt(noise, 2))
            if( shares(j) <= noise ) cycle
            do iend = 1, 2
                do ir = 1, nrun
                    c1 = matmul(transpose(runs(ir)%d), pcs(:,j))
                    if( iend == 1 )then
                        ext(ir) = minloc(c1, 1)
                    else
                        ext(ir) = maxloc(c1, 1)
                    endif
                    xc(ir) = c1(ext(ir))
                end do
                msg = '  PC'//int2str(j)//' '//trim(merge('low ', 'high', iend == 1))//':'
                do ir = 1, nrun
                    if( ir > 1 ) msg = msg//' |'
                    msg = msg//' '//runs(ir)%name//' s'//pad2(runs(ir)%ks(ext(ir)))//' '//rfmt(xc(ir), 3)//' ('//&
                        &rfmt(runs(ir)%p(ext(ir)), 2)//')'
                end do
                call emit(msg)
                do ia = 1, nrun
                    do ib = 1, nrun
                        if( ia == ib ) cycle
                        kb = maxloc([(xmatch(runs(ia), ext(ia), runs(ib), k), k=1,runs(ib)%k)], dim=1)
                        c1 = matmul(transpose(runs(ib)%d), pcs(:,j))
                        if( iend == 1 )then
                            nrank = count(c1 < c1(kb)) + 1
                        else
                            nrank = count(c1 > c1(kb)) + 1
                        endif
                        call emit('    '//runs(ia)%name//'''s extreme in '//runs(ib)%name//': s'//pad2(runs(ib)%ks(kb))//&
                            &' (rank '//int2str(nrank)//', '//rfmt(c1(kb), 3)//'), '//rfmt(match(runs(ia), ext(ia), runs(ib), kb), 2)//&
                            &', cross-half '//rfmt(xmatch(runs(ia), ext(ia), runs(ib), kb), 2)//', Jaccard '//&
                            &member_jaccard(runs(ia), ext(ia), runs(ib), kb)//', FSC 0.143 '//&
                            &resfmt(fsc_res(runs(ia)%fe(ext(ia)), runs(ia)%fo(ext(ia)), runs(ib)%fe(kb), runs(ib)%fo(kb)))//&
                            &' (own '//resfmt(runs(ia)%res_own(ext(ia)))//' / '//resfmt(runs(ib)%res_own(kb))//')')
                    end do
                end do
            end do
        end do
        if( present(polefile) ) call poles
        do ir = 1, nrun
            do k = 1, runs(ir)%k
                call runs(ir)%fe(k)%kill
                call runs(ir)%fo(k)%kill
            end do
        end do
        call vol%kill
        call fclose(fout)

    contains

        subroutine load_run( spec, r )
            character(len=*), intent(in)    :: spec
            type(lrun),       intent(inout) :: r
            real,     allocatable :: e(:,:,:,:), o(:,:,:,:), me(:,:,:), mo(:,:,:)
            real(dp), allocatable :: el(:,:), ol(:,:), f(:,:), de(:,:), do_(:,:), m(:)
            character(len=:), allocatable :: mdir
            integer :: kk, nmax
            r%name = word(spec, 1)
            mdir   = word(spec, 2)
            r%lab  = read_ints(word(spec, 3))
            nmax   = maxval(r%lab)
            r%ks   = pack([(kk, kk=1,nmax)], [(count(r%lab == kk) > 0, kk=1,nmax)])
            r%k    = size(r%ks)
            if( r%k < 2 ) THROW_HARD('landscape_compare: run '//r%name//' has fewer than two states')
            r%p    = [(real(count(r%lab == r%ks(kk)),dp), kk=1,r%k)]
            r%p    = r%p / sum(r%p)
            if( .not. allocated(sup) ) call setup(mapname(mdir, r%ks(1), 'even'))
            allocate(e(n,n,n,r%k), o(n,n,n,r%k), me(n,n,n), mo(n,n,n))
            me = 0.
            mo = 0.
            do kk = 1, r%k
                e(:,:,:,kk) = read_map(mapname(mdir, r%ks(kk), 'even'))
                o(:,:,:,kk) = read_map(mapname(mdir, r%ks(kk), 'odd'))
                me = me + real(r%p(kk)) * e(:,:,:,kk)
                mo = mo + real(r%p(kk)) * o(:,:,:,kk)
            end do
            allocate(r%fe(r%k), r%fo(r%k), r%res_own(r%k), el(nsup,r%k), ol(nsup,r%k))
            do kk = 1, r%k
                call ft_image(msk * (e(:,:,:,kk) - me), r%fe(kk))
                call ft_image(msk * (o(:,:,:,kk) - mo), r%fo(kk))
                r%res_own(kk) = fsc_res(r%fe(kk), r%fo(kk), r%fe(kk), r%fo(kk))
                el(:,kk) = pack(real(lowpass(e(:,:,:,kk)),dp), sup)
                ol(:,kk) = pack(real(lowpass(o(:,:,:,kk)),dp), sup)
            end do
            deallocate(e, o)
            f   = 0.5_dp * (el + ol)
            de  = el
            do_ = ol
            m   = matmul(f, r%p)
            do kk = 1, r%k
                f(:,kk) = f(:,kk) - m
            end do
            m   = matmul(el, r%p)
            do kk = 1, r%k
                de(:,kk) = el(:,kk) - m
            end do
            m   = matmul(ol, r%p)
            do kk = 1, r%k
                do_(:,kk) = ol(:,kk) - m
            end do
            allocate(r%rel(r%k), r%rh(r%k), r%d(nsup,r%k), r%he(nsup,r%k), r%ho(nsup,r%k))
            do kk = 1, r%k
                r%rh(kk)   = sum(w * de(:,kk) * do_(:,kk)) / max(sqrt(sum(w * de(:,kk)**2) * sum(w * do_(:,kk)**2)), 1.e-300_dp)
                r%rel(kk)  = 2._dp * r%rh(kk) / (1._dp + r%rh(kk))
                r%d(:,kk)  = f(:,kk)   * sqrt(w)
                r%he(:,kk) = de(:,kk)  * sqrt(w)
                r%ho(:,kk) = do_(:,kk) * sqrt(w)
            end do
        end subroutine load_run

        subroutine setup( fname )
            character(len=*), intent(in) :: fname
            integer :: ldim(3), ifoo
            call find_ldim_nptcls(string(fname), ldim, ifoo)
            n    = ldim(1)
            smpd = find_img_smpd(string(fname))
            call vol%read_and_crop(maskfile, find_img_smpd(maskfile), n, smpd)
            msk  = min(max(vol%get_rmat(), 0.), 1.)
            call vol%kill
            call vol%new([n,n,n], smpd)
            res  = vol%get_res()
            sup  = msk > 0.
            nsup = count(sup)
            if( nsup == 0 ) THROW_HARD('landscape_compare: empty mask '//maskfile%to_char())
            w    = pack(real(msk,dp), sup)
            w    = w / sum(w)
        end subroutine setup

        function read_map( fname ) result( x )
            character(len=*), intent(in) :: fname
            real, allocatable :: x(:,:,:)
            integer :: ld(3), ifo
            call find_ldim_nptcls(string(fname), ld, ifo)
            if( ld(1) /= n ) THROW_HARD('landscape_compare: all maps must share one box; '//fname)
            call vol%set_ft(.false.)
            call vol%read(string(fname))
            x = vol%get_rmat()
        end function read_map

        subroutine ft_image( x, img )
            real,        intent(in)    :: x(:,:,:)
            type(image), intent(inout) :: img
            call img%new([n,n,n], smpd)
            call img%set_rmat(x, .false.)
            call img%fft
        end subroutine ft_image

        function lowpass( x ) result( y )
            real, intent(in) :: x(:,:,:)
            real, allocatable :: y(:,:,:)
            complex, allocatable :: cm(:,:,:)
            integer :: hh, kk, ll, ky, lz
            real    :: f2
            call vol%set_ft(.false.)
            call vol%set_rmat(x, .false.)
            call vol%fft
            cm = vol%get_cmat()
            do ll = 1, size(cm,3)
                lz = ll - 1 - merge(n, 0, ll > n/2)
                do kk = 1, size(cm,2)
                    ky = kk - 1 - merge(n, 0, kk > n/2)
                    do hh = 1, size(cm,1)
                        f2 = real((hh-1)**2 + ky*ky + lz*lz) / real(n*n) / smpd**2
                        cm(hh,kk,ll) = cm(hh,kk,ll) * exp(-0.5 * f2 * (lp / 1.2)**2)
                    end do
                end do
            end do
            call vol%set_cmat(cm)
            call vol%ifft
            y = vol%get_rmat()
        end function lowpass

        real function fsc_res( e1, o1, e2, o2 )
            type(image), intent(in) :: e1, o1, e2, o2
            real, allocatable :: ca(:), cb(:)
            real :: f05
            allocate(ca(size(res)), cb(size(res)))
            call e1%fsc(o2, ca)
            call o1%fsc(e2, cb)
            call get_resolution(0.5 * (ca + cb), res, f05, fsc_res)
        end function fsc_res

        function explained( ra, kk, rb ) result( str )
            type(lrun), intent(in) :: ra, rb
            integer,    intent(in) :: kk
            character(len=:), allocatable :: str
            real(dp) :: gb(rb%k,rb%k), gy(rb%k), yy, best
            integer  :: jj
            gb   = matmul(transpose(rb%d), rb%d)
            gy   = matmul(transpose(rb%d), ra%d(:,kk))
            yy   = sum(ra%d(:,kk)**2)
            best = -huge(best)
            do jj = 1, rb%k
                best = max(best, match(ra, kk, rb, jj))
            end do
            str = rfmt(sum(gy * pinv_solve(gb, gy)) / max(yy, 1.e-300_dp) / max(ra%rel(kk), 1.e-6_dp), 2)//'/'//rfmt(best, 2)
        end function explained

        integer function best_match( ra, kk, rb )
            type(lrun), intent(in) :: ra, rb
            integer,    intent(in) :: kk
            real(dp) :: cc(rb%k)
            integer  :: jj
            do jj = 1, rb%k
                cc(jj) = match(ra, kk, rb, jj)
            end do
            best_match = maxloc(cc, 1)
        end function best_match

        real(dp) function match( ra, ka, rb, kb_ )
            type(lrun), intent(in) :: ra, rb
            integer,    intent(in) :: ka, kb_
            match = sum(ra%d(:,ka) * rb%d(:,kb_)) / max(sqrt(sum(ra%d(:,ka)**2) * sum(rb%d(:,kb_)**2)), 1.e-300_dp) / &
                &sqrt(max(ra%rel(ka) * rb%rel(kb_), 1.e-6_dp))
        end function match

        real(dp) function xmatch( ra, ka, rb, kb_ )
            type(lrun), intent(in) :: ra, rb
            integer,    intent(in) :: ka, kb_
            real(dp) :: c1, c2
            c1 = sum(ra%he(:,ka) * rb%ho(:,kb_)) / max(sqrt(sum(ra%he(:,ka)**2) * sum(rb%ho(:,kb_)**2)), 1.e-300_dp)
            c2 = sum(ra%ho(:,ka) * rb%he(:,kb_)) / max(sqrt(sum(ra%ho(:,ka)**2) * sum(rb%he(:,kb_)**2)), 1.e-300_dp)
            xmatch = 0.5_dp * (c1 + c2) / sqrt(max(ra%rh(ka) * rb%rh(kb_), 1.e-6_dp))
        end function xmatch

        function member_jaccard( ra, ka, rb, kb_ ) result( str )
            type(lrun), intent(in) :: ra, rb
            integer,    intent(in) :: ka, kb_
            character(len=:), allocatable :: str
            integer :: nboth, neither
            if( size(ra%lab) /= size(rb%lab) )then
                str = 'n/a'
                return
            endif
            nboth   = count(ra%lab == ra%ks(ka) .and. rb%lab == rb%ks(kb_))
            neither = count(ra%lab == ra%ks(ka) .or.  rb%lab == rb%ks(kb_))
            str = rfmt(real(nboth,dp) / real(max(1, neither),dp), 2)
        end function member_jaccard

        subroutine subspace( r, v, ev )
            type(lrun),            intent(in)  :: r
            real(dp), allocatable, intent(out) :: v(:,:), ev(:)
            real(dp), allocatable :: dw(:,:)
            integer :: kk
            allocate(dw(nsup,r%k))
            do kk = 1, r%k
                dw(:,kk) = r%d(:,kk) * sqrt(r%p(kk))
            end do
            call top_vectors(dw, NSUB, v, ev)
        end subroutine subspace

        function occ( x, p ) result( str )
            real(dp), intent(in) :: x(:), p(:)
            character(len=:), allocatable :: str
            real(dp), allocatable :: q(:)
            real(dp) :: mu
            mu  = sum(p * x)
            q   = quantiles(x, p)
            str = rfmt(mu, 3)//' / '//rfmt(sqrt(sum(p * (x - mu)**2)), 3)//' / '//rfmt(q(5), 3)//'..'//rfmt(q(95), 3)
        end function occ

        subroutine poles
            character(len=4096), allocatable :: pspecs(:)
            real,     allocatable :: post(:,:), v(:)
            real(dp), allocatable :: rk(:,:), sp(:,:)
            logical,  allocatable :: m(:)
            integer,  allocatable :: grp(:)
            integer :: np_, ip, jp, gi, gj, it
            call read_specs(polefile, pspecs)
            np_ = size(pspecs)
            if( np_ < 2 ) return
            do ip = 1, np_
                v = read_reals(word(pspecs(ip), 2))
                if( ip == 1 ) allocate(post(size(v),np_))
                if( size(v) /= size(post,1) ) THROW_HARD('landscape_compare: the pole files differ in length')
                post(:,ip) = v
            end do
            m = all(post >= 0., dim=2)
            allocate(rk(count(m),np_), sp(np_,np_))
            do ip = 1, np_
                rk(:,ip) = avg_ranks(pack(post(:,ip), m))
            end do
            do ip = 1, np_
                do jp = 1, np_
                    sp(ip,jp) = pearson(rk(:,ip), rk(:,jp))
                end do
            end do
            grp = [(ip, ip=1,np_)]
            do ip = 1, np_
                do jp = ip + 1, np_
                    if( sp(ip,jp) >= SPEAR_GROUP )then
                        gi = grp(ip)
                        gj = grp(jp)
                        where( grp == gj ) grp = gi
                    endif
                end do
            end do
            call emit('== (5) pole rankings over the '//int2str(count(m))//' particles ranked in every file, grouped at '//&
                &'Spearman >= 0.99; each member against its group''s first: Spearman / top-1000 / top-2500 Jaccard')
            do ip = 1, np_
                if( any(grp(1:ip-1) == grp(ip)) ) cycle
                msg = '  group:'
                do jp = ip, np_
                    if( grp(jp) == grp(ip) ) msg = msg//' '//word(pspecs(jp), 1)
                end do
                call emit(msg)
                do jp = ip + 1, np_
                    if( grp(jp) /= grp(ip) ) cycle
                    msg = '    '//word(pspecs(jp), 1)//' vs '//word(pspecs(ip), 1)//': '//rfmt(sp(ip,jp), 3)
                    do it = 1, size(NTOP)
                        msg = msg//' / '//rfmt(top_jaccard(pack(post(:,jp), m), pack(post(:,ip), m), NTOP(it)), 2)
                    end do
                    call emit(msg)
                end do
            end do
        end subroutine poles

        subroutine emit( str )
            character(len=*), intent(in) :: str
            write(logfhandle,'(A)') str
            write(fout,'(A)') str
        end subroutine emit

    end subroutine landscape_compare

    subroutine top_vectors( x, m, v, shares )
        real(dp),              intent(in)  :: x(:,:)
        integer,               intent(in)  :: m
        real(dp), allocatable, intent(out) :: v(:,:), shares(:)
        real(dp), allocatable :: gm(:,:), u(:,:), lam(:)
        integer,  allocatable :: ord(:)
        integer :: nrot, ii, mm, ng, nv
        gm = matmul(transpose(x), x)
        ng = size(gm,1)
        allocate(u(ng,ng), lam(ng))
        call jacobi(gm, ng, ng, lam, u, nrot)
        lam = max(lam, 0._dp)
        ord = order_desc(lam)
        mm  = min(m, ng)
        nv  = count(lam(ord(1:mm)) > EIG_TOL * lam(ord(1)))
        allocate(v(size(x,1),nv), shares(mm))
        do ii = 1, mm
            shares(ii) = lam(ord(ii)) / max(sum(lam), 1.e-300_dp)
            if( ii <= nv ) v(:,ii) = matmul(x, u(:,ord(ii))) / sqrt(lam(ord(ii)))
        end do
    end subroutine top_vectors

    function cosines( a, b ) result( cs )
        real(dp), intent(in)  :: a(:,:), b(:,:)
        real(dp), allocatable :: cs(:)
        real(dp), allocatable :: mm(:,:), u(:,:), lam(:)
        integer,  allocatable :: ord(:)
        integer :: nrot, ng
        mm = matmul(transpose(a), b)
        mm = matmul(mm, transpose(mm))
        ng = size(mm,1)
        allocate(u(ng,ng), lam(ng))
        call jacobi(mm, ng, ng, lam, u, nrot)
        ord = order_desc(lam)
        cs  = sqrt(max(lam(ord(1:min(size(a,2),size(b,2)))), 0._dp))
    end function cosines

    function pinv_solve( a, b ) result( x )
        real(dp), intent(in)  :: a(:,:), b(:)
        real(dp), allocatable :: x(:)
        real(dp) :: am(size(a,1),size(a,1)), u(size(a,1),size(a,1)), lam(size(a,1)), t(size(a,1))
        integer  :: nrot, ii
        am = a
        call jacobi(am, size(a,1), size(a,1), lam, u, nrot)
        t  = matmul(transpose(u), b)
        do ii = 1, size(lam)
            if( lam(ii) > EIG_TOL * maxval(lam) )then
                t(ii) = t(ii) / lam(ii)
            else
                t(ii) = 0._dp
            endif
        end do
        x = matmul(u, t)
    end function pinv_solve

    function quantiles( x, p ) result( q )
        real(dp), intent(in)  :: x(:), p(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: xs(:), cpos(:)
        integer,  allocatable :: ord(:)
        real(dp) :: lev, acc
        integer  :: ii, jj, nx
        nx  = size(x)
        ord = order_desc(-x)
        xs  = x(ord)
        allocate(cpos(nx), q(NQ))
        acc = 0._dp
        do ii = 1, nx
            acc      = acc + p(ord(ii))
            cpos(ii) = acc - 0.5_dp * p(ord(ii))
        end do
        do jj = 1, NQ
            lev = real(jj,dp) / real(NQ + 1,dp)
            if( lev <= cpos(1) )then
                q(jj) = xs(1)
            else if( lev >= cpos(nx) )then
                q(jj) = xs(nx)
            else
                ii = 1
                do while( cpos(ii+1) < lev )
                    ii = ii + 1
                end do
                q(jj) = xs(ii) + (xs(ii+1) - xs(ii)) * (lev - cpos(ii)) / max(cpos(ii+1) - cpos(ii), 1.e-300_dp)
            endif
        end do
    end function quantiles

    function order_desc( v ) result( ord )
        real(dp), intent(in) :: v(:)
        integer, allocatable :: ord(:)
        real,    allocatable :: key(:)
        integer :: ii
        key = real(-v)
        ord = [(ii, ii=1,size(v))]
        call hpsort(key, ord)
    end function order_desc

    real(dp) function pearson( a, b )
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: ma, mb
        ma = sum(a) / real(size(a),dp)
        mb = sum(b) / real(size(b),dp)
        pearson = sum((a - ma) * (b - mb)) / max(sqrt(sum((a - ma)**2) * sum((b - mb)**2)), 1.e-300_dp)
    end function pearson

    real(dp) function top_jaccard( a, b, ntop )
        real,    intent(in) :: a(:), b(:)
        integer, intent(in) :: ntop
        logical, allocatable :: ma(:), mb(:)
        ma = top_mask(a, ntop)
        mb = top_mask(b, ntop)
        top_jaccard = real(count(ma .and. mb),dp) / real(max(1, count(ma .or. mb)),dp)
    end function top_jaccard

    function top_mask( a, ntop ) result( mk )
        real,    intent(in) :: a(:)
        integer, intent(in) :: ntop
        logical, allocatable :: mk(:)
        real,    allocatable :: key(:)
        integer, allocatable :: ord(:)
        integer :: ii
        key = -a
        ord = [(ii, ii=1,size(a))]
        call hpsort(key, ord)
        allocate(mk(size(a)), source=.false.)
        mk(ord(1:min(ntop,size(a)))) = .true.
    end function top_mask

    subroutine read_specs( fname, specs )
        class(string),                    intent(in)  :: fname
        character(len=4096), allocatable, intent(out) :: specs(:)
        character(len=4096) :: ln
        integer :: fu, io
        call fopen(fu, fname, status='old', action='read', iostat=io)
        call fileiochk('landscape_compare: cannot read '//fname%to_char(), io)
        allocate(specs(0))
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            ln = adjustl(ln)
            if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
            specs = [specs, ln]
        end do
        call fclose(fu)
    end subroutine read_specs

    function read_ints( fname ) result( v )
        character(len=*), intent(in) :: fname
        integer, allocatable :: v(:)
        integer :: fu, io, i
        allocate(v(nlines(string(fname))))
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        call fileiochk('landscape_compare: cannot read '//fname, io)
        do i = 1, size(v)
            read(fu,*) v(i)
        end do
        call fclose(fu)
    end function read_ints

    function read_reals( fname ) result( v )
        character(len=*), intent(in) :: fname
        real, allocatable :: v(:)
        integer :: fu, io, i
        allocate(v(nlines(string(fname))))
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        call fileiochk('landscape_compare: cannot read '//fname, io)
        do i = 1, size(v)
            read(fu,*) v(i)
        end do
        call fclose(fu)
    end function read_reals

    function word( str, iw ) result( wd )
        character(len=*), intent(in) :: str
        integer,          intent(in) :: iw
        character(len=:), allocatable :: wd
        integer :: i0, i1, cnt, pos
        i0 = 1
        i1 = 0
        do cnt = 1, iw
            pos = 0
            if( i1 < len(str) ) pos = verify(str(i1+1:), ' ')
            if( pos == 0 ) THROW_HARD('landscape_compare: missing field '//int2str(iw)//' in: '//trim(str))
            i0  = i1 + pos
            pos = scan(str(i0:), ' ')
            i1  = len(str)
            if( pos > 0 ) i1 = i0 + pos - 2
        end do
        wd = str(i0:i1)
    end function word

    function mapname( mdir, s, half ) result( fname )
        character(len=*), intent(in) :: mdir, half
        integer,          intent(in) :: s
        character(len=:), allocatable :: fname
        fname = mdir//'/recvol_state'//pad2(s)//'_'//half//'_unfil.mrc'
    end function mapname

    function padname( name ) result( str )
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: str
        str = name//repeat(' ', max(0, LNAME - len(name)))
    end function padname

    function fmt_list( v, nd ) result( str )
        real(dp), intent(in) :: v(:)
        integer,  intent(in) :: nd
        character(len=:), allocatable :: str
        integer :: ii
        str = ''
        do ii = 1, size(v)
            str = str//' '//rfmt(v(ii), nd)
        end do
    end function fmt_list

    function rfmt( x, nd ) result( str )
        real(dp), intent(in) :: x
        integer,  intent(in) :: nd
        character(len=:), allocatable :: str
        character(len=32) :: s, fm
        write(fm,'(A,I0,A)') '(F24.', nd, ')'
        write(s,fm) x
        str = trim(adjustl(s))
    end function rfmt

    function resfmt( r ) result( str )
        real, intent(in) :: r
        character(len=:), allocatable :: str
        if( r > 0. )then
            str = rfmt(real(r,dp), 1)
        else
            str = 'none'
        endif
    end function resfmt

    function pad2( i ) result( str )
        integer, intent(in) :: i
        character(len=:), allocatable :: str
        str = int2str(i)
        if( i < 10 ) str = '0'//str
    end function pad2

end module simple_landscape_compare
