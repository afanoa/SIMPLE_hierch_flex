!@descr: the kernel-latent EM tree in Fortran: what the bash/Python chain did around the board. A finished board
!! (klem_latent.txt, klem_states.txt, klem_resp.bin, klem_halves.bin) is read back; families by spectral bisection of
!! the cell-deviation graph; data-derived regions from the cross-half reproducible variability (spatial components and
!! motion-pattern groups) written as full-box soft masks; states = average-linkage clusters of surviving cells with the
!! top particles by responsibility mass; ranking along the node's own main deviation axis; identification of cells
!! against deposited maps (readout only); the landscape count over the leaves. Numerically the twin of
!! expand_families.py, derive_regions.py, derive_regions2.py, node_states.py, rank_by_target.py, deposited_match.py
!! and count_states.py (real-space work at the full box instead of the coarse grid; masks therefore slightly smoother).
module simple_klem_tree
use simple_core_module_api
use simple_image,  only: image
use simple_string, only: string
use simple_linalg, only: jacobi, eigsrt
!$ use omp_lib
implicit none
public :: board_t, read_board, kill_board, tree_families, tree_regions, tree_node_states, tree_rank_pc1, tree_rank_ppca, tree_deposited_match, tree_count_states, tree_kmeans
private
#include "simple_local_flags.inc"

type board_t
    integer :: n = 0, ncl = 0, ngrid = 0, kmax = 0, n3 = 0, n1 = 0, n2 = 0
    integer, allocatable :: pinds(:), am(:), i1(:), i2(:)
    real,    allocatable :: z1(:), z2(:), gm(:), maxr(:), r(:,:), crep(:), w(:)
    complex, allocatable :: D0(:,:,:,:), D1(:,:,:,:)   ! (-kmax:kmax,-kmax:kmax,-kmax:kmax, ngrid) cross-half deviations
end type board_t

contains

    ! ================================================================== board I/O
    subroutine read_board( dir, B, with_halves )
        character(len=*),  intent(in)    :: dir
        type(board_t),     intent(inout) :: B
        logical, optional, intent(in)    :: with_halves
        character(len=1024) :: line
        integer :: u, io, i, ic, ih, isg, n1, n2, m
        real    :: a1, a2, a5, a6, a7, pop, rep
        logical :: lh
        lh = .true.; if( present(with_halves) ) lh = with_halves
        call kill_board(B)
        ! responsibilities (stream): n ncl ngrid n1 n2, pinds, i1, i2, r rows
        open(newunit=u, file=trim(dir)//'/klem_resp.bin', access='stream', form='unformatted', status='old', action='read')
        read(u) B%n, B%ncl, B%ngrid, n1, n2; B%n1 = n1; B%n2 = n2
        allocate(B%pinds(B%n), B%i1(B%ngrid), B%i2(B%ngrid), B%r(B%ncl,B%n))
        read(u) B%pinds; read(u) B%i1; read(u) B%i2
        do i = 1,B%n
            read(u) B%r(:,i)
        enddo
        close(u)
        ! latent: pind z1 z2 argmax maxr gridmass layer2mass
        allocate(B%am(B%n), B%z1(B%n), B%z2(B%n), B%gm(B%n), B%maxr(B%n)); i = 0
        open(newunit=u, file=trim(dir)//'/klem_latent.txt', status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1
            if( i > B%n ) exit
            read(line,*) m, a1, a2, ic, a5, a6, a7
            B%z1(i) = a1; B%z2(i) = a2; B%am(i) = ic; B%maxr(i) = a5; B%gm(i) = a6
            if( m /= B%pinds(i) ) THROW_HARD('klem_latent.txt and klem_resp.bin disagree on the particle order')
        enddo
        close(u)
        if( i /= B%n ) THROW_HARD('klem_latent.txt row count differs from klem_resp.bin')
        ! states: cell population reproducibility isgrid
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
        ! halves (stream): ncl ngrid kmax n3, w(ngrid), then for ih = 1,2; ic = 1..ngrid: the deviation cube (h fastest)
        open(newunit=u, file=trim(dir)//'/klem_halves.bin', access='stream', form='unformatted', status='old', action='read')
        read(u) ic, m, B%kmax, B%n3
        if( ic /= B%ncl .or. m /= B%ngrid ) THROW_HARD('klem_halves.bin header disagrees with klem_resp.bin')
        allocate(B%w(B%ngrid)); read(u) B%w
        allocate(B%D0(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,B%ngrid), B%D1(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,B%ngrid))
        do ih = 1,2
            do ic = 1,B%ngrid
                if( ih == 1 )then
                    read(u) B%D0(:,:,:,ic)
                else
                    read(u) B%D1(:,:,:,ic)
                endif
            enddo
        enddo
        close(u)
    end subroutine read_board

    subroutine kill_board( B )
        type(board_t), intent(inout) :: B
        if( allocated(B%pinds) ) deallocate(B%pinds, B%am, B%i1, B%i2, B%r, B%z1, B%z2, B%gm, B%maxr, B%crep)
        if( allocated(B%w) ) deallocate(B%w, B%D0, B%D1)
        B%n = 0; B%ncl = 0; B%ngrid = 0; B%kmax = 0; B%n3 = 0
    end subroutine kill_board

    !> correlation of two deviation cubes over the shells in [slo, shi] (shell = nint(|k|)), real part
    function cube_corr( A, Bc, kmax, slo, shi ) result( c )
        integer, intent(in) :: kmax, slo, shi
        complex, intent(in) :: A(-kmax:kmax,-kmax:kmax,-kmax:kmax), Bc(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real :: c
        real(dp) :: num, na, nb
        integer  :: h, k, l, s
        num = 0.d0; na = 0.d0; nb = 0.d0
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = -kmax,kmax
                    s = nint(sqrt(real(h*h + k*k + l*l)))
                    if( s < slo .or. s > shi ) cycle
                    num = num + real(conjg(A(h,k,l)) * Bc(h,k,l), dp); na = na + real(abs(A(h,k,l)),dp)**2; nb = nb + real(abs(Bc(h,k,l)),dp)**2
                enddo
            enddo
        enddo
        c = real(num / max(sqrt(na * nb), 1.d-30))
    end function cube_corr

    !> real-space map (full box) of a coarse-grid Fourier cube: the cube is embedded in the box's Fourier space
    subroutine cube_to_realspace( C, kmax, box, smpd, rm )
        integer, intent(in)  :: kmax, box
        complex, intent(in)  :: C(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real,    intent(in)  :: smpd
        real,    intent(out) :: rm(box,box,box)
        type(image) :: img
        integer :: h, k, l, phys(3)
        call img%new([box,box,box], smpd); call img%zero_and_flag_ft()
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = 0,kmax
                    phys = img%comp_addr_phys(h,k,l); call img%set_fcomp([h,k,l], phys, C(h,k,l) / real(box)**3)
                enddo
            enddo
        enddo
        call img%ifft(); call img%get_rmat_sub(rm); call img%kill
    end subroutine cube_to_realspace

    !> coarse-grid Fourier cube of a full-box real-space map (Fourier crop), with an optional soft spherical mask (radius in A)
    subroutine realspace_to_cube( rm, box, smpd, kmax, C, rmask )
        integer, intent(in)  :: box, kmax
        real,    intent(in)  :: rm(box,box,box), smpd
        complex, intent(out) :: C(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real, optional, intent(in) :: rmask
        type(image) :: img
        real, allocatable :: tmp(:,:,:)
        integer :: h, k, l, phys(3), i, j, m
        real    :: rr
        allocate(tmp(box,box,box)); tmp = rm
        if( present(rmask) )then
            do m = 1,box
                do j = 1,box
                    do i = 1,box
                        rr = sqrt(real((i - box/2 - 1)**2 + (j - box/2 - 1)**2 + (m - box/2 - 1)**2)); tmp(i,j,m) = tmp(i,j,m) / (1. + exp((rr - rmask / smpd) / 2.))
                    enddo
                enddo
            enddo
        endif
        call img%new([box,box,box], smpd); call img%set_rmat(tmp, .false.); call img%fft()
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = -kmax,kmax
                    phys = img%comp_addr_phys(h,k,l); C(h,k,l) = img%get_fcomp([h,k,l], phys) * real(box)**3
                enddo
            enddo
        enddo
        call img%kill; deallocate(tmp)
    end subroutine realspace_to_cube

    ! ================================================================== families
    !> spectral bisection of the cell graph (correlation of the half-averaged deviations); fam(ic) in 1..nfam; particles
    !! follow their argmax cell; subsets written as family%d_subset.txt (1-based project rows) + families.txt
    subroutine tree_families( B, nfam, outdir, nmem )
        type(board_t),    intent(in)  :: B
        integer,          intent(in)  :: nfam
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: nmem(nfam)
        complex, allocatable :: Dm(:,:,:)
        real,    allocatable :: S(:,:), W(:,:), Lsym(:,:), V(:,:), d(:), dg(:), emb(:,:), cent(:,:)
        integer, allocatable :: fam(:), lab(:)
        integer :: K, i, j, nrot, u, ic, kk, it, ngrid
        real    :: med
        ngrid = B%ngrid; K = ngrid; allocate(S(K,K), source=0.); allocate(Dm(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax))
        ! correlation of the half-averaged deviations over the coarse grid (all shells >= 1)
        do i = 1,K
            do j = i,K
                Dm = 0.5 * (B%D0(:,:,:,i) + B%D1(:,:,:,i))
                S(i,j) = cube_corr(Dm, 0.5 * (B%D0(:,:,:,j) + B%D1(:,:,:,j)), B%kmax, 1, 2 * B%kmax); S(j,i) = S(i,j)
            enddo
        enddo
        allocate(W(K,K), Lsym(K,K), V(K,K), d(K), dg(K)); W = max(S, 0.)
        do i = 1,K
            W(i,i) = 0.
        enddo
        do i = 1,K
            dg(i) = sum(W(i,:)) + 1.e-9
        enddo
        do i = 1,K
            do j = 1,K
                Lsym(i,j) = merge(1., 0., i == j) - W(i,j) / sqrt(dg(i) * dg(j))
            enddo
        enddo
        nrot = 0; call jacobi(Lsym, K, K, d, V, nrot); call eigsrt(d, V, K, K)   ! descending: the Fiedler vector is column K-1
        allocate(fam(K))
        if( nfam == 2 )then
            med = median_copy(V(:,K-1))
            do i = 1,K
                fam(i) = merge(2, 1, V(i,K-1) > med)
            enddo
        else   ! k-means on the nfam-1 smallest non-trivial eigenvectors
            allocate(emb(K,nfam-1))
            do j = 1,nfam-1
                emb(:,j) = V(:,K-j)
            enddo
            allocate(lab(K), cent(nfam,nfam-1)); call tree_kmeans(emb, K, nfam-1, nfam, 10, lab, cent); fam = lab; deallocate(emb, lab, cent)
        endif
        ! particles follow their argmax cell
        nmem = 0
        do kk = 1,nfam
            open(newunit=u, file=trim(outdir)//'/family'//int2str(kk)//'_subset.txt', status='replace', action='write')
            do i = 1,B%n
                ic = B%am(i)
                if( B%gm(i) > 0.5 .and. ic >= 1 .and. ic <= K )then
                    if( fam(ic) == kk )then
                        write(u,'(I0)') B%pinds(i); nmem(kk) = nmem(kk) + 1
                    endif
                endif
            enddo
            close(u)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F5.2)') '>>> family ', kk, ': ', count(fam == kk), ' cells, ', nmem(kk), ' particles; within-family mean cell similarity ', fam_mean_sim(kk, kk)
        enddo
        if( nfam == 2 ) write(logfhandle,'(A,F5.2)') '>>> between-family mean cell similarity ', fam_mean_sim(1, 2)
        open(newunit=u, file=trim(outdir)//'/families.txt', status='replace', action='write'); write(u,'(A)') '# cell family'
        do i = 1,K
            write(u,'(I0,1X,I0)') i, fam(i)
        enddo
        close(u)
        deallocate(S, W, Lsym, V, d, dg, fam, Dm)

      contains

        function fam_mean_sim( a, b ) result( ms )
            integer, intent(in) :: a, b
            real :: ms
            real(dp) :: acc
            integer  :: cnt, i_, j_
            acc = 0.d0; cnt = 0
            do i_ = 1,K
                do j_ = 1,K
                    if( fam(i_) == a .and. fam(j_) == b )then
                        acc = acc + S(i_,j_); cnt = cnt + 1
                    endif
                enddo
            enddo
            ms = real(acc / max(cnt, 1))
        end function fam_mean_sim

    end subroutine tree_families

    ! ================================================================== regions
    !> data-derived focal regions of a board: variability v = sum_k dev0_k x dev1_k (real space, full box) over the
    !! surviving cells; mode 1 = spatial connected components of (fm > thr), mode 2 = motion-pattern k-means (kgroups)
    !! on the moving voxels' profiles over the cells, merged above cmerge, split into components. Masks region_N.mrc
    !! (largest variability mass first, up to maxreg), all_moving.mrc; dilation/smoothing scaled from coarse voxels.
    subroutine tree_regions( B, box, smpd, mode, thr, minvox_c, dilate_c, maxreg, kgroups, cmerge, minrep, outdir, nreg, within )
        type(board_t),    intent(in)  :: B
        integer,          intent(in)  :: box, mode, minvox_c, dilate_c, maxreg, kgroups
        real,             intent(in)  :: smpd, thr, cmerge, minrep
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: nreg
        character(len=*), optional, intent(in) :: within   ! inherited full-box mask: only moving voxels inside it
        real,    allocatable :: v(:,:,:), fm(:,:,:), r0(:,:,:), r1(:,:,:), prof(:,:), cent(:,:), soft(:,:,:), union(:,:,:), devm(:,:,:,:), rw(:,:,:)
        logical, allocatable :: mov(:,:,:), m(:,:,:), comp(:,:,:)
        integer, allocatable :: labc(:,:,:), lab(:), parent(:), idx(:,:), order(:)
        real(dp), allocatable :: mass(:)
        integer :: nkeep, k, i, j, l, nl, c, nmov, g, ng, u, nv, minvox, dilate_it, cnt, ngs
        real    :: scalef, p99, sig_c, smpd_c, rad, radmax, cx, cy, cz, sm
        type(image) :: img
        type(string) :: fn
        scalef = real(box) / real(B%n3); smpd_c = real(box) * smpd / real(B%n3)
        minvox = nint(real(minvox_c) * scalef**3); dilate_it = max(1, nint(real(dilate_c) * scalef)); sig_c = 1.5 * scalef
        allocate(v(box,box,box), source=0.); allocate(r0(box,box,box), r1(box,box,box)); nkeep = 0
        if( mode == 2 ) allocate(devm(box,box,box,B%ngrid))
        do k = 1,B%ngrid
            if( B%crep(k) < minrep ) cycle
            nkeep = nkeep + 1
            call cube_to_realspace(B%D0(:,:,:,k), B%kmax, box, smpd, r0); call cube_to_realspace(B%D1(:,:,:,k), B%kmax, box, smpd, r1)
            v = v + r0 * r1
            if( mode == 2 ) devm(:,:,:,nkeep) = 0.5 * (r0 + r1)
        enddo
        deallocate(r0, r1)
        v = max(v, 0.); call gauss3d(v, box, sig_c); p99 = percentile3d(v, box, 99.); allocate(fm(box,box,box)); fm = sqrt(min(max(v / max(p99, 1.e-30), 0.), 1.))
        allocate(mov(box,box,box)); mov = fm > thr
        if( present(within) )then
            allocate(rw(box,box,box)); fn = trim(within); call img%new([box,box,box], smpd); call img%read(fn); call img%get_rmat_sub(rw); call img%kill
            mov = mov .and. (rw > 0.5); deallocate(rw)
            write(logfhandle,'(A,I0,A)') '>>> restricted to the inherited region: ', count(mov), ' moving voxels inside it'
        endif
        nmov = count(mov); allocate(labc(box,box,box), source=0); allocate(union(box,box,box), source=0.); call simple_mkdir(trim(outdir))
        if( mode == 1 )then
            allocate(m(box,box,box)); m = mov; call binary_open(m, box); call label_cc(m, box, labc, nl)
        else
            ! motion-pattern k-means on the moving voxels' profiles over the cells
            allocate(idx(3,nmov)); cnt = 0
            do l = 1,box
                do j = 1,box
                    do i = 1,box
                        if( mov(i,j,l) )then
                            cnt = cnt + 1; idx(:,cnt) = [i,j,l]
                        endif
                    enddo
                enddo
            enddo
            allocate(prof(nmov,nkeep))
            do c = 1,nmov
                prof(c,:) = devm(idx(1,c),idx(2,c),idx(3,c),1:nkeep); prof(c,:) = prof(c,:) - sum(prof(c,:)) / real(nkeep)
                sm = sqrt(sum(prof(c,:)**2)); prof(c,:) = prof(c,:) / max(sm, 1.e-30)
            enddo
            allocate(lab(nmov), cent(kgroups,nkeep)); call tree_kmeans(prof, nmov, nkeep, kgroups, 12, lab, cent)
            ! merge groups whose centroids correlate above cmerge
            allocate(parent(kgroups))
            do g = 1,kgroups
                parent(g) = g; sm = sqrt(sum(cent(g,:)**2)); cent(g,:) = cent(g,:) / max(sm, 1.e-30)
            enddo
            do i = 1,kgroups
                do j = i+1,kgroups
                    if( sum(cent(i,:) * cent(j,:)) > cmerge ) parent(j) = parent(i)
                enddo
            enddo
            ngs = 0
            do g = 1,kgroups
                if( parent(g) == g ) ngs = ngs + 1
            enddo
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> motion-pattern k-means: ', nmov, ' moving voxels, k = ', kgroups, ' -> ', ngs, ' groups after merging (', nkeep, ' cells)'
            ! connected components within each group, labelled consecutively
            allocate(m(box,box,box)); nl = 0
            do g = 1,kgroups
                if( parent(g) /= g ) cycle
                m = .false.
                do c = 1,nmov
                    if( parent(lab(c)) == g ) m(idx(1,c),idx(2,c),idx(3,c)) = .true.
                enddo
                call label_cc_offset(m, box, labc, nl)
            enddo
            deallocate(prof, lab, cent, parent, idx, devm)
        endif
        ! rank components by variability mass, keep >= minvox
        allocate(mass(max(nl,1)), source=0.d0)
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    if( labc(i,j,l) > 0 ) mass(labc(i,j,l)) = mass(labc(i,j,l)) + v(i,j,l)
                enddo
            enddo
        enddo
        allocate(order(max(nl,1))); call argsort_desc_dp(mass, nl, order); nreg = 0
        allocate(comp(box,box,box), soft(box,box,box))
        write(logfhandle,'(A,I0,A,I0,A,F4.2,A,I0,A,I0,A)') '>>> variability: ', nkeep, ' of ', B%ngrid, ' cells survive; above ', thr, ': ', nmov, ' voxels in ', nl, ' components'
        do c = 1,nl
            k = order(c); nv = count(labc == k)
            if( nv < minvox ) cycle
            nreg = nreg + 1
            if( nreg > maxreg )then
                nreg = maxreg; exit
            endif
            comp = labc == k; cnt = 0; cx = 0.; cy = 0.; cz = 0.; radmax = 0.
            do l = 1,box
                do j = 1,box
                    do i = 1,box
                        if( comp(i,j,l) )then
                            cnt = cnt + 1; cx = cx + i; cy = cy + j; cz = cz + l
                            rad = sqrt(real((i - box/2 - 1)**2 + (j - box/2 - 1)**2 + (l - box/2 - 1)**2)) * smpd; radmax = max(radmax, rad)
                        endif
                    enddo
                enddo
            enddo
            cx = cx / cnt - box/2 - 1; cy = cy / cnt - box/2 - 1; cz = cz / cnt - box/2 - 1; rad = sqrt(cx*cx + cy*cy + cz*cz) * smpd
            call binary_dilate(comp, box, dilate_it); soft = merge(1., 0., comp); call gauss3d(soft, box, 1.0 * scalef); sm = maxval(soft); soft = min(max(soft / max(sm, 1.e-30), 0.), 1.)
            union = max(union, soft)
            fn = trim(outdir)//'/region_'//int2str(nreg)//'.mrc'; call img%new([box,box,box], smpd); call img%set_rmat(soft, .false.); call img%write(fn); call img%kill
            write(logfhandle,'(A,I0,A,I0,A,F5.1,A,F5.0,A,F5.0,A)') '>>>   region ', nreg, ': ', nv, ' voxels (', 100. * real(mass(k) / max(sum(mass), 1.d-30)), ' % of the variability mass), centroid at ', rad, ' A radius, max radius ', radmax, ' A'
        enddo
        fn = trim(outdir)//'/all_moving.mrc'; call img%new([box,box,box], smpd); call img%set_rmat(union, .false.); call img%write(fn); call img%kill
        deallocate(v, fm, mov, labc, union, mass, order, comp, soft)
        if( allocated(m) ) deallocate(m)
    end subroutine tree_regions

    ! ================================================================== states of a node
    !> surviving cells clustered by cross-half similarity (average linkage at sim); particles ranked by their
    !! responsibility mass in their cluster; the top N (min(nmax, frac x mass)) per cluster get state = rank;
    !! states_clusters.txt (nrows rows) and clusters.txt written; returns the number of states
    subroutine tree_node_states( B, nrows, sim, nmax, frac, minmass, maxk, minrep, outdir, nstates )
        type(board_t),    intent(in)  :: B
        integer,          intent(in)  :: nrows, nmax, maxk
        real,             intent(in)  :: sim, frac, minmass, minrep
        character(len=*), intent(in)  :: outdir
        integer,          intent(out) :: nstates
        integer, allocatable :: keep(:), lab(:), order(:), states(:), cnt(:), amc(:), idx(:)
        real,    allocatable :: S(:,:), M(:,:), col(:)
        real(dp), allocatable :: massc(:)
        integer :: nk, i, j, ii, jj, K, c, u, nsel, rank, nmine
        real    :: c01, c10
        nk = count(B%crep(1:B%ngrid) >= minrep); allocate(keep(nk)); j = 0
        do i = 1,B%ngrid
            if( B%crep(i) >= minrep )then
                j = j + 1; keep(j) = i
            endif
        enddo
        allocate(S(nk,nk))
        do ii = 1,nk
            do jj = ii,nk
                i = keep(ii); j = keep(jj)
                c01 = cube_corr(B%D0(:,:,:,i), B%D1(:,:,:,j), B%kmax, 1, 2 * B%kmax); c10 = cube_corr(B%D1(:,:,:,i), B%D0(:,:,:,j), B%kmax, 1, 2 * B%kmax)
                S(ii,jj) = 0.5 * (c01 + c10); S(jj,ii) = S(ii,jj)
            enddo
        enddo
        allocate(lab(nk)); call average_linkage(S, nk, 1. - sim, lab, K)
        allocate(massc(K), source=0.d0); allocate(M(B%n,K), source=0.)
        do c = 1,K
            do ii = 1,nk
                if( lab(ii) == c ) M(:,c) = M(:,c) + B%r(keep(ii),:)
            enddo
            massc(c) = sum(real(M(:,c), dp))
        enddo
        allocate(order(K)); call argsort_desc_dp(massc, K, order)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F3.1,A,*(I0,1X))') '>>> node states: ', nk, ' of ', B%ngrid, ' cells survive, ', K, ' clusters at similarity ', sim, '; masses: ', (nint(massc(order(c))), c=1,min(K,12))
        allocate(amc(B%n))
        do i = 1,B%n
            amc(i) = maxloc(M(i,:), 1)
        enddo
        allocate(states(nrows), source=0); nstates = 0; call simple_mkdir(trim(outdir))
        open(newunit=u, file=trim(outdir)//'/clusters.txt', status='replace', action='write'); write(u,'(A)') '# state ncells mass nselected cells'
        do rank = 1,min(K, maxk)
            c = order(rank)
            if( massc(c) < real(minmass,dp) ) cycle
            nmine = count(amc == c); nsel = min(nmax, int(frac * massc(c)), nmine)
            if( nsel < 1 ) cycle
            nstates = nstates + 1
            allocate(idx(nmine), col(nmine)); j = 0
            do i = 1,B%n
                if( amc(i) == c )then
                    j = j + 1; idx(j) = i; col(j) = -M(i,c)
                endif
            enddo
            call argsort_real_idx(col, nmine)   ! idx reordered by ascending -M (descending mass)
            do j = 1,nsel
                states(B%pinds(idx(j))) = nstates
            enddo
            write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,*(I0,1X))') nstates, count(lab == c), nint(massc(c)), nsel, pack(keep, lab == c)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>>   state ', nstates, ': ', count(lab == c), ' cells, mass ', nint(massc(c)), ', selected ', nsel, ' particles'
            deallocate(idx, col)
        enddo
        close(u)
        open(newunit=u, file=trim(outdir)//'/states_clusters.txt', status='replace', action='write')
        do i = 1,nrows
            write(u,'(I0)') states(i)
        enddo
        close(u)
        deallocate(keep, S, lab, massc, M, order, amc, states)

      contains

        subroutine argsort_real_idx( key, nn )   ! sorts key ascending and permutes idx along (insertion sort; small)
            integer, intent(in)    :: nn
            real,    intent(inout) :: key(nn)
            integer :: a_, b_, t_
            real    :: kt
            do a_ = 2,nn
                kt = key(a_); t_ = idx(a_); b_ = a_ - 1
                do while( b_ >= 1 )
                    if( key(b_) <= kt ) exit
                    key(b_+1) = key(b_); idx(b_+1) = idx(b_); b_ = b_ - 1
                enddo
                key(b_+1) = kt; idx(b_+1) = t_
            enddo
        end subroutine argsort_real_idx

    end subroutine tree_node_states

    ! ================================================================== ranking along the node's own main axis
    !> cells scored by the first principal axis of their mean-centred (half-averaged) deviations; particles by their
    !! responsibility-weighted cell score; top n -> state 1, bottom n -> state 2 (states_top.txt, nrows rows)
    subroutine tree_rank_pc1( B, n, nrows, outdir )
        type(board_t),    intent(in) :: B
        integer,          intent(in) :: n, nrows
        character(len=*), intent(in) :: outdir
        complex, allocatable :: Dc(:,:,:,:), Dmean(:,:,:)
        real,    allocatable :: G(:,:), V(:,:), d(:), cs(:), s(:)
        integer, allocatable :: order(:), st(:)
        real(dp), allocatable :: sd(:)
        integer :: K, i, j, nrot, u
        K = B%ngrid; allocate(Dc(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,K), Dmean(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax)); Dmean = cmplx(0.,0.)
        do i = 1,K
            Dc(:,:,:,i) = 0.5 * (B%D0(:,:,:,i) + B%D1(:,:,:,i)); Dmean = Dmean + Dc(:,:,:,i)
        enddo
        Dmean = Dmean / real(K)
        do i = 1,K
            Dc(:,:,:,i) = Dc(:,:,:,i) - Dmean
        enddo
        allocate(G(K,K), V(K,K), d(K))
        do i = 1,K
            do j = i,K
                G(i,j) = real(sum(conjg(Dc(:,:,:,i)) * Dc(:,:,:,j))); G(j,i) = G(i,j)
            enddo
        enddo
        nrot = 0; call jacobi(G, K, K, d, V, nrot); call eigsrt(d, V, K, K)
        allocate(cs(K)); cs = sqrt(max(d(1), 0.)) * V(:,1); cs = cs - sum(cs) / real(K)
        write(logfhandle,'(A,4F6.2)') '>>> own main axis: eigenvalues (first 4, relative) ', (max(d(i),0.) / max(d(1),1.e-30), i=1,min(4,K))
        allocate(s(B%n), sd(B%n))
        do i = 1,B%n
            s(i) = sum(B%r(1:K,i) * cs); sd(i) = -real(s(i), dp)
        enddo
        allocate(order(B%n)); call argsort_desc_dp(-sd, B%n, order)   ! descending s
        allocate(st(nrows), source=0)
        do i = 1,min(n, B%n)
            st(B%pinds(order(i))) = 1
        enddo
        do i = max(1, B%n - n + 1),B%n
            st(B%pinds(order(i))) = 2
        enddo
        call simple_mkdir(trim(outdir)); open(newunit=u, file=trim(outdir)//'/states_top.txt', status='replace', action='write')
        do i = 1,nrows
            write(u,'(I0)') st(i)
        enddo
        close(u)
        write(logfhandle,'(A,F6.3,A,F6.3,A,I0,A,F6.3)') '>>> cell score range ', minval(cs), ' .. ', maxval(cs), '; top ', n, ' particles >= ', s(order(min(n,B%n)))
        deallocate(Dc, Dmean, G, V, d, cs, s, sd, order, st)
    end subroutine tree_rank_pc1

    !> level-2 ranking by the in-house probabilistic PCA (simple_ppca) fitted to the node's cell deviations, every cell
    !! weighted by its responsibility mass (flex_pca's model on the cells instead of the particles); rank by BIC up to qmax
    subroutine tree_rank_ppca( B, n, nrows, outdir, qmax )
        use simple_ppca, only: ppca
        type(board_t),    intent(in) :: B
        integer,          intent(in) :: n, nrows, qmax
        character(len=*), intent(in) :: outdir
        type(ppca) :: P
        complex, allocatable :: Dc(:,:,:,:), Dmean(:,:,:)
        real,    allocatable :: X(:,:), w(:), cs(:), s(:), feat(:), G(:,:), V(:,:), d(:), c1(:)
        real(dp), allocatable :: sd(:), bic(:)
        integer, allocatable :: order(:), st(:)
        integer :: K, i, j, nc, nd, u, q, qbest, nrot
        real    :: agree
        K = B%ngrid; nc = (2 * B%kmax + 1)**3; nd = 2 * nc
        allocate(w(K), Dc(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,K), Dmean(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax)); Dmean = cmplx(0.,0.)
        do i = 1,K
            w(i) = sum(B%r(i,:)); Dc(:,:,:,i) = 0.5 * (B%D0(:,:,:,i) + B%D1(:,:,:,i))
        enddo
        w = w / max(sum(w), 1.e-30)
        do i = 1,K
            Dmean = Dmean + w(i) * Dc(:,:,:,i)
        enddo
        allocate(X(nd,K))
        do i = 1,K
            Dc(:,:,:,i) = Dc(:,:,:,i) - Dmean
            X(1:nc,i) = sqrt(w(i) * real(K)) * real(reshape(Dc(:,:,:,i), [nc])); X(nc+1:nd,i) = sqrt(w(i) * real(K)) * aimag(reshape(Dc(:,:,:,i), [nc]))
        enddo
        allocate(bic(qmax)); qbest = 1
        do q = 1,qmax   ! BIC rank selection with the class's own criterion
            call P%new(K, nd, q); call P%master(X, 100); bic(q) = P%calc_bic(X)
            if( bic(q) < bic(qbest) ) qbest = q
        enddo
        call P%new(K, nd, qbest); call P%master(X, 200)
        allocate(cs(K))
        do i = 1,K
            feat = P%get_feat(i); cs(i) = feat(1) / sqrt(max(w(i) * real(K), 1.e-30)); deallocate(feat)   ! the cell's coordinate (row weight undone)
        enddo
        cs = cs - sum(w * cs)
        ! agreement with the unweighted Gram-matrix axis (tree_rank_pc1) as a readout
        allocate(G(K,K), V(K,K), d(K), c1(K))
        do i = 1,K
            do j = i,K
                G(i,j) = real(sum(conjg(Dc(:,:,:,i)) * Dc(:,:,:,j))); G(j,i) = G(i,j)
            enddo
        enddo
        nrot = 0; call jacobi(G, K, K, d, V, nrot); call eigsrt(d, V, K, K); c1 = V(:,1) - sum(V(:,1)) / real(K)
        agree = abs(sum(cs * c1)) / max(sqrt(sum(cs * cs) * sum(c1 * c1)), 1.e-30)
        write(logfhandle,'(A,I0,A,*(ES10.3,1X))') '>>> cell PPCA (simple_ppca, mass-weighted cells): BIC rank ', qbest, ' of candidates; BIC ', (bic(q), q=1,qmax)
        write(logfhandle,'(A,F5.2)') '>>> cell PPCA axis vs the Gram-matrix axis (|corr| over cells) ', agree
        allocate(s(B%n), sd(B%n))
        do i = 1,B%n
            s(i) = sum(B%r(1:K,i) * cs); sd(i) = -real(s(i), dp)
        enddo
        allocate(order(B%n)); call argsort_desc_dp(-sd, B%n, order)
        allocate(st(nrows), source=0)
        do i = 1,min(n, B%n)
            st(B%pinds(order(i))) = 1
        enddo
        do i = max(1, B%n - n + 1),B%n
            st(B%pinds(order(i))) = 2
        enddo
        call simple_mkdir(trim(outdir)); open(newunit=u, file=trim(outdir)//'/states_top.txt', status='replace', action='write')
        do i = 1,nrows
            write(u,'(I0)') st(i)
        enddo
        close(u)
        write(logfhandle,'(A,F6.3,A,F6.3,A,I0,A,F6.3)') '>>> cell score range ', minval(cs), ' .. ', maxval(cs), '; top ', n, ' particles >= ', s(order(min(n,B%n)))
        call P%kill; deallocate(Dc, Dmean, X, w, cs, s, sd, order, st, G, V, d, c1, bic)
    end subroutine tree_rank_ppca

    ! ================================================================== identification against deposited maps (readout)
    !> cells' deviations (from the mean cell) vs deposited maps' deviations from the consensus, correlated at lo-hi A
    !! on the coarse grid inside a 100 A soft mask; list file: lines "id path"; writes per-cell best match when save given
    subroutine tree_deposited_match( B, box, smpd, listfile, consfile, lo, hi, save )
        type(board_t),    intent(in) :: B
        integer,          intent(in) :: box
        real,             intent(in) :: smpd, lo, hi
        character(len=*), intent(in) :: listfile, consfile
        character(len=*), optional, intent(in) :: save
        type(image) :: img
        type(string) :: fn
        character(len=1024) :: line, path
        integer, allocatable :: ids(:), best(:), nbest(:)
        real,    allocatable :: cons(:,:,:), rm(:,:,:), S(:,:)
        complex, allocatable :: Dep(:,:,:,:), Dc(:,:,:,:), Dmean(:,:,:)
        integer :: nd, i, j, u, io, K, slo, shi, nid
        real    :: bc
        ! deposited list
        nd = 0
        open(newunit=u, file=trim(listfile), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) > 0 .and. line(1:1) /= '#' ) nd = nd + 1
        enddo
        rewind(u); allocate(ids(nd)); allocate(Dep(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,nd)); allocate(cons(box,box,box), rm(box,box,box))
        fn = trim(consfile); call img%new([box,box,box], smpd); call img%read(fn); call img%get_rmat_sub(cons); call img%kill
        i = 0
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1; line = adjustl(line); j = index(line, ' '); read(line(1:j-1),*) ids(i); path = adjustl(line(j+1:))   ! not list-directed: a slash in the path ends a list-directed record
            fn = trim(path); call img%new([box,box,box], smpd); call img%read(fn); call img%get_rmat_sub(rm); call img%kill
            rm = rm - cons; call realspace_to_cube(rm, box, smpd, B%kmax, Dep(:,:,:,i), rmask=100.)
        enddo
        close(u)
        K = B%ngrid; allocate(Dc(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax,K), Dmean(-B%kmax:B%kmax,-B%kmax:B%kmax,-B%kmax:B%kmax)); Dmean = cmplx(0.,0.)
        do i = 1,K
            Dc(:,:,:,i) = 0.5 * (B%D0(:,:,:,i) + B%D1(:,:,:,i)); Dmean = Dmean + Dc(:,:,:,i)
        enddo
        Dmean = Dmean / real(K)
        slo = max(1, nint(real(box) * smpd / lo)); shi = min(B%kmax, nint(real(box) * smpd / hi))
        allocate(S(K,nd)); S = 0.; allocate(best(K), nbest(nd), source=0)
        do i = 1,K
            Dc(:,:,:,i) = Dc(:,:,:,i) - Dmean
            do j = 1,nd
                S(i,j) = cube_corr(Dc(:,:,:,i), Dep(:,:,:,j), B%kmax, slo, shi)
            enddo
            best(i) = maxloc(S(i,:), 1); nbest(best(i)) = nbest(best(i)) + 1
        enddo
        write(logfhandle,'(A,I0,A,F4.0,A,F4.0,A,I0,A)') '>>> identification vs the ', nd, ' deposited conformers (deviation corr at ', lo, '-', hi, ' A on the coarse grid), ', K, ' cells:'
        write(logfhandle,'(A)') '>>>   conformer  best cell corr  cells whose best match it is'
        nid = 0
        do j = 1,nd
            bc = maxval(S(:,j))
            if( nbest(j) >= 1 .and. bc >= 0.3 ) nid = nid + 1
            write(logfhandle,'(A,I5,A,F6.2,A,I3,A)') '>>>   ', ids(j), '      ', bc, '        ', nbest(j), merge('   identified', '   -         ', nbest(j) >= 1 .and. bc >= 0.3)
        enddo
        write(logfhandle,'(A,I0,A,I0,A,F5.2)') '>>>   identified ', nid, ' / ', nd, '; mean best-match corr over cells ', sum(maxval(S, 2)) / real(K)
        if( present(save) )then
            open(newunit=u, file=trim(save), status='replace', action='write'); write(u,'(A,*(I0,1X))') '# cell best_conformer best_corr | corr per conformer ', ids
            do i = 1,K
                write(u,'(I0,1X,I0,1X,F8.4,1X,*(F7.3,1X))') i, ids(best(i)), S(i,best(i)), S(i,:)
            enddo
            close(u)
        endif
        deallocate(ids, best, nbest, cons, rm, S, Dep, Dc, Dmean)
    end subroutine tree_deposited_match

    ! ================================================================== the landscape over the leaves
    !> surviving cells of several boards clustered by cross-half similarity (average linkage) at several similarity cuts
    subroutine tree_count_states( dirs, ndirs, minrep, sims, nsims, box )
        character(len=*), intent(in) :: dirs(ndirs)
        integer,          intent(in) :: ndirs, nsims, box
        real,             intent(in) :: minrep, sims(nsims)
        type(board_t) :: Bd
        complex, allocatable :: C0(:,:,:,:), C1(:,:,:,:)
        real,    allocatable :: S(:,:), W(:), wsz(:)
        integer, allocatable :: lab(:), src(:)
        integer :: j, k, ntot, c, kmax, i, nc, nK, is
        real    :: c01, c10
        ntot = 0; kmax = -1
        do j = 1,ndirs
            call read_board(dirs(j), Bd); ntot = ntot + count(Bd%crep(1:Bd%ngrid) >= minrep)
            if( kmax < 0 ) kmax = Bd%kmax
            if( Bd%kmax /= kmax ) THROW_HARD('count_states: the boards have different coarse grids')
            call kill_board(Bd)
        enddo
        allocate(C0(-kmax:kmax,-kmax:kmax,-kmax:kmax,ntot), C1(-kmax:kmax,-kmax:kmax,-kmax:kmax,ntot), W(ntot), src(ntot)); c = 0
        do j = 1,ndirs
            call read_board(dirs(j), Bd)
            do k = 1,Bd%ngrid
                if( Bd%crep(k) < minrep ) cycle
                c = c + 1; C0(:,:,:,c) = Bd%D0(:,:,:,k); C1(:,:,:,c) = Bd%D1(:,:,:,k); W(c) = Bd%w(k); src(c) = j
            enddo
            call kill_board(Bd)
        enddo
        write(logfhandle,'(A,I0,A,F4.2,A,I0,A)') '>>> LANDSCAPE: ', ntot, ' surviving cells (reproducibility >= ', minrep, ') from ', ndirs, ' boards'
        allocate(S(ntot,ntot))
        do i = 1,ntot
            do k = i,ntot
                c01 = cube_corr(C0(:,:,:,i), C1(:,:,:,k), kmax, 1, 2 * kmax); c10 = cube_corr(C1(:,:,:,i), C0(:,:,:,k), kmax, 1, 2 * kmax)
                S(i,k) = 0.5 * (c01 + c10); S(k,i) = S(i,k)
            enddo
        enddo
        allocate(lab(ntot))
        do is = 1,nsims
            call average_linkage(S, ntot, 1. - sims(is), lab, nK); allocate(wsz(nK), source=0.)
            do i = 1,ntot
                wsz(lab(i)) = wsz(lab(i)) + W(i)
            enddo
            write(logfhandle,'(A,F3.1,A,I0,A,I0,A,I0,A,I0,A)') '>>>   similarity ', sims(is), ': ', nK, ' states; particles per state median ', nint(median_copy(wsz)), ' (min ', nint(minval(wsz)), ', max ', nint(maxval(wsz)), ')'
            deallocate(wsz)
        enddo
        deallocate(C0, C1, W, src, S, lab)
    end subroutine tree_count_states

    ! ================================================================== numerical helpers
    !> average-linkage agglomeration on a similarity matrix, merging while the closest pair is below distance t
    subroutine average_linkage( S, n, t, lab, K )
        integer, intent(in)  :: n
        real,    intent(in)  :: S(n,n), t
        integer, intent(out) :: lab(n), K
        real,    allocatable :: D(:,:)
        integer, allocatable :: sz(:)
        logical, allocatable :: act(:)
        integer :: i, j, a, b, c
        real    :: dmin, dnew
        allocate(D(n,n), sz(n), act(n)); D = min(max(1. - S, 0.), 2.); sz = 1; act = .true.
        do i = 1,n
            lab(i) = i; D(i,i) = 0.
        enddo
        do
            dmin = huge(dmin); a = 0; b = 0
            do i = 1,n
                if( .not. act(i) ) cycle
                do j = i+1,n
                    if( .not. act(j) ) cycle
                    if( D(i,j) < dmin )then
                        dmin = D(i,j); a = i; b = j
                    endif
                enddo
            enddo
            if( a == 0 .or. dmin > t ) exit
            do c = 1,n
                if( .not. act(c) .or. c == a .or. c == b ) cycle
                dnew = (real(sz(a)) * D(a,c) + real(sz(b)) * D(b,c)) / real(sz(a) + sz(b)); D(a,c) = dnew; D(c,a) = dnew
            enddo
            act(b) = .false.; sz(a) = sz(a) + sz(b)
            do i = 1,n
                if( lab(i) == b ) lab(i) = a
            enddo
        enddo
        ! relabel 1..K
        K = 0
        do i = 1,n
            if( act(i) )then
                K = K + 1; sz(i) = K
            endif
        enddo
        do i = 1,n
            lab(i) = sz(lab(i))
        enddo
        deallocate(D, sz, act)
    end subroutine average_linkage

    !> k-means (Lloyd) with nrestart random starts; X(n,p), labels 1..k, centroids cent(k,p)
    subroutine tree_kmeans( X, n, p, k, nrestart, lab, cent )
        integer, intent(in)  :: n, p, k, nrestart
        real,    intent(in)  :: X(n,p)
        integer, intent(out) :: lab(n)
        real,    intent(out) :: cent(k,p)
        real,    allocatable :: C(:,:), Cn(:,:), rnd(:)
        integer, allocatable :: l(:), cnt(:)
        real(dp) :: inertia, best
        integer  :: rs, it, i, j, jb, pick
        real     :: d, dmin
        allocate(C(k,p), Cn(k,p), l(n), cnt(k), rnd(n)); best = huge(best)
        do rs = 1,nrestart
            call random_number(rnd)
            do j = 1,k
                pick = 1 + int(rnd(min(j,n)) * real(n - 1)); C(j,:) = X(pick,:)
            enddo
            do it = 1,60
                do i = 1,n
                    dmin = huge(dmin); jb = 1
                    do j = 1,k
                        d = sum((X(i,:) - C(j,:))**2)
                        if( d < dmin )then
                            dmin = d; jb = j
                        endif
                    enddo
                    l(i) = jb
                enddo
                Cn = 0.; cnt = 0
                do i = 1,n
                    Cn(l(i),:) = Cn(l(i),:) + X(i,:); cnt(l(i)) = cnt(l(i)) + 1
                enddo
                do j = 1,k
                    if( cnt(j) > 0 )then
                        Cn(j,:) = Cn(j,:) / real(cnt(j))
                    else
                        Cn(j,:) = C(j,:)
                    endif
                enddo
                if( maxval(abs(Cn - C)) < 1.e-6 ) exit
                C = Cn
            enddo
            inertia = 0.d0
            do i = 1,n
                inertia = inertia + real(sum((X(i,:) - C(l(i),:))**2), dp)
            enddo
            if( inertia < best )then
                best = inertia; lab = l; cent = C
            endif
        enddo
        deallocate(C, Cn, l, cnt, rnd)
    end subroutine tree_kmeans

    !> 6-connected components of a logical volume; labels 1..nl
    subroutine label_cc( m, box, lab, nl )
        integer, intent(in)  :: box
        logical, intent(in)  :: m(box,box,box)
        integer, intent(out) :: lab(box,box,box), nl
        lab = 0; nl = 0; call label_cc_offset(m, box, lab, nl)
    end subroutine label_cc

    !> as label_cc but adds components to an existing labelling (labels continue from nl)
    subroutine label_cc_offset( m, box, lab, nl )
        integer, intent(in)    :: box
        logical, intent(in)    :: m(box,box,box)
        integer, intent(inout) :: lab(box,box,box), nl
        integer, allocatable :: stack(:,:)
        integer :: i, j, l, top, x, y, z, dx(6), dy(6), dz(6), q, nx, ny, nz
        dx = [1,-1,0,0,0,0]; dy = [0,0,1,-1,0,0]; dz = [0,0,0,0,1,-1]; allocate(stack(3, box*box*box))
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    if( .not. m(i,j,l) .or. lab(i,j,l) /= 0 ) cycle
                    nl = nl + 1; top = 1; stack(:,1) = [i,j,l]; lab(i,j,l) = nl
                    do while( top > 0 )
                        x = stack(1,top); y = stack(2,top); z = stack(3,top); top = top - 1
                        do q = 1,6
                            nx = x + dx(q); ny = y + dy(q); nz = z + dz(q)
                            if( nx < 1 .or. nx > box .or. ny < 1 .or. ny > box .or. nz < 1 .or. nz > box ) cycle
                            if( m(nx,ny,nz) .and. lab(nx,ny,nz) == 0 )then
                                lab(nx,ny,nz) = nl; top = top + 1; stack(:,top) = [nx,ny,nz]
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(stack)
    end subroutine label_cc_offset

    subroutine binary_dilate( m, box, iters )
        integer, intent(in)    :: box, iters
        logical, intent(inout) :: m(box,box,box)
        logical, allocatable :: t(:,:,:)
        integer :: it, i, j, l
        allocate(t(box,box,box))
        do it = 1,iters
            t = m
            do l = 1,box
                do j = 1,box
                    do i = 1,box
                        if( m(i,j,l) ) cycle
                        if( i > 1 ) t(i,j,l) = t(i,j,l) .or. m(i-1,j,l)
                        if( i < box ) t(i,j,l) = t(i,j,l) .or. m(i+1,j,l)
                        if( j > 1 ) t(i,j,l) = t(i,j,l) .or. m(i,j-1,l)
                        if( j < box ) t(i,j,l) = t(i,j,l) .or. m(i,j+1,l)
                        if( l > 1 ) t(i,j,l) = t(i,j,l) .or. m(i,j,l-1)
                        if( l < box ) t(i,j,l) = t(i,j,l) .or. m(i,j,l+1)
                    enddo
                enddo
            enddo
            m = t
        enddo
        deallocate(t)
    end subroutine binary_dilate

    subroutine binary_erode( m, box )
        integer, intent(in)    :: box
        logical, intent(inout) :: m(box,box,box)
        logical, allocatable :: t(:,:,:)
        integer :: i, j, l
        allocate(t(box,box,box)); t = m
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    if( .not. m(i,j,l) ) cycle
                    if( i == 1 .or. i == box .or. j == 1 .or. j == box .or. l == 1 .or. l == box )then
                        t(i,j,l) = .false.
                    else
                        t(i,j,l) = m(i-1,j,l) .and. m(i+1,j,l) .and. m(i,j-1,l) .and. m(i,j+1,l) .and. m(i,j,l-1) .and. m(i,j,l+1)
                    endif
                enddo
            enddo
        enddo
        m = t; deallocate(t)
    end subroutine binary_erode

    subroutine binary_open( m, box )
        integer, intent(in)    :: box
        logical, intent(inout) :: m(box,box,box)
        call binary_erode(m, box); call binary_dilate(m, box, 1)
    end subroutine binary_open

    !> separable Gaussian smoothing (sigma in voxels)
    subroutine gauss3d( v, box, sigma )
        integer, intent(in)    :: box
        real,    intent(inout) :: v(box,box,box)
        real,    intent(in)    :: sigma
        real, allocatable :: ker(:), t(:,:,:)
        integer :: r, q, i, j, l, a, idx
        real    :: s
        r = max(1, nint(3. * sigma)); allocate(ker(-r:r)); s = 0.
        do q = -r,r
            ker(q) = exp(-0.5 * real(q*q) / sigma**2); s = s + ker(q)
        enddo
        ker = ker / s; allocate(t(box,box,box))
        t = 0.
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    do q = -r,r
                        idx = i + q
                        if( idx >= 1 .and. idx <= box ) t(i,j,l) = t(i,j,l) + ker(q) * v(idx,j,l)
                    enddo
                enddo
            enddo
        enddo
        v = 0.
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    do q = -r,r
                        idx = j + q
                        if( idx >= 1 .and. idx <= box ) v(i,j,l) = v(i,j,l) + ker(q) * t(i,idx,l)
                    enddo
                enddo
            enddo
        enddo
        t = 0.
        do l = 1,box
            do j = 1,box
                do i = 1,box
                    do q = -r,r
                        idx = l + q
                        if( idx >= 1 .and. idx <= box ) t(i,j,l) = t(i,j,l) + ker(q) * v(i,j,idx)
                    enddo
                enddo
            enddo
        enddo
        v = t; deallocate(ker, t)
    end subroutine gauss3d

    function percentile3d( v, box, pc ) result( val )
        integer, intent(in) :: box
        real,    intent(in) :: v(box,box,box), pc
        real :: val
        real(dp), allocatable :: w(:)
        integer :: k
        allocate(w(box*box*box)); w = real(reshape(v, [box*box*box]), dp); call hpsort_local(w); k = max(1, min(size(w), nint(pc / 100. * real(size(w))))); val = real(w(k)); deallocate(w)
    end function percentile3d

    function median_copy( v ) result( med )
        real, intent(in) :: v(:)
        real :: med
        real(dp), allocatable :: w(:)
        integer :: n
        allocate(w(size(v))); w = real(v, dp); call hpsort_local(w); n = size(w)
        if( mod(n,2) == 1 )then
            med = real(w((n+1)/2))
        else
            med = real(0.5d0 * (w(n/2) + w(n/2+1)))
        endif
        deallocate(w)
    end function median_copy

    !> indices of v sorted descending
    subroutine argsort_desc_dp( v, n, order )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: v(n)
        integer,  intent(out) :: order(n)
        integer :: i, j, t
        do i = 1,n
            order(i) = i
        enddo
        do i = 2,n
            t = order(i); j = i - 1
            do while( j >= 1 )
                if( v(order(j)) >= v(t) ) exit
                order(j+1) = order(j); j = j - 1
            enddo
            order(j+1) = t
        enddo
    end subroutine argsort_desc_dp

    subroutine hpsort_local( ra )   ! heapsort, Numerical Recipes form
        real(dp), intent(inout) :: ra(:)
        integer  :: n, i, ir, j, l
        real(dp) :: rra
        n = size(ra)
        if( n < 2 ) return
        l = n/2 + 1; ir = n
        do
            if( l > 1 )then
                l = l - 1; rra = ra(l)
            else
                rra = ra(ir); ra(ir) = ra(1); ir = ir - 1
                if( ir == 1 )then
                    ra(1) = rra; return
                endif
            endif
            i = l; j = l + l
            do while( j <= ir )
                if( j < ir )then
                    if( ra(j) < ra(j+1) ) j = j + 1
                endif
                if( rra < ra(j) )then
                    ra(i) = ra(j); i = j; j = j + j
                else
                    j = ir + 1
                endif
            enddo
            ra(i) = rra
        enddo
    end subroutine hpsort_local

end module simple_klem_tree
