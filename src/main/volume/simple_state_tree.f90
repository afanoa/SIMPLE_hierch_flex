!@descr: bottom-up tree of states from the unbiased cross-half difference energy of their even/odd maps
module simple_state_tree
use simple_defs,         only: dp, logfhandle
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: fopen, fclose, fileiochk
use simple_syslib,       only: simple_mkdir
use simple_imghead,      only: find_ldim_nptcls
use simple_image,        only: image
implicit none

public :: state_tree
private
#include "simple_local_flags.inc"

real(dp), parameter :: ENV_LP = 15._dp, ENV_FRAC = 0.1_dp, ENV_SIGMA = 2._dp
integer,  parameter :: ENV_GROW = 4, ENV_REACH = 16, MAXLINE = 262144

type :: str_slot
    character(len=:), allocatable :: s
end type str_slot

type :: state_tree
    private
    integer               :: nleaves = 0, nsteps = 0, nbudget = -1
    integer,  allocatable :: names(:), steps(:,:)
    real(dp), allocatable :: counts(:), gram(:,:), noisev(:,:), nw(:), zs(:), d2s(:), wsteps(:,:)
    real(dp)              :: infl = 1._dp
    real                  :: budget = 0.
    character(len=8)      :: linkage = 'complete'
    logical               :: l_exists = .false.
  contains
    procedure :: new
    procedure :: build
    procedure :: write
    procedure :: read
    procedure :: cut
    procedure :: get_nleaves
    procedure :: get_names
    procedure :: get_steps
    procedure :: kill
end type state_tree

contains

    subroutine new( self, evens, odds, names, counts, smpd, hp, lp )
        class(state_tree), intent(inout) :: self
        class(string),     intent(in)    :: evens(:), odds(:)
        integer,           intent(in)    :: names(:)
        real(dp),          intent(in)    :: counts(:)
        real,              intent(in)    :: smpd, hp, lp
        type(image)           :: vol
        real(dp), allocatable :: avg(:,:,:), msk(:,:,:), wgt(:)
        complex,  allocatable :: ev(:,:), od(:,:)
        integer,  allocatable :: phys(:,:), shell(:), shmap(:)
        real(dp) :: res, s2, acc
        integer  :: ldim(3), ifoo, n, box, i, j, a, b, h, k, l, nb, nsh, ipass, lims(3,2)
        call self%kill
        n = size(names)
        if( n < 2 ) THROW_HARD('the state tree needs at least two states')
        if( size(evens) /= n .or. size(odds) /= n .or. size(counts) /= n ) THROW_HARD('inconsistent state map lists')
        self%nleaves = n
        self%names   = names
        self%counts  = counts
        call find_ldim_nptcls(evens(1), ldim, ifoo)
        box = ldim(1)
        write(logfhandle,'(A)') '>>> '//int2str(n)//' states, box '//int2str(box)//', band '//fmtr(real(hp,dp),0)//&
            &' - '//fmtr(real(lp,dp),0)//' A'
        call vol%new(ldim, smpd)
        allocate(avg(box,box,box), source=0._dp)
        do i = 1, n
            call add_map(evens(i), 0.5_dp*counts(i))
            call add_map(odds(i),  0.5_dp*counts(i))
        end do
        avg = avg / sum(counts)
        call calc_envelope(avg, smpd, msk)
        deallocate(avg)
        self%infl = real(box,dp)**3 * sum(msk**4) / sum(msk**2)**2
        lims = vol%loop_lims(2)
        do ipass = 1, 2
            nb = 0
            do l = lims(3,1), lims(3,2)
                do k = lims(2,1), lims(2,2)
                    do h = lims(1,1), lims(1,2)
                        s2 = real(h*h + k*k + l*l, dp)
                        if( s2 < 0.5_dp ) cycle
                        res = real(box,dp) * real(smpd,dp) / sqrt(s2)
                        if( res > real(hp,dp) .or. res < real(lp,dp) ) cycle
                        nb = nb + 1
                        if( ipass == 1 ) cycle
                        phys(:,nb) = vol%comp_addr_phys(h,k,l)
                        shell(nb)  = nint(sqrt(real(s2)))
                        wgt(nb)    = merge(1._dp, 2._dp, h == 0 .or. h == box/2)
                    end do
                end do
            end do
            if( ipass == 1 ) allocate(phys(3,nb), shell(nb), wgt(nb))
        end do
        allocate(shmap(0:maxval(shell)), source=0)
        do j = 1, nb
            shmap(shell(j)) = 1
        end do
        nsh = 0
        do i = 0, ubound(shmap,1)
            if( shmap(i) == 0 ) cycle
            nsh      = nsh + 1
            shmap(i) = nsh
        end do
        shell = shmap(shell)
        write(logfhandle,'(A)') '>>> envelope volume fraction '//fmtr(sum(msk)/real(box,dp)**3,3)//&
            &' (correlation inflation '//fmtr(self%infl,1)//'); '//int2str(nb)//' band coefficients in '//&
            &int2str(nsh)//' shells'
        allocate(ev(nb,n), od(nb,n))
        do i = 1, n
            call masked_coefs(evens(i), ev(:,i))
            call masked_coefs(odds(i),  od(:,i))
        end do
        call vol%kill
        allocate(self%nw(nsh), self%noisev(nsh,n), self%gram(n,n), source=0._dp)
        do j = 1, nb
            self%nw(shell(j)) = self%nw(shell(j)) + wgt(j)
        end do
        do i = 1, n
            do j = 1, nb
                self%noisev(shell(j),i) = self%noisev(shell(j),i) + wgt(j) * 0.5_dp * real(abs(ev(j,i) - od(j,i)),dp)**2
            end do
            self%noisev(:,i) = self%noisev(:,i) / self%nw
        end do
        !$omp parallel do collapse(2) default(shared) private(a,b,j,acc) schedule(static) proc_bind(close)
        do b = 1, n
            do a = 1, n
                acc = 0._dp
                do j = 1, nb
                    acc = acc + wgt(j) * (real(real(ev(j,a)),dp) * real(real(od(j,b)),dp) + &
                        &real(aimag(ev(j,a)),dp) * real(aimag(od(j,b)),dp))
                end do
                self%gram(a,b) = acc
            end do
        end do
        !$omp end parallel do
        deallocate(ev, od, phys, shell, wgt, shmap, msk)
        self%l_exists = .true.

    contains

        subroutine add_map( fname, w )
            class(string), intent(in) :: fname
            real(dp),      intent(in) :: w
            real, pointer :: rmat(:,:,:)
            call vol%set_ft(.false.)
            call vol%read(fname)
            call vol%get_rmat_ptr(rmat)
            avg = avg + w * real(rmat, dp)
        end subroutine add_map

        subroutine masked_coefs( fname, coefs )
            class(string), intent(in)  :: fname
            complex,       intent(out) :: coefs(:)
            complex, allocatable :: cmat(:,:,:)
            real, pointer :: rmat(:,:,:)
            integer :: jj
            call vol%set_ft(.false.)
            call vol%read(fname)
            call vol%get_rmat_ptr(rmat)
            rmat = rmat * real(msk)
            call vol%fft
            cmat = vol%get_cmat()
            do jj = 1, size(coefs)
                coefs(jj) = cmat(phys(1,jj),phys(2,jj),phys(3,jj))
            end do
        end subroutine masked_coefs

    end subroutine new

    subroutine build( self, linkage, budget )
        class(state_tree), intent(inout) :: self
        character(len=*),  intent(in)    :: linkage
        real,              intent(in)    :: budget
        real(dp), allocatable :: cg(:,:), vg(:,:), dg(:,:), w(:)
        integer,  allocatable :: alive(:)
        real(dp) :: d, bestd, wp, wq, s, z
        integer  :: n, nn, nalive, ii, jj, kk, p, q, r, k, bp, bq
        logical  :: l_budget, l_found
        if( .not. self%l_exists ) THROW_HARD('state tree has no states; build')
        select case(trim(linkage))
            case('complete','average','centroid')
                self%linkage = trim(linkage)
            case DEFAULT
                THROW_HARD('unsupported linkage: '//trim(linkage))
        end select
        self%budget = budget
        n  = self%nleaves
        nn = 2*n - 1
        allocate(cg(0:nn-1,0:nn-1), dg(0:nn-1,0:nn-1), w(0:nn-1), source=0._dp)
        allocate(vg(size(self%nw),0:nn-1), source=0._dp)
        cg(0:n-1,0:n-1) = self%gram
        vg(:,0:n-1)     = self%noisev
        w(0:n-1)        = self%counts
        do p = 0, n-1
            do q = p+1, n-1
                dg(p,q) = d2val(cg(p,p), cg(q,q), cg(p,q), cg(q,p))
                dg(q,p) = dg(p,q)
            end do
        end do
        if( allocated(self%steps) ) deallocate(self%steps, self%zs, self%d2s, self%wsteps)
        allocate(self%steps(3,n-1), self%wsteps(2,n-1), self%zs(n-1), self%d2s(n-1))
        allocate(alive(n))
        alive  = [(p, p=0,n-1)]
        nalive = n
        self%nsteps  = 0
        self%nbudget = -1
        l_budget = budget > 0.
        do while( nalive > 1 )
            l_found = .false.
            bestd   = 0._dp
            bp      = -1
            bq      = -1
            do ii = 1, nalive
                p = alive(ii)
                do jj = ii+1, nalive
                    q = alive(jj)
                    if( l_budget .and. w(p) >= real(budget,dp) .and. w(q) >= real(budget,dp) ) cycle
                    if( self%linkage == 'centroid' )then
                        d = d2val(cg(p,p), cg(q,q), cg(p,q), cg(q,p))
                    else
                        d = dg(p,q)
                    endif
                    if( .not. l_found .or. d < bestd )then
                        bestd   = d
                        bp      = p
                        bq      = q
                        l_found = .true.
                    endif
                end do
            end do
            if( .not. l_found )then
                self%nbudget = self%nsteps
                l_budget     = .false.
                write(logfhandle,'(A)') '>>> budget '//int2str(int(budget))//': '//int2str(nalive)//&
                    &' states, particles smallest '//int2str(int(minval(w(alive(1:nalive)))))
                cycle
            endif
            p  = bp
            q  = bq
            z  = zval(cg(p,p), cg(q,q), cg(p,q), cg(q,p), vg(:,p), vg(:,q), self%nw, self%infl)
            wp = w(p)
            wq = w(q)
            s  = wp + wq
            r  = n + self%nsteps
            do kk = 1, nalive
                k = alive(kk)
                cg(r,k) = (wp*cg(p,k) + wq*cg(q,k)) / s
                cg(k,r) = (wp*cg(k,p) + wq*cg(k,q)) / s
            end do
            do kk = 1, nalive
                k = alive(kk)
                if( self%linkage == 'average' )then
                    dg(r,k) = (wp*dg(p,k) + wq*dg(q,k)) / s
                else
                    dg(r,k) = max(dg(p,k), dg(q,k))
                endif
                dg(k,r) = dg(r,k)
            end do
            cg(r,r) = (wp*wp*cg(p,p) + wp*wq*(cg(p,q) + cg(q,p)) + wq*wq*cg(q,q)) / s**2
            vg(:,r) = (wp*wp*vg(:,p) + wq*wq*vg(:,q)) / s**2
            w(r)    = s
            self%nsteps = self%nsteps + 1
            self%steps(:,self%nsteps)  = [p, q, r]
            self%zs(self%nsteps)       = z
            self%d2s(self%nsteps)      = bestd
            self%wsteps(:,self%nsteps) = [wp, wq]
            write(logfhandle,'(A)') 'merge '//int2str(self%nsteps)//': d2 '//fmtr(bestd,4)//' z '//fmtr(z,1)//&
                &'  ['//int2str(p)//'] ('//int2str(int(wp))//') + ['//int2str(q)//'] ('//int2str(int(wq))//')'
            alive(1:nalive-2) = pack(alive(1:nalive), alive(1:nalive) /= p .and. alive(1:nalive) /= q)
            nalive = nalive - 1
            alive(nalive) = r
        end do
        deallocate(cg, dg, w, vg, alive)
    end subroutine build

    subroutine write( self, labels )
        class(state_tree), intent(in) :: self
        integer, optional, intent(in) :: labels(:)
        type(str_slot), allocatable :: nwk(:)
        real(dp), allocatable :: zleaf(:,:), dleaf(:,:), ht(:)
        type(string) :: fname
        integer :: funit, ios, n, i, j, p, q, r, g
        n = self%nleaves
        allocate(zleaf(n,n), dleaf(n,n), source=0._dp)
        do i = 1, n
            do j = i+1, n
                zleaf(i,j) = zval(self%gram(i,i), self%gram(j,j), self%gram(i,j), self%gram(j,i), &
                    &self%noisev(:,i), self%noisev(:,j), self%nw, self%infl)
                dleaf(i,j) = d2val(self%gram(i,i), self%gram(j,j), self%gram(i,j), self%gram(j,i))
                zleaf(j,i) = zleaf(i,j)
                dleaf(j,i) = dleaf(i,j)
            end do
        end do
        call write_table(string('leaf_z.txt'),  zleaf, 1)
        call write_table(string('leaf_d2.txt'), dleaf, 5)
        write(logfhandle,'(A)') '>>> leaf pair z (cross-half difference energy / sd): min '//&
            &fmtr(minval([((zleaf(i,j), j=i+1,n), i=1,n-1)]),1)//', max '//fmtr(maxval(zleaf),1)
        fname = 'merge_tree.txt'
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call fileiochk('state_tree: cannot write '//fname%to_char(), ios)
        write(funit,'(A)') '# step node_a node_b new_node z d2 count_a count_b (merged by d2; nodes < '//&
            &int2str(n)//' are the leaves '//joined(self%names, ' ')//')'
        do i = 1, self%nsteps
            write(funit,'(A)') int2str(i)//' '//int2str(self%steps(1,i))//' '//int2str(self%steps(2,i))//' '//&
                &int2str(self%steps(3,i))//' '//fmtr(self%zs(i),2)//' '//fmtr(self%d2s(i),5)//' '//&
                &int2str(int(self%wsteps(1,i)))//' '//int2str(int(self%wsteps(2,i)))
        end do
        call fclose(funit)
        allocate(nwk(0:2*n-2), ht(0:2*n-2))
        ht = 0._dp
        do i = 0, n-1
            nwk(i)%s = int2str(self%names(i+1))
        end do
        do i = 1, self%nsteps
            p = self%steps(1,i)
            q = self%steps(2,i)
            r = self%steps(3,i)
            ht(r) = self%d2s(i)
            nwk(r)%s = '('//nwk(p)%s//':'//fmtr(max(ht(r)-ht(p),0._dp),4)//','//nwk(q)%s//':'//&
                &fmtr(max(ht(r)-ht(q),0._dp),4)//')'
        end do
        fname = 'tree.nwk'
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call fileiochk('state_tree: cannot write '//fname%to_char(), ios)
        if( self%nsteps > 0 ) write(funit,'(A)') nwk(self%steps(3,self%nsteps))%s//';'
        call fclose(funit)
        g = n
        do while( g > 1 )
            g = (g + 1) / 2
            call write_cut('cut'//int2str(g), 'cut '//int2str(g), n - g)
        end do
        if( self%nbudget >= 0 ) call write_cut('budget'//int2str(int(self%budget)), &
            &'budget'//int2str(int(self%budget)), self%nbudget)
        deallocate(zleaf, dleaf, nwk, ht)

    contains

        subroutine write_cut( dname, tag, nmerged )
            character(len=*), intent(in) :: dname, tag
            integer,          intent(in) :: nmerged
            integer,  allocatable :: alive(:), head(:), tail(:), nxt(:), order(:), lut(:)
            real(dp), allocatable :: gcount(:)
            type(string) :: cname
            character(len=:), allocatable :: line
            integer :: nalive, s, a, b, c, m, gi, cunit, ii, jj
            allocate(alive(n), head(0:2*n-2), tail(0:2*n-2), nxt(n), source=0)
            do a = 0, n-1
                alive(a+1) = a
                head(a)    = a + 1
                tail(a)    = a + 1
            end do
            nalive = n
            do s = 1, nmerged
                a = self%steps(1,s)
                b = self%steps(2,s)
                c = self%steps(3,s)
                head(c)       = head(a)
                nxt(tail(a))  = head(b)
                tail(c)       = tail(b)
                alive(1:nalive-2) = pack(alive(1:nalive), alive(1:nalive) /= a .and. alive(1:nalive) /= b)
                nalive = nalive - 1
                alive(nalive) = c
            end do
            allocate(gcount(nalive), order(nalive))
            do ii = 1, nalive
                gcount(ii) = 0._dp
                m = head(alive(ii))
                do while( m > 0 )
                    gcount(ii) = gcount(ii) + self%counts(m)
                    m = nxt(m)
                end do
                order(ii) = ii
            end do
            do ii = 2, nalive
                jj = ii
                do while( jj > 1 )
                    if( gcount(order(jj-1)) >= gcount(order(jj)) ) exit
                    order(jj-1:jj) = [order(jj), order(jj-1)]
                    jj = jj - 1
                end do
            end do
            call simple_mkdir(dname, verbose=.false.)
            cname = dname//'/groups.txt'
            call fopen(cunit, cname, status='replace', action='write', iostat=ios)
            call fileiochk('state_tree: cannot write '//cname%to_char(), ios)
            write(cunit,'(A)') '# bottom-up state ('//trim(self%linkage)//' linkage on d2, '//tag//&
                &'): particles, members (leaf states)'
            allocate(lut(0:max(1,maxval(self%names))), source=0)
            do gi = 1, nalive
                line = ''
                m = head(alive(order(gi)))
                do while( m > 0 )
                    if( len(line) > 0 ) line = line//'+'
                    line = line//int2str(self%names(m))
                    lut(self%names(m)) = gi
                    m = nxt(m)
                end do
                write(cunit,'(A)') int2str(gi)//' '//int2str(int(gcount(order(gi))))//' '//line
            end do
            call fclose(cunit)
            if( present(labels) )then
                cname = dname//'/states.txt'
                call fopen(cunit, cname, status='replace', action='write', iostat=ios)
                call fileiochk('state_tree: cannot write '//cname%to_char(), ios)
                do ii = 1, size(labels)
                    if( labels(ii) > 0 .and. labels(ii) <= ubound(lut,1) )then
                        write(cunit,'(I0)') lut(labels(ii))
                    else
                        write(cunit,'(I0)') 0
                    endif
                end do
                call fclose(cunit)
            endif
            deallocate(alive, head, tail, nxt, gcount, order, lut)
        end subroutine write_cut

        subroutine write_table( tname, vals, ndec )
            class(string), intent(in) :: tname
            real(dp),      intent(in) :: vals(:,:)
            integer,       intent(in) :: ndec
            character(len=:), allocatable :: line
            integer :: tunit, ii, jj
            call fopen(tunit, tname, status='replace', action='write', iostat=ios)
            call fileiochk('state_tree: cannot write '//tname%to_char(), ios)
            write(tunit,'(A)') '# '//joined(self%names, ' ')
            do ii = 1, size(vals,1)
                line = fmtr(vals(ii,1),ndec)
                do jj = 2, size(vals,2)
                    line = line//' '//fmtr(vals(ii,jj),ndec)
                end do
                write(tunit,'(A)') line
            end do
            call fclose(tunit)
        end subroutine write_table

    end subroutine write

    subroutine read( self, fname )
        class(state_tree), intent(inout) :: self
        class(string),     intent(in)    :: fname
        character(len=MAXLINE) :: line
        integer :: funit, ios, ipos, jpos, n, nrows, i, step
        real(dp) :: z, d, wa, wb
        call self%kill
        call fopen(funit, fname, status='old', action='read', iostat=ios)
        call fileiochk('state_tree: cannot read '//fname%to_char(), ios)
        read(funit,'(A)',iostat=ios) line
        ipos = index(line, 'are the leaves ')
        if( ios /= 0 .or. ipos == 0 ) THROW_HARD('no leaf names in the header of '//fname%to_char())
        ipos = ipos + len('are the leaves ')
        jpos = index(line(ipos:), ')') + ipos - 2
        n    = nwords(line(ipos:jpos))
        allocate(self%names(n), self%counts(n))
        read(line(ipos:jpos),*) self%names
        self%counts  = 0._dp
        self%nleaves = n
        nrows = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            nrows = nrows + 1
        end do
        if( nrows /= n - 1 ) THROW_HARD('incomplete merge tree in '//fname%to_char())
        allocate(self%steps(3,nrows), self%wsteps(2,nrows), self%zs(nrows), self%d2s(nrows))
        rewind(funit)
        read(funit,'(A)') line
        i = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            i = i + 1
            read(line,*) step, self%steps(:,i), z, d, wa, wb
            self%zs(i)       = z
            self%d2s(i)      = d
            self%wsteps(:,i) = [wa, wb]
        end do
        call fclose(funit)
        self%nsteps   = nrows
        self%l_exists = .true.
    end subroutine read

    subroutine cut( self, ngroups, grp )
        class(state_tree),    intent(in)  :: self
        integer,              intent(in)  :: ngroups
        integer, allocatable, intent(out) :: grp(:)
        integer, allocatable :: owner(:), nodes(:)
        logical, allocatable :: l_alive(:)
        integer :: n, s, i, gi, a, b, c
        n = self%nleaves
        if( ngroups < 1 .or. ngroups > n ) THROW_HARD('tree cut out of range')
        allocate(owner(n), l_alive(0:2*n-2))
        owner   = [(i, i=0,n-1)]
        l_alive = .false.
        l_alive(0:n-1) = .true.
        do s = 1, n - ngroups
            a = self%steps(1,s)
            b = self%steps(2,s)
            c = self%steps(3,s)
            where( owner == a .or. owner == b ) owner = c
            l_alive(a) = .false.
            l_alive(b) = .false.
            l_alive(c) = .true.
        end do
        nodes = pack([(i, i=0,2*n-2)], l_alive)
        allocate(grp(n))
        do gi = 1, size(nodes)
            where( owner == nodes(gi) ) grp = gi
        end do
        deallocate(owner, l_alive, nodes)
    end subroutine cut

    integer function get_nleaves( self )
        class(state_tree), intent(in) :: self
        get_nleaves = self%nleaves
    end function get_nleaves

    subroutine get_names( self, names )
        class(state_tree),    intent(in)  :: self
        integer, allocatable, intent(out) :: names(:)
        names = self%names
    end subroutine get_names

    subroutine get_steps( self, steps )
        class(state_tree),    intent(in)  :: self
        integer, allocatable, intent(out) :: steps(:,:)
        steps = self%steps(:,1:self%nsteps)
    end subroutine get_steps

    subroutine kill( self )
        class(state_tree), intent(inout) :: self
        if( allocated(self%names)  ) deallocate(self%names)
        if( allocated(self%steps)  ) deallocate(self%steps)
        if( allocated(self%counts) ) deallocate(self%counts)
        if( allocated(self%gram)   ) deallocate(self%gram)
        if( allocated(self%noisev) ) deallocate(self%noisev)
        if( allocated(self%nw)     ) deallocate(self%nw)
        if( allocated(self%zs)     ) deallocate(self%zs)
        if( allocated(self%d2s)    ) deallocate(self%d2s)
        if( allocated(self%wsteps) ) deallocate(self%wsteps)
        self%nleaves  = 0
        self%nsteps   = 0
        self%nbudget  = -1
        self%infl     = 1._dp
        self%l_exists = .false.
    end subroutine kill

    pure real(dp) function zval( caa, cbb, cab, cba, va, vb, nw, infl )
        real(dp), intent(in) :: caa, cbb, cab, cba, va(:), vb(:), nw(:), infl
        zval = (caa + cbb - cab - cba) / sqrt(sum(nw * (va + vb)**2) * infl)
    end function zval

    pure real(dp) function d2val( caa, cbb, cab, cba )
        real(dp), intent(in) :: caa, cbb, cab, cba
        d2val = (caa + cbb - cab - cba) / sqrt(max(caa, 1.e-30_dp) * max(cbb, 1.e-30_dp))
    end function d2val

    subroutine calc_envelope( avg, smpd, msk )
        real(dp),              intent(in)  :: avg(:,:,:)
        real,                  intent(in)  :: smpd
        real(dp), allocatable, intent(out) :: msk(:,:,:)
        type(image) :: img
        complex,  allocatable :: cmat(:,:,:)
        real,     allocatable :: sm(:,:,:)
        logical,  allocatable :: bin(:,:,:), grown(:,:,:)
        integer,  allocatable :: offs(:,:)
        real(dp) :: f2
        integer  :: n, lims(3,2), h, k, l, phys(3), i, j, m, o, noffs
        n = size(avg,1)
        call img%new([n,n,n], smpd)
        call img%set_rmat(real(avg), .false.)
        call img%fft
        cmat = img%get_cmat()
        lims = img%loop_lims(2)
        do l = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    f2   = real(h*h + k*k + l*l, dp) * (real(smpd,dp) / real(n,dp))**2
                    phys = img%comp_addr_phys(h,k,l)
                    cmat(phys(1),phys(2),phys(3)) = cmat(phys(1),phys(2),phys(3)) * real(exp(-2._dp * ENV_LP**2 * f2))
                end do
            end do
        end do
        call img%set_cmat(cmat)
        call img%ifft
        sm  = img%get_rmat()
        bin = sm > real(ENV_FRAC) * maxval(sm)
        call img%kill
        allocate(offs(3,(2*ENV_GROW+1)**3))
        noffs = 0
        do k = -ENV_GROW, ENV_GROW
            do j = -ENV_GROW, ENV_GROW
                do i = -ENV_GROW, ENV_GROW
                    if( i*i + j*j + k*k > ENV_GROW**2 ) cycle
                    noffs = noffs + 1
                    offs(:,noffs) = [i, j, k]
                end do
            end do
        end do
        allocate(grown(n,n,n), source=.false.)
        do m = 1, n
            do j = 1, n
                do i = 1, n
                    if( .not. bin(i,j,m) ) cycle
                    do o = 1, noffs
                        grown(wrap(i+offs(1,o)), wrap(j+offs(2,o)), wrap(m+offs(3,o))) = .true.
                    end do
                end do
            end do
        end do
        allocate(msk(n,n,n))
        msk = merge(1._dp, 0._dp, grown)
        do i = 1, 3
            call gauss_axis(msk, i)
        end do
        msk = min(max(msk, 0._dp), 1._dp)
        deallocate(cmat, sm, bin, grown, offs)

    contains

        pure integer function wrap( x )
            integer, intent(in) :: x
            wrap = modulo(x - 1, n) + 1
        end function wrap

        subroutine gauss_axis( v, ax )
            real(dp), intent(inout) :: v(:,:,:)
            integer,  intent(in)    :: ax
            real(dp), allocatable :: g(:), u(:,:,:)
            integer :: ii, jj, kk, d, reach
            reach = min(ENV_REACH, n/2)
            allocate(g(-reach:reach))
            do d = -reach, reach
                g(d) = exp(-real(d*d,dp) / (2._dp * ENV_SIGMA**2))
            end do
            g = g / sum([(exp(-real(min(d,n-d)**2,dp) / (2._dp * ENV_SIGMA**2)), d=0,n-1)])
            allocate(u(n,n,n), source=0._dp)
            !$omp parallel do collapse(2) default(shared) private(ii,jj,kk,d) schedule(static) proc_bind(close)
            do kk = 1, n
                do jj = 1, n
                    do ii = 1, n
                        do d = -reach, reach
                            select case(ax)
                                case(1)
                                    u(ii,jj,kk) = u(ii,jj,kk) + g(d) * v(wrap(ii+d),jj,kk)
                                case(2)
                                    u(ii,jj,kk) = u(ii,jj,kk) + g(d) * v(ii,wrap(jj+d),kk)
                                case DEFAULT
                                    u(ii,jj,kk) = u(ii,jj,kk) + g(d) * v(ii,jj,wrap(kk+d))
                            end select
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            v = u
            deallocate(g, u)
        end subroutine gauss_axis

    end subroutine calc_envelope

    pure integer function nwords( str )
        character(len=*), intent(in) :: str
        integer :: i
        logical :: l_in
        nwords = 0
        l_in   = .false.
        do i = 1, len_trim(str)
            if( str(i:i) /= ' ' )then
                if( .not. l_in ) nwords = nwords + 1
                l_in = .true.
            else
                l_in = .false.
            endif
        end do
    end function nwords

    function joined( vals, sep ) result( str )
        integer,          intent(in)  :: vals(:)
        character(len=*), intent(in)  :: sep
        character(len=:), allocatable :: str
        integer :: i
        str = ''
        do i = 1, size(vals)
            if( i > 1 ) str = str//sep
            str = str//int2str(vals(i))
        end do
    end function joined

    function fmtr( x, ndec ) result( str )
        real(dp), intent(in)          :: x
        integer,  intent(in)          :: ndec
        character(len=:), allocatable :: str
        character(len=40) :: buf
        character(len=16) :: fmt
        write(fmt,'(A,I0,A)') '(F40.', ndec, ')'
        write(buf,fmt) x
        str = trim(adjustl(buf))
        if( ndec == 0 ) str = str(1:len(str)-1)
    end function fmtr

end module simple_state_tree
