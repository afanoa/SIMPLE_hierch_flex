!@descr: normalized same-factor ancestral partition initialization on fixed leaf rows
module simple_flex_ancestral_init
use, intrinsic :: iso_fortran_env, only: real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: flex_ancestral_mix

contains

    !> No project I/O. Caller verifies physical image and factor/mask provenance before entry.
    !! Output columns sum over particles; memberships sum over targets. Labels are initial only.
    subroutine flex_ancestral_mix( leaf, ancestor, rows, ancestor_rows, eo, ancestor_eo, targets, ancestor_targets, &
        &training, uniforms, blend, factor, ancestor_factor, verified, memberships, columns, transition, prior, &
        &labels, status, message )
        real(real64),     intent(in)  :: leaf(:,:), ancestor(:,:), uniforms(:), blend
        integer,          intent(in)  :: rows(:), ancestor_rows(:), eo(:), ancestor_eo(:)
        integer,          intent(in)  :: targets(:), ancestor_targets(:)
        logical,          intent(in)  :: training(:), verified
        character(len=*), intent(in)  :: factor, ancestor_factor
        real(real64), allocatable, intent(out) :: memberships(:,:), columns(:,:), transition(:,:), prior(:)
        integer,      allocatable, intent(out) :: labels(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        real(real64), allocatable :: lr(:,:), ar(:,:), mapped(:,:), wl(:,:), wa(:,:), mix(:,:), prob(:,:), tr(:,:), pr(:)
        integer,      allocatable :: lab(:)
        real(real64) :: mass, cumulative
        integer :: n, k, ka, nt, i, s, t, label_index
        status = 1
        message = ''
        n = size(rows)
        k = size(leaf,2)
        ka = size(ancestor,2)
        if( n < 1 .or. k < 1 .or. ka < 1 )then
            message = 'empty ancestral initializer'; return
        endif
        if( size(leaf,1) /= n .or. size(ancestor,1) /= n .or. size(ancestor_rows) /= n .or. &
            &size(eo) /= n .or. size(ancestor_eo) /= n .or. size(training) /= n .or. size(uniforms) /= n .or. &
            &size(targets) /= k .or. size(ancestor_targets) /= ka )then
            message = 'ancestral initializer dimensions differ'; return
        endif
        if( .not. verified .or. len_trim(factor) == 0 .or. trim(factor) /= trim(ancestor_factor) )then
            message = 'same-factor provenance is not verified'; return
        endif
        if( any(rows < 1) .or. any(rows /= ancestor_rows) .or. any(eo /= ancestor_eo) .or. &
            &any(eo < 0) .or. any(eo > 1) )then
            message = 'original row or half identity differs'; return
        endif
        do i = 2, n
            if( rows(i) <= rows(i-1) )then
                message = 'original rows must be unique and ordered'; return
            endif
        enddo
        if( any(targets < 1) .or. any(targets > k) .or. any(ancestor_targets < 1) .or. any(ancestor_targets > ka) )then
            message = 'target identity outside contiguous state set'; return
        endif
        do s = 1, k
            if( count(targets == targets(s)) /= 1 )then
                message = 'duplicate leaf target identity'; return
            endif
        enddo
        do s = 1, ka
            if( count(ancestor_targets == ancestor_targets(s)) /= 1 )then
                message = 'duplicate ancestor target identity'; return
            endif
        enddo
        if( .not. all(ieee_is_finite(leaf)) .or. .not. all(ieee_is_finite(ancestor)) .or. &
            &.not. all(ieee_is_finite(uniforms)) .or. .not. ieee_is_finite(blend) )then
            message = 'non-finite ancestral initializer input'; return
        endif
        if( any(leaf < 0._real64) .or. any(ancestor < 0._real64) .or. &
            &any(leaf > 1.000001_real64) .or. any(ancestor > 1.000001_real64) .or. any(uniforms < 0._real64) .or. &
            &any(uniforms >= 1._real64) .or. blend < 0._real64 .or. blend > 1._real64 )then
            message = 'invalid partition, uniform draw or blend fraction'; return
        endif
        nt = count(training)
        if( nt == 0 )then
            message = 'no correspondence training rows'; return
        endif
        allocate(lr(n,k), ar(n,ka), tr(ka,k), pr(k), wl(n,k), wa(n,k), mix(n,k), prob(n,k), lab(n))
        do i = 1, n
            mass = sum(leaf(i,:))
            if( mass <= 0._real64 .or. abs(mass - 1._real64) > 1.e-6_real64 )then
                message = 'leaf row has zero or non-partition mass'; return
            endif
            lr(i,:) = leaf(i,:) / mass
            mass = sum(ancestor(i,:))
            if( mass <= 0._real64 .or. abs(mass - 1._real64) > 1.e-6_real64 )then
                message = 'ancestor row has zero or non-partition mass'; return
            endif
            ar(i,:) = ancestor(i,:) / mass
        enddo
        do t = 1, k
            pr(t) = sum(lr(:,t), mask=training) / real(nt,real64)
        enddo
        if( any(pr <= 0._real64) )then
            message = 'leaf target absent from correspondence training rows'; return
        endif
        do s = 1, ka
            do t = 1, k
                tr(s,t) = sum(ar(:,s) * lr(:,t), mask=training)
            enddo
            mass = sum(tr(s,:))
            if( mass > 0._real64 )then
                tr(s,:) = tr(s,:) / mass
            else
                tr(s,:) = pr
            endif
        enddo
        mapped = matmul(ar, tr)
        do t = 1, k
            mass = sum(lr(:,t))
            if( mass <= 0._real64 )then
                message = 'zero-support leaf target'; return
            endif
            wl(:,t) = lr(:,t) / mass
            mass = sum(mapped(:,t))
            if( mass <= 0._real64 )then
                message = 'zero-support mapped ancestor target'; return
            endif
            wa(:,t) = mapped(:,t) / mass
            mix(:,t) = (1._real64 - blend) * wl(:,t) + blend * wa(:,t)
            mass = sum(mix(:,t))
            if( mass <= 0._real64 )then
                message = 'zero-support mixed target'; return
            endif
            mix(:,t) = mix(:,t) / mass
        enddo
        do i = 1, n
            prob(i,:) = mix(i,:) * pr
            mass = sum(prob(i,:))
            if( mass <= 0._real64 )then
                message = 'zero-support mixed row'; return
            endif
            prob(i,:) = prob(i,:) / mass
            cumulative = 0._real64
            label_index = 1
            do t = 1, k - 1
                cumulative = cumulative + prob(i,t)
                if( uniforms(i) >= cumulative ) label_index = label_index + 1
            enddo
            lab(i) = targets(label_index)
        enddo
        call move_alloc(prob, memberships)
        call move_alloc(mix, columns)
        call move_alloc(tr, transition)
        call move_alloc(pr, prior)
        call move_alloc(lab, labels)
        status = 0
    end subroutine flex_ancestral_mix

end module simple_flex_ancestral_init
