!
SUBROUTINE xc_cider_x( length, nfeat, ns, np, rho, grho, tau, feat, ex, v1x, v2x, v3x, vfeat )
    !
    USE kind_l,        ONLY: DP
    !
    IMPLICIT NONE
    !
    integer, intent(in) :: length
    integer, intent(in) :: nfeat
    integer, intent(in) :: ns
    integer, intent(in) :: np
    real(dp), intent(in) :: rho(length,ns)
    real(dp), intent(in) :: grho(3,length,ns)
    real(dp), intent(in) :: tau(length,ns)
    real(dp), intent(in) :: feat(length,nfeat,ns)
    real(dp), intent(inout) :: ex(length)
    real(dp), intent(inout) :: v1x(length,ns)
    real(dp), intent(inout) :: v2x(length,ns)
    real(dp), intent(inout) :: v3x(length,ns)
    real(dp), intent(out) :: vfeat(length,nfeat,ns)
    !
    integer :: is,k,sgn
    real(dp) :: x43,x13
    real(dp) :: xf,rval
    real(dp), allocatable :: grho2(:,:)
    real(dp), allocatable :: rho_tot(:)
    !
    !polarized = .false.
    !if (ns == 2) then
    !    polarized = .true.
    !endif
    !
    !pol_unpol = ns
    !
    allocate( rho_tot(length) )
    x43 = 1.33333333333333333333333_DP
    x13 = 0.33333333333333333333333_DP
    xf = 0.1_DP !* SQRT(4 * 3.141592653)
    vfeat = 0.0_DP
    rho_tot = 0.0_DP
    !
    ! TODO spin polarization requires different factors for nsp/sp
    do is = 1, ns
        rho_tot = rho_tot + rho(:,is)
    enddo
    if (ns==1) then
        do k = 1, length
            sgn = SIGN(1.0_dp, rho(k,1))
            rval = abs(rho(k,1))
            ex(k) = ex(k) + 2 * xf * feat(k,1,1)/2.0_dp * (rval/2.0_dp)**x43 * sgn
            v1x(k,1) = v1x(k,1) + xf * feat(k,1,1)/2.0_dp * x43 * (rval/2.0_dp)**x13
            vfeat(k,1,1) = vfeat(k,1,1) + xf * (rval/2.0_dp)**x43
        enddo
    else
        do k = 1, length
            do is=1,ns
                sgn = SIGN(1.0_dp, rho(k,1))
                ex(k) = ex(k) + xf * feat(k,1,is) * abs(rho(k,is))**x43 * sgn
                v1x(k,is) = v1x(k,is) + xf * feat(k,is,1) * x43 * abs(rho(k,is))**x13
                vfeat(k,1,is) = vfeat(k,1,is) + xf * abs(rho(k,is))**x43
            enddo
        enddo
    endif
    !
    return
    !
END SUBROUTINE xc_cider_x
