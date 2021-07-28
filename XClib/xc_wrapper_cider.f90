!
SUBROUTINE xc_cider_x( length, ns, np, rho, grho, tau, feat, ex, v1x, v2x, v3x, vfeat )
    !
    USE kind_l,        ONLY: DP
    !
    IMPLICIT NONE
    !
    integer, intent(in) :: length
    integer, intent(in) :: ns
    integer, intent(in) :: np
    real(dp), intent(in) :: rho(length,ns)
    real(dp), intent(in) :: grho(3,length,ns)
    real(dp), intent(in) :: tau(length,ns)
    real(dp), intent(in) :: feat(length,ns,1)
    real(dp), intent(inout) :: ex(length)
    real(dp), intent(inout) :: v1x(length,ns)
    real(dp), intent(inout) :: v2x(length,ns)
    real(dp), intent(inout) :: v3x(length,ns)
    real(dp), intent(out) :: vfeat(length,ns,1)
    !
    integer :: is,k
    real(dp) :: x43,x13
    real(dp) :: xf
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
    xf = 0.01_DP
    vfeat = 0.0_DP
    rho_tot = 0.0_DP
    !
    ! TODO spin polarization requires different factors for nsp/sp
    do is = 1, ns
        rho_tot = rho_tot + rho(:,is)
    enddo
    do k = 1, length
        do is = 1, ns
            ex(k) = ex(k) + xf * feat(k,is,1) * rho(k,is)**x43 / (ABS(rho_tot(k)) + 1.0e-8_DP)
            v1x(k,is) = v1x(k,is) + xf * feat(k,is,1) * x43 * Abs(rho(k,is))**x13
            vfeat(k,is,1) = vfeat(k,is,1) + xf * rho(k,is)**x43
            if (isnan(v1x(k,is))) then
                print *, k, is
            endif
        enddo
    enddo
    !
    return
    !
END SUBROUTINE xc_cider_x
