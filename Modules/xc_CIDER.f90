
subroutine add_cider_gbas( rhog, feat, ibas1, ibas2)
    !
    !! CIDER features from density in reciprocal space
    !
    USE constants,         ONLY : fpi, e2, pi
    USE input_parameters,  ONLY : cider_params
    USE kinds,             ONLY : dp
    USE fft_base,          ONLY : dfftp
    USE gvect,             ONLY : ngm, gg
    USE cell_base,         ONLY : omega, tpiba2
    !
    implicit none
    !
    complex(dp), intent(in) :: rhog(ngm)
    real(dp), intent(inout) :: feat(ngm) ! TODO: should just be out?
    integer, intent(in) :: ibas1, ibas2
    !
    ! ... local variables
    !
    real(dp)                    :: fac 
    real(dp)                    :: rgtot_re, rgtot_im
    integer                     :: ig
    real(dp)                    :: aexp
    !
    call start_clock( 'add_cider_gbas' )
    !
    aexp = cider_params(1) / cider_params(2)**(ibas1-1) + cider_params(1) / cider_params(2)**(ibas2-1)
    !
!$omp parallel do private( fac, rgtot_re, rgtot_im )
    do ig = 1, ngm
        !
        fac = (pi/aexp)**1.5 * exp(-gg(ig) * tpiba2 / (4.D0 * aexp))
        !rgtot_re = real(  rhog(ig) )
        !rgtot_im = aimag( rhog(ig) )
        feat(ig) = feat(ig) + fac * rhog(ig)
        !
    enddo
!$omp end parallel do
    !
    call stop_clock( 'add_cider_gbas' )
    !
    return
    !
end subroutine add_cider_gbas
!
subroutine get_cider_coefs( length, cider_exp, fc, dfc )
    !
    ! Takes CIDER auxiliary basis features as input, along with
    ! CIDER exponents, and returns the coefficnets for constructing
    ! CIDER features from the basis functions
    !
    use kind_l,             only: dp
    use input_parameters,   only: cider_params, cider_nbas
    use constants,          only: pi
    !
    implicit none
    !
    integer, intent(in)     :: length
    real(dp), intent(in)    :: cider_exp(length)
    real(dp), intent(in)    :: fc(length,cider_nbas)
    real(dp), intent(out)   :: dfc(length,cider_nbas)
    !
    integer :: i,j
    real(dp), allocatable :: bas_exp(:)
    real(dp), allocatable :: sinv(:,:)
    real(dp), allocatable :: cvec(:,:)
    real(dp), allocatable :: dcvec(:,:)
    real(dp)              :: cider_maxexp, cider_xexp, power
    !
    allocate ( bas_exp(cider_nbas) )
    allocate ( sinv(cider_nbas, cider_nbas) )
    allocate ( cvec(length, cider_nbas) )
    allocate ( dcvec(length, cider_nbas) )
    cider_maxexp = cider_params(1)
    cider_xexp = cider_params(2)
    !
    power = -1.5_dp
    !
    do i=1,cider_nbas
        bas_exp(i) = cider_maxexp / cider_xexp**(i-1)
    enddo
    do i=1,cider_nbas
        do j=1,cider_nbas
            sinv(i,j) = (bas_exp(i) + bas_exp(j))**power
        enddo
    enddo
    !
    CALL MatInv('G', cider_nbas, sinv)
    ! if not G, need to fill matrix into square after calling above subroutine
    !
    do i=1,length
        do j=1,cider_nbas
            cvec(i,j) = (cider_exp(i) + bas_exp(j))**power
            dcvec(i,j) = power * cvec(i,j) / (cider_exp(i) + bas_exp(j))
        enddo
    enddo
    !
    fc = matmul(cvec, sinv) ! d_{i,sigma}
    dfc = matmul(dcvec, sinv) ! deriv wrt exponent
    !
    deallocate(bas_exp, sinv, cvec, dcvec)
    !
    return
    !
end subroutine get_cider_coefs
!
!
SUBROUTINE get_cider_alpha( length, rho, kin, cider_consts, cider_exp )
  !) rho%of_r, grho, rho%kin_r/e2, cider_exp(:,ialpha,1) )
  USE kind_l, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE constants, ONLY : pi
  IMPLICIT NONE
  
  integer,  intent(in) :: length
  real(dp), intent(in) :: rho(length)
  real(dp), intent(in) :: grho(length)
  real(dp), intent(in) :: kin(length)
  real(dp), intent(in) :: cider_consts(4)
  real(dp), intent(out) :: cider_exp(length)

  real(dp) :: const1, const2, fac
  integer :: i
  integer :: tot

  fac = cider_consts(3) * 1.2 * (6 * pi**2)**(2.0/3) / pi
  const1 = cider_consts(2) - fac
  const2 = fac / (0.3 * (3 * pi**2)**(2.0_DP/3))
  IF (nspin==1) THEN
    const1 = const1 * pi / 2.0_DP**(2.0_DP/3)
  ELSE
    const1 = const1 * pi
  ENDIF
  const2 = const2 * pi / 2**(2.0_DP/3)

  cider_exp = cider_consts(1) + const1 * (abs(rho))**(2.0_DP/3) &
              + const2 * abs(kin) / (abs(rho)+1e-16)

  do i=1,length
    if (cider_exp(i) < cider_consts(4)) then
      cider_exp(i) = cider_consts(4)*exp(cider_exp(i)/cider_consts(4)-1.0_DP)
    endif
  enddo

  return
END SUBROUTINE get_cider_alpha
!
SUBROUTINE get_cider_lpot_exp ( length, rho, kin, cider_consts, &
                                cider_exp, vexp, v1x, v2x, v3x )
  USE kind_l, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE constants, ONLY : pi
  IMPLICIT NONE

  integer,  intent(in) :: length
  real(dp), intent(in) :: rho(length)
  real(dp), intent(in) :: grho(length)
  real(dp), intent(in) :: kin(length)
  real(dp), intent(in) :: cider_consts(4)
  real(dp), intent(in) :: cider_exp(length)
  real(dp), intent(inout) :: vexp(length)
  real(dp), intent(inout) :: v1x(length)
  real(dp), intent(inout) :: v2x(length)
  real(dp), intent(inout) :: v3x(length)

  real(dp) :: const1, const2, fac
  integer :: i
  real(dp), allocatable :: cider_tmp(:)
  integer :: tot

  fac = cider_consts(3) * 1.2 * (6 * pi**2)**(2.0/3) / pi
  const1 = cider_consts(2) - fac
  const2 = fac / (0.3 * (3 * pi**2)**(2.0_DP/3))
  IF (nspin==1) THEN
    const1 = const1 * pi / 2.0_DP**(2.0_DP/3)
  ELSE
    const1 = const1 * pi
  ENDIF
  const2 = const2 * pi / 2**(2.0_DP/3)

  allocate(cider_tmp(length))
  cider_tmp = cider_consts(1) + const1 * (abs(rho))**(2.0_DP/3) &
              + const2 * abs(kin) / (abs(rho)+1e-16)

  do i=1,length
    if (cider_tmp(i) < cider_consts(4)) then
      vexp(i) = vexp(i) / cider_consts(4) * cider_exp(i)
    endif
  enddo

  v1x = v1x + 2.0_DP/3.0_DP * vexp &
              * const1 / (abs(rho)**(1.0_DP/3)+1e-16)
  v1x = v1x - vexp * const2 * abs(kin) / (abs(rho)**2+1e-16)
  v3x = v3x + vexp * const2 / (abs(rho)+1e-16)

  deallocate(cider_tmp)

  return
END SUBROUTINE get_cider_lpot_exp
!
subroutine v_xc_cider( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
    !
    !! Non-local XC potential using CIDER descriptors
    !
    use kinds,              only : dp
    use constants,          only : e2, eps8, pi
    use input_parameters,   only : cider_consts, cider_nalpha, cider_nfeat, &
                                   cider_nbas, a_list, cider_params, &
                                   cider_nset, cider_ls, ialphas, isets
    use io_global,          only : stdout
    use fft_base,           only : dfftp
    use gvect,              only : g, gg, ngm
    use lsda_mod,           only : nspin
    use cell_base,          only : omega
    use funct,              only : dft_is_nonlocc, nlc
    use scf,                only : scf_type, rhoz_or_updw
    use mp,                 only : mp_sum
    use mp_bands,           only : intra_bgrp_comm
    !
    implicit none
    !
    type (scf_type), intent(inout) :: rho
    !! the valence charge
    REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    !! the core charge in real space
    COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    !! the core charge in reciprocal space
    REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
    !! V_xc potential
    REAL(DP), INTENT(INOUT) :: kedtaur(dfftp%nnr,nspin)
    !! local K energy density                     
    REAL(DP), INTENT(INOUT) :: vtxc
    !! integral V_xc * rho
    REAL(DP), INTENT(INOUT) :: etxc
    !! E_xc energy
    LOGICAL, INTENT(IN) :: is_meta
    !
    ! ... local variables
    !
    real(dp) :: zeta, rh, sgn(2)
    integer  :: k, ipol, is, np
    !
    real(dp), allocatable : ex(:) ec(:)
    real(dp), allocatable : v1x(:,:), v2x(:,:), v3x(:,:)
    real(dp), allocatable : v1c(:,:), v2c(:,:,:), v3c(:,:)
    !
    ! CIDER arrays
    !
    !real(dp), allocatable :: const(:,:,:)
    real(dp), allocatable :: cider_exp(:,:,:) ! nnr, nalpha, nspin
    real(dp), allocatable :: bas(:,:,:)   ! nnr, nbas, nspin
    real(dp), allocatable :: gbas(:,:,:)  ! 2, ngm, nbas, nspin
    real(dp), allocatable :: dbas(:,:,:)  ! nnr, 3, nbas, nspin
    real(dp), allocatable :: vbas(:,:,:)  ! nnr, nbas, nspin
    real(dp), allocatable :: feat(:,:,:)  ! nnr, nfeat, nspin
    real(dp), allocatable :: dfeat(:,:,:) ! nnr, nfeat, nspin; deriv wrt cider_exp
    real(dp), allocatable :: vfeat(:,:,:) ! nnr, nfeat, nspin
    real(dp), allocatable :: vexp(:,:,:)  ! nnr, nalpha, nspin
    real(dp), allocatable :: fc(:,:,:,:)  ! nnr, nbas, napha, nspin
    real(dp), allocatable :: dfc(:,:,:,:) ! nnr, nbas, napha, nspin
    real(dp), allocatable :: wqrho(:,:,:) ! nnr, nbas, nspin
    ! nalpha should be 4: 3 for l=0 scales and 1 for core damper
    ! nfeat should be 3 or 6, depending on how l=1 feat is stored, 6 probably better
    ! nbas should be for 1 set of exponents??? TODO
    integer :: ibas, ibas1, ibas2
    !
    ! initialize timer and meta-gga stuff
    !
    CALL start_clock( 'v_xc_cider' )
    !
    etxc = zero
    vtxc = zero
    v(:,:) = zero
    rhoneg(:) = zero
    sgn(1) = 1._dp  ;   sgn(2) = -1._dp
    fac = 1.D0 / DBLE( nspin )
    np = 1
    IF (nspin==2) np=3
    !
    ALLOCATE( grho(3,dfftp%nnr,nspin) )
    ALLOCATE( h(3,dfftp%nnr,nspin) )
    ALLOCATE( rhogsum(ngm) )
    !
    ALLOCATE( ex(dfftp%nnr), ec(dfftp%nnr) )
    ALLOCATE( v1x(dfftp%nnr,nspin), v2x(dfftp%nnr,nspin)   , v3x(dfftp%nnr,nspin) )
    ALLOCATE( v1c(dfftp%nnr,nspin), v2c(np,dfftp%nnr,nspin), v3c(dfftp%nnr,nspin) )
    !
    ! initialize CIDER components
    !
    allocate( feat(dfftp%nnr,cider_nfeat,nspin) )
    allocate( dfeat(dfftp%nnr,cider_nfeat,nspin) )
    !
    vexp(:,:,:) = zero
    vbas(:,:,:) = zero
    vfeat(:,:,:) = zero
    h(:,:,:) = zero
    !
    cider_nfeat0 = cider_nalpha - 1
    !
    if (nspin == 2) then
        call rhoz_or_updw( rho, 'both', '->updw' )
    endif
    !
    do is=1, nspin
        rhogsum(:) = fac*rhog_core(:) + rho%of_g(:,is) ! rhoz is transformed to updw above
        call fft_gradient_g2r(dfftp, rhogsum, g, grho(1,1,is))
        do ialpha=1,cider_nalpha
            call get_cider_alpha( dfftp%nnr, rho%of_r(:,is), grho(:,:,is), rho%kin_r(:,is)/e2, &
                                  cider_consts(:,ialpha), cider_exp(:,ialpha,is) )
            call get_cider_coefs( dfftp%nnr, cider_exp(:,ialpha,is), 
                                  fc(:,:,ialpha,is), dfc(:,:,ialpha,is) )
        enddo
        do ibas1=1,cider_nbas
            rhogsum(:) = zero
            psi = rho%of_r(:,is) * fc(:,ibas1,cider_nalpha,is) ! wq * rho
            call fwfft( 'Rho', psi, dfftp )
            call fftx_threed2oned( dfftp, psi, rhogsum )
            do ibas2=1,cider_nbas
                call add_cider_gbas( rhogsum, gbas(:,ibas2,is), ibas1, ibas2 )
            enddo
        enddo
        do ibas1=1,cider_nbas
            call rho_g2r( dfftp, gbas(:,ibas1,is), bas(:,ibas1,is) )
            if ( cider_uses_grad ) then
                call fft_gradient_g2r( dfftp, gbas(:,ibas1,is), dbas(:,:,ibas1,is) )
            enddo
        enddo
        const = cider_consts(2,ialpha)**1.5_dp
        const = const * cider_consts(1,ialpha) / (cider_consts(1,ialpha) + const)
        do ialpha=1,cider_nfeat0
            feat(:,ialpha,is)  = const * sum(  fc(:,:,ialpha,is) * bas(:,:,is), 2 )
            dfeat(:,ialpha,is) = const * sum( dfc(:,:,ialpha,is) * bas(:,:,is), 2 )
        enddo
        l1c = sqrt(3/pi) * (8*pi/3)**(1.0_dp/3)
        if ( cider_uses_grad ) then
            const = cider_consts(2,ialpha)**1.5_dp * l1c
            const = const * cider_consts(1,ialpha) / (cider_consts(1,ialpha) + const)
            ialpha = ialpha_grad
            do rindex=1,3
                ifeat = cider_nfeat0 + rindex
                feat(:,ifeat,is)  = const * sum(  fc(:,:,ialpha,is) * dbas(:,rindex,:,is), 2 )
                dfeat(:,ifeat,is) = const * sum( dfc(:,:,ialpha,is) * dbas(:,rindex,:,is), 2 )
            enddo
        endif
    enddo
    !
    ! ... call the semi-local DFT part and scale exchange by mixing factor, then call CIDER
    !
    call xc_metagcx( dfftp%nnr, nspin, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c )
    ex = ex * (1 - cider_params(3))
    v1x = v1x * (1 - cider_params(3))
    v2x = v2x * (1 - cider_params(3))
    v3x = v3x * (1 - cider_params(3))
    call xc_cider_x_py( dfftp%nnr, cider_nfeat, nspin, np, rho%of_r, grho, &
                        rho%kin_r/e2, feat, ex, v1x, v2x, v3x, vfeat, h )
    !
    ! ... backpropagate nonlocal functional derivatives
    !
    do is=1, nspin
        !
        ! ... contract vfeat into vexp and vbas, vdbas
        !
        const = cider_consts(2,ialpha)**1.5_dp
        const = const * cider_consts(1,ialpha) / (cider_consts(1,ialpha) + const)
        do ialpha=1,cider_nfeat0
            vexp(:,ialpha,is) = vexp(:,ialpha,is) + vfeat(:,ialpha,is) * dfeat(:,ialpha,is)
            do ibas=1,cider_nbas
                vbas(:,ibas,is) = vbas(:,ibas,is) + vfeat(:,ialpha,is) * const * fc(:,ibas,ialpha,is)
            enddo
        enddo
        if (cider_uses_grad) then
            const = cider_consts(2,ialpha)**1.5_dp * l1c
            const = const * cider_consts(1,ialpha) / (cider_consts(1,ialpha) + const)
            do rindex=1,3
                ifeat = cider_nfeat0 + rindex
                ialpha = ialpha_grad
                vexp(:,ialpha,is) = vexp(:,ialpha,is) + vfeat(:,ifeat,is) * dfeat(:,ifeat,is)
                do ibas=1,cider_nbas
                    vdbas(:,:,ibas,is) = vdbas(:,:,ibas,is) + vfeat(:,ifeat,is) * const * fc(:,ibas,ialpha,is)
                enddo
            enddo
            ! use fft_graddot to transform vdbas into vbas
            do ibas=1,cider_nbas
                vbas_tmp(:) = zero
                call fft_graddot( dfftp, vdbas(:,:,ibas,is), g, vbas_tmp )
                vbas(:,ibas,is) = vbas(:,ibas,is) - vbas_tmp
            enddo
        endif
        ! transform rho back to deriv wrt xq
        gbas(:,:,:) = zero
        do ibas1=1,cider_nbas
            psi(:) = vbas(:,ibas1,is) ! wq * rho
            call fwfft( 'Rho', psi, dfftp )
            call fftx_threed2oned( dfftp, psi, rhogsum )
            do ibas2=1,cider_nbas
                call add_cider_gbas( rhogsum, gbas(:,ibas2,is), ibas1, ibas2 )
            enddo
        enddo
        bas(:,:,:) = zero
        do ibas1=1,cider_nbas
            vbas_tmp(:) = zero
            ! psi = rho%of_r(:,is) * fc(:,ibas1,nalpha,is)
            call rho_g2r( dfftp, gbas(:,ibas1,is), vbas_tmp )
            vexp(:,cider_nalpha,is) = vexp(:,cider_nalpha,is) + vbas_tmp * rho%of_r(:,is) * dfc(:,ibas1,cider_nalpha,is)
            v1x(:,is) = v1x(:,is) + vbas_tmp * fc(:,ibas1,cider_nalpha,is)
        enddo
        !
        ! ... Finally, now that all exponent constribs are computed, we need to 
        ! ... backpropagate into the v1x, v2x, v3x.
        !
        do ialpha=1,cider_nalpha
            call get_cider_lpot_exp( dfftp%nnr, rho%of_r(:,is), rho%kin_r(:,is)/e2, &
                                     cider_consts(:,ialpha), cider_exp(:,ialpha,is), &
                                     vexp(:,ialpha,is), v1x(:,is), v2x(:,is), v3x(:,is) )
        enddo
    enddo
    !
    ! ... Now do the rest of the meta-GGA stuff.
    ! ... TODO: The code below is the same as for meta-GGA, maybe it can
    ! ... be moved to another function somewhere for reproducibility?
    !
    IF (nspin == 1) THEN
        DO k = 1, dfftp%nnr
            !
            v(k,1) = (v1x(k,1)+v1c(k,1)) * e2
            !
            ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
            (:,k,1) = h(:,k,1) * e2 + (v2x(k,1)+v2c(1,k,1)) * grho(:,k,1) * e2 
            !
            kedtaur(k,1) = (v3x(k,1)+v3c(k,1)) * 0.5d0 * e2
            !
            etxc = etxc + (ex(k)+ec(k)) * e2
            vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * e2 * ABS(rho%of_r(k,1))
            !
            IF (rho%of_r(k,1) < zero) rhoneg(1) = rhoneg(1)-rho%of_r(k,1)
            !
        ENDDO
        !
    ELSE
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        DO k = 1, dfftp%nnr
            !
            v(k,1) = (v1x(k,1) + v1c(k,1)) * e2
            v(k,2) = (v1x(k,2) + v1c(k,2)) * e2
            !
            ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
            !
            h(:,k,1) = h(:,k,1) * e2 + (v2x(k,1) * grho(:,k,1) + v2c(:,k,1)) * e2
            h(:,k,2) = h(:,k,2) * e2 + (v2x(k,2) * grho(:,k,2) + v2c(:,k,2)) * e2
            !
            kedtaur(k,1) = (v3x(k,1) + v3c(k,1)) * 0.5d0 * e2
            kedtaur(k,2) = (v3x(k,2) + v3c(k,2)) * 0.5d0 * e2
            !
            etxc = etxc + (ex(k)+ec(k)) * e2
            vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * ABS(rho%of_r(k,1)) * e2
            vtxc = vtxc + (v1x(k,2)+v1c(k,2)) * ABS(rho%of_r(k,2)) * e2
            !
            IF ( rho%of_r(k,1) < 0.d0 ) rhoneg(1) = rhoneg(1) - rho%of_r(k,1)
            IF ( rho%of_r(k,2) < 0.d0 ) rhoneg(2) = rhoneg(2) - rho%of_r(k,2)
            !
        ENDDO
        !
        CALL rhoz_or_updw( rho, 'both', '->rhoz' ) ! transform rho back
        !
    ENDIF
    !
    DEALLOCATE( ex, ec )
    DEALLOCATE( v1x, v2x, v3x )
    DEALLOCATE( v1c, v2c, v3c )
    ! TODO fix deallocation statements for new version
    DEALLOCATE( bas, vbas, vexp, fc, dfc, ylm, cider_exp, const )
    DEALLOCATE( feat, dfeat, vfeat )
    !
    ALLOCATE( dh( dfftp%nnr ) )    
    !
    ! ... second term of the gradient correction :
    ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
    !
    ALLOCATE (rhoout(dfftp%nnr))
    DO is = 1, nspin
        !
        CALL fft_graddot( dfftp, h(1,1,is), g, dh )
        !
        v(:,is) = v(:,is) - dh(:)
        !
        ! ... rhoout is in (up,down) format 
        !
        rhoout(:) = ( rho%of_r(:,1) + sgn(is)*rho%of_r(:,nspin) )*0.5D0
        vtxc = vtxc - SUM( dh(:) * rhoout(:) )
        !
    END DO
    DEALLOCATE(rhoout)
    DEALLOCATE(dh)
    !
    CALL mp_sum( rhoneg, intra_bgrp_comm )
    !
    rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    !
    IF ((rhoneg(1) > eps8) .OR. (rhoneg(2) > eps8)) THEN
        write (stdout, '(/,5x, "negative rho (up,down): ", 2es10.3)') rhoneg(:)
    ENDIF
    !
    vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) 
    etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    !
    IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
    !
    CALL mp_sum(  vtxc , intra_bgrp_comm )
    CALL mp_sum(  etxc , intra_bgrp_comm )
    !
    DEALLOCATE(grho)
    DEALLOCATE(h)
    !
    CALL stop_clock( 'v_xc_cider' )
    !
    RETURN
    !
end subroutine v_xc_cider



