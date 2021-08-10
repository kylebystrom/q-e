!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
  !----------------------------------------------------------------------------
  !! This routine computes the Hartree and Exchange and Correlation
  !! potential and energies which corresponds to a given charge density
  !! The XC potential is computed in real space, while the
  !! Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE ions_base,        ONLY : nat, tau
  USE ldaU,             ONLY : lda_plus_u, lda_plus_u_kind, ldmx_b, &
                               nsg, v_nsg 
  USE xc_lib,           ONLY : xclib_dft_is
  USE dft_mod,          ONLY : xclib_get_id
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : ts_vdw, mbd_vdw
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
  USE libmbd_interface, ONLY : mbd_interface
  !
  USE lsda_mod, ONLY : nspin
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v
  !! the scf (Hxc) potential 
  !=================> NB: NOTE that in F90 derived data type must be INOUT and 
  !=================> not just OUT because otherwise their allocatable or pointer
  !=================> components are NOT defined 
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc
  !! the integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! the E_xc energy
  REAL(DP), INTENT(OUT) :: ehart
  !! the hartree energy
  REAL(DP), INTENT(OUT) :: eth
  !! the hubbard energy
  REAL(DP), INTENT(OUT) :: charge
  !! the integral of the charge
  REAL(DP) :: eth1
  !! the hubbard energy coming from the background states
  REAL(DP), INTENT(INOUT) :: etotefield
  !! electric field energy - inout due to the screwed logic of add_efield
  !
  INTEGER :: is, ir
  !
  !TEST variables
  real(dp), allocatable :: test_vr(:,:)
  real(dp), allocatable :: test_vk(:,:)
  real(dp) :: etxc2
  allocate(test_vr(dfftp%nnr,nspin))
  allocate(test_vk(dfftp%nnr,nspin))
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  IF (xclib_get_id('LDA', 'EXCH') == -1) then
     CALL v_xc_cider( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
  ELSEIF (xclib_dft_is('meta')) then
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
  ELSE
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  ENDIF
  !
  ! ... add a magnetic field  (if any)
  !
  CALL add_bfield( v%of_r, rho%of_r )
  !
  ! ... calculate hartree potential
  !
  CALL v_h( rho%of_g(:,1), ehart, charge, v%of_r )
  !
  ! ... DFT+U(+V): build up (extended) Hubbard potential 
  !
  IF (lda_plus_u) THEN
     !
     IF (lda_plus_u_kind == 0) THEN
        !
        ! DFT+U (simplified)
        !
        CALL v_hubbard (rho%ns, v%ns, eth)
        !
        ! Background
        IF (ldmx_b.GT.0) THEN
           CALL v_hubbard_b (rho%nsb, v%nsb, eth1)
           eth = eth + eth1
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 1) THEN
        !
        ! DFT+U (full)
        !
        IF (noncolin) THEN
           CALL v_hubbard_full_nc (rho%ns_nc, v%ns_nc, eth)
        ELSE
           CALL v_hubbard_full (rho%ns, v%ns, eth)
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 2) THEN
        !
        ! DFT+U+V (simplified)
        !
        CALL v_hubbard_extended (nsg, v_nsg, eth)
        !
     ELSE
        !
        CALL errore('v_of_rho', 'Not allowed value of lda_plus_u_kind',1)
        !
     ENDIF
     !
  ENDIF
  !
  ! ... add an electric field
  ! 
  DO is = 1, nspin_lsda
     CALL add_efield(v%of_r(1,is), etotefield, rho%of_r(:,1), .false. )
  END DO
  !
  ! ... add Tkatchenko-Scheffler potential (factor 2: Ha -> Ry)
  !
  IF (ts_vdw .or. mbd_vdw) THEN
     CALL tsvdw_calculate(tau*alat,rho%of_r(:,1))
     DO is = 1, nspin_lsda
        DO ir=1,dfftp%nnr
           v%of_r(ir,is)=v%of_r(ir,is)+2.0d0*UtsvdW(ir)
        END DO
     END DO
  END IF
  !
  IF (mbd_vdw) THEN
    call mbd_interface() ! self-consistent but only up to TS level
  END IF
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
!
!----------------------------------------------------------------------------
SUBROUTINE v_xc_cider( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur)
  !--------------------------------------------------------------------------
  !! Non-local XC potential using CIDER descriptors
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8, pi
  USE input_parameters, ONLY : cider_consts, cider_nalpha, cider_nfeat, &
                               cider_lmax, cider_nl, cider_nbas, &
                               a_list, l_list, lm_list
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g, gg, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE funct,            ONLY : dft_is_nonlocc, nlc
  USE scf,              ONLY : scf_type, rhoz_or_updw
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
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
  !
  ! ... local variables
  !
  REAL(DP) :: zeta, rh, sgn(2)
  INTEGER  :: k, ipol, is, np
  !
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:,:), v3c(:,:)
  !
  ! CIDER arrays
  REAL(DP), ALLOCATABLE :: bas(:,:,:,:)
  REAL(DP), ALLOCATABLE :: feat(:,:,:)
  REAL(DP), ALLOCATABLE :: dfeat(:,:,:)
  REAL(DP), ALLOCATABLE :: vfeat(:,:,:)
  REAL(DP), ALLOCATABLE :: vbas(:,:,:,:)
  REAL(DP), ALLOCATABLE :: vexp(:,:,:)
  REAL(DP), ALLOCATABLE :: ylm(:,:)
  REAL(DP), ALLOCATABLE :: fc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: fc2(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dfc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: cider_exp(:,:,:)
  REAL(DP), ALLOCATABLE :: const(:,:,:)
  !
  REAL(DP) :: fac
       
  REAL(DP), DIMENSION(2) :: grho2, rhoneg
  REAL(DP), DIMENSION(3) :: grhoup, grhodw
  !
  REAL(DP), ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: rhoout(:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:)
  REAL(DP), PARAMETER :: eps12 = 1.0d-12, zero=0._dp
  !
  INTEGER :: ialpha, ifeat, l, lmind, ibas, m
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
  ALLOCATE( feat(dfftp%nnr,cider_nfeat,nspin) )
  ALLOCATE( dfeat(dfftp%nnr,cider_nfeat,nspin) )
  ALLOCATE( vfeat(dfftp%nnr,cider_nfeat,nspin) )
  ALLOCATE( bas(dfftp%nnr,cider_nbas,cider_nl,nspin) )
  ALLOCATE( vbas(dfftp%nnr,cider_nbas,cider_nl,nspin) )
  ALLOCATE( vexp(dfftp%nnr,cider_nalpha,nspin) )
  ALLOCATE( fc(dfftp%nnr,cider_nbas,cider_nalpha,nspin) )
  ALLOCATE( fc2(dfftp%nnr,cider_nbas,cider_nalpha,nspin) )
  ALLOCATE( dfc(dfftp%nnr,cider_nbas,cider_nalpha,nspin) )
  ALLOCATE( ylm(ngm,cider_nl) )
  ALLOCATE( cider_exp(dfftp%nnr,cider_nalpha,nspin) )
  ALLOCATE( const(dfftp%nnr,cider_nfeat,nspin) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  ! ... in LSDA case rhogsum is in (up,down) format
  !
  CALL ylmr2(cider_nl, ngm, g, gg, ylm)
  !
  vexp(:,:,:) = zero
  vbas(:,:,:,:) = zero
  vfeat(:,:,:) = zero
  !
  IF (nspin == 2) THEN
    CALL rhoz_or_updw( rho, 'both', '->updw' )
  ENDIF
  !
  DO is = 1, nspin
    !
    rhogsum(:) = fac*rhog_core(:) + rho%of_g(:,is) ! rhoz is transformed to updw above
    !
    CALL fft_gradient_g2r( dfftp, rhogsum, g, grho(1,1,is) )
    !
    do ibas=1,cider_nbas
      lmind=1
      do l=0,cider_lmax
        do m=-l,l
          CALL get_cider_bas( rhogsum, bas(:,ibas,lmind,is), ylm(:,lmind), l, ibas )
          lmind=lmind+1
        enddo
      enddo
    enddo
    !
    do ialpha=1,cider_nalpha
      CALL get_cider_alpha( dfftp%nnr, rho%of_r(:,is), grho(:,:,is), rho%kin_r(:,is)/e2, &
                            cider_consts(:,ialpha), cider_exp(:,ialpha,is) )
      CALL get_cider_coefs( dfftp%nnr, cider_exp(:,ialpha,is), &
                            fc(:,:,ialpha,is), dfc(:,:,ialpha,is) )
      !CALL get_cider_coefs( dfftp%nnr, cider_exp(:,ialpha,is)+1e-8, &
      !                      fc2(:,:,ialpha,is), dfc(:,:,ialpha,is) )
    enddo
    do ifeat=1,cider_nfeat
      l = l_list(ifeat)
      lmind = lm_list(ifeat)
      ialpha = a_list(ifeat)
      const(:,ifeat,is) = cider_consts(2,ialpha)**1.5_DP &
                          * SQRT(4*pi**(1-l))  &
                          * (8*pi/3)**(l/3.0) * cider_exp(:,ialpha,is)**(l/2.0)
      const(:,ifeat,is) = cider_consts(1,ialpha) / (cider_consts(1,ialpha) + const(:,ifeat,is)) &
                          + const(:,ifeat,is)
      feat(:,ifeat,is) = SUM(fc(:,:,ialpha,is) * bas(:,:,lmind,is), 2)
      feat(:,ifeat,is) = feat(:,ifeat,is) * const(:,ifeat,is)
      dfeat(:,ifeat,is) = SUM(dfc(:,:,ialpha,is) * bas(:,:,lmind,is), 2) * const(:,ifeat,is)
      dfeat(:,ifeat,is) = dfeat(:,ifeat,is) + feat(:,ifeat,is) &
                          * DBLE(l) / (2.0 * cider_exp(:,ialpha,is))

      !dfeat(:,ifeat,is) = SUM(fc2(:,:,ialpha,is) * bas(:,:,lmind,is), 2)
      !dfeat(:,ifeat,is) = dfeat(:,ifeat,is) * const(:,ifeat,is)
      !dfeat(:,ifeat,is) = (dfeat(:,ifeat,is) - feat(:,ifeat,is)) / 1e-8
      ! TODO numerical derivative seems more stable for some reason
      ! and so might want to change back to it
    enddo
    !
  ENDDO
  !
  DEALLOCATE(rhogsum)
  !
  ! TODO problems for nspin>2
  CALL xc_metagcx( dfftp%nnr, nspin, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c )
  CALL xc_cider_x_py( dfftp%nnr, cider_nfeat, nspin, np, rho%of_r, grho, &
                      rho%kin_r/e2, feat, ex, &
                      v1x, v2x, v3x, vfeat )
  !
  DO is=1,nspin
    !
    do ifeat=1,cider_nfeat
      lmind = lm_list(ifeat)
      l = l_list(ifeat)
      ialpha = a_list(ifeat)
      vexp(:,ialpha,is) = vexp(:,ialpha,is) + vfeat(:,ifeat,is) * dfeat(:,ifeat,is)
      do ibas=1,cider_nbas
        vbas(:,ibas,lmind,is) = vbas(:,ibas,lmind,is) + const(:,ifeat,is) &
                                * vfeat(:,ifeat,is) * fc(:,ibas,ialpha,is)
      enddo
    enddo
    !
    do ialpha=1,cider_nalpha
      CALL get_cider_lpot_exp( dfftp%nnr, rho%of_r(:,is), &
                               grho(:,:,is), rho%kin_r(:,is)/e2, &
                               cider_consts(:,ialpha), cider_exp(:,ialpha,is), &
                               vexp(:,ialpha,is), v1x(:,is), &
                               v2x(:,is), v3x(:,is) )
    enddo
    do ibas=1,cider_nbas
      lmind=1
      do l=0,cider_lmax
        do m=-l,l
          CALL get_cider_lpot_bas( vbas(:,ibas,lmind,is), ylm(:,lmind), &
                                   l, ibas, v1x(:,is) )
          lmind=lmind+1
        enddo
      enddo
    enddo
  ENDDO
  !
  IF (nspin == 1) THEN
    DO k = 1, dfftp%nnr
       !
       v(k,1) = (v1x(k,1)+v1c(k,1)) * e2
       !
       ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       h(:,k,1) = (v2x(k,1)+v2c(1,k,1)) * grho(:,k,1) * e2 
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
       h(:,k,1) = (v2x(k,1) * grho(:,k,1) + v2c(:,k,1)) * e2
       h(:,k,2) = (v2x(k,2) * grho(:,k,2) + v2c(:,k,2)) * e2
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
  DEALLOCATE( bas, vbas, vexp, fc, dfc, ylm, cider_exp, const )
  DEALLOCATE( feat, dfeat, vfeat )
  !
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
END SUBROUTINE v_xc_cider
!
!
!----------------------------------------------------------------------------
SUBROUTINE v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
  !----------------------------------------------------------------------------
  !! Exchange-Correlation potential (meta) Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE funct,            ONLY : dft_is_nonlocc, nlc
  USE scf,              ONLY : scf_type, rhoz_or_updw
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
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
  !
  ! ... local variables
  !
  REAL(DP) :: zeta, rh, sgn(2)
  INTEGER  :: k, ipol, is, np
  !
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:,:), v3c(:,:)
  !
  REAL(DP) :: fac
       
  REAL(DP), DIMENSION(2) :: grho2, rhoneg
  REAL(DP), DIMENSION(3) :: grhoup, grhodw
  !
  REAL(DP), ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: rhoout(:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:)
  REAL(DP), PARAMETER :: eps12 = 1.0d-12, zero=0._dp
  !
  CALL start_clock( 'v_xc_meta' )
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
  ! ... calculate the gradient of rho + rho_core in real space
  ! ... in LSDA case rhogsum is in (up,down) format
  !
  DO is = 1, nspin
     !
     rhogsum(:) = fac*rhog_core(:) + ( rho%of_g(:,1) + sgn(is)*rho%of_g(:,nspin) )*0.5D0
     !
     CALL fft_gradient_g2r( dfftp, rhogsum, g, grho(1,1,is) )
     !
  ENDDO
  DEALLOCATE(rhogsum)
  !
  !
  IF (nspin == 1) THEN
    !
    CALL xc_metagcx( dfftp%nnr, 1, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                      v1x, v2x, v3x, v1c, v2c, v3c )
    !
    DO k = 1, dfftp%nnr
       !
       v(k,1) = (v1x(k,1)+v1c(k,1)) * e2
       !
       ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       h(:,k,1) = (v2x(k,1)+v2c(1,k,1)) * grho(:,k,1) * e2 
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
    CALL rhoz_or_updw( rho, 'only_r', '->updw' )
    !
    CALL xc_metagcx( dfftp%nnr, 2, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c )
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
       h(:,k,1) = (v2x(k,1) * grho(:,k,1) + v2c(:,k,1)) * e2
       h(:,k,2) = (v2x(k,2) * grho(:,k,2) + v2c(:,k,2)) * e2
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
    CALL rhoz_or_updw( rho, 'only_r', '->rhoz' )
    !
  ENDIF
  !
  DEALLOCATE( ex, ec )
  DEALLOCATE( v1x, v2x, v3x )
  DEALLOCATE( v1c, v2c, v3c )
  !
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
  CALL stop_clock( 'v_xc_meta' )
  !
  RETURN
  !
END SUBROUTINE v_xc_meta
!
!
SUBROUTINE v_xc( rho, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !! Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : nlc, dft_is_nonlocc
  USE scf,              ONLY : scf_type
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin)
  !! V_xc potential
  REAL(DP), INTENT(OUT) :: vtxc
  !! integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg(2), vs
  !
  !REAL(DP), ALLOCATABLE :: arhox(:), amag(:), zeta(:)
  REAL(DP) :: arho, amag
  REAL(DP) :: rhoup2, rhodw2
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
  ! In order:
    ! the absolute value of the total charge
    ! the absolute value of the magnetization
    ! zeta = amag / arhox
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  !
  !
  CALL start_clock( 'v_xc' )
  !
  ALLOCATE( ex(dfftp%nnr) )
  ALLOCATE( ec(dfftp%nnr) )
  ALLOCATE( vx(dfftp%nnr,nspin) )
  ALLOCATE( vc(dfftp%nnr,nspin) )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  !
  rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     ! ... spin-unpolarized case
     !
     CALL xc( dfftp%nnr, 1, 1, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr
        v(ir,1) = e2*( vx(ir,1) + vc(ir,1) )
        etxc = etxc + e2*( ex(ir) + ec(ir) )*rho%of_r(ir,1)
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtxc = vtxc + v(ir,1)*rho%of_r(ir,1)
        IF (rho%of_r(ir,1) < 0.D0) rhoneg(1) = rhoneg(1)-rho%of_r(ir,1)
     ENDDO
     !
     !
  ELSEIF ( nspin == 2 ) THEN
     ! ... spin-polarized case
     !
     CALL xc( dfftp%nnr, 2, 2, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr   !OMP ?
        v(ir,:) = e2*( vx(ir,:) + vc(ir,:) )
        etxc = etxc + e2*( (ex(ir) + ec(ir))*rho%of_r(ir,1) )
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtxc = vtxc + ( ( v(ir,1) + v(ir,2) )*rho%of_r(ir,1) + &
                        ( v(ir,1) - v(ir,2) )*rho%of_r(ir,2) )
        !
        rhoup2 = rho%of_r(ir,1)+rho%of_r(ir,2)
        rhodw2 = rho%of_r(ir,1)-rho%of_r(ir,2)
        IF (rhoup2 < 0.d0) rhoneg(1) = rhoneg(1) - rhoup2
        IF (rhodw2 < 0.d0) rhoneg(2) = rhoneg(2) - rhodw2
     ENDDO
     !
     vtxc   = 0.5d0 * vtxc
     rhoneg = 0.5d0 * rhoneg
     !
     !
  ELSE IF ( nspin == 4 ) THEN
     ! ... noncolinear case
     !
     CALL xc( dfftp%nnr, 4, 2, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr  !OMP ?
        arho = ABS( rho%of_r(ir,1) )
        IF ( arho < vanishing_charge ) CYCLE
        vs = 0.5D0*( vx(ir,1) + vc(ir,1) - vx(ir,2) - vc(ir,2) )
        v(ir,1) = e2*( 0.5D0*( vx(ir,1) + vc(ir,1) + vx(ir,2) + vc(ir,2) ) )
        !
        amag = SQRT( SUM( rho%of_r(ir,2:4)**2 ) )
        IF ( amag > vanishing_mag ) THEN
           v(ir,2:4) = e2 * vs * rho%of_r(ir,2:4) / amag
           vtxc = vtxc + SUM( v(ir,2:4) * rho%of_r(ir,2:4) )
        ENDIF
        etxc = etxc + e2*( ex(ir) + ec(ir) ) * arho
        !
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        IF ( amag / arho > 1.D0 )  rhoneg(2) = rhoneg(2) + 1.D0/omega
        vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
     ENDDO
     !
     !
  ENDIF
  !
  DEALLOCATE( ex, ec )
  DEALLOCATE( vx, vc )
  !
  CALL mp_sum(  rhoneg , intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2ES10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  ! ... add gradient corrections (if any)
  !
  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, v )
  !
  ! ... add non local corrections (if any)
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )
  !
  CALL stop_clock( 'v_xc' )
  !
  RETURN
  !
END SUBROUTINE v_xc
!
!----------------------------------------------------------------------------
SUBROUTINE v_h( rhog, ehart, charge, v )
  !----------------------------------------------------------------------------
  !! Hartree potential VH(r) from n(G)
  !
  USE constants,         ONLY : fpi, e2
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : invfft
  USE gvect,             ONLY : ngm, gg, gstart
  USE lsda_mod,          ONLY : nspin
  USE cell_base,         ONLY : omega, tpiba2
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE esm,               ONLY : do_comp_esm, esm_hartree, esm_bc
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_2D, cutoff_hartree  
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: rhog(ngm)
  !! the charge density in reciprocal space
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !! Hartree potential
  REAL(DP), INTENT(OUT) :: ehart
  !! Hartree energy
  REAL(DP), INTENT(OUT) :: charge
  !
  !  ... local variables
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  !
  CALL start_clock( 'v_h' )
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  charge = 0.D0
  !
  IF ( gstart == 2 ) THEN
     !
     charge = omega*REAL( rhog(1) )
     !
  END IF
  !
  CALL mp_sum(  charge , intra_bgrp_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... calculate modified Hartree potential for ESM
     !
     CALL esm_hartree (rhog, ehart, aux)
     !
  ELSE
     !
     ehart     = 0.D0
     aux1(:,:) = 0.D0
     !
     IF (do_cutoff_2D) THEN  !TS
        CALL cutoff_hartree(rhog(:), aux1, ehart)
     ELSE
!$omp parallel do private( fac, rgtot_re, rgtot_im ), reduction(+:ehart)
        DO ig = gstart, ngm
           !
           fac = 1.D0 / gg(ig) 
           !
           rgtot_re = REAL(  rhog(ig) )
           rgtot_im = AIMAG( rhog(ig) )
           !
           ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
           !
           aux1(1,ig) = rgtot_re * fac
           aux1(2,ig) = rgtot_im * fac
           !
        ENDDO
!$omp end parallel do
     ENDIF
     !
     fac = e2 * fpi / tpiba2
     !
     ehart = ehart * fac
     !
     aux1 = aux1 * fac
     !
     IF ( gamma_only ) THEN
        !
        ehart = ehart * omega
        !
     ELSE 
        !
        ehart = ehart * 0.5D0 * omega
        !
     END IF
     ! 
     if (do_comp_mt) then
        ALLOCATE( vaux( ngm ), rgtot(ngm) )
        rgtot(:) = rhog(:)
        CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
        aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
        aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
        ehart = ehart + eh_corr
        DEALLOCATE( rgtot, vaux )
     end if
     !
     CALL mp_sum(  ehart , intra_bgrp_comm )
     ! 
     aux(:) = 0.D0
     !
     aux(dfftp%nl(1:ngm)) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )
     !
     IF ( gamma_only ) THEN
        !
        aux(dfftp%nlm(1:ngm)) = CMPLX( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
        !
     END IF
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL invfft('Rho', aux, dfftp)
  !
  ! ... add hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     v(:,1) = v(:,1) + DBLE(aux(:))
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + DBLE(aux(:))
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( aux, aux1 )
  !
  CALL stop_clock( 'v_h' )
  !
  RETURN
  !
END SUBROUTINE v_h
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !! DFT+U+J0: B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_alpha, Hubbard_J0, Hubbard_beta
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity, dfpt_hub
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  REAL(DP) :: effU, sgn(2) 
  INTEGER  :: is, isop, na, nt, m1, m2
  !
  eth    = 0.d0
  sgn(1) =  1.d0  
  sgn(2) = -1.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
        !
        IF (Hubbard_J0(nt) /= 0.d0) THEN
           effU = Hubbard_U(nt) - Hubbard_J0(nt)
        ELSE
           effU = Hubbard_U(nt)
        ENDIF 
        ! 
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + ( Hubbard_alpha(nt) + 0.5D0*effU )*ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                                   Hubbard_alpha(nt)  + 0.5D0*effU
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 eth = eth - 0.5D0 * effU * ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                        effU * ns(m2,m1,is,na)
              ENDDO
           ENDDO
        ENDDO
        !
     ENDIF
     !
     IF (Hubbard_J0(nt) /= 0.d0 .OR. Hubbard_beta(nt) /= 0.d0) THEN
        !
        DO is = 1, nspin
           isop = 1
           IF ( nspin == 2 .AND. is == 1) isop = 2
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + sgn(is)*Hubbard_beta(nt) * ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + sgn(is)*Hubbard_beta(nt)
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 eth = eth + 0.5D0*Hubbard_J0(nt)*ns(m2,m1,is,na)*ns(m1,m2,isop,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + Hubbard_J0(nt) * &
                                                            ns(m2,m1,isop,na)
              ENDDO
           ENDDO
        ENDDO
        !
     END IF
     !   
  ENDDO
  !
  IF (nspin==1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,*) '--- in v_hubbard ---'
     WRITE(stdout,'("Hubbard energy ",f9.4)') eth
     WRITE(stdout,*) '-------'
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard
!-----------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE v_hubbard_b (ns, v_hub, eth)
  !-------------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy for background states.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_J0, Hubbard_beta, Hubbard_U_back, &
                                   ldim_back, ldmx_b, Hubbard_alpha_back, &
                                   is_hubbard_back
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity, dfpt_hub
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: v_hub(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: eth
  REAL(DP) :: effU
  INTEGER :: is, is1, na, nt, m1, m2, m3, m4
  !
  eth = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_J0(nt).NE.0.d0 .OR. Hubbard_beta(nt).NE.0.d0) &
     CALL errore('v_hubbard_b', 'Hubbard_J0 and Hubbard_beta are not supported',1) 
     !
     IF (is_hubbard_back(nt)) THEN
        !
        effU = Hubbard_U_back(nt)
        !
        DO is = 1, nspin
           !
           DO m1 = 1, ldim_back(nt) 
              !
              eth = eth + ( Hubbard_alpha_back(nt) + 0.5D0 * effU ) * &
                              ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                            ( Hubbard_alpha_back(nt) + 0.5D0 * effU )
              !
              DO m2 = 1, ldim_back(nt)
                 !
                 eth = eth - 0.5D0 * effU * &
                            ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                       effU * ns(m2,m1,is,na)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF (nspin.EQ.1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,*) '--- in v_hubbard_b ---'
     WRITE(stdout,'(''Hubbard_back energy '',f9.4)') eth
     WRITE(stdout,*) '-------'
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_b
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_full( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  REAL(DP) :: n_tot, n_spin, eth_dc, eth_u, mag2
  INTEGER  :: is, isop, is1, na, nt, m1, m2, m3, m4
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth    = 0.d0
  eth_dc = 0.d0
  eth_u  = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt)/=0.d0) THEN
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                               Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        !
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_tot = n_tot + ns(m1,m1,is,na)
           ENDDO
        ENDDO
        !
        IF (nspin==1) n_tot = 2.d0 * n_tot
        !
        mag2  = 0.d0
        ! 
        IF (nspin==2) THEN
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              mag2 = mag2 + ns(m1,m1,1,na) - ns(m1,m1,2,na)
           ENDDO
        ENDIF
        mag2  = mag2**2
        !
        ! Hubbard energy: DC term
        !
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*n_tot*(n_tot-1.d0) - &
                                  Hubbard_J(1,nt)*n_tot*(0.5d0*n_tot-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )
        !
        DO is = 1, nspin
           !
           ! n_spin = up/down N
           !
           n_spin = 0.d0
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_spin = n_spin + ns(m1,m1,is,na)
           ENDDO
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              !
              ! Hubbard potential: DC contribution  
              !
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_spin + &
                         0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot
              !
              ! +U contributions 
              !
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 DO m3 = 1, 2 * Hubbard_l(nt) + 1
                    DO m4 = 1, 2 * Hubbard_l(nt) + 1
                       !
                       DO is1 = 1, nspin
                          v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + (MOD(nspin,2)+1) * &
                                               u_matrix(m1,m3,m2,m4) * ns(m3,m4,is1,na)
                       ENDDO
                       !
                       v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                            u_matrix(m1,m3,m4,m2) * ns(m3,m4,is,na)
                       !
                       eth_u = eth_u + 0.5d0*( ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) ) * &
                                       ns(m1,m3,is,na)*ns(m2,m4,is,na)+u_matrix(m1,m2,m3,m4)   * &
                                       ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                    ENDDO ! m4
                 ENDDO ! m3
              ENDDO ! m2
              !
           ENDDO ! m1
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  ! 
  IF (nspin==1) eth_u = 2.d0 * eth_u
  eth = eth_u - eth_dc
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,*) '--- in v_hubbard ---'
     WRITE(stdout,'("Hubbard energies (dc, U, total) ",3f9.4)') eth_dc, eth_u, eth
     WRITE(stdout,*) '-------'
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full
!---------------------------------------------------------------

!---------------------------------------------------------------
SUBROUTINE v_hubbard_full_nc( ns, v_hub, eth )
  !-------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy (noncollinear case).
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE noncollin_module,     ONLY : noncolin
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP) :: eth
  !! Hubbard matrix
  !
  !  ... local variables
  !
  REAL(DP) :: eth_dc, eth_noflip, eth_flip, mx, my, mz, mag2
  INTEGER :: is, is1, js, i, j, na, nt, m1, m2, m3, m4
  COMPLEX(DP) :: n_tot, n_aux
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth        = 0.d0
  eth_dc     = 0.d0  
  eth_noflip = 0.d0  
  eth_flip   = 0.d0  
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat  
     !
     nt = ityp (na)  
     !
     IF (Hubbard_U(nt) /= 0.d0) THEN  
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                             Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        mx    = 0.d0
        my    = 0.d0
        mz    = 0.d0
        DO m1 = 1, 2*Hubbard_l(nt)+1
          n_tot = n_tot + ns(m1,m1,1,na) + ns(m1,m1,4,na)
          mx = mx + DBLE( ns(m1, m1, 2, na) + ns(m1, m1, 3, na) )
          my = my + 2.d0 * AIMAG( ns(m1, m1, 2, na) )
          mz = mz + DBLE( ns(m1, m1, 1, na) - ns(m1, m1, 4, na) )
        ENDDO  
        mag2 = mx**2 + my**2 + mz**2  
        !
        ! Hubbard energy: DC term
        !
        mx = REAL(n_tot)
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*mx*(mx-1.d0) - &
                                  Hubbard_J(1,nt)*mx*(0.5d0*mx-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )   
        !
        DO is = 1, nspin  
           !
           IF (is == 2) THEN
            is1 = 3
           ELSEIF (is == 3) THEN
            is1 = 2
           ELSE
            is1 = is
           ENDIF
           !
           ! Hubbard energy:
           !
           IF (is1 == is) THEN
             !
             ! Non spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                DO m3 = 1, 2*Hubbard_l(nt)+1
                 DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_noflip = eth_noflip + 0.5d0*(                            &
                              ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) )*  & 
                              ns(m1,m3,is,na)*ns(m2,m4,is,na) +                 &
                      u_matrix(m1,m2,m3,m4)*ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                 ENDDO
                ENDDO
              ENDDO
             ENDDO
             !
           ELSE
             ! 
             ! Spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_flip = eth_flip - 0.5d0*u_matrix(m1,m2,m4,m3)* &
                                     ns(m1,m3,is,na)*ns(m2,m4,is1,na) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! Hubbard potential: non spin-flip contribution 
           !
           IF (is1 == is) THEN
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                    u_matrix(m1,m3,m2,m4)*( ns(m3,m4,1,na)+ns(m3,m4,4,na) ) 
                ENDDO 
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! n_aux = /sum_{i} n_{i,i}^{sigma2, sigma1} for DC term
           !
           n_aux = 0.d0
           DO m1 = 1, 2*Hubbard_l(nt)+1
              n_aux = n_aux + ns(m1,m1,is1,na)  
           ENDDO
           !
           DO m1 = 1, 2*Hubbard_l(nt)+1
             ! 
             ! Hubbard potential: DC contribution  
             !
             v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_aux
             !
             IF (is1 == is) THEN
                v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                     0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot  
             ENDIF
             !
             ! Hubbard potential: spin-flip contribution
             !
             DO m2 = 1, 2*Hubbard_l(nt)+1  
              DO m3 = 1, 2*Hubbard_l(nt)+1
               DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                             u_matrix(m1,m3,m4,m2) * ns(m3,m4,is1,na) 
               ENDDO 
              ENDDO
             ENDDO
             !
           ENDDO
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  !
  eth = eth_noflip + eth_flip - eth_dc
  !
  ! Hubbard energies
  !
  IF ( iverbosity > 0 ) THEN
    WRITE(stdout,*) '--- in v_hubbard ---'
    WRITE(stdout,'("Hub. E (dc, noflip, flip, total) ",4f9.4)') &
                                 eth_dc, eth_noflip, eth_flip, eth 
    WRITE(stdout,*) '-------'
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full_nc
!----------------------------------------------------------------------------

!------------------------------------------------------------------------------------
SUBROUTINE v_hubbard_extended (nsg, v_hub, eth)
  !-----------------------------------------------------------------------------------
  !
  !! Computes extended Hubbard potential and Hubbard energy.
  !! DFT+U+V: Simplified rotationally-invariant formulation by
  !! V.L. Campo Jr and M. Cococcioni, J. Phys.: Condens. Matter 22, 055602 (2010).
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                Hubbard_V, Hubbard_alpha_back, is_hubbard, is_hubbard_back
  USE lsda_mod,          ONLY : nspin
  USE control_flags,     ONLY : iverbosity, dfpt_hub
  USE io_global,         ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  REAL(DP),    INTENT(OUT) :: eth
  ! 
  ! Local variables 
  !
  INTEGER :: is, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2, i_type
  INTEGER, EXTERNAL :: type_interaction, find_viz
  !
  eth  = 0.d0
  v_hub(:,:,:,:,:) = (0.d0, 0.d0)
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        !
        DO is = 1, nspin
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              equiv_na2 = at_sc(na2)%at
              nt2 = ityp(equiv_na2)
              !
              IF ((is_hubbard(nt2).OR.is_hubbard_back(nt2)) .AND. &
                  (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,2).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,3).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,4).NE.0.d0) ) THEN
                  !
                  ! For both standard and background states of a center atom
                  DO m1 = 1, ldim_u(nt1)
                     ! For both standard and background states of the neighbor atom
                     DO m2 = 1, ldim_u(nt2)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m2)
                        !
                        v_hub(m2,m1,viz,na1,is) = &
                           - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,i_type)
                        !
                        eth = eth - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                    * Hubbard_V(na1,na2,i_type) * 0.5d0
                        !
                     ENDDO
                  ENDDO
                  !
                  IF ( na1.EQ.na2 ) THEN
                     !
                     na = find_viz(na1,na1)
                     !
                     ! This is the diagonal term (like in the DFT+U only case)
                     ! 
                     DO m1 = 1, ldim_u(nt1)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m1)
                        !
                        v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                       + Hubbard_V(na1,na1,i_type) * 0.5d0
                        ! 
                        eth = eth + nsg(m1,m1,na,na1,is) &
                                       * Hubbard_V(na1,na1,i_type) * 0.5d0
                        !
                     ENDDO
                     !
                     ! Hubbard_J0 (only on-site)
                     !
                     IF ( nspin.EQ.2 .AND. &
                          (Hubbard_J0(nt1).NE.0.d0 .OR. Hubbard_beta(nt1).NE.0.d0) ) THEN
                          !
                          IF (is.EQ.1) THEN
                             isop = 2
                          ELSE
                             isop = 1
                          ENDIF
                          !
                          DO m1 = 1, 2*Hubbard_l(nt1)+1
                             !
                             v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                      Hubbard_J0(nt1) * 0.5d0
                             eth = eth - 0.5d0 * Hubbard_J0(nt1) * nsg(m1,m1,na,na1,is)
                             !
                             IF (is.EQ.1) THEN
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + &
                                                         Hubbard_beta(nt1)
                                eth = eth + Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ELSE
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                         Hubbard_beta(nt1)
                                eth = eth - Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ENDIF
                             !
                             DO m2 = 1, 2*Hubbard_l(nt1)+1
                                v_hub(m2,m1,na,na1,is) = v_hub(m2,m1,na,na1,is) + &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1)
                                eth = eth + nsg(m2,m1,na,na1,is) *    &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1) * 0.5d0
                             ENDDO
                             !
                          ENDDO ! m1
                          !
                     ENDIF
                     !
                  ENDIF
                  !
              ENDIF
              !
           ENDDO ! viz
           !
        ENDDO ! is
        !  
     ENDIF
     !
     ! Hubbard_alpha or Hubbard_alpha_back
     !
     IF ( ldim_u(nt1).GT.0 .AND. &
          (Hubbard_alpha(nt1).NE.0.0d0 .OR. Hubbard_alpha_back(nt1).NE.0.0d0) ) THEN
          !
          na = find_viz(na1,na1)
          !
          DO is = 1, nspin
             !
             DO m1 = 1, 2*Hubbard_l(nt1)+1
                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha(nt1)
                eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha(nt1)
             ENDDO
             !
             IF ( ldim_u(nt1).GT.2*Hubbard_l(nt1)+1 .AND. &
                  Hubbard_alpha_back(nt1).NE.0.0d0 ) THEN
                ! Background states
                DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                   v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha_back(nt1)
                   eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha_back(nt1)
                ENDDO
             ENDIF
             !
          ENDDO
          !
     ENDIF
     !
  ENDDO ! na1
  !
  IF (nspin.EQ.1) eth = eth * 2.d0
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,*) '--- in v_hubbard_extended ---'
     WRITE(stdout,'("Hubbard energy ",f9.4)') eth
     WRITE(stdout,*) '-----------------------------'
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_extended
!---------------------------------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE v_h_of_rho_r( rhor, ehart, charge, v )
  !----------------------------------------------------------------------------
  !! Hartree potential VH(r) from a density in R space n(r) 
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rhor( dfftp%nnr )
  REAL( DP ), INTENT(INOUT)  :: v( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: ehart, charge
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhog( : )
  COMPLEX( DP ), ALLOCATABLE :: aux( : )
  REAL( DP ), ALLOCATABLE :: vaux(:,:)
  INTEGER :: is
  !
  ! ... bring the (unsymmetrized) rho(r) to G-space (use aux as work array)
  !
  ALLOCATE( rhog( dfftp%ngm ) )
  ALLOCATE( aux( dfftp%nnr ) )
  aux = CMPLX(rhor,0.D0,kind=dp)
  CALL fwfft ('Rho', aux, dfftp)
  rhog(:) = aux(dfftp%nl(:))
  DEALLOCATE( aux )
  !
  ! ... compute VH(r) from n(G)
  !
  ALLOCATE( vaux( dfftp%nnr, nspin ) )
  vaux = 0.D0
  CALL v_h( rhog, ehart, charge, vaux )
  v(:) = v(:) + vaux(:,1)
  !
  DEALLOCATE( rhog )
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE v_h_of_rho_r
!
!----------------------------------------------------------------------------
SUBROUTINE gradv_h_of_rho_r( rho, gradv )
  !----------------------------------------------------------------------------
  !! Gradient of Hartree potential in R space from a total 
  !! (spinless) density in R space n(r)
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE constants,       ONLY : fpi, e2
  USE control_flags,   ONLY : gamma_only
  USE cell_base,       ONLY : tpiba, omega
  USE gvect,           ONLY : ngm, gg, gstart, g
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rho(dfftp%nnr)
  REAL( DP ), INTENT(OUT)    :: gradv(3, dfftp%nnr)
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhoaux(:)
  COMPLEX( DP ), ALLOCATABLE :: gaux(:)
  COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
  REAL( DP )                 :: fac, eh_corr
  INTEGER                    :: ig, ipol
  !
  ! ... Bring rho to G space
  !
  ALLOCATE( rhoaux( dfftp%nnr ) )
  rhoaux( : ) = CMPLX( rho( : ), 0.D0, KIND=dp ) 
  !
  CALL fwfft( 'Rho', rhoaux, dfftp )
  !
  ! ... Compute total potential in G space
  !
  ALLOCATE( gaux( dfftp%nnr ) )
  !
  DO ipol = 1, 3
    !
    gaux(:) = (0.0_dp,0.0_dp)
    !
    DO ig = gstart, ngm
      !
      fac = g(ipol,ig) / gg(ig)
      gaux(dfftp%nl(ig)) = CMPLX(-AIMAG(rhoaux(dfftp%nl(ig))),REAL(rhoaux(dfftp%nl(ig))),KIND=dp)*fac 
      !
    ENDDO
    !
    ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of 
    !  V = e2 * fpi divided by the 2\pi/a factor missing in G  
    !
    fac = e2 * fpi / tpiba
    gaux = gaux * fac 
    !
    ! ...add martyna-tuckerman correction, if needed
    ! 
    IF (do_comp_mt) THEN
       ALLOCATE( vaux( ngm ), rgtot(ngm) )
       rgtot(1:ngm) = rhoaux(dfftp%nl(1:ngm))
       CALL wg_corr_h( omega, ngm, rgtot, vaux, eh_corr )
       DO ig = gstart, ngm
         fac = g(ipol,ig) * tpiba
         gaux(dfftp%nl(ig)) = gaux(dfftp%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac 
       END DO
       DEALLOCATE( rgtot, vaux )
    ENDIF
    !
    IF ( gamma_only ) THEN
      !
      gaux(dfftp%nlm(:)) = &
        CMPLX( REAL( gaux(dfftp%nl(:)) ), -AIMAG( gaux(dfftp%nl(:)) ) ,kind=DP)
       !
    END IF
    !
    ! ... bring back to R-space, (\grad_ipol a)(r) ...
    !
    CALL invfft( 'Rho', gaux, dfftp )
    !
    gradv(ipol,:) = REAL( gaux(:) )
    !
  ENDDO
  !
  DEALLOCATE(gaux)
  !
  DEALLOCATE(rhoaux)
  !
  RETURN
  !
END SUBROUTINE gradv_h_of_rho_r
!
!----------------------------------------------------------------------------
SUBROUTINE get_cider_bas( rhog, feat, ylm, l, ibas )
  !----------------------------------------------------------------------------
  !! Cider features from density in reciprocal space
  !
  USE constants,         ONLY : fpi, e2, pi
  USE input_parameters,  ONLY : cider_params
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : invfft
  USE gvect,             ONLY : ngm, gg
  USE lsda_mod,          ONLY : nspin
  USE cell_base,         ONLY : omega, tpiba2
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: rhog(ngm)
  REAL(DP), INTENT(IN) :: ylm(ngm)
  INTEGER, INTENT(IN) :: l
  INTEGER, INTENT(IN) :: ibas
  !! the charge density in reciprocal space
  REAL(DP), INTENT(INOUT) :: feat(dfftp%nnr)
  !
  !  ... local variables
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  REAL(DP)              :: aexp
  !
  CALL start_clock( 'get_cider_bas' )
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  !
  ! ... calculate CIDER features in reciprocal space
  !
  aux1(:,:) = 0.D0
  !aexp = 0.2
  aexp = cider_params(1) / cider_params(2)**(ibas-1)
  !
!$omp parallel do private( fac, rgtot_re, rgtot_im )
  DO ig = 1, ngm
     !
     fac = (pi/aexp)**1.5 * EXP(-gg(ig) * tpiba2 / (4.D0 * aexp))
     fac = fac * (gg(ig) / (4 * aexp * aexp))**(DBLE(l)/2.0_DP) * ylm(ig) ! real sph_harm
     !
     rgtot_re = REAL(  rhog(ig) )
     rgtot_im = AIMAG( rhog(ig) )
     !
     aux1(1,ig) = rgtot_re * fac
     aux1(2,ig) = rgtot_im * fac
     !
  ENDDO
!$omp end parallel do
  !
  ! 
  aux(:) = 0.D0
  !
  aux(dfftp%nl(1:ngm)) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )
  !
  IF ( gamma_only ) THEN
     !
     aux(dfftp%nlm(1:ngm)) = CMPLX( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
     !
  END IF
  !
  ! ... transform CIDER potential to real space
  !
  CALL invfft('Rho', aux, dfftp)
  !
  ! ... set to features
  !
  feat = DBLE(aux(:))
  !
  DEALLOCATE( aux, aux1 )
  !
  CALL stop_clock( 'get_cider_bas' )
  !
  RETURN
  !
END SUBROUTINE get_cider_bas
!
SUBROUTINE get_cider_lpot_bas( vr, ylm, l, ibas, v )
  !----------------------------------------------------------------------------
  !! CIDER potential wrt density from potential wrt features
  !
  USE constants,         ONLY : fpi, e2, pi
  USE input_parameters,  ONLY : cider_params
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : invfft, fwfft
  USE fft_helper_subroutines, ONLY: fftx_threed2oned
  USE gvect,             ONLY : ngm, gg, gstart
  USE lsda_mod,          ONLY : nspin
  USE cell_base,         ONLY : omega, tpiba2
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: vr(dfftp%nnr)
  REAL(DP), INTENT(IN) :: ylm(ngm)
  INTEGER, INTENT(IN) :: l
  INTEGER, INTENT(IN) :: ibas
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr)
  !
  !  ... local variables
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  COMPLEX(dp), ALLOCATABLE :: psi(:)
  COMPLEX(DP), ALLOCATABLE :: vg(:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  REAL(DP)              :: aexp
  !
  CALL start_clock( 'get_cider_lpot_bas' )
  !
  ALLOCATE(vg(ngm),psi(dfftp%nnr))
  !
  ! ... Transform the real-space de/dfeat(r) to de/dfeat(k)
  !
  vg(:) = 0.0_dp
  psi(:) = 0.0_dp
  psi = CMPLX(vr,0.D0,kind=dp)
  CALL fwfft ('Rho', psi, dfftp)
  CALL fftx_threed2oned( dfftp, psi, vg )
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  !
  ! ... Convolve potential in reciprocal space by multiplying by
  !     the feature kernel, converting de/dfeat(k) to de/drho(k).
  !
  aux1(:,:) = 0.D0
  !aexp = 0.2
  aexp = cider_params(1) / cider_params(2)**(ibas-1)
  !
!$omp parallel do private( fac, rgtot_re, rgtot_im )
  DO ig = 1, ngm
     !
     fac = (pi/aexp)**1.5 * EXP(-gg(ig) * tpiba2 / (4.D0 * aexp))
     fac = fac * (gg(ig) / (4 * aexp * aexp))**(DBLE(l)/2.0_DP) * ylm(ig) ! real sph_harm
     !
     rgtot_re = REAL(  vg(ig) )
     rgtot_im = AIMAG( vg(ig) )
     !
     aux1(1,ig) = rgtot_re * fac
     aux1(2,ig) = rgtot_im * fac
     !
  ENDDO
!$omp end parallel do
  ! 
  aux(:) = 0.D0
  !
  aux(dfftp%nl(1:ngm)) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )
  !
  IF ( gamma_only ) THEN
     !
     aux(dfftp%nlm(1:ngm)) = CMPLX( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
     !
  END IF
  !
  ! ... transform CIDER xc potential to real space
  !
  CALL invfft('Rho', aux, dfftp)
  !
  ! ... add contribution to CIDER potential
  !
  v = v + DBLE(aux(:))
  !
  DEALLOCATE( aux, aux1, vg, psi )
  !
  CALL stop_clock( 'get_cider_lpot_bas' )
  !
  RETURN
  !
END SUBROUTINE get_cider_lpot_bas
!
!
SUBROUTINE get_cider_coefs( length, cider_exp, fc, dfc )
    !
    ! Takes CIDER auxiliary basis features as input, along with
    ! CIDER exponents, and returns the basis-transformed version
    !
    USE kind_l,            ONLY: DP
    USE input_parameters,  ONLY: cider_params, cider_nbas
    USE constants,         ONLY: pi
    !
    IMPLICIT NONE
    !
    integer, intent(in)   :: length
    real(dp), intent(in)  :: cider_exp(length)
    real(dp), intent(out) :: fc(length,cider_nbas)
    real(dp), intent(out) :: dfc(length,cider_nbas)
    !
    integer :: i,j
    real(dp), allocatable :: bas_exp(:)
    real(dp), allocatable :: sinv(:,:)
    real(dp), allocatable :: cvec(:,:)
    real(dp), allocatable :: dcvec(:,:)
    real(dp)              :: cider_maxexp, cider_xexp
    !
    allocate ( bas_exp(cider_nbas) )
    allocate ( sinv(cider_nbas, cider_nbas) )
    allocate ( cvec(length, cider_nbas) )
    allocate ( dcvec(length, cider_nbas) )
    cider_maxexp = cider_params(1)
    cider_xexp = cider_params(2)
    !
    do i=1,cider_nbas
        bas_exp(i) = cider_maxexp / cider_xexp**(i-1)
    enddo
    do i=1,cider_nbas
        do j=1,cider_nbas
            sinv(i,j) = sqrt( pi / (bas_exp(i) + bas_exp(j)) )**3.0_DP / (4 * pi)
        enddo
    enddo
    !
    CALL MatInv('G', cider_nbas, sinv)
    ! if not G, need to fill matrix into square after calling above subroutine
    !
    do i=1,length
        do j=1,cider_nbas
            cvec(i,j) = sqrt( pi / (cider_exp(i) + bas_exp(j)) )**3.0_DP / (4 * pi)
            dcvec(i,j) = -1.5_dp * cvec(i,j) / (cider_exp(i) + bas_exp(j))
        enddo
    enddo
    !
    fc = matmul(cvec, sinv) ! d_{i,sigma}
    dfc = matmul(dcvec, sinv) ! deriv wrt exponent
    !
    DEALLOCATE(bas_exp)
    DEALLOCATE(sinv)
    DEALLOCATE(cvec)
    DEALLOCATE(dcvec)
    !
    return
    !
END SUBROUTINE get_cider_coefs
!
!
SUBROUTINE get_cider_alpha( length, rho, grho, kin, cider_consts, cider_exp )
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
  real(dp) :: ns
  integer :: tot
  ns = DBLE(nspin)

  fac = cider_consts(3) * 1.2 * (6 * pi**2)**(2.0/3) / pi
  const1 = cider_consts(2) - fac
  const2 = fac / (0.3 * (3 * pi**2)**(2.0_DP/3))
  const1 = const1 * pi / 2**(2.0_DP/3)
  const2 = const2 * pi / 2**(2.0_DP/3)

  cider_exp = cider_consts(1) + const1 * (ns*rho)**(2.0_DP/3) + const2 * kin / rho

  do i=1,length
    if (cider_exp(i) < cider_consts(4)) then
      cider_exp(i) = cider_consts(4)*exp(cider_exp(i)/cider_consts(4)-1.0_DP)
    endif
  enddo

  return
END SUBROUTINE get_cider_alpha
!
SUBROUTINE get_cider_lpot_exp ( length, rho, grho, kin, cider_consts, &
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
  real(dp) :: ns
  integer :: tot
  ns = DBLE(nspin)

  fac = cider_consts(3) * 1.2 * (6 * pi**2)**(2.0/3) / pi
  const1 = cider_consts(2) - fac
  const2 = fac / (0.3 * (3 * pi**2)**(2.0_DP/3))
  const1 = const1 * pi / 2**(2.0_DP/3)
  const2 = const2 * pi / 2**(2.0_DP/3)

  allocate(cider_tmp(length))
  cider_tmp = cider_consts(1) + const1 * (ns*rho)**(2.0_DP/3) + const2 * kin / rho

  do i=1,length
    if (cider_tmp(i) < cider_consts(4)) then
      vexp(i) = vexp(i) / cider_consts(4) * cider_exp(i)
    endif
  enddo

  v1x = v1x + 2.0_DP/3.0_DP * DBLE(ns) * vexp &
              * const1 / (ns*abs(rho)+1e-8)**(1.0_DP/3)
  v1x = v1x - vexp * const2 * kin / (abs(rho)+1e-8)**2
  v3x = v3x + vexp * const2 / (abs(rho)+1e-8)

  DEALLOCATE(cider_tmp)

  return
END SUBROUTINE get_cider_lpot_exp
!
SUBROUTINE xc_cider_x_py(length, nfeat, ns, np, rho, grho, tau, feat, ex, v1x, v2x, v3x, vfeat )
    !
    USE forpy_mod
    !
    USE kind_l,        ONLY: DP
    USE input_parameters, ONLY: cider_py_obj
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
    integer :: is,k,sgn,ierror
    real(dp) :: x43,x13
    real(dp) :: xf,rval
    real(dp), allocatable :: grho2(:,:)
    real(dp), allocatable :: rho_tot(:)

    type(tuple) :: args
    type(ndarray) :: py_rho, py_grho, py_tau, py_feat, py_ex, py_v1x, py_v2x, py_v3x, py_vfeat

    ierror = forpy_initialize()
    print *, "init", ierror

    ierror = ndarray_create(py_rho, rho)
    print *, "numpy", ierror
    ierror = ndarray_create(py_grho, grho)
    ierror = ndarray_create(py_tau, tau)
    ierror = ndarray_create(py_feat, feat)
    ierror = ndarray_create_nocopy(py_ex, ex)
    ierror = ndarray_create_nocopy(py_v1x, v1x)
    ierror = ndarray_create_nocopy(py_v2x, v2x)
    ierror = ndarray_create_nocopy(py_v3x, v3x)
    ierror = ndarray_create_nocopy(py_vfeat, vfeat)

    ierror = tuple_create(args,9)
    ierror = args%setitem(0, py_rho)
    ierror = args%setitem(1, py_grho)
    ierror = args%setitem(2, py_tau)
    ierror = args%setitem(3, py_feat)
    ierror = args%setitem(4, py_ex)
    ierror = args%setitem(5, py_v1x)
    ierror = args%setitem(6, py_v2x)
    ierror = args%setitem(7, py_v3x)
    ierror = args%setitem(8, py_vfeat)
    print *, "args", ierror

    ierror = call_py_noret(cider_py_obj, "get_xc_fortran", args)
    print *, "call", ierror

    call py_rho%destroy
    call py_grho%destroy
    call py_tau%destroy
    call py_feat%destroy
    call py_ex%destroy
    call py_v1x%destroy
    call py_v2x%destroy
    call py_v3x%destroy
    call py_vfeat%destroy

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
END SUBROUTINE xc_cider_x_py