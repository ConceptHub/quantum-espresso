!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
!
!----------------------------------------------------------------------------
SUBROUTINE rrmmdiagg_gpu( h_psi_ptr, s_psi_ptr, npwx, npw, nbnd, psi, hpsi, spsi, e, &
                      g2kin_d, btype, ethr, ndiis, uspp, do_hpsi, is_exx, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE util_param,     ONLY : DP, eps14, eps16
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_bands_util, ONLY : gstart, inter_bgrp_comm, intra_bgrp_comm, me_bgrp, root_bgrp, &
                            root_bgrp_id, use_bgrp_in_hpsi
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: ethr
  INTEGER,     INTENT(IN)    :: ndiis
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: do_hpsi
  LOGICAL,     INTENT(IN)    :: is_exx
  INTEGER,     INTENT(OUT)   :: notconv
  REAL(DP),    INTENT(OUT)   :: rmm_iter
  !
  ! ... local variables
  !
  INTEGER                  :: ierr
  INTEGER                  :: idiis
  INTEGER                  :: motconv
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: ibnd, ibnd_start, ibnd_end, ibnd_size
  INTEGER,     ALLOCATABLE :: ibnd_index(:)
  INTEGER,     ALLOCATABLE :: jbnd_index(:)
  REAL(DP)                 :: empty_ethr
  REAL(DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  !
  REAL(DP),    PARAMETER   :: SREF = 0.50_DP
  REAL(DP),    PARAMETER   :: SMIN = 0.05_DP
  REAL(DP),    PARAMETER   :: SMAX = 1.00_DP
  !
  !
  ! ... device variables
  !
  INTEGER :: ii, jj, kk
  COMPLEX(DP), INTENT(INOUT)  :: psi (npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT)  :: hpsi(npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT)  :: spsi(npwx,nbnd)
  REAL(DP), INTENT(IN)        :: g2kin_d(npwx)
  COMPLEX(DP), ALLOCATABLE    :: phi_d(:,:,:), hphi_d(:,:,:), sphi_d(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  !$acc declare device_resident(kpsi, hkpsi, skpsi)
#if defined(__CUDA)
  attributes(device) :: g2kin_d
  attributes(device) :: phi_d, hphi_d, sphi_d
#endif
  EXTERNAL :: h_psi_ptr, s_psi_ptr
    ! h_psi_ptr(npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nbnd,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nbnd)
  !
  CALL start_clock( 'rrmmdiagg' )
  !
  IF ( gstart == -1 ) CALL errore( ' rrmmdiagg ', 'gstart variable not initialized', 1 )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx
  !
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  !
  ibnd_size = MAX( ibnd_end - ibnd_start + 1, 0 )
  !
  IF( ibnd_size == 0 ) CALL errore( ' rrmmdiagg ', ' ibnd_size == 0 ', 1 )
  !
  ALLOCATE( phi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate phi_d ', ABS(ierr) )
  !
  ALLOCATE( hphi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hphi_d ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sphi_d ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( kpsi( npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( npwx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !
  END IF
  !
!
  ALLOCATE( hr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hr ', ABS(ierr) )
  !
  ALLOCATE( sr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sr ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
  !
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( sw( nbnd ) )
  !
  ALLOCATE( conv( nbnd ) )
  ALLOCATE( ibnd_index( nbnd ) )
  ALLOCATE( jbnd_index( ibnd_start:ibnd_end ) )
  !
!$cuf kernel do(3)
  DO ii = 1, npwx 
    DO jj = ibnd_start, ibnd_end
      DO kk = 1, ndiis 
        phi_d(ii,jj,kk) = ZERO
        hphi_d(ii,jj,kk) = ZERO
      END DO 
    END DO 
  END DO 
  IF ( uspp ) THEN  
!$cuf kernel do(3)
    DO ii = 1, npwx 
      DO jj = ibnd_start, ibnd_end
        DO kk = 1, ndiis 
          sphi_d(ii,jj,kk) = ZERO
        END DO 
      END DO 
    END DO 
  END IF
  !
  !$acc kernels
  kpsi = ZERO
  hkpsi = ZERO
  IF ( uspp ) skpsi = ZERO 
  !$acc end kernels 
  !
  hr = 0._DP
  sr = 0._DP
  !
  php = 0._DP
  psp = 0._DP
  !
  ew = e
  hw = e
  sw = 1._DP
  !
  conv = .FALSE.
  ibnd_index = 0
  jbnd_index = 0
  !
  FORALL( ibnd = 1:nbnd )              ibnd_index(ibnd) = ibnd
  FORALL( ibnd = ibnd_start:ibnd_end ) jbnd_index(ibnd) = ibnd - ibnd_start + 1
  !
  rmm_iter = 0._DP
  notconv  = nbnd
  motconv  = ibnd_size
  !
  ! ... Calculate H |psi> and S |psi>, if required
  !
  IF ( do_hpsi ) THEN
     !
     CALL calc_hpsi_gpu( )
     !
  END IF
  !
  ! ... Set Im[ psi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !
     !$acc kernels
     DO jj  = 1, nbnd 
       psi (1,jj) = CMPLX( DBLE( psi (1,jj) ), 0._DP, kind=DP )
       hpsi(1,jj) = CMPLX( DBLE( hpsi(1,jj) ), 0._DP, kind=DP )
     END DO 
     IF ( uspp ) THEN 
       DO jj  = 1, nbnd 
         spsi(1,jj) = CMPLX( DBLE( spsi(1,jj) ), 0._DP, kind=DP )
       END DO 
     END IF 
     !$acc end kernels
     !
  END IF
  !
  ! ... RMM-DIIS's loop
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + DBLE( notconv )
     !
     ! ... Perform DIIS
     !
     CALL do_diis_gpu( idiis )
     !
     ! ... Line searching
     !
     CALL rr_line_search_gpu( )
     !
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( )
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  !
  rmm_iter = rmm_iter / DBLE( nbnd )
  !
  ! ... Merge wave functions
  !
  IF ( ibnd_start > 1 ) THEN
     !
     !$acc parallel loop collapse(2)
     DO ii = 1, npwx 
       DO jj = 1,ibnd_start-1  
           psi(ii,jj)  = ZERO
           hpsi(ii,jj) = ZERO
       END DO
     END DO
     !
     IF ( uspp ) THEN
        !
        !$acc parallel loop collapse(2)
        DO ii = 1, npwx 
          DO jj = 1,ibnd_start-1  
            spsi(ii,jj) = ZERO
          END DO
        END DO
        !
     END IF
     !
  END IF
  !
  IF ( ibnd_end < nbnd ) THEN
     !
     !$acc parallel loop collapse(2)
     DO ii = 1, npwx 
       DO jj = ibnd_end+1,nbnd  
           psi (ii,jj) = ZERO
           hpsi(ii,jj) = ZERO
       END DO
     END DO
     !
     IF ( uspp ) THEN
        !$acc parallel loop collapse(2)
        DO ii = 1, npwx 
          DO jj = ibnd_end+1,nbnd
              spsi(ii,jj) = ZERO
          END DO
        END DO
        !
     END IF
     !
  END IF
  !
  !$acc host_data use_device(psi, hpsi, spsi)
  CALL mp_sum( psi,  inter_bgrp_comm )
  CALL mp_sum( hpsi, inter_bgrp_comm )
  IF ( uspp ) &
  CALL mp_sum( spsi, inter_bgrp_comm )
  !$acc end host_data
  !
  DEALLOCATE( phi_d )
  DEALLOCATE( hphi_d )
  IF ( uspp ) DEALLOCATE( sphi_d )
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
!
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  DEALLOCATE( ibnd_index )
  DEALLOCATE( jbnd_index )
  !
  CALL stop_clock( 'rrmmdiagg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE calc_hpsi_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd
    !
    REAL(DP), EXTERNAL :: gpu_DDOT
    !
    ! ... Operate the Hamiltonian : H |psi>
    !
    !$acc kernels
    hpsi = ZERO 
    !$acc end kernels
    !
    !$acc host_data use_device(psi, hpsi)
    CALL h_psi_ptr( npwx, npw, nbnd, psi, hpsi )
    !$acc end host_data
    !
    ! ... Operate the Overlap : S |psi>
    !
    IF ( uspp ) THEN
       !
       !$acc kernels
       spsi = ZERO
       !$acc end kernels
       !
       !$acc host_data use_device(psi, spsi)
       CALL s_psi_ptr( npwx, npw, nbnd, psi, spsi )
       !$acc end host_data
       !
    END IF
    !
    ! ... Matrix element : <psi| H |psi>
    !
    !$acc host_data use_device(psi, hpsi)
    DO ibnd = ibnd_start, ibnd_end
       !
       hw(ibnd) = 2._DP * gpu_DDOT( npw2, psi(1,ibnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) &
          hw(ibnd)= hw(ibnd) - gpu_DDOT( 1, psi(1,ibnd), 1, hpsi(1,ibnd),1)
       !
    END DO
    !$acc end host_data
    !
    CALL mp_sum( hw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Matrix element : <psi| S |psi>
    !
    !$acc host_data use_device(psi, spsi)
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( uspp ) THEN
          !
          sw(ibnd) = 2._DP * gpu_DDOT( npw2, psi(1,ibnd), 1, spsi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) &
             sw(ibnd)= sw(ibnd) - gpu_DDOT( 1, psi(1,ibnd), 1, spsi(1,ibnd),1)
          !
       ELSE
          !
          sw(ibnd) = 2._DP * gpu_DDOT( npw2, psi(1,ibnd), 1, psi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) &
             sw(ibnd )= sw(ibnd) - gpu_DDOT( 1, psi(1,ibnd), 1, psi(1,ibnd),1)
          !
       END IF
       !
    END DO
    !$acc end host_data
    !
    CALL mp_sum( sw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' rrmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE calc_hpsi_gpu
  !
  SUBROUTINE do_diis_gpu( idiis )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: idiis
    !
    INTEGER                  :: ibnd, jbnd, kbnd
    INTEGER                  :: kdiis
    REAL(DP)                 :: norm
    REAL(DP)                 :: er
    REAL(DP),    ALLOCATABLE :: vr(:)
    REAL(DP),    ALLOCATABLE :: tr(:,:)
    !
    ! device variables
    !
    REAL(DP) :: kvr  ! vr(kdiis)
    COMPLEX(DP) :: ctmp , cctmp
    COMPLEX(DP), ALLOCATABLE :: vec1_d(:)
    COMPLEX(DP), ALLOCATABLE :: vec2_d(:,:)
    REAL(DP),    ALLOCATABLE :: tr_d(:,:)
#if defined (__CUDA) 
    attributes(device) :: vec1_d, vec2_d, tr_d
#endif
    !
    ALLOCATE( vec1_d( npwx ) )
    ALLOCATE( vec2_d( npwx, idiis ) )
    IF ( idiis > 1 )   ALLOCATE( vr( idiis ) )
    IF ( motconv > 0 ) ALLOCATE( tr( idiis, motconv ) )
    IF ( motconv > 0 ) ALLOCATE( tr_d( idiis, motconv ) )
    !
    ! ... Save current wave functions and matrix elements
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       !$acc host_data use_device(psi, hpsi, spsi)
       CALL DCOPY_gpu( npw2, psi (1,ibnd), 1, phi_d (1,ibnd,idiis), 1 )
       CALL DCOPY_gpu( npw2, hpsi(1,ibnd), 1, hphi_d(1,ibnd,idiis), 1 )
       IF ( uspp ) &
       CALL DCOPY_gpu( npw2, spsi(1,ibnd), 1, sphi_d(1,ibnd,idiis), 1 )
       !$acc end host_data
       !
       php(ibnd,idiis) = hw(ibnd)
       psp(ibnd,idiis) = sw(ibnd)
       !
    END DO
    !
    ! ... <R_i|R_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       ! ... Residual vectors : |R> = (H - e S) |psi>
       !
       DO kdiis = 1, idiis
          !
          er = php(ibnd,kdiis)
          !
          CALL DCOPY_gpu( npw2, hphi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL DAXPY_gpu( npw2, -er, sphi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
             !
          ELSE
             !
             CALL DAXPY_gpu( npw2, -er, phi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
             !
          END IF
          !
       END DO
       !
       er = php(ibnd,idiis)
       !
       CALL DCOPY_gpu( npw2, hphi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL DAXPY_gpu( npw2, -er, sphi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       ELSE
          !
          CALL DAXPY_gpu( npw2, -er, phi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       END IF
       !
       CALL DGEMV_gpu( 'T', npw2, idiis, 2._DP, vec2_d(1,1), npwx2, &
                   vec1_d(1), 1, 0._DP, tr_d(1,jbnd), 1 )
       !
       IF ( gstart == 2 ) THEN  
!$cuf kernel do(1)
         DO ii = 1, idiis
           tr_d(ii,jbnd) = tr_d(ii,jbnd) - DBLE( vec2_d(1,ii) ) * DBLE( vec1_d(1) )
         END DO 
       END IF 
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( tr_d, intra_bgrp_comm )
       !
    END IF
    !
    tr = tr_d
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       hr(1:idiis,idiis,ibnd) = tr(1:idiis,jbnd)
       hr(idiis,1:idiis,ibnd) = tr(1:idiis,jbnd)
       !
    END DO
    !
    ! ... <phi_i| S |phi_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       DO kdiis = 1, idiis
          !
          CALL DCOPY_gpu( npw2, phi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
          !
       END DO
       !
       IF ( uspp ) THEN
          !
          CALL DCOPY_gpu( npw2, sphi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       ELSE
          !
          CALL DCOPY_gpu( npw2, phi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       END IF
       !
       CALL DGEMV_gpu( 'T', npw2, idiis, 2._DP, vec2_d(1,1), npwx2, &
                   vec1_d(1), 1, 0._DP, tr_d(1,jbnd), 1 )
       !
       IF ( gstart == 2 ) THEN 
!$cuf kernel do(1)
         DO ii = 1, idiis
           tr_d(ii,jbnd) = tr_d(ii,jbnd) - DBLE( vec2_d(1,ii) ) * DBLE( vec1_d(1) )
         END DO 
       END IF 
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( tr_d, intra_bgrp_comm )
       !
    END IF
    !
    tr = tr_d
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       sr(1:idiis,idiis,ibnd) = tr(1:idiis,jbnd)
       sr(idiis,1:idiis,ibnd) = tr(1:idiis,jbnd)
       !
    END DO
    !
    ! ... Update current wave functions and residual vectors
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       kbnd = ibnd_index(ibnd)
       !
       IF ( idiis > 1 ) THEN
          !
          ! ... solve Rv = eSv
          !
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vr(:) )
          CALL mp_bcast( vr, root_bgrp, intra_bgrp_comm )
          !
          !$acc kernels
          psi(:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          kpsi(:,kbnd) = ZERO
          IF ( uspp ) spsi(:,ibnd) = ZERO
          !$acc end kernels
          !
          DO kdiis = 1, idiis
             !
             ! ... Wave functions
             !
             kvr = vr(kdiis) 
             !$acc host_data use_device(psi, hpsi, spsi)
             CALL DAXPY_gpu( npw2, kvr, phi_d (1,ibnd,kdiis), 1, psi (1,ibnd), 1 )
             CALL DAXPY_gpu( npw2, kvr, hphi_d(1,ibnd,kdiis), 1, hpsi(1,ibnd), 1 )
             IF ( uspp ) &
             CALL DAXPY_gpu( npw2, kvr, sphi_d(1,ibnd,kdiis), 1, spsi(1,ibnd), 1 )
             !$acc end host_data
             !
             ! ... Residual vectors
             !
             er = php(ibnd,kdiis)
             !
             CALL DCOPY_gpu( npw2, hphi_d(1,ibnd,kdiis), 1, vec1_d(1), 1 )
             !
             IF ( uspp ) THEN
                !
                CALL DAXPY_gpu( npw2, -er, sphi_d(1,ibnd,kdiis), 1, vec1_d(1), 1 )
                !
             ELSE
                !
                CALL DAXPY_gpu( npw2, -er, phi_d(1,ibnd,kdiis), 1, vec1_d(1), 1 )
                !
             END IF
             !
             !$acc host_data use_device(kpsi)
             CALL DAXPY_gpu( npw2, kvr, vec1_d(1), 1, kpsi(1,kbnd), 1 )
             !$acc end host_data
             !
          END DO
          !
       ELSE
          !
          ! ... Wave functions
          !
          norm = SQRT( sw(ibnd) )
          !$acc host_data use_device(psi, hpsi, spsi)
          CALL DSCAL_gpu( npw2, 1._DP / norm, psi (1,ibnd), 1 )
          CALL DSCAL_gpu( npw2, 1._DP / norm, hpsi(1,ibnd), 1 )
          IF ( uspp ) &
          CALL DSCAL_gpu( npw2, 1._DP / norm, spsi(1,ibnd), 1 )
          !$acc end host_data
          !
          ! ... Residual vectors
          !
          er = hw(ibnd)
          !
          !$acc host_data use_device(hpsi, kpsi)
          CALL DCOPY_gpu( npw2, hpsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
          !$acc end host_data
          !
          IF ( uspp ) THEN
             !
             !$acc host_data use_device(spsi, kpsi)
             CALL DAXPY_gpu( npw2, -er, spsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !$acc end host_data
             !
          ELSE
             !
!civn: here changed spsi --> psi due to a possible typo in the original version 
             !CALL DAXPY_gpu( npw2, -er, spsi_d(1,ibnd), 1, kpsi_d(1,kbnd), 1 )
             !$acc host_data use_device(psi, kpsi)
             CALL DAXPY_gpu( npw2, -er, psi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !$acc end host_data
             !
          END IF
          !
       END IF
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
          !
          !$acc kernels
          DO ii = 1, 1
            psi (1,ibnd) = CMPLX( DBLE( psi (1,ibnd) ), 0._DP, kind=DP )
            hpsi(1,ibnd) = CMPLX( DBLE( hpsi(1,ibnd) ), 0._DP, kind=DP )
            kpsi(1,kbnd) = CMPLX( DBLE( kpsi(1,kbnd) ), 0._DP, kind=DP )
          END DO 
          IF ( uspp ) THEN 
            DO ii = 1, 1
              spsi(1,ibnd) = CMPLX( DBLE( spsi(1,ibnd) ), 0._DP, kind=DP )
            END DO 
          END IF 
          !$acc end kernels
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( vec1_d )
    DEALLOCATE( vec2_d )
    IF ( idiis > 1 )   DEALLOCATE( vr )
    IF ( motconv > 0 ) DEALLOCATE( tr )
    IF ( motconv > 0 ) DEALLOCATE( tr_d )
    !
    RETURN
    !
  END SUBROUTINE do_diis_gpu
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vr )
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: ibnd
    INTEGER,  INTENT(IN)  :: idiis
    REAL(DP), INTENT(OUT) :: vr(idiis)
    !
    INTEGER               :: info
    INTEGER               :: ndim, kdim
    INTEGER               :: i, imin
    REAL(DP)              :: emin
    REAL(DP)              :: vnrm
    REAL(DP), ALLOCATABLE :: h1(:,:)
    REAL(DP), ALLOCATABLE :: h2(:,:)
    REAL(DP), ALLOCATABLE :: h3(:,:)
    REAL(DP), ALLOCATABLE :: s1(:,:)
    REAL(DP), ALLOCATABLE :: x1(:,:)
    REAL(DP), ALLOCATABLE :: u1(:)
    REAL(DP), ALLOCATABLE :: e1(:)
    INTEGER               :: nwork
    REAL(DP), ALLOCATABLE :: work(:)
    !
    REAL(DP), EXTERNAL    :: DDOT
    !
    ndim  = idiis
    nwork = 3 * ndim
    !
    ALLOCATE( h1( ndim, ndim ) )
    ALLOCATE( h2( ndim, ndim ) )
    ALLOCATE( h3( ndim, ndim ) )
    ALLOCATE( s1( ndim, ndim ) )
    ALLOCATE( x1( ndim, ndim ) )
    ALLOCATE( u1( ndim ) )
    ALLOCATE( e1( ndim ) )
    ALLOCATE( work( nwork ) )
    !
    h1(1:ndim,1:ndim) = hr(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DSYEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    kdim = 0
    !
    x1 = 0._DP
    !
    DO i = 1, ndim
       !
       IF ( e1(i) > eps14 ) THEN
          !
          kdim = kdim + 1
          !
          x1(:,kdim) = s1(:,i) / SQRT(e1(i))
          !
       END IF
       !
    END DO
    !
    IF ( kdim <= 1 ) THEN
       !
       vr        = 0._DP
       vr(idiis) = 1._DP
       !
       GOTO 10
       !
    END IF
    !
    h2 = 0._DP
    !
    CALL DGEMM( 'N', 'N', ndim, kdim, ndim, 1._DP, h1, ndim, x1, ndim, 0._DP, h2, ndim )
    !
    h3 = 0._DP
    !
    CALL DGEMM( 'T', 'N', kdim, kdim, ndim, 1._DP, x1, ndim, h2, ndim, 0._DP, h3, ndim )
    !
    e1 = 0._DP
    !
    CALL DSYEV( 'V', 'U', kdim, h3, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    imin = 1
    emin = e1(1)
    !
    DO i = 2, kdim
       !
       IF ( ABS( e1(i) ) < ABS( emin ) ) THEN
          !
          imin = i
          emin = e1(i)
          !
       END IF
       !
    END DO
    !
    CALL DGEMV( 'N', ndim, kdim, 1._DP, x1, ndim, h3(:,imin), 1, 0._DP, vr, 1 )
    !
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DGEMV( 'N', ndim, ndim, 1._DP, s1, ndim, vr, 1, 0._DP, u1, 1 )
    !
    vnrm = SQRT( DDOT( ndim, vr, 1, u1, 1 ) )
    !
    vr = vr / vnrm
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( x1 )
    DEALLOCATE( u1 )
    DEALLOCATE( e1 )
    DEALLOCATE( work )
    !
    RETURN
    !
  END SUBROUTINE diag_diis
  !
  !
  SUBROUTINE rr_line_search_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER               :: ig
    INTEGER               :: ibnd, jbnd, kbnd
    LOGICAL               :: para_hpsi
    REAL(DP)              :: psir, psii, psi2
    REAL(DP)              :: kdiag, k1, k2
    REAL(DP)              :: x, x2, x3, x4
    REAL(DP), ALLOCATABLE :: ekin(:)
    REAL(DP)              :: a, b
    REAL(DP)              :: ene0, ene1
    REAL(DP)              :: step, norm
    REAL(DP)              :: php, khp, khk
    REAL(DP)              :: psp, ksp, ksk
    REAL(DP), ALLOCATABLE :: hmat(:,:), smat(:,:)
    REAL(DP), ALLOCATABLE :: heig(:), seig(:)
    REAL(DP)              :: c1, c2
    REAL(DP), ALLOCATABLE :: coef(:,:)
    !
    REAL(DP), EXTERNAL    :: gpu_DDOT
!civn 
    INTEGER :: idx
    REAL(DP) :: ekinj
    REAL(DP) :: rtmp
!
    !
    IF ( motconv > 0 ) THEN
       !
       ALLOCATE( ekin( motconv ) )
       ALLOCATE( hmat( 3, motconv ) )
       ALLOCATE( smat( 3, motconv ) )
       ALLOCATE( heig( motconv ) )
       ALLOCATE( seig( motconv ) )
       ALLOCATE( coef( 2, motconv ) )
       !
    END IF
    !
    ! ... Kinetic energy
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       ekinj = 0._DP 
       !
       !$acc kernels
       DO ig = gstart, npw
          ekinj = ekinj + 2._DP * g2kin_d(ig) * & 
              (DBLE ( psi(ig,ibnd) ) * DBLE ( psi(ig,ibnd) ) + AIMAG( psi(ig,ibnd) ) * AIMAG( psi(ig,ibnd) ))
       END DO
       !$acc end kernels
       !
       IF ( gstart == 2 ) THEN
          !
          !$acc kernels
          DO ii = 1, 1
            ekinj = ekinj + g2kin_d(1) * DBLE ( psi(1,ibnd) ) * DBLE ( psi(1,ibnd) ) 
          END DO 
          !$acc end kernels
          !
       END IF
       !
       ekin(jbnd) = ekinj
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( ekin, intra_bgrp_comm )
       !
    END IF
    !
    ! ... Preconditioning vectors : K (H - e S) |psi>
    !
    ! ... G.Kresse and J.Furthmuller, PRB 54, 11169 (1996)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       ekinj = ekin(jbnd)
       !
       kbnd = ibnd_index(ibnd)
       !
       !$acc kernels
       DO ig = 1, npw
          x  = g2kin_d(ig) / ( 1.5_DP * ekinj )
          x2 = x * x
          x3 = x * x2
          x4 = x * x3
          k1 = 27._DP + 18._DP * x + 12._DP * x2 + 8._DP * x3
          k2 = k1 + 16._DP * x4
          kdiag = ( -4._DP / 3._DP / ekinj ) * k1 / k2
          kpsi(ig,kbnd) = kdiag * kpsi(ig,kbnd)
       END DO
       !$acc end kernels
       !
    END DO
    !
    ! ... Share kpsi for all band-groups
    !
    IF ( use_bgrp_in_hpsi .AND. ( .NOT. is_exx ) .AND. ( notconv > 1 ) ) THEN
       !
       para_hpsi = .TRUE.
       !
    ELSE
       !
       para_hpsi = .FALSE.
       !
    END IF
    !
    IF ( ( .NOT. para_hpsi ) .OR. ( notconv /= nbnd ) ) THEN
       !
       DO ibnd = 1, ( ibnd_start - 1)
          !
          idx = ibnd_index(ibnd)
          !
          IF ( .NOT. conv(ibnd) ) THEN 
            !$acc kernels
            kpsi(:,idx) = ZERO
            !$acc end kernels
          END IF 
          !
       END DO
       !
       DO ibnd = ( ibnd_end + 1 ), nbnd
          !
          idx = ibnd_index(ibnd)
          !
          IF ( .NOT. conv(ibnd) ) THEN 
            !$acc kernels
            kpsi(:,idx) = ZERO
            !$acc end kernels
          END IF 
          !
       END DO
       !
       !$acc host_data use_device(kpsi)
       CALL mp_sum( kpsi(:,1:notconv), inter_bgrp_comm )
       !$acc end host_data
       !
    END IF
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN 
      !$acc kernels
      DO ii = 1, notconv 
        kpsi (1,ii) = CMPLX( DBLE( kpsi (1,ii) ), 0._DP, kind=DP )
      END DO 
      !$acc end kernels
    END IF 
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    !$acc host_data use_device(kpsi, hkpsi)
    CALL h_psi_ptr( npwx, npw, notconv, kpsi, hkpsi )
    !$acc end host_data
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    !$acc host_data use_device(kpsi, skpsi)
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, notconv, kpsi, skpsi )
    !$acc end host_data
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
       !
       !$acc kernels
       DO ii = 1, nbnd
         hkpsi(1,ii) = CMPLX( DBLE( hkpsi(1,ii) ), 0._DP, kind=DP )
       END DO 
       IF ( uspp ) THEN 
         DO ii = 1, nbnd
           skpsi(1,ii) = CMPLX( DBLE( skpsi(1,ii) ), 0._DP, kind=DP )
         END DO 
       END IF 
       !$acc end kernels
       !
    END IF
    !
    ! ... Create 2 x 2 matrix
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       kbnd = ibnd_index(ibnd)
       !
       !$acc host_data use_device(psi, hpsi, kpsi, hkpsi)
       php = 2._DP * gpu_DDOT( npw2, psi (1,ibnd), 1, hpsi (1,ibnd), 1 )
       khp = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, hpsi (1,ibnd), 1 )
       khk = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, hkpsi(1,kbnd), 1 )
       !$acc end host_data
       !
       IF ( gstart == 2 ) THEN
          !
          !$acc kernels copyout(rtmp)
          DO ii = 1, 1 
            rtmp = DBLE( psi (1,ibnd) ) * DBLE ( hpsi (1,ibnd) )
          END DO 
          !$acc end kernels
          php = php - rtmp  
          !$acc kernels copyout(rtmp)
          DO ii = 1, 1 
            rtmp = DBLE( kpsi(1,kbnd) ) * DBLE ( hpsi (1,ibnd) ) 
          END DO 
          !$acc end kernels
          khp = khp - rtmp 
          !$acc kernels copyout(rtmp)
          DO ii = 1, 1 
            rtmp = DBLE( kpsi(1,kbnd) ) * DBLE ( hkpsi(1,kbnd) ) 
          END DO 
          !$acc end kernels
          khk = khk - rtmp  
          !
       END IF
       !
       IF ( uspp ) THEN
          !
          !$acc host_data use_device(psi, spsi, kpsi, skpsi)
          psp = 2._DP * gpu_DDOT( npw2, psi (1,ibnd), 1, spsi (1,ibnd), 1 )
          ksp = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, spsi (1,ibnd), 1 )
          ksk = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, skpsi(1,kbnd), 1 )
          !$acc end host_data
          !
          IF ( gstart == 2 ) THEN
             !
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1 
               rtmp = DBLE( psi (1,ibnd) ) * DBLE ( spsi (1,ibnd) )
             END DO 
             !$acc end kernels
             psp = psp - rtmp 
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1 
               rtmp = DBLE( kpsi(1,kbnd) ) * DBLE ( spsi (1,ibnd) )
             END DO 
             !$acc end kernels
             ksp = ksp - rtmp 
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1 
               rtmp = DBLE( kpsi(1,kbnd) ) * DBLE ( skpsi(1,kbnd) ) 
             END DO 
             !$acc end kernels
             ksk = ksk - rtmp 
             !
          END IF
          !
       ELSE
          !
          !$acc host_data use_device(psi, kpsi)
          psp = 2._DP * gpu_DDOT( npw2, psi (1,ibnd), 1, psi (1,ibnd), 1 )
          ksp = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, psi (1,ibnd), 1 )
          ksk = 2._DP * gpu_DDOT( npw2, kpsi(1,kbnd), 1, kpsi(1,kbnd), 1 )
          !$acc end host_data
          !
          IF ( gstart == 2 ) THEN
             !
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1
               rtmp =DBLE( psi (1,ibnd) ) * DBLE ( psi (1,ibnd) ) 
             END DO 
             !$acc end kernels
             psp = psp - rtmp
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1
               rtmp =DBLE( kpsi(1,kbnd) ) * DBLE ( psi (1,ibnd) ) 
             END DO 
             !$acc end kernels
             ksp = ksp - rtmp 
             !$acc kernels copyout(rtmp)
             DO ii = 1, 1
               rtmp =DBLE( kpsi(1,kbnd) ) * DBLE ( kpsi(1,kbnd) ) 
             END DO 
             !$acc end kernels
             ksk = ksk -rtmp 
             !
          END IF
          !
       END IF
       !
       hmat(1,jbnd) = php
       hmat(2,jbnd) = khp
       hmat(3,jbnd) = khk
       !
       smat(1,jbnd) = psp
       smat(2,jbnd) = ksp
       smat(3,jbnd) = ksk
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( hmat, intra_bgrp_comm )
       CALL mp_sum( smat, intra_bgrp_comm )
       !
    END IF
    !
    ! ... Line searching for each band
    !
    IF ( me_bgrp == root_bgrp ) THEN
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          IF ( conv(ibnd) ) CYCLE
          !
          jbnd = jbnd_index(ibnd)
          !
          php = hmat(1,jbnd)
          khp = hmat(2,jbnd)
          khk = hmat(3,jbnd)
          !
          psp = smat(1,jbnd)
          ksp = smat(2,jbnd)
          ksk = smat(3,jbnd)
          IF( psp <= eps16 ) CALL errore( ' rrmmdiagg ', ' psp <= 0 ', 1 )
          !
          norm = psp + 2._DP * ksp * SREF + ksk * SREF * SREF
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 1 )
          !
          ene0 = php / psp
          ene1 = ( php + 2._DP * khp * SREF + khk * SREF * SREF ) / norm
          !
          a = 2._DP * ( khp * psp - php * ksp ) / psp / psp
          b = ( ene1 - ene0 - a * SREF ) / SREF / SREF
          !
          IF( ABS( b ) > eps16 ) THEN
             step  = -0.5_DP * a / b
          ELSE
             IF ( a < 0._DP ) THEN
                step = SMAX
             ELSE
                step = SMIN
             END IF
          END IF
          !
          step  = MAX( SMIN, step )
          step  = MIN( SMAX, step )
          !
          norm  = psp + 2._DP * ksp * step + ksk * step * step
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 2 )
          norm  = SQRT( norm )
          !
          coef(1,jbnd) = 1._DP / norm
          coef(2,jbnd) = step  / norm
          !
          ! ... Update current matrix elements
          !
          c1 = coef(1,jbnd)
          c2 = coef(2,jbnd)
          !
          heig(jbnd) = php * c1 * c1 + 2._DP * khp * c1 * c2 + khk * c2 * c2
          seig(jbnd) = psp * c1 * c1 + 2._DP * ksp * c1 * c2 + ksk * c2 * c2
          !
       END DO
       !
    END IF
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_bcast( coef, root_bgrp, intra_bgrp_comm )
       CALL mp_bcast( heig, root_bgrp, intra_bgrp_comm )
       CALL mp_bcast( seig, root_bgrp, intra_bgrp_comm )
       !
    END IF
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       hw(ibnd) = heig(jbnd)
       sw(ibnd) = seig(jbnd)
       !
    END DO
    !
    ! ... Update current wave functions
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       kbnd = ibnd_index(ibnd)
       !
       c1 = coef(1,jbnd)
       c2 = coef(2,jbnd)
       !
       !$acc host_data use_device(psi, hpsi, spsi, kpsi, hkpsi, skpsi)
       CALL DSCAL_gpu( npw2, c1, psi (1,ibnd), 1 )
       CALL DAXPY_gpu( npw2, c2, kpsi(1,kbnd), 1, psi(1,ibnd), 1 )
       !
       CALL DSCAL_gpu( npw2, c1, hpsi (1,ibnd), 1 )
       CALL DAXPY_gpu( npw2, c2, hkpsi(1,kbnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL DSCAL_gpu( npw2, c1, spsi (1,ibnd), 1 )
          CALL DAXPY_gpu( npw2, c2, skpsi(1,kbnd), 1, spsi(1,ibnd), 1 )
          !
       END IF
       !$acc end host_data
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
          !
          !$acc kernels
          DO ii = 1, 1
             psi (1,ibnd) = CMPLX( DBLE( psi (1,ibnd) ), 0._DP, kind=DP )
             hpsi(1,ibnd) = CMPLX( DBLE( hpsi(1,ibnd) ), 0._DP, kind=DP )
          END DO 
          IF ( uspp ) THEN  
            DO ii = 1, 1
              spsi(1,ibnd) = CMPLX( DBLE( spsi(1,ibnd) ), 0._DP, kind=DP )
            END DO 
          END IF 
          !$acc end kernels
          !
       END IF
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       DEALLOCATE( ekin )
       DEALLOCATE( hmat )
       DEALLOCATE( smat )
       DEALLOCATE( heig )
       DEALLOCATE( seig )
       DEALLOCATE( coef )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE rr_line_search_gpu
  !
  !
  SUBROUTINE eigenvalues( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' rrmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    ! ... Check convergence
    !
    WHERE( btype(1:nbnd) == 1 )
       !
       conv(1:nbnd) = conv(1:nbnd) .OR. ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < ethr ) )
       !
    ELSEWHERE
       !
       conv(1:nbnd) = conv(1:nbnd) .OR. ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < empty_ethr ) )
       !
    END WHERE
    !
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm )
    !
    ! ... Count not converged bands
    !
    notconv = 0
    !
    DO ibnd = 1, nbnd
       !
       IF ( conv(ibnd) ) THEN
          !
          ibnd_index(ibnd) = 0
          !
       ELSE
          !
          notconv = notconv + 1
          !
          ibnd_index(ibnd) = notconv
          !
       END IF
       !
    END DO
    !
    motconv = 0
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) THEN
          !
          jbnd_index(ibnd) = 0
          !
       ELSE
          !
          motconv = motconv + 1
          !
          jbnd_index(ibnd) = motconv
          !
       END IF
       !
    END DO
    !
    ! ... Save current eigenvalues
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE eigenvalues
  !
  !
END SUBROUTINE rrmmdiagg_gpu
