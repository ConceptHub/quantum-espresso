!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE ctb09
  !--------------------------------------------------------------------------
  !! Mixing parameter (c) of the Tran-Blaha (TB09) meta-GGA exchange
  !! potential [F. Tran and P. Blaha, PRL 102, 226401 (2009)]
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: tb09_c
  !
CONTAINS
!
!----------------------------------------------------------------------------
FUNCTION tb09_c( nnr, nspin, rho, rho_core, grho ) RESULT( c )
  !--------------------------------------------------------------------------
  !! Compute the mixing parameter as a functional of the density.
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE fft_base,  ONLY : dfftp
  USE mp,        ONLY : mp_sum
  USE mp_bands,  ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nnr
  !! number of points of the dense FFT grid, local to this process
  INTEGER, INTENT(IN) :: nspin
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho(nnr)
  !! total valence charge density, in real space
  REAL(DP), INTENT(IN) :: rho_core(nnr)
  !! core charge of the nonlinear core correction, in real space
  REAL(DP), INTENT(IN) :: grho(3,nnr,nspin)
  !! gradient of the spin densities, core charge included
  REAL(DP) :: c
  !! the mixing parameter
  !
  ! ... local variables
  !
  INTEGER  :: k
  REAL(DP) :: gx, gy, gz, grho2, gbar
  REAL(DP), PARAMETER :: alpha = -0.012_DP, beta = 1.023_DP
  !! intercept and slope, in bohr^(1/2), of the c(g) fit
  REAL(DP), PARAMETER :: eps12 = 1.0E-12_DP
  !
  gbar = 0._DP
  !
  ! ... numerator and denominator must contain the same density the
  ! ... functional is evaluated with, that is, the valence density plus
  ! ... the core charge of the nonlinear core correction
  !
  IF ( nspin == 1 ) THEN
    !$acc parallel loop reduction(+:gbar) present(rho,rho_core,grho) &
    !$acc& private(gx,gy,gz,grho2)
    DO k = 1, nnr
      gx = grho(1,k,1)
      gy = grho(2,k,1)
      gz = grho(3,k,1)
      grho2 = gx*gx + gy*gy + gz*gz
      gbar = gbar + SQRT(grho2) / MAX( ABS(rho(k)+rho_core(k)), eps12 )
    ENDDO
  ELSE
    !$acc parallel loop reduction(+:gbar) present(rho,rho_core,grho) &
    !$acc& private(gx,gy,gz,grho2)
    DO k = 1, nnr
      gx = grho(1,k,1) + grho(1,k,2)
      gy = grho(2,k,1) + grho(2,k,2)
      gz = grho(3,k,1) + grho(3,k,2)
      grho2 = gx*gx + gy*gy + gz*gz
      gbar = gbar + SQRT(grho2) / MAX( ABS(rho(k)+rho_core(k)), eps12 )
    ENDDO
  ENDIF
  !
  CALL mp_sum( gbar, intra_bgrp_comm )
  !
  gbar = gbar / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  c = alpha + beta * SQRT( MAX(gbar,0._DP) )
  !
  WRITE( stdout, '(/5X,"TB09: <|grad rho|/rho> = ",F12.8,", c = ",F12.8)' ) &
         gbar, c
  !
  RETURN
  !
END FUNCTION tb09_c
!
END MODULE ctb09
