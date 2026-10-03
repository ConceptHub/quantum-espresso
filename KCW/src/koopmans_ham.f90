!!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!-----------------------------------------------------------------------
SUBROUTINE koopmans_ham ()
  !---------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout, ionode
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nkstot, nks
  USE lsda_mod,              ONLY : nspin, lsda, isk
  USE control_kcw,           ONLY : num_wann, Hamlt, nqstot, l_alpha_corr, evc0, &
                                    alpha_final, num_wann_occ, on_site_only, iuwfc_wann, &
                                    nkstot_eff, x_q, spin_component
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : lrwfc, iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE noncollin_module,      ONLY : npol, nspin_lsda, nspin_gga, nspin_mag
  USE mp,                    ONLY : mp_sum, mp_max, mp_min
  USE mp_pools,              ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !! The global index of a local (pool) k-point
  !
  INTEGER :: ik, ik_eff
  !
  REAL(DP), ALLOCATABLE :: et_eff(:,:)
  !
  REAL(DP), ALLOCATABLE :: et_ks_eff(:,:)
  !
  ! The scalar part (independent on k) <rho_0i|v_0i|rho_0i>delta_ij
  COMPLEX(DP) :: deltah_scal (num_wann, num_wann)
  !
  ! The Hamiltonian due to the real contribution to the potential \int f_Hxc(r,r') rho_0i(r')
  COMPLEX(DP) deltah_real (num_wann, num_wann)
  !
  ! the KI hamiltonian, the KI contribution, and the new eigenvectors at a given k-point
  COMPLEX(DP) :: ham(num_wann,num_wann), deltaH(num_wann,num_wann), eigvc(npwx*npol,num_wann)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(num_wann)
  !
  ! The correction to the diagonal term beyond second order
  REAL(DP) :: ddH(num_wann)
  !
  INTEGER :: i, iwann, jwann
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  INTEGER  :: lrwannfc
  !
  !
  IF (on_site_only) WRITE(stdout, '(/,5X, "INFO: Skipping off-diag: only R=0 and i=j")') 
  !
  ! The scalar term R=0 i=j 
  deltah_scal=CMPLX(0.D0,0.D0,kind=DP)
  CALL ham_scalar (deltah_scal)
  ! 
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian:")')
  DO iwann=1, num_wann
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_scal(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  ! ... The correction beyond 2nd order 
  IF (l_alpha_corr) CALL beyond_2nd (deltah_scal, ddH)
  !
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian plus correction:")')
  ham = deltah_scal
  DO iwann=1, num_wann
    ham(iwann,iwann) = ham(iwann,iwann)-ddH(iwann)
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (ham(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  !
  ! The new (KI) eigenvalues, filled only for the k-points owned by this pool
  ! (zero elsewhere) then summed across pools right after the loop.
  ALLOCATE ( et_eff(num_wann, nkstot_eff) )
  ALLOCATE ( et_ks_eff(num_wann, nkstot_eff) )
  et_eff = 0.D0
  et_ks_eff = 0.D0
  !
  DO ik = 1, nks
    !
    IF ( lsda .AND. isk(ik) /= spin_component ) CYCLE
    !
    ik_eff = global_kpoint_index (nkstot, ik) - (spin_component-1)*nkstot_eff
    !
    deltah_real = CMPLX(0.D0, 0.D0, kind = DP)
    !
    IF (.NOT. on_site_only) THEN
       ! General routine for empty state hamiltonian at k
       CALL koopmans_ham_real_k ( ik, ik_eff, deltah_real )
    ELSE
      ! SKIP off-diagonal elements in REAL SPACE (i.e. R/=0 i/=j) 
      ! If only R=0 and i=j, then the "real" contribution for empty states
      ! is <wi| {\int f_Hxc wi^2} |wi> (which is simply -2.D0 times
      ! the scalar controbution -0.5*<wi^2|f_Hxc|wi^2>) 
      DO iwann = num_wann_occ+1, num_wann
        deltah_real(iwann,iwann) = -2.D0*deltah_scal(iwann,iwann)
      ENDDO
    ENDIF
    !
    !
#ifdef DEBUG
    WRITE( stdout,'(/,5X," Real term Hamiltonian ik =:", i5)') ik
    DO iwann=1, num_wann
      WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_real(iwann,jwann), jwann=1,num_wann)
    ENDDO
#endif
    !
    ! The KS hamiltonian in the Wannier Gauge (just to check)
    ham(:,:) = Hamlt(ik_eff,:,:)
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    !
    ! The KI contribution at k
    deltaH = deltah_scal + deltah_real
    !
    ! Apply the screening coefficient
!    DO iwann = 1, num_wann; deltaH(:,iwann) = alpha_final (iwann)* deltaH(:,iwann) ; ENDDO
    DO iwann = 1, num_wann
      DO jwann = iwann , num_wann
        deltaH(iwann,jwann) = deltaH(iwann,jwann)*alpha_final(jwann)
        deltaH(jwann,iwann) = CONJG(deltaH(iwann,jwann))
      ENDDO
    ENDDO

    !
    ! Add to the KS Hamiltonian
    Hamlt(ik_eff,:,:) = Hamlt(ik_eff,:,:) + deltaH(:,:)
    !
#ifdef DEBUG
    WRITE(stdout, '("KI Contribution to the Hamiltonian at k = ", i4)') ik_eff
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(deltaH(iwann,jwann)),AIMAG(deltaH(iwann,jwann)), jwann=1,num_wann)
    ENDDO
    !
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik_eff
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(Hamlt(ik_eff, iwann,jwann)),AIMAG(Hamlt(ik_eff,iwann,jwann)), jwann=1,num_wann)
    ENDDO
#endif
    !
    ! Stash the KS-in-Wannier-gauge eigenvalues for this k-point; the actual
    ! "KS"/"KI" report is printed after the pool loop (see below), once every
    ! process has the full, gathered table - printing here would only ever
    ! reach the log for the k-points owned by ionode's own pool.
    DO iwann = 1, num_wann
      et_ks_eff(iwann,ik_eff) = eigvl(iwann)
    ENDDO
    !
    ehomo_ks = MAX ( ehomo_ks, eigvl(num_wann_occ ) )
    IF (num_wann .gt. num_wann_occ) elumo_ks = MIN ( elumo_ks, eigvl(num_wann_occ+1 ) )
    !
    ham(:,:) = Hamlt(ik_eff,:,:)
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    !
    nbnd = num_wann
    ! Canonical wfc at each k point (overwrite the evc from DFT)
    lrwannfc = num_wann*npwx*npol
    !
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    ! Retrive the ks function at k (in the Wannier Gauge) - LOCAL ik: this pool's own buffer
    CALL ZGEMM( 'N','N', npw*npol, num_wann, num_wann, ONE, evc0, npwx*npol, eigvc, num_wann, &
                 ZERO, evc, npwx*npol )
    lrwfc = nbnd * npwx*npol
    !
    CALL save_buffer ( evc, lrwfc, iuwfc, ik )
    ! LOCAL ik: iuwfc follows the standard (per-pool) PW buffer convention
    !
    !nbnd = num_wann
    DO iwann = 1, nbnd
      et_eff(iwann,ik_eff) = eigvl(iwann)
    ENDDO
    !
    ehomo = MAX ( ehomo, eigvl(num_wann_occ ) )
    IF (num_wann > num_wann_occ) elumo = MIN ( elumo, eigvl(num_wann_occ+1 ) )
    !
  ENDDO
  !
  ! ... Gather across pools the k-point-resolved quantities computed above: each
  ! pool has filled only the entries for the k-points it owns (zero elsewhere, since
  ! Hamlt/et_eff start at zero and each effective k-point index is owned by exactly
  ! one pool), so a sum over inter_pool_comm reconstructs the full array everywhere.
  !
  CALL mp_sum ( Hamlt(1:nkstot_eff,:,:), inter_pool_comm )
  CALL mp_sum ( et_eff, inter_pool_comm )
  CALL mp_sum ( et_ks_eff, inter_pool_comm )
  et(1:num_wann,1:nkstot_eff) = et_eff(1:num_wann,1:nkstot_eff)
  CALL mp_max ( ehomo_ks, inter_pool_comm )
  CALL mp_max ( ehomo, inter_pool_comm )
  CALL mp_min ( elumo_ks, inter_pool_comm )
  CALL mp_min ( elumo, inter_pool_comm )
  !
  ! ... The per-k "KS"/"KI" report, now that every process holds the full,
  ! gathered table: print once, in k-point order, from ionode only (x_q is
  ! already replicated on every process - see compute_map_ikq_single.f90).
  !
  IF (ionode) THEN
    DO ik_eff = 1, nkstot_eff
      WRITE( stdout, 9020 ) ( x_q(i,ik_eff), i = 1, 3 )
      WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (et_ks_eff(iwann,ik_eff)*rytoev, iwann=1,num_wann)
      WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (et_eff(iwann,ik_eff)*rytoev, iwann=1,num_wann)
    ENDDO
  ENDIF
  DEALLOCATE ( et_eff )
  DEALLOCATE ( et_ks_eff )
  !
  IF ( elumo < 1d+6) THEN
     WRITE( stdout, 9042 ) ehomo_ks*rytoev, elumo_ks*rytoev
     WRITE( stdout, 9044 ) ehomo*rytoev, elumo*rytoev
  ELSE
     WRITE( stdout, 9043 ) ehomo_ks*rytoev
     WRITE( stdout, 9045 ) ehomo*rytoev
  END IF
  !
  ! For an isolated molecule the full Hamiltonian is also available
  ik=1
  IF (nkstot/nspin == 1 .AND. .NOT. nspin_mag==4 ) CALL full_ham ( ik )
  ! ... formats
  !
9043 FORMAT(/,8x,'KS       highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS       highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI[2nd]  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI[2nd]  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
  !
  RETURN
  !
  CONTAINS 
  !
  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (deltaH, ddH)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, alpha_final, alpha_final_full
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE  
    !
    COMPLEX(DP), INTENT(INOUT) :: deltaH(num_wann,num_wann)
    REAL(DP), INTENT(OUT) :: ddH(num_wann)
    !
    REAL(DP) :: alpha_(num_wann), alpha_fd
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    REAL(DP) second_der(num_wann), delta 
    !
    INTEGER :: iwann
    !
    ddH = 0.D0
    !
    WRITE( stdout, '(/,5X, "INFO: Correction beyond 2nd order ...",/)')
    IF (lrpa) THEN 
      !
      WRITE(*, '(8X,"INFO: l_alpha_corr and lrpa are NOT consistent.At RPA")')
      WRITE(*, '(8X,"      level there is no contribution beyond 2nd order.")')
      WRITE(*, '(8X,"      Nothing to do here. RETURN")')
      !
      RETURN
      !
    ENDIF 
    !
    DO iwann = 1, num_wann
      second_der(iwann) = -REAL(deltaH(iwann,iwann))
    ENDDO
    !
    !DO iwann = 1, num_wann_occ
    DO iwann = 1, num_wann
      delta =0.D0 
      alpha_(iwann) = alpha_final(iwann) 
      !
      !
       ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
       ! ... value of the energy in the frozen orbital approximation ...
       CALL alpha_corr (iwann, delta)
       ddH(iwann) = delta
       !deltaH(iwann,iwann) = deltaH(iwann,iwann)-ddH(iwann)
       !
       ! ... The new alpha that should be closer to the Finite-difference one ...
       ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
       alpha_fd = (alpha_final(iwann)*second_der(iwann) + delta)/ (second_der(iwann)+delta)
       IF(nkstot/nspin == 1) alpha_final_full(iwann) = alpha_fd
       !
       ! ... Since the ham in the frozen approximation is approximated to second
       ! ... order only, this is the alpha we want to use. Only the
       ! ... numerator matters.
       alpha_(iwann) = (alpha_final(iwann)*second_der(iwann) + delta)/second_der(iwann)
       !
       ! ... Write it just to compare with the FD one from CP ... 
       WRITE(stdout,'(5X, "INFO: iwann, LR-alpha, FD-alpha, alpha", i3, 3f12.8)') iwann, alpha_final(iwann),alpha_fd,  alpha_(iwann)
       !
       !WRITE(stdout,'("Nicola", i3, 6f12.8)') iwann, deltaH(iwann,iwann)
       !
       ! Re-define the corrected screening parameter. 
       alpha_final(iwann) = alpha_(iwann) 
       WRITE( stdout, '(5X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      !
    ENDDO
    !
  END SUBROUTINE beyond_2nd
  !
  !----------------------------------------------------------------
  SUBROUTINE ham_scalar (deltah_scal)
    !----------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : num_wann, nqstot, iurho_wann, &
                                     spin_component, nrho, x_q, kcw_iverbosity
    USE fft_base,             ONLY : dffts
    USE cell_base,            ONLY : omega
    USE gvecs,                ONLY : ngms
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE buffers,              ONLY : get_buffer
    USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
    !
    IMPLICIT NONE
    !
    ! The scalar contribution to the hamiltonian 
    COMPLEX(DP), INTENT(INOUT) :: deltah_scal (num_wann, num_wann)
    !
    ! Couters for the q point, wannier index. record length for the wannier density
    INTEGER :: iq, iwann, lrrho, ip, ii
    !
    ! The periodic part of the wannier orbital density
    !!COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr,nrho)
    COMPLEX(DP), ALLOCATABLE :: rhowann(:, :, :), rhor(:,:)
    !!COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
    COMPLEX(DP), ALLOCATABLE :: delta_vr(:,:), delta_vr_(:,:)
    !
    ! The self Hartree
    COMPLEX(DP) :: sh(num_wann)
    !
    ! Auxiliary variables 
    COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:), wann_c(:,:,:), rho_c(:,:,:)
    !
    ! The weight of each q point
    REAL(DP) :: weight(nqstot)
    !
    COMPLEX(DP) :: struct_fact,int_rho,int_wann, zpom
   !!JA  COMPLEX(DP) :: phase(dffts%nnr)
    REAL(DP) :: xq_(3)
    !
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... START")')
    !
    ALLOCATE ( rhog (ngms,nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
    ALLOCATE(rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr,nrho) )
    ALLOCATE( delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag))
    !
    ALLOCATE ( wann_c(dffts%nnr,num_wann,nrho),rho_c(dffts%nnr,num_wann,nrho) )
    !$acc enter data create(rhor, rhog, vh_rhog, delta_vr, delta_vr_, delta_vg, delta_vg_)
    !
    wann_c = ZERO
    rho_c = ZERO
    sh = ZERO
    !
    DO iq = 1, nqstot
      !
      !!!CALL structure_factor(iq, struct_fact)
      !
      lrrho=num_wann*dffts%nnr*nrho
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      !! Retrive the rho_wann_q(r) from buffer in REAL space
      !
      weight(iq) = 1.D0/nqstot ! No SYMM 
      !
      DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
         !
         !$acc kernels present(rhog, delta_vg, vh_rhog, rhor)
         rhog(:,:)       = (0.0_dp,0.0_dp)
         delta_vg(:,:)   = (0.0_dp,0.0_dp)
         vh_rhog(:)      = (0.0_dp,0.0_dp)
         rhor(:,:)       = (0.0_dp,0.0_dp)
         !$acc end kernels
         !
         rhor(:,:) = rhowann(:,iwann,:)
         !$acc update device(rhor)
         !! The periodic part of the orbital density in real space
         !
         CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
         !! The periodic part of the perturbation DeltaV_q(G)
         ! 
       !!  sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:))*weight(iq)*omega
         zpom = (0.0_dp, 0.0_dp)
         !$acc parallel loop reduction(+:zpom) present(rhog, vh_rhog)
         DO ii = 1, ngms
            zpom = zpom + CONJG(rhog (ii,1)) * vh_rhog(ii)
         END DO   
         sh(iwann) = sh(iwann) + 0.5D0 * zpom*weight(iq)*omega


         IF (nspin_mag==2) THEN
           zpom = (0.0_dp, 0.0_dp)
           !$acc parallel loop reduction(+:zpom) present(rhog, delta_vg)
           DO ii = 1, ngms
               zpom = zpom + CONJG(rhog (ii,1)) * delta_vg(ii,spin_component)
           END DO
           !!deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component)) &
           deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * zpom &
                                       * weight(iq) * omega
         ELSE
           zpom = (0.0_dp, 0.0_dp)
           !$acc parallel loop collapse(2) reduction(+:zpom) present(rhog, delta_vg)
           DO ip =1,nspin_mag
              DO ii = 1, ngms
                 zpom = zpom + CONJG(rhog (ii,ip)) * delta_vg(ii,ip) 
              END DO
           END DO
           deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * zpom &
                                          * weight(iq) * omega
        END IF
      ENDDO
      ! 
    ENDDO ! qpoints

    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... END")')
    !
    !$acc exit data delete(rhor, rhog, vh_rhog, delta_vr, delta_vr_, delta_vg, delta_vg_)
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    DEALLOCATE(rhowann, rhor, delta_vr, delta_vr_)
    !
    CALL mp_sum (deltah_scal, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
    !
  END SUBROUTINE ham_scalar
  !
  !-----------------------------------------------------------------------
  SUBROUTINE koopmans_ham_real_k (ik, ik_eff, deltaH)
    !---------------------------------------------------------------------
    !
    ! This routine compute the KI real term correction to second order to the KS
    ! Hamiltonian at a given k point. ik is the LOCAL (pool) index of this k-point
    ! (used for the per-pool iuwfc_wann buffer and the local igk_k/ngk arrays);
    ! ik_eff is the "effective" index (1:nkstot_eff) used for the k+q/k-q map and
    ! for the iuwfc_wann_allk/igk_k_all/ngk_all buffers, which hold ALL k-points
    ! (not just the ones owned by this pool - needed since k+q, k-q in general
    ! belong to a different pool, see bcast_wfc.f90).
    !
    USE kinds,                ONLY : DP
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE fft_wave,             ONLY : invfft_wave
    USE klist,                ONLY : igk_k, ngk
    USE mp,                   ONLY : mp_sum
    USE control_kcw,          ONLY : spin_component, num_wann, x_q, &
                                     num_wann_occ, evc0, iuwfc_wann, iuwfc_wann_allk, &
                                     igk_k_all, ngk_all, &
                                     map_ikq, shift_1bz, nrho, &
                                     map_ikq_minus, shift_1bz_minus
    USE buffers,              ONLY : get_buffer, save_buffer
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : iurho_wann
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE dv_of_drho_lr,        ONLY : dv_of_drho
    USE control_lr,           ONLY : lgamma
    USE lsda_mod,             ONLY : nspin
    USE gvecs,                ONLY : ngms
    USE solve_linter_koop_mod 
    USE qpoint,               ONLY : xq
    USE wvfct,                ONLY : npwx 
    !
    !USE mp_world,             ONLY : mpime
    !
    USE cell_base,            ONLY : omega
    USE noncollin_module,     ONLY : npol
    !
    USE constants,            ONLY : rytoev
    USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
    ! 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, ik_eff
    ! ik: local (pool) k-point index; ik_eff: effective (1:nkstot_eff) k-point index
    !
    ! The KI real term contribution to the Hamiltonian
    COMPLEX(DP), INTENT (INOUT) :: deltaH(num_wann, num_wann)
    !
    INTEGER :: iq, nqs
    ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
    !
    INTEGER :: iwann, jwann, lrrho, lrwfc
    ! Band counters, leght of the rho record
    !
  !!  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann,nrho), rhor(dffts%nnr,nrho), delta_vr(dffts%nnr,nspin_mag), sh(num_wann), &
  !!                 delta_vr_(dffts%nnr,nspin_mag)
    COMPLEX(DP), ALLOCATABLE :: rhowann(:, :,:), rhor(:,:), delta_vr(:,:), sh(:), &
                   delta_vr_(:,:)

    ! The periodic part of the wannier orbital density in r space
    ! The perturbig potential in real space
    ! The self-hartree 
    ! The perturbig potential in real space (without the g=0 contribution) 
    !
    COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    ! The periodic parte of the wannier density in G-space 
    ! The perturbig potential in G space
    ! The hartree potential due to the wannier in G-space 
    ! The perturbig potential in G space without the g=0 contribution 
    !
    REAL(DP) :: weight(nqstot)
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    INTEGER :: ikq, npw_k, npw_kq
    ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
    INTEGER :: ikq_m, npw_kq_m
    !
    !!COMPLEX(DP) ::  evc0_kq(npwx*npol, num_wann)
    COMPLEX(DP), ALLOCATABLE ::  evc0_kq(:, :)
    ! Auxiliary vector to store the wfc at k+q
    !
    !!COMPLEX(DP) :: rho_r_nm(dffts%nnr,nrho), rho_g_nm(ngms,nrho), aux(dffts%nnr)
    COMPLEX(DP), ALLOCATABLE :: rho_r_nm(:,:), rho_g_nm(:,:), aux(:)
    REAL(DP) :: g_vect(3),g_vect_m(3)
    ! G vector that shift the k+q inside the 1BZ
    !
    !!COMPLEX(DP) :: evc_k_g (npwx*npol), evc_k_r (dffts%nnr,npol), phase(dffts%nnr)
    COMPLEX(DP), ALLOCATABLE :: evc_k_g (:), evc_k_r (:,:), phase(:)
    ! Auxiliary wfc in reciprocal and real space, the phase associated to the hift k+q-> k', and to k-q
    !
    !!COMPLEX(DP) :: evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr,npol)
    COMPLEX(DP), ALLOCATABLE :: evc_kq_g (:), evc_kq_r (:,:)
    ! Auxiliary wfc in reciprocal and real space
    !
    LOGICAL :: off_diag = .TRUE., debug_nc = .true.
    ! compute Off-diagonal elements. NsC: not sure off_diag=.false. here makes sense: DO NOT CHANGE!!!!
    INTEGER :: ip, is, iss, ii
    COMPLEX(dp) :: zpom
    !
    IF (nspin_mag==2 .AND. debug_nc) &
      WRITE (stdout, '(/,5X,"INFO: debug_nc = ", L2, " Note: the k-q formula will be used.")') debug_nc
    WRITE( stdout, '(/,/,5X,"INFO: KI[2nd] HAMILTONIAN CALCULATION ik= ", i4, " ...", /)') ik_eff
    !
    nqs = nqstot
    !
    ALLOCATE(rho_r_nm(dffts%nnr,nrho), rho_g_nm(ngms,nrho), aux(dffts%nnr))
    ALLOCATE(evc_k_g (npwx*npol), evc_k_r (dffts%nnr,npol), phase(dffts%nnr))
    ALLOCATE(evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr,npol) )
    ALLOCATE( evc0_kq(npwx*npol, num_wann) )
    !
    !
    lrwfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
    ! Retrive the ks function at k (in the Wannier Gauge)
    ! IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)') 
    !
    CALL compute_map_ikq_single (ik_eff,.true.)
    ! find the map k+q --> k'+G and store the res 
    ! find also the map k-q --> k'+G and store the res 

    ALLOCATE(rhowann(dffts%nnr, num_wann,nrho), rhor(dffts%nnr,nrho), delta_vr(dffts%nnr,nspin_mag), sh(num_wann), &
                   delta_vr_(dffts%nnr,nspin_mag) )
    !
    ALLOCATE ( rhog (ngms,nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
    !$acc enter data create(rhor, rhog, vh_rhog, delta_vr, delta_vr_, delta_vg, delta_vg_)
    !$acc enter data copyin(igk_k_all)
    !
    deltaH = CMPLX(0.D0,0.D0,kind=DP)
    rho_r_nm = CMPLX(0.D0,0.D0,kind=DP)
    sh     = CMPLX(0.D0,0.D0,kind=DP)
    !
    DO iq = 1, nqs
      !! Sum over the BZ 
      !
      xq(1:3)  = x_q(1:3,iq)
      !
      lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
      !
      lrrho=num_wann* dffts%nnr *nrho
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      ! Retrive the rho_wann_q(r) from buffer in REAL space
      !IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
      !
      !
      weight(iq) = 1.D0/nqs ! No SYMM 
      !
      !WRITE(stdout, '("weight =", i5, f12.8)') iq, weight(iq)
      !
      IF (nspin==4 .or. debug_nc ) THEN 
        ikq_m = map_ikq_minus(iq)
        g_vect_m(:) = shift_1bz_minus(:,iq)
        phase(:) = 0.D0
        CALL calculate_phase(g_vect_m, phase) 
        ! Calculate the phase associated to the k-q-> ikq map: exp[ -i(G_bar * r) ]
      ELSE
        ikq = map_ikq(iq)
        g_vect(:) = shift_1bz(:,iq)
        ! The index ikq in the 1BZ corresponding at the k+q, and the vector G_bar defining the shift 
        ! xk(:,ikq)+G_bar = xk(:,k+q)
        ! see compute_map_ikq_single.f90  
        phase(:) = 0.D0
        CALL calculate_phase(g_vect, phase) 
        ! Calculate the phase associated to the k+q-> ikq map: exp[ -i(G_bar * r) ]
      END IF
      !
      evc0_kq = CMPLX(0.D0,0.D0,kind=DP)
      lrwfc = num_wann * npwx * npol
      IF (nspin==4 .or. debug_nc ) THEN
        CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann_allk, ikq_m )
        !Retrive the ks function at k-q (in the Wannier Gauge): ikq_m may belong to
        ! a different pool, hence the ALL-k buffer (see bcast_wfc.f90)
      ELSE
        CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann_allk, ikq )
        ! Retrive the ks function at k+q (in the Wannier Gauge): ikq may belong to
        ! a different pool, hence the ALL-k buffer (see bcast_wfc.f90)
      ENDIF
      !
      !IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_kq(g) RETRIEVED")') 
      !
      DO iwann = num_wann_occ+1, num_wann
!      DO iwann = 1, num_wann
         !         
         npw_k = ngk(ik)
         evc_k_g(:) =  evc0(:,iwann)
         !$acc enter data copyin(evc_k_g) create(evc_k_r)
         CALL invfft_wave (npwx, npw_k, igk_k (1,ik), evc_k_g , evc_k_r )
        !! !$acc exit data copyout(evc_k_r) delete(evc_k_g)
         !$acc exit data delete(evc_k_g)
         !! The wfc R=0 n=iwann in R-space at k
         !
         !DO jwann = iwann+1, num_wann 
         !$acc enter data create(evc_kq_r, rho_r_nm, rho_g_nm, aux) copyin(phase)
         DO jwann = iwann, num_wann 
            !
            IF (.NOT. off_diag .AND. jwann /= iwann) CYCLE 
            !
            !$acc kernels present(rhog, delta_vg, vh_rhog, rhor)
            rhog(:,:)       = (0.0_dp,0.0_dp)
            delta_vg(:,:)   = (0.0_dp,0.0_dp)
            vh_rhog(:)      = (0.0_dp,0.0_dp)
            rhor(:,:)       = (0.0_dp,0.0_dp)
            !$acc end kernels
            !
            rhor(:,:) = rhowann(:,jwann,:)
            !$acc update device(rhor)
            ! The periodic part of the orbital density R=0, n=iwann in real space
            !
            CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
            ! The periodic part of the perturbation DeltaV_q(G)
            !
            evc_kq_g = evc0_kq(:,jwann)
            !$acc enter data copyin(evc_kq_g) 
            IF(nspin==4 .or. debug_nc ) THEN
              npw_kq_m = ngk_all(ikq_m)
              CALL invfft_wave (npwx, npw_kq_m, igk_k_all (1,ikq_m), evc_kq_g , evc_kq_r )
            ELSE
              npw_kq = ngk_all(ikq)
              CALL invfft_wave (npwx, npw_kq, igk_k_all (1,ikq), evc_kq_g , evc_kq_r )
            END IF
            !! !$acc exit data copyout(evc_kq_r) delete(evc_kq_g)
            !$acc exit data delete(evc_kq_g)

            ! The wfc in R-space at k' <-- k+q where k' = (k+q)-G_bar
            ! evc_k+q(r) = sum_G exp[iG r] c_(k+q+G) = sum_G exp[iG r] c_k'+G_bar+G 
            !            = exp[-iG_bar r] sum_G' exp[iG'r] c_k'+G' = exp[-iG_bar r] *evc_k'(r)
            !

            IF (nspin==4) THEN
              !rho_r_nm(:,1) = rho_r_nm(:,1) + conjg(evc_k_r(:,2))*evc_kq_r(:,2)*phase(:)/nqs 
              !rho_r_nm(:,1) = CONJG(rho_r_nm(:,1))
              !rho_r_nm(:,2) = (conjg(evc_k_r(:,1))*evc_kq_r(:,2)+conjg(evc_k_r(:,2))*evc_kq_r(:,1))*phase(:)/nqs 
              !rho_r_nm(:,3) = CMPLX(0.D0,1.D0, kind=DP)*(conjg(evc_k_r(:,2))*evc_kq_r(:,1)-&
              !                                                          conjg(evc_k_r(:,1))*evc_kq_r(:,2))*phase(:)/nqs 
              !rho_r_nm(:,4) = (conjg(evc_k_r(:,1))*evc_kq_r(:,1)-conjg(evc_k_r(:,2))*evc_kq_r(:,2))*phase(:)/nqs
              !
              !phase = conjg(phase)
              ! \sum_{s1,s2} [u^{s1}_{k-q}(r)]^* \sigma_i^{s1,s2} [u^{}s2_{k}(r)] = \rho^i_{k-q,k}(r)
              ! i is the pauli matrices index i=0,x,y,z
              !$acc kernels present(rho_r_nm, evc_k_r, evc_kq_r, phase)
              rho_r_nm(:,1) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:) + conjg(evc_kq_r(:,2))*evc_k_r(:,2)*phase(:) )/nqs 
              rho_r_nm(:,2) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,2)*phase(:) + conjg(evc_kq_r(:,2))*evc_k_r(:,1)*phase(:) )/nqs 
              rho_r_nm(:,3) = (-CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_kq_r(:,1))*evc_k_r(:,2)*phase(:) & 
                               +CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_kq_r(:,2))*evc_k_r(:,1)*phase(:) )/nqs
              rho_r_nm(:,4) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:) - conjg(evc_kq_r(:,2))*evc_k_r(:,2)*phase(:) )/nqs
              !$acc end kernels
            ELSE IF (nspin==2 .and. debug_nc ) THEN
              !phase = conjg(phase)
              !$acc kernels present(rho_r_nm, evc_k_r, evc_kq_r, phase)
              rho_r_nm(:,1) = conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:)/nqs 
              !$acc end kernels
              ! rho_{k-q,k}^{ji}
            ELSE IF (nspin==2 .and. .not. debug_nc) THEN
              !$acc kernels present(rho_r_nm, evc_k_r, evc_kq_r, phase)
              rho_r_nm(:,1) = conjg(evc_k_r(:,1))*evc_kq_r(:,1)*phase(:)/nqs 
              !$acc end kernels
              ! rho_{k,k+q}^{ij}
            ENDIF

            ! generalized density in r-space 
  !          IF (jwann == iwann) THEN
  !            WRITE(*,'("NICOLA evc i", 6F20.15)') evc_k_r(1:3)
  !            WRITE(*,'("NICOLA evc j", 6F20.15)') evc_kq_r(1:3)
  !            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_r_nm(1:3)
  !            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhor(1:3)
  !          ENDIF
            !WRITE(*,'("NICOLA R", 2i5, 2F20.15)'), iwann, jwann, SUM( delta_vr(:,spin_component)*rho_r_nm(:) )/( dffts%nr1*dffts%nr2*dffts%nr3 )
            !
            DO is=1,nrho
                !$acc kernels present(aux, rho_r_nm)
                aux(:) = rho_r_nm(:,is)/omega    
                !$acc end kernels
                !$acc host_data use_device(aux)
                CALL fwfft ('Rho', aux, dffts)
                !$acc end host_data
                !$acc kernels present(rho_g_nm, aux, dffts, dffts%nl) 
               !! rho_g_nm(:,is) = aux(dffts%nl(:))
                rho_g_nm(:,is) = aux(dffts%nl(:))
                !$acc end kernels
            END DO
            ! generalized density in G-spage
  !          IF (jwann == iwann) THEN
  !            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_g_nm(1:3)
  !            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhog(1:3)
  !            WRITE(*,'("NICOLA vi_g", 6F20.15)') delta_vg(1:3,spin_component)
  !            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
  !            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
  !          ENDIF
            !
            IF (nspin==2 .and. .not. debug_nc) THEN
              !!deltaH(iwann, jwann) = deltaH(iwann,jwann) + SUM((rho_g_nm(:,1))*CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
              zpom = (0.0_dp, 0.0_dp)
              !$acc parallel loop reduction(+:zpom) present(delta_vg, rho_g_nm)
              DO ii =1, ngms
                 zpom = zpom + rho_g_nm(ii,1)*CONJG(delta_vg(ii,spin_component))
              END DO   
              deltaH(iwann, jwann) = deltaH(iwann,jwann) + zpom*weight(iq)*omega

            ELSE IF (nspin==2 .and. debug_nc ) THEN 
             !! deltaH(iwann, jwann) = deltaH(iwann,jwann) + SUM(delta_vg(:,spin_component)*conjg(rho_g_nm(:,1)))*weight(iq)*omega
              zpom = (0.0_dp, 0.0_dp)
              !$acc parallel loop reduction(+:zpom) present(delta_vg, rho_g_nm)
              DO ii =1, ngms
                 zpom = zpom+ delta_vg(ii,spin_component)*conjg(rho_g_nm(ii,1))
              END DO   
              deltaH(iwann, jwann) = deltaH(iwann,jwann) + zpom*weight(iq)*omega

            ELSE
              zpom = (0.0_dp, 0.0_dp)
              !$acc parallel loop collapse(2) reduction(+:zpom) present(delta_vg, rho_g_nm)
              DO is=1,nspin_mag
              DO ii = 1, ngms
                  zpom = zpom+ CONJG(rho_g_nm(ii,is))*(delta_vg(ii,is))
              END DO    
              END DO
              !!deltaH(iwann, jwann) = deltaH(iwann,jwann) + SUM(CONJG(rho_g_nm(:,is))*(delta_vg(:,is)))*weight(iq)*omega
              deltaH(iwann, jwann) = deltaH(iwann,jwann) + zpom*weight(iq)*omega
            END IF
            !deltaH(jwann, iwann) = deltaH(jwann,iwann) + SUM(rho_g_nm(:)*CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
            !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
            !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (rho_g_nm (:) * CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
            deltaH(jwann,iwann) = CONJG(deltaH(iwann,jwann))
            !
         ENDDO ! jwann
         !$acc exit data delete(evc_kq_r, rho_r_nm, rho_g_nm, aux, phase)
         !$acc exit data delete(evc_k_r) 
         ! 
         !
      ENDDO ! iwann
      !
      !
      !    
    ENDDO ! qpoints
    !$acc exit data delete(rhor, rhog, delta_vg, vh_rhog, delta_vg_, delta_vr, delta_vr_)
    !$acc exit data delete(igk_k_all)

    !WRITE( stdout, '(5X,"INFO: KC HAMILTONIAN CALCULATION ik= ", i4, " ... DONE")') ik
    !
    deltaH = nqstot*deltaH
    CALL mp_sum (deltaH, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
    ! Sum over different processes (G vectors) 
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    DEALLOCATE(rhowann, rhor, delta_vr, sh, delta_vr_ )

    DEALLOCATE(rho_r_nm, rho_g_nm, aux)
    DEALLOCATE(evc_k_g, evc_k_r, phase)
    DEALLOCATE(evc_kq_g , evc_kq_r, evc0_kq )

    RETURN 
    !
  END subroutine 
  ! 
END SUBROUTINE koopmans_ham
