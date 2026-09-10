  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !--------------------------------------------------------------------------
  MODULE polaron_hamiltonian
  !--------------------------------------------------------------------------
  !!
  !! Builds B and the polaron Hamiltonian, and applies H, S and the preconditioner
  !! for the iterative solvers. Solving lives in polaron_diagonalization.
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : test_tags_plrn, nbnd_plrn, nbnd_g_plrn, lword_h, lword_g, &
                             hblocksize, select_bands_plrn, kpg_map, etf_all, Hamil,   &
                             eigvec, Bmat, gq_model, epf, epfall
  USE polaron_grid,   ONLY : ikq_all, ikqLocal2Global, isGVec
  USE buffers,        ONLY : get_buffer, save_buffer
  USE io_var,         ONLY : iepfall, ihamil

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: build_plrn_bmat, build_plrn_hamil, h_psi_plrn, s_psi_plrn, g_psi_plrn
  PUBLIC :: get_cfac, norm_plrn_wf

  CONTAINS

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_bmat(bqv)
    !-----------------------------------------------------------------------
    !! Create the Bmat
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nqtotf, xqf, nktotf
    USE input,         ONLY : model_vertex_plrn, io_lvl_plrn,               &
                              g_start_energy_plrn, g_end_energy_plrn,       &
                              g_start_band_plrn, model_vertex_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, eps8
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(out) :: bqv(:, :)
    !! Polaron displacement coefficients in phonon basis, Bqv
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: iqpg
    !! Mirror q-point index
    REAL(KIND = DP) :: eig
    !! KS eigenvalue
    !
    bqv = czero
    DO iq = 1, nqtotf
      IF (model_vertex_plrn) THEN
        epf = czero
        epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
      ELSE
        IF (io_lvl_plrn == 0) THEN
          epf(:, :, :, :) = epfall(:, :, :, :, iq)
        ELSE IF (io_lvl_plrn == 1) THEN
          CALL get_buffer(epf, lword_g, iepfall, iq)
        ENDIF
      ENDIF
      IF (test_tags_plrn(1)) THEN
        epf = 1E-6
        epf(:, :, nmodes, :) = 2E-2
      ELSE IF (test_tags_plrn(2)) THEN
        epf = ABS(epf)
      ELSE IF (test_tags_plrn(3)) THEN
        IF (NORM2(xqf(:,iq)) > 1E-5) epf(:, :, nmodes, :) = 0.005 / NORM2(xqf(:, iq))
      ENDIF
      !
      ! energy cutoff for g
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_g_plrn
          eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
          IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
            epf(ibnd, :, :, ik) = czero
            epf(:, ibnd, :, ik) = czero
          ENDIF
        ENDDO
      ENDDO
      !
      iqpg = kpg_map(iq)
      ! if iq is the gamma point, the first three modes should be
      ! dropped because wf(q=0, 1:3) will be zero
      IF (isGVec(xqf(1:3, iq))) THEN
        start_mode = 4
      ELSE
        start_mode = 1
      ENDIF
      !
      IF ( iq > iqpg ) THEN
        DO inu = start_mode, nmodes!
          ! Enforce the relation B_q = B*_{G-q}
          bqv(iq, inu) = CONJG(bqv(iqpg, inu))
        ENDDO
      ELSE
        DO inu = start_mode, nmodes
          bqv(iq, inu)   = cal_Bmat(iq, inu)
        ENDDO
      ENDIF
      !TODO: to be consistent with Denny, whether this is correct?
      ! IF ( ABS(wf(1, iq)) < eps2 ) bqv(iq, 1) = czero
    ENDDO
    ! cal_bqv only sum over local k, so we have to do mp_sum
    CALL mp_sum(bqv, inter_pool_comm )
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE build_plrn_bmat
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_hamil(bqv)
    !-----------------------------------------------------------------------
    !! Build the effective polaron Hamiltonian,
    !! Eq.(61) of PRB 99, 235139 (2019)
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nqtotf, xqf, nktotf
    USE input,         ONLY : model_vertex_plrn, io_lvl_plrn,  type_plrn,           &
                              nhblock_plrn, g_start_energy_plrn, g_end_energy_plrn, &
                              g_start_band_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, eps8, twopi, ci
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: bqv(:,:)
    !! FIXME
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Globar k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: index_blk
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    REAL(KIND = DP) :: eig
    !! KS eigenvalue
    COMPLEX(KIND = DP) :: ctemp
    !! Prefactor
    !
    test_tags_plrn(1) = .FALSE.
    test_tags_plrn(2) = .FALSE.
    test_tags_plrn(3) = .FALSE.
    !
    ! Calculate the Hamiltonian with Bq $$H_{n\bk,n'\bk'} = \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} -\frac{2}{N_p} \sum_{\nu} B^*_{\bq,\nu}g_{nn'\nu}(\bk',\bq)$$
    ! H_{n\bk,n'\bk'} -> Hamil(ik, ibnd, ikq, jbnd)
    ! B^*_{\bq,\nu} -> conj(bqv(iq, inu))
    ! g_{nn'\nu}(\bk',\bq) -> epf(ibnd, jbnd, inu, ikq, iq)
    ! if q == 0, \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} is diagonal matrix with \varepsilon_{n\bk}
    !
    ! G == (0,0,0) means this is the diagonal term with k=k'
    ! ikq is the global index, i.e. the second index
    ! ik is the local index, i.e. the first index
    Hamil = czero
    DO iq = 1, nqtotf
      IF (model_vertex_plrn) THEN
        epf = czero
        epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
      ELSE
        CALL start_clock('read_gmat') !nbndsub*nbndsub*nmodes*nkf
        IF(io_lvl_plrn == 0) THEN
          epf(:, :, :, :) = epfall(:, :, :, :, iq)
        ELSE IF (io_lvl_plrn == 1) THEN
          CALL get_buffer(epf, lword_g, iepfall, iq)
        ENDIF
        CALL stop_clock('read_gmat')
      ENDIF
      ! energy cutoff of g
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_g_plrn
          eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
          IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
            epf(ibnd, :, :, ik) = czero
            epf(:, ibnd, :, ik) = czero
          ENDIF
        ENDDO
      ENDDO
      !
      DO ik = 1, nkf
        ikq = ikq_all(ik, iq)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          !
          IF (nhblock_plrn == 1) THEN
            index_loc = indexkn1
            index_blk = 1
          ELSE
            index_loc = MOD(indexkn1 - 1, hblocksize) + 1
            index_blk = INT((indexkn1 - 1) / hblocksize) + 1
          ENDIF
          !
          IF (iq /= 1 .AND. index_loc == 1 .AND. nhblock_plrn /= 1) THEN
            CALL start_clock('read_Hmat')
            CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
            CALL stop_clock('read_Hmat')
          ENDIF
          !
          IF (isGVec(xqf(1:3, iq))) THEN
            ! Note that, ik is local index while ikq is global index,
            ! so even when q=0, ik \= ikq, but ik_global == ikq
            ! delta_{nn' kk'} epsilon_{nk}
            ctemp = etf_all(select_bands_plrn(ibnd), ikq)
            indexkn2 = (ikq - 1) * nbnd_plrn + ibnd
            Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
          ENDIF
          !
          DO jbnd = 1, nbnd_plrn
            indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
            DO inu = 1, nmodes
              ctemp = type_plrn * two / REAL(nqtotf, KIND = DP) * (bqv(iq, inu)) * &
                 CONJG(epf(select_bands_plrn(jbnd) - g_start_band_plrn + 1, &
                       select_bands_plrn(ibnd)- g_start_band_plrn + 1, inu, ik))
              Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
            ENDDO
          ENDDO
          IF (nhblock_plrn /= 1) THEN
            IF ((index_loc == hblocksize .OR. indexkn1 == nkf * nbnd_plrn)) THEN
              CALL start_clock('Write_Hmat')
              CALL save_buffer(Hamil, lword_h, ihamil, index_blk)
              Hamil = czero
              CALL stop_clock('Write_Hmat')
            ENDIF
          ENDIF
        ENDDO !ibnd
      ENDDO !ik
    ENDDO ! iq
    !-----------------------------------------------------------------------
    END SUBROUTINE build_plrn_hamil
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE h_psi_plrn(lda, n, m, psi, hpsi)
    !-----------------------------------------------------------------------
    ! Calculate Hpsi with psi as input to use the diagon solver in KS_solver cegterg
    ! cegterg take two external subroutine to calculate Hpsi and Spsi to calculate
    ! ( H - e S ) * evc = 0, since H and S is not saved due to their sizes
    ! Hamil need to be passed to h_psi because the parameter space is fixed
    ! to meet the requirement of Davidson diagonalization.
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nktotf
    USE input,         ONLY : type_plrn, nhblock_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, ci
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE ep_constants,  ONLY : czero
    USE global_var,    ONLY : nkf, nktotf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    !! leading dimension of arrays psi, spsi, hpsi, which is nkf * nbnd_plrn
    INTEGER, INTENT(in) :: n
    !! true dimension of psi, spsi, hpsi
    INTEGER, INTENT(in) :: m
    !! number of states psi
    COMPLEX(KIND = DP), INTENT(inout) :: psi(lda, m)
    !! the wavefunction
    COMPLEX(KIND = DP), INTENT(out) :: hpsi(lda, m)
    !! FIXME
    !
    ! Local variables
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: jbnd
    !! Electron band index
    INTEGER :: indexkn1
    !! Combined electron and k-point index
    INTEGER :: indexkn2
    !! Combined electron and k-point index
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    !
    ! Gather psi (dimension nkf) to form eigvec (dimension nktotf)
    CALL start_clock('cal_hpsi')
    IF (lda < nkf * nbnd_plrn) CALL errore('h_psi_plrn', 'leading dimension of arrays psi is not correct', 1)
    eigvec = czero
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigvec(indexkn2, 1:m) = psi(indexkn1, 1:m)
      ENDDO
    ENDDO
    CALL mp_sum(eigvec, inter_pool_comm)
    !
    ! Iterative diagonalization only get the lowest eigenvalues,
    ! however, we will need the highest eigenvalues if we are calculating hole polaron
    ! so, we multiply hpsi by -1, and get eigenvalues in diagonalization
    ! and then multiply eigenvalues by -1.
    hpsi(1:lda, 1:m) = czero
    DO ik = 1, nkf
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        IF(nhblock_plrn == 1) THEN
          index_loc = indexkn1
          index_blk = 1
        ELSE
          index_loc = MOD(indexkn1 - 1, hblocksize) + 1
          index_blk = INT((indexkn1 - 1) / hblocksize) + 1
        ENDIF
        IF (index_loc == 1 .AND. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
        DO ikq = 1, nktotf
          DO jbnd = 1, nbnd_plrn
            indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
            hpsi(indexkn1, 1:m) = hpsi(indexkn1, 1:m) - &
               type_plrn * Hamil(indexkn2, index_loc) * eigvec(indexkn2, 1:m)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL stop_clock('cal_hpsi')
    !-----------------------------------------------------------------------
    END SUBROUTINE h_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE s_psi_plrn(lda, n, m, psi, spsi)
    !-----------------------------------------------------------------------
    !! FIXME ???
    !! SP - what is the point of this subroutine ? Can it be removed ?
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    !! FIXME
    INTEGER, INTENT(in) :: n
    !! FIXME
    INTEGER, INTENT(in) :: m
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: psi(lda,m)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: spsi(lda,m)
    !! FIXME
    CALL errore('s_psi_plrn', "WARNING: This function should not be called at all!", 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE s_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE g_psi_plrn(lda, n, m, npol, psi, e)
    !-----------------------------------------------------------------------
    !! This routine computes an estimate of the inverse Hamiltonian
    !! and applies it to m wavefunctions.
    !! SP: This routines does nothing !!
    !! FIXME
    !
    IMPLICIT NONE
    !
    INTEGER     :: lda, n, m, npol
    COMPLEX(KIND = DP) :: psi(lda, npol, m)
    REAL(KIND = DP)    :: e(m)
    !-----------------------------------------------------------------------
    END SUBROUTINE g_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE get_cfac(xk, nrr_k, irvec_r, cfac)
    !-----------------------------------------------------------------------
    !! Compute the exponential factor.
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : twopi, ci, czero
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    REAL(KIND = DP), INTENT(in) :: xk(3)
    !! k-point coordinates
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: cfac(nrr_k)
    !! Exponential prefactor
    !
    ! Local Variables
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! Dot product between k-point and R WS vector
    !
    cfac = czero
    rdotk = czero
    !
    CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xk, 1, 0.0_dp, rdotk, 1 )
    cfac(:) = EXP(ci * rdotk(:))
    !-----------------------------------------------------------------------
    END SUBROUTINE
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    FUNCTION cal_Bmat(iq, inu)
    !-----------------------------------------------------------------------
    !!
    !! This function calculates the Bq matrix:
    !! B_{qu} = 1/N_p \sum_{mnk} A^*_{mk+q}A_{nk} [g_{mnu}(k,q)/\hbar\omega_{qu}]
    !! Eq.(38) of PRB 99, 235139 (2019)
    !!
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nktotf, wf, nqtotf
    USE input,         ONLY : eps_acoustic, &
                              g_start_band_plrn, istate_relax_plrn
    USE ep_constants,  ONLY : czero, one, eps2, cone, eps8
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q-point counter
    INTEGER, INTENT(in) :: inu
    !! Phonon mode counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: jbnd
    !! Electron-band index
    !INTEGER :: iplrn
    !!! Polaron state index
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    COMPLEX(KIND = DP) :: cal_Bmat
    !! Polaron displacement coefficients in phonon basis, Bqv
    COMPLEX(KIND = DP) :: prefac
    !! Prefactor variable
    !
    ! sum k = sum 1 to nkf + mp_sum (inter_pool)
    ! mp_sum is in polaron_scf
    cal_Bmat = czero
    DO ik = 1, nkf
      ikq = ikq_all(ik, iq)
      ik_global = ikqLocal2Global(ik, nktotf)
      ! TODO : what should do for iplrn?
      ! KL: we should not run the loop over iplrn and
      !     istate_relax_plrn is used to select an excited-state polaron
      !DO iplrn = 1, 1
        DO ibnd = 1, nbnd_plrn
          DO jbnd = 1, nbnd_plrn
            indexkn1 = (ikq - 1) * nbnd_plrn + ibnd
            indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd
            IF (wf(inu, iq) > eps_acoustic ) THEN
              prefac = cone / (wf(inu, iq) * REAL(nqtotf, DP))
            ELSE
              prefac = czero
            ENDIF
            ! B_{q\nu} = \frac{1}{N_p}\sum_{nn'k}A^*_{n'k+q}\frac{g_{n'n\nu}(k, q)}{\hbar \omega_{q\nu}} A_{nk}
            ! cal_Bmat = cal_Bmat + prefac * (eigvec(indexkn2, iplrn)) * CONJG(eigvec(indexkn1, iplrn)) * &
            cal_Bmat = cal_Bmat + prefac * (eigvec(indexkn2, istate_relax_plrn)) * CONJG(eigvec(indexkn1, istate_relax_plrn)) * &
               (epf(select_bands_plrn(ibnd) - g_start_band_plrn + 1, &
               select_bands_plrn(jbnd) - g_start_band_plrn + 1, inu, ik)) !conjg
          ENDDO
        ENDDO
      !ENDDO
    ENDDO
    ! JLB - discard zero or imaginary frequency modes
    IF (wf(inu, iq) < eps_acoustic) THEN
      cal_Bmat = czero
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION cal_Bmat
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    SUBROUTINE norm_plrn_wf(eigvec_coef, norm_new)
    !----------------------------------------------------------------------
    !! Computes the norm of the polaron wavefunction
    !----------------------------------------------------------------------
    USE global_var,    ONLY : nktotf
    USE input,         ONLY : nstate_plrn
    USE ep_constants,  ONLY : czero, one, two, cone
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: norm_new
    !! Wave function normalization
    COMPLEX(KIND = DP), INTENT(inout) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    !
    ! Local variable
    INTEGER :: iplrn
    !! Polaron state counter
    REAL(KIND = DP) :: norm
    !! Normalization
    !
    DO iplrn = 1, nstate_plrn
      norm = REAL(DOT_PRODUCT(eigvec_coef(1:nbnd_plrn * nktotf, iplrn), eigvec_coef(1:nbnd_plrn * nktotf, iplrn)))
      eigvec_coef(:, iplrn) = eigvec_coef(:, iplrn) / DSQRT(norm) * SQRT(norm_new)
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE norm_plrn_wf

  END MODULE polaron_hamiltonian
