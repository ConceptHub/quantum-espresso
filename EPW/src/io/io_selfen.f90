  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE io_selfen
  !----------------------------------------------------------------------
  !!
  !! This module contains various writing or reading routines related to self-energies.
  !! Most of them are for restart purposes.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_ph_write()
    !----------------------------------------------------------------------------
    !!
    !! SP: Added lambda and phonon lifetime writing to file.
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : gtemp, nqtotf, lambda_all, wf, gamma_all
    USE ep_constants,  ONLY : ryd2ev, kelvin2eV, ryd2mev
    USE input,         ONLY : nstemp, lsda
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_var,        ONLY : lambda_phself, linewidth_phself
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    ! Local variable
    CHARACTER(LEN = 20) :: tp
    !! string for temperature
    CHARACTER(LEN = 256) :: filephselfen
    !! file name of phonon selfenergy
    CHARACTER(LEN = 30)  :: myfmt
    !! Variable used for formatting output
    CHARACTER(LEN = 256) :: filephlinewid
    !! file name of phonon linewidth
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name

    INTEGER :: itempphen
    !! Temperature counter for writing phonon selfen
    INTEGER :: iqq
    !! Counter on coarse q-point grid
    INTEGER :: imode
    !! Counter on mode
    !
    fnm = ''
    IF (TRIM(lsda) == 'down') fnm = '.down'
    IF (mpime == ionode_id) THEN
      !
      DO itempphen = 1, nstemp
        WRITE(tp, "(f8.3)") gtemp(itempphen) * ryd2ev / kelvin2eV
        filephselfen = 'lambda.phself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
        OPEN(UNIT = lambda_phself, FILE = filephselfen)
        WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
        WRITE(lambda_phself, *) '#Modes     ',(imode, imode = 1, nmodes)
        DO iqq = 1, nqtotf
          !
          myfmt = "(1000(3x,E15.5))"
          WRITE(lambda_phself,'(i9,4x)', ADVANCE = 'no') iqq
          WRITE(lambda_phself, FMT = myfmt) (REAL(lambda_all(imode, iqq, 1, itempphen)), imode = 1, nmodes)
          !
        ENDDO
        CLOSE(lambda_phself)
        !
        ! SP - 03/2019
        ! \Gamma = 1/\tau = phonon lifetime
        ! \Gamma = - 2 * Im \Pi^R where \Pi^R is the retarted phonon self-energy.
        ! Im \Pi^R = pi*k-point weight*[f(E_k+q) - f(E_k)]*delta[E_k+q - E_k - w_q]
        ! Since gamma_all = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q] we have
        ! \Gamma = 2 * gamma_all
        filephlinewid = 'linewidth.phself.' // trim(adjustl(tp)) // 'K'// TRIM(fnm)
        OPEN(UNIT = linewidth_phself, FILE = filephlinewid)
        WRITE(linewidth_phself, '(a)') '# Phonon frequency and phonon lifetime in meV '
        WRITE(linewidth_phself, '(a)') '# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)'
        DO iqq = 1, nqtotf
          DO imode = 1, nmodes
            WRITE(linewidth_phself, '(i9,i6,E20.8,E22.10)') iqq, imode, &
                                   ryd2mev * wf(imode, iqq), 2.0d0 * ryd2mev * REAL(gamma_all(imode, iqq, 1, itempphen))
          ENDDO
        ENDDO
        CLOSE(linewidth_phself)
      ENDDO ! itempphen
    ENDIF ! mpime
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_ph_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE input,         ONLY : nstemp, lsda
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    IF (meta_ionode) THEN
      !
      fnm = 'sigma_restart'
      IF (TRIM(lsda) == 'down') fnm = 'down.sigma_restart'
      lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmar_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmai_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = zi_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufilsigma_all, TRIM(fnm), lsigma_all, exst)
      CALL davcio(aux, lsigma_all, iufilsigma_all, 1, +1)
      CLOSE(iufilsigma_all)
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN
      sigmar_all(:, 1:lower_bnd - 1, :) = zero
      sigmai_all(:, 1:lower_bnd - 1, :) = zero
      zi_all(:, 1:lower_bnd - 1, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
      zi_all(:, upper_bnd + 1:nktotf, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_read(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY :  zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : mpime, world_comm
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE input,         ONLY : nstemp, lsda
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !! File name
    CHARACTER(LEN = 256) :: fnm
    !! Buffer variables for file name
    !
    IF (meta_ionode) THEN
      !
      ! First inquire if the file exists
      fnm = TRIM(prefix)
      IF (TRIM(lsda) == 'down') fnm = TRIM(prefix) // '.down'
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.sigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.sigma_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        fnm = 'sigma_restart'
        IF (TRIM(lsda) == 'down') fnm = 'down.sigma_restart'
        lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
        CALL diropn(iufilsigma_all, TRIM(fnm), lsigma_all, exst)
        CALL davcio(aux, lsigma_all, iufilsigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('selfen_el_read', &
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmar_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmai_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilsigma_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, meta_ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, meta_ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, meta_ionode_id, world_comm)
      CALL mp_bcast(zi_all, meta_ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
        zi_all(:, upper_bnd + 1:nktotf, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_write_wfpt(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy for WFPT. The only difference is that we write sigmar_dw_all.
    !! TODO (JML) : merge with selfen_el_write
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE input,         ONLY : nstemp, lsda
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmar_dw_all(nbndfst, nktotf, nstemp)
    !! Debye-Waller electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    IF (mpime == ionode_id) THEN
      !
      lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmar_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmai_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = zi_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmar_dw_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      fnm = 'sigma_restart'
      IF (TRIM(lsda) == 'down')  fnm = 'down.sigma_restart'
      CALL diropn(iufilsigma_all, TRIM(fnm), lsigma_all, exst)
      CALL davcio(aux, lsigma_all, iufilsigma_all, 1, +1)
      CLOSE(iufilsigma_all)
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN
      sigmar_all(:, 1:lower_bnd - 1, :) = zero
      sigmai_all(:, 1:lower_bnd - 1, :) = zero
      zi_all(:, 1:lower_bnd - 1, :) = zero
      sigmar_dw_all(:, 1:lower_bnd - 1, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
      zi_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmar_dw_all(:, upper_bnd + 1:nktotf, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_write_wfpt
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_read_wfpt(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading for WFPT.  The only difference is that we read sigmar_dw_all.
    !! TODO (JML) : merge with selfen_el_read
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY :  zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : mpime, world_comm
    USE io_global,     ONLY : ionode_id
    USE input,         ONLY : nstemp, lsda
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmar_dw_all(nbndfst, nktotf, nstemp)
    !! Debyw-Waller electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !! File name
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
      fnm = TRIM(prefix)
      IF (TRIM(lsda) == 'down') fnm = TRIM(prefix) // '.down'
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.sigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.sigma_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
        fnm = 'sigma_restart'
        IF (TRIM(lsda) == 'down') fnm = 'down.sigma_restart'
        CALL diropn(iufilsigma_all, TRIM(fnm), lsigma_all, exst)
        CALL davcio(aux, lsigma_all, iufilsigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('selfen_el_read', &
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmar_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmai_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmar_dw_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilsigma_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, ionode_id, world_comm)
      CALL mp_bcast(zi_all, ionode_id, world_comm)
      CALL mp_bcast(sigmar_dw_all, ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
        sigmar_dw_all(:, 1:lower_bnd - 1, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
        zi_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmar_dw_all(:, upper_bnd + 1:nktotf, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_read_wfpt
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,     ONLY : DP
    USE global_var,ONLY : lower_bnd, upper_bnd, nbndfst
    USE input,     ONLY : nstemp, wmin_specfun, wmax_specfun, nw_specfun, lsda
    USE io_var,    ONLY : iufilesigma_all
    USE io_files,  ONLY : diropn
    USE ep_constants,      ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: esigmar_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: esigmai_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperature
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2)
    !! Vector to store the array
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! energy range and spacing for spectral function
      !
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1.d0)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      !
      lesigma_all = 2 * nbndfst * nktotf * nw_specfun * nstemp + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmar_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmai_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      fnm = 'esigma_restart'
      IF (TRIM(lsda) == 'down') fnm = 'down.esigma_restart'
      CALL diropn(iufilesigma_all, TRIM(fnm), lesigma_all, exst)
      CALL davcio(aux, lesigma_all, iufilesigma_all, 1, +1)
      CLOSE(iufilesigma_all)
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN
      esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
      esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      esigmar_all(:, upper_bnd + 1:nktotf, :, :) = zero
      esigmai_all(:, upper_bnd + 1:nktotf, :, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_write_scgd0(totq, nktotf, esigmar_all, esigmai_all, nelec_w)
    !----------------------------------------------------------------------------
    !!
    !! This subroutine is used for the scGD0 calculation. Here we write down the 
    !! self-energy components in the binary files after each iteration, as well as 
    !! the iteration number, the number of bare electrons and the fermi level.
    !! Number of bare electrons is also computed here before writing it down.  
    !!
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : lower_bnd, upper_bnd, nbndfst, iter_scgd0, nkqf,&
                                 wkf, etf, nkqtotf, efnew, ibndmin, gtemp, mu_t
    USE input,            ONLY : nstemp, wmin_specfun, wmax_specfun, nw_specfun, &
                                 nbndsub, degaussw, efermi_read, fermi_energy, &
                                 fsthick, ahc_win_min, ahc_win_max, lwfpt
    USE io_var,           ONLY : iufilesigmasc_all
    USE io_files,         ONLY : diropn
    USE ep_constants,	  ONLY : zero, ryd2mev, pi, ryd2ev
    USE mp,               ONLY : mp_barrier, mp_sum
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : ionode_id
    USE mp_global,        ONLY : my_pool_id, inter_pool_comm
    USE parallelism,      ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: nelec_w
    !! number of electrons inside the active window
    REAL(KIND = DP), INTENT(inout) :: esigmar_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: esigmai_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperature
    INTEGER :: ierr
    !! Integer to check (de)allocation error
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: ef0
    !! Fermi level
    REAL(KIND = DP) :: ekk
    !! electron energy w.r.t. the fermi level
    REAL(KIND = DP), ALLOCATABLE :: aux(:)
    !! Vector to store the array
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! Collect k point weights from all pools in parallel case
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    REAL(KIND = DP), EXTERNAL :: wgauss, w0gauss
    !! Fermi-Dirac distribution function (when -99)
    !
    ALLOCATE(aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2 + nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_write_scgd0', 'Error allocating aux', 1)
    !
    ! energy range and spacing for spectral function
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1.d0)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    IF (efermi_read) THEN
      ef0 = fermi_energy
    ELSE
      ef0 = efnew
    ENDIF
    ! Here we compute the number of bare  electrons - at the first temperature only
    !
    IF (iter_scgd0 == 0) THEN
      nelec_w = 0.0
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_write_scgd0', 'Error allocating etf_all', 1)
      ALLOCATE(wkf_all(nkqtotf), STAT = ierr)
      wkf_all(:) = zero
      etf_all(:, :) = zero
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL poolgather2(1, nkqtotf, nkqf, wkf, wkf_all)
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst 
          IF (ABS(etf_all(ibndmin - 1 + ibnd, ik * 2 - 1) - ef0) > fsthick) CYCLE
          IF (lwfpt) THEN
            !
            ! Skip active states outside the ahc window
            IF (etf_all(ibnd -1+ibndmin, ik*2-1) < ahc_win_min .OR. etf_all(ibndmin-1+ibnd, ik*2-1) > ahc_win_max) CYCLE
          ENDIF
          ekk = etf_all(ibndmin - 1 + ibnd, ik * 2 - 1) - mu_t(1)
          nelec_w = nelec_w + wgauss(-ekk/gtemp(1), -99) * wkf_all(2*ik-1)  
        ENDDO
      ENDDO
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_write_scgd0', 'Error deallocating etf_all', 1)
      DEALLOCATE(wkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_write_scgd0', 'Error deallocating wkf_all', 1)
    ENDIF
    IF (my_pool_id == ionode_id) THEN
      !
      lesigma_all = 2 * nbndfst * nktotf * nw_specfun * nstemp + 2 + nstemp
      ! First element is the iteration
      aux(1) = INT(iter_scgd0) + 1
      ! Second element is the total number electrons
      aux(2) = REAL(nelec_w, KIND = DP)
      ! third element is the fermi level
      i = 2
      DO itemp = 1, nstemp
        i = i + 1
        aux(i) = REAL(mu_t(itemp), KIND = DP)
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmar_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmai_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      fnm = 'esigmasc_restart'
      CALL diropn(iufilesigmasc_all, TRIM(fnm), lesigma_all, exst)
      CALL davcio(aux, lesigma_all, iufilesigmasc_all, 1, +1)
      CLOSE(iufilesigmasc_all)
    ENDIF
    DEALLOCATE(aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_write_scgd0', 'Error deallocating aux', 1)    
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_write_scgd0
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_read(iqq, totq, nktotf, esigmar_all, esigmai_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE input,         ONLY : nstemp, nw_specfun, lsda
    USE io_var,        ONLY : iufilesigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: esigmar_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: esigmai_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !! File name
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! First inquire if the file exists
      fnm = TRIM(prefix)
      IF (TRIM(lsda) == 'down') fnm = TRIM(prefix) // '.down'
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigma_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        lesigma_all = 2 * nbndfst * nktotf * nw_specfun * nstemp + 2
        fnm = 'esigma_restart'
        IF (TRIM(lsda) == 'down') fnm = 'down.esigma_restart'
        CALL diropn(iufilesigma_all, TRIM(fnm), lesigma_all, exst)
        CALL davcio(aux, lesigma_all, iufilesigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('electron_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO iw = 1, nw_specfun
                i = i + 1
                esigmar_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO iw = 1, nw_specfun
                i = i + 1
                esigmai_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilesigma_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(esigmar_all, ionode_id, world_comm)
      CALL mp_bcast(esigmai_all, ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        esigmar_all(:, upper_bnd + 1:nktotf, :, :) = zero
        esigmai_all(:, upper_bnd + 1:nktotf, :, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_read_scgd0_check(iter_rest)
    !----------------------------------------------------------------------------
    !!
    !! This subroutine is used for the scGD0 calculation with a restart option. 
    !! It is used just to quickly check whether iteration = 0 is done or not.
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iufilesigmasc_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : world_comm, mpime
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : my_pool_id, inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: iter_rest
    !! Current iteration
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    !
    REAL(KIND = DP) :: first_val
    !! first value in the file
    CHARACTER(LEN = 256) :: name1
    !! File name
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    exst =.FALSE.
    iter_rest = 0
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! First inquire if the file exists
      fnm = TRIM(prefix)
      !
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigmasc_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigmasc_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
	fnm = 'esigmasc_restart'
        !
        CALL diropn(iufilesigmasc_all, TRIM(fnm), 1, exst)
        CALL davcio(first_val, 1, iufilesigmasc_all, 1, -1)
        iter_rest = INT(first_val)
        CLOSE(iufilesigmasc_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    CALL mp_bcast(iter_rest, ionode_id, world_comm)
    !
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_read_scgd0_check
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_read_scgd0(nktotf, esigmar_all, esigmai_all, iter_rest, nelec_w, ef)
    !----------------------------------------------------------------------------
    !!
    !! This subroutine is used for scGD0 calculations. Here we read the self-energy components
    !! from the previous run. THe quantities that are read consist of the real and imaginary
    !! energy dependendetn Fan-Migdal terms computed on the grid of nw_specfun points and, if
    !! WFPT is used, also the static DW and FM terms from the active and rest space.
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, gtemp,        &
                              sigmar_dw_all, sigma_ahc_uf, sigma_ahc_hdw
    USE input,         ONLY : nstemp, nw_specfun, lsda, lwfpt
    USE io_var,        ONLY : iufilesigmasc_all, iuelself_wfpt
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY : zero, ryd2ev, kelvin2eV, ryd2mev
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : world_comm, mpime
    USE io_global,     ONLY : ionode_id
    USE supercond_common, ONLY : nbndfs
    USE mp_global,     ONLY : my_pool_id, inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: iter_rest
    !! Current iteration
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: esigmar_all(nbndfs, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: esigmai_all(nbndfs, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: nelec_w
    !! Number of electrons inside the frequency window
    REAL(KIND = DP), INTENT(out) :: ef(nstemp)
    !! Fermi energy from the previous run 
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: ierr
    !! Error status
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ios
    !! integer to check if the file is opened correctly
    !
    REAL(KIND = DP), ALLOCATABLE :: aux(:)
    !! Vector to store the array
    REAL(KIND = DP) :: ik_
    !! This and the following 7 variables are needed to real the elself_wfpt_sup file.
    !! ik_ stands for the k point.
    REAL(KIND = DP) :: ibnd_
    !! band from the elself_wfpt_sup file
    REAL(KIND = DP) :: eks_
    !! energy from the elself_wfpt_sup file
    REAL(KIND = DP) :: resig
    !! active FM contribution to the real part of the self-energy from the elself_wfpt_sup file
    REAL(KIND = DP) :: dw
    !! active DW contribution to the real part of the self-energy from the elself_wfpt_sup file
    REAL(KIND = DP) :: uf
    !! rest FM contribution to the real part of the self-energy from the elself_wfpt_sup file
    REAL(KIND = DP) :: hdw
    !! rest DW contribution to the real part of the self-energy from the elself_wfpt_sup file
    REAL(KIND = DP) :: imsig
    !! imaginary part of the self-energy from the elself_wfpt_sup file
    !
    CHARACTER(LEN = 256) :: name1
    !! File name
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    CHARACTER(LEN = 20) :: tp
    !! string for temperature
    CHARACTER(LEN = 256) :: line
    !! lines to skip inside the elself_wfpt_sup file
    !
    exst =.FALSE.
    nelec_w = 0.0
    ALLOCATE(aux(2 * nbndfs * nktotf * nw_specfun * nstemp + 2 + nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_read_scgd0', 'Error allocating aux', 1)
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! First inquire if the file exists
      fnm = TRIM(prefix)
      !
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigmasc_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(fnm) // '.esigmasc_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        fnm = 'esigmasc_restart'
	lesigma_all = 2 * nbndfs * nktotf * nw_specfun * nstemp + 2 + nstemp
        CALL diropn(iufilesigmasc_all, TRIM(fnm), lesigma_all, exst)
        CALL davcio(aux, lesigma_all, iufilesigmasc_all, 1, -1)
        !
	! First element is the iteration number
        iter_rest = INT(aux(1))
        nelec_w = REAL(aux(2))
        !
	i = 2
	DO itemp = 1, nstemp
          i = i + 1
          ef(itemp) = REAL(aux(i))
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfs
              DO iw = 1, nw_specfun
                i = i + 1
                esigmar_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
	DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfs
              DO iw = 1, nw_specfun
                i = i + 1
                esigmai_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
	CLOSE(iufilesigmasc_all)
      ELSE
        esigmai_all(:, :, :, :) = zero
        esigmar_all(:, :, :, :) = zero
      ENDIF
      IF (lwfpt) THEN
        DO itemp = 1, nstemp
          ! Read AHC decomposition from elself_wfpt_sup.* files
          !
          fnm = ''
          IF (TRIM(lsda) == 'down') fnm = '.down'
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          name1 = 'elself_wfpt_sup.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
          !
          INQUIRE(FILE = name1, EXIST = exst)
          !
          IF (exst) THEN
            !
            OPEN(unit=iuelself_wfpt, FILE = name1, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
            IF (ios /= 0) CALL errore('spectral_read_scgd0', 'opening file ' // name1, ABS(ios))
            !
            ! Skip header (2 lines)
            READ(iuelself_wfpt, '(A)', iostat=ios) line
            READ(iuelself_wfpt, '(A)', iostat=ios) line
            !
            DO ibnd = 1, nbndfs
              DO ik = 1, nktotf
                !
                ! the elself_wfpt_sup file consists of 8 columns.
                ! we are interested in reading the active DW, rest FM, and rest DW terms.
                ! these are saved in columns 5, 6, 7
                READ(iuelself_wfpt, *, iostat=ios) ik_, ibnd_, eks_, resig, dw, uf, hdw, imsig
                !
                IF (ios /= 0) CALL errore('spectral_read_scgd0', 'Read error in ' // name1, ABS(ios))
                !
                sigmar_dw_all(ibnd, ik, itemp) = dw / ryd2mev
                sigma_ahc_uf(ibnd, ik, itemp)  = uf / ryd2mev
                sigma_ahc_hdw(ibnd, ik, itemp) = hdw / ryd2mev
                !
              ENDDO
              !
              ! skip blank line between bands
              READ(iuelself_wfpt, '(A)', iostat=ios) line
              !
            ENDDO
            CLOSE(iuelself_wfpt)
            !
          ELSE
            sigmar_dw_all(:, :, itemp) = zero
            sigma_ahc_uf(:, :, itemp)  = zero
            sigma_ahc_hdw(:, :, itemp) = zero
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    CALL mp_bcast(ef, ionode_id, world_comm)
    CALL mp_bcast(iter_rest, ionode_id, world_comm)
    CALL mp_bcast(esigmar_all, ionode_id, world_comm)
    CALL mp_bcast(esigmai_all, ionode_id, world_comm)
    CALL mp_bcast(nelec_w, ionode_id, world_comm)
    IF (lwfpt) THEN
      CALL mp_bcast(sigmar_dw_all, ionode_id, world_comm)
      CALL mp_bcast(sigma_ahc_uf, ionode_id, world_comm)
      CALL mp_bcast(sigma_ahc_hdw, ionode_id, world_comm)
    ENDIF
    !
    DEALLOCATE(aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_read_scgd0', 'Error deallocating aux', 1)
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_read_scgd0
    !----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_selfen
  !------------------------------------------------------------------------------
