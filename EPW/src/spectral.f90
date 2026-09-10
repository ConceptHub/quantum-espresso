  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE spectral
  !----------------------------------------------------------------------
  !!
  !! This module contains the various spectral function routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_el_print()
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation.
    !!
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode, ionode_id
    USE io_var,        ONLY : iospectral_sup, iospectral, iunkf
    USE input,         ONLY : nbndsub, wmin_specfun, wmax_specfun, nw_specfun, &
                              efermi_read, fermi_energy, nstemp, lsda,         &
                              specfun_el_scgd0, lwfpt, ncarrier, system_2d,    &
                              carrier, opt_cond, specfun_el, fsthick
    USE global_var,    ONLY : gtemp, etf, ibndmin, nkqf, nktotf, efnew, &
                              xkf, nkqtotf, esigmar_all, esigmai_all, a_all,   &
                              nbndfst, iter_scgd0, lower_bnd, upper_bnd, nkf,  &
                              a_all_ibnd, mu_t, wkf
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm, my_pool_id
    USE parallelism,   ONLY : poolgather2
    USE supercond_common, ONLY : xkfs_all, ekfs, ixkff, ekfs, nkfs, nbndfs, nkfs_all, ekfs_all
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    CHARACTER(LEN = 10) :: itr
    !! string for the number of scGD0 iterations
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ie
    !! Counter on the frequency
    INTEGER :: iw_plus
    !! Index of w + e on the ww grid
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! Spectral function
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: ww_plus   
    !! Value of w + e
    REAL(KIND = DP) :: temp_a_all
    !! temporary varibale fo rthe spectral function
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! Collect k-point weight from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    ! The k points are distributed among pools: here we collect them
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error allocating etf_all', 1)
    ALLOCATE(wkf_all(nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error allocating wkf_all', 1)
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    wkf_all(:) = zero
    !
    IF (opt_cond) THEN
      a_all_ibnd(:, :, :, :) = zero
    ENDIF
    IF (specfun_el_scgd0 .AND. iter_scgd0 > 0) THEN
      a_all(:, :, :) = zero
      DO ik = 1, nkfs_all
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        xkf_all(:, ikk) = xkfs_all(:, ik)
        DO ibnd = 1, nbndfs
          etf_all(ibndmin - 1 + ibnd, ikk) = ekfs_all(ibnd, ik)
        ENDDO
      ENDDO
      !
    ELSE
      CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL poolgather2(1, nkqtotf, nkqf, wkf, wkf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(esigmar_all, inter_image_comm)
      CALL mp_sum(esigmai_all, inter_image_comm)
    ENDIF
    !
    DO itemp = 1, nstemp ! second temperature loop to write data
      !
      ! Fermi level
      ef0 = mu_t(itemp)
      !
      inv_eptemp = one / gtemp(itemp)
      !
      ! Output electron spectral function here after looping over all q-points
      ! (with their contributions summed in a etc.)
      !
      WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      ! construct the trace of the spectral function (assume diagonal selfenergy
      ! and constant matrix elements for dipole transitions)
      !
      IF (ionode) THEN
        fnm = ''
        IF (TRIM(lsda) == 'down') fnm = '.down'
        WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
        IF (specfun_el_scgd0 .AND. iter_scgd0 /= 0) THEN
           filespec = 'specfun.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0' 
           filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0' 
            OPEN(UNIT = iospectral, FILE = TRIM(filespec) )
            OPEN(UNIT = iospectral_sup, FILE = TRIM(filespecsup) )
        ELSE
          filespec = 'specfun.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
          filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
          OPEN(UNIT = iospectral, FILE = TRIM(filespec) )
          OPEN(UNIT = iospectral_sup, FILE = TRIM(filespecsup) )
        ENDIF
        WRITE(iospectral, '(/2x, a/)') '#Electronic spectral function (meV)'
        WRITE(iospectral_sup, '(/2x, a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
        WRITE(iospectral, '(/2x, a/)') '#K-point     Energy[eV]     A(k,w)[meV^-1]'
        WRITE(iospectral_sup, '(/2x, a/)') '#K-point    Band   e_nk[eV]   w[eV]      Real Sigma[meV]  Im Sigma[meV]'
        !
      ENDIF
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        IF (.NOT. (specfun_el_scgd0 .OR. opt_cond)) THEN
          WRITE(stdout, '(/5x, "ik = ", i5, " coord.: ", 3f12.7, " Temp. : ", f8.3)') ik, xkf_all(:, ikk), &
                                                                                    gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
        ENDIF
        !
        DO iw = 1, nw_specfun
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k
            !
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
               ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
            IF (opt_cond) THEN
              a_all_ibnd(iw, ik, ibnd, itemp) = ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
               ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
            ENDIF
            !
          ENDDO
          !
          IF (.NOT. (specfun_el_scgd0 .OR. opt_cond)) THEN
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), a_all(iw, ik, itemp) / ryd2mev
          ENDIF
          !
        ENDDO
        !
        IF (.NOT. (specfun_el_scgd0 .OR. opt_cond))  WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
        IF (ik == 1 .AND. (specfun_el_scgd0 .OR. opt_cond)) THEN
          WRITE(stdout, '(/5x, a)') REPEAT('-', 67)
          WRITE(stdout, '(5x, "K-point        Energy (eV)         Spectral function (1/meV)")')
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          DO iw = 1, nw_specfun
            temp_a_all = a_all(iw, ik, itemp) / ryd2mev
            IF (ABS(temp_a_all) < 1.0d-8) temp_a_all = 0.0d0
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), temp_a_all
          ENDDO
          WRITE(stdout, '(/5x, a)')
          WRITE(stdout, '(5x, "Additional k-points are provided in the specfun.elself file")')
          WRITE(stdout, '(5x, a/)')
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
        ENDIF
        !
      ENDDO ! k pts
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik -1
        ! The spectral function should integrate to 1 for each k-point
        specfun_sum = 0.0
        !
        DO iw = 1, nw_specfun
          !
          fermi(iw) = wgauss(-ww(iw) * inv_eptemp, -99)
          !
          specfun_sum = specfun_sum + a_all(iw, ik, itemp) * fermi(iw) * dw
          !
          IF (ionode) WRITE(iospectral, '(2x, i7, 2x, f10.5, 2x, E12.5)') ik, ryd2ev * ww(iw), &
                                                                            a_all(iw, ik, itemp) / ryd2mev
          !
        ENDDO
        !
        IF (ionode) WRITE(iospectral, '(a)') ' '
        IF (ionode) WRITE(iospectral, '(2x, a, 2x, E12.5)') '# Integrated spectral function ', specfun_sum
        !
      ENDDO
      !
      IF (ionode) CLOSE(iospectral)
      !
      DO ibnd = 1, nbndfst
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          !  the energy of the electron at k
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
          !
          DO iw = 1, nw_specfun
            !
            IF (.NOT. (specfun_el_scgd0 .OR. opt_cond)) THEN
              WRITE(stdout, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
            ENDIF        
            !
            IF (ionode) &
              WRITE(iospectral_sup, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                  ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
            !
          ENDDO
          !
        ENDDO
        !
        IF (.NOT. (specfun_el_scgd0 .OR. opt_cond)) THEN
          WRITE(stdout, *) ' '
        ENDIF
        !
      ENDDO
      !
      IF (ionode) CLOSE(iospectral_sup)
      !
      IF (opt_cond) CALL optical_conductivity()
      !
    ENDDO ! itemp
    DEALLOCATE(wkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error deallocating wkf_all', 1)
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error deallocating etf_all', 1)
    !
    101 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', ES13.5, ' meV^-1')
    102 FORMAT(2i9, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_el_print
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_el_interpolate()
    !-----------------------------------------------------------------------
    !!
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode, ionode_id
    USE io_var,        ONLY : iospectral_sup, iospectral, iunkf
    USE input,         ONLY : nbndsub, wmin_specfun, wmax_specfun, nw_specfun, &
                              efermi_read, fermi_energy, nstemp, lsda,         &
                              specfun_el_scgd0, lwfpt, ncarrier, system_2d,    &
                              carrier, opt_cond, filkf
    USE global_var,    ONLY : gtemp, etf, ibndmin, nkqf, nktotf, efnew, &
                              xkf, nkqtotf, esigmar_all, esigmai_all, a_all,   &
                              nbndfst, iter_scgd0, lower_bnd, upper_bnd, nkf,  &
                              a_all_ibnd, mu_t
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm, my_pool_id
    USE parallelism,   ONLY : poolgather2
    USE supercond_common, ONLY : xkfs_all, ekfs, ekfs, nkfs, nbndfs, ekfs_all, nkfs_all
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    CHARACTER(LEN = 10) :: itr
    !! string for the number of scGD0 iterations
    CHARACTER(LEN = 10) :: coordinate_type
    !! filkf coordinate type (crystal or cartesian)
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ie
    !! Counter on the frequency
    INTEGER :: iw_plus
    !! Index of w + e on the ww grid
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    INTEGER :: nk_path
    !! Counter on xkf_path
    !
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! Spectral function
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: ww_plus   
    !! Value of w + e
    REAL(KIND = DP) :: a_path(nw_specfun)
    !! scgd0 spectral function to be interpolated to a desired k pt.
    REAL(KIND = DP) :: esigmaisc_interp(nw_specfun)
    !! imaginary part of the scgd0 self-energy to be interpolated to a desired k pt.
    REAL(KIND = DP) :: esigmarsc_interp(nw_specfun)
    !! real part of the scgd0 self-energy to be interpolated to a desired k pt.
    REAL(KIND = DP) :: temp_a_all
    !! temporary varibale fo rthe spectral function
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: xkf_path(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    ! The k points are distributed among pools: here we collect them
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error allocating etf_all', 1)
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    !
    IF (iter_scgd0 > 0) THEN
      a_all(:, :, :) = zero
      IF (opt_cond) THEN
        a_all_ibnd(:, :, :, :) = zero
      ENDIF
      DO ik = 1, nkfs_all
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        xkf_all(:, ikk) = xkfs_all(:, ik)
        DO ibnd = 1, nbndfs
          etf_all(ibndmin - 1 + ibnd, ikk) = ekfs_all(ibnd, ik)
        ENDDO
      ENDDO
      !
    ELSE
      CALL poolgather2(3,	nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(esigmar_all, inter_image_comm)
      CALL mp_sum(esigmai_all, inter_image_comm)
    ENDIF
    !
    IF (my_pool_id == ionode_id) THEN
      !
      WRITE(stdout, '(5x, "A path file is given and the spectral function will be interpolated.")')
      OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'opening file ' // filkf, ABS(ierr))
      READ(iunkf, *) nk_path, coordinate_type
      ALLOCATE(xkf_path(3, nk_path), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error allocating xkf_path', 1)
        DO ik = 1, nk_path
          READ(iunkf, *) xkf_path(:, ik )
        ENDDO
      CLOSE(iunkf)
    ENDIF
    CALL mp_bcast(nk_path,  ionode_id, inter_pool_comm)
    IF (my_pool_id /= ionode_id) THEN
      ALLOCATE(xkf_path(3, nk_path), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error allocating xkf_path', 1)
    ENDIF
    CALL mp_bcast(xkf_path, ionode_id, inter_pool_comm)
    !
    DO itemp = 1, nstemp ! temperature loop to write data
      !
      inv_eptemp = one / gtemp(itemp)
      !
      ! Output electron spectral function here after looping over all q-points
      ! (with their contributions summed in a etc.)
      !
      WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      ! construct the trace of the spectral function (assume diagonal selfenergy
      ! and constant matrix elements for dipole transitions)
      !
      IF (ionode) THEN
        fnm = ''
        IF (TRIM(lsda) == 'down') fnm = '.down'
        WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
        IF (iter_scgd0 == 0) THEN
          filespec = 'specfun.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) 
          filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
        ELSE
          filespec = 'specfun.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0'
          filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0'
        ENDIF
        OPEN(UNIT = iospectral, FILE = TRIM(filespec) )
        OPEN(UNIT = iospectral_sup, FILE = TRIM(filespecsup) )
        WRITE(iospectral, '(/2x, a/)') '#Electronic spectral function (meV)'
        WRITE(iospectral_sup, '(/2x, a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
        WRITE(iospectral, '(/2x, a/)') '#K-point     Energy[eV]     A(k,w)[meV^-1]'
        WRITE(iospectral_sup, '(/2x, a/)') '#K-point    Band   e_nk[eV]   w[eV]      Real Sigma[meV]  Im Sigma[meV]'
        !
      ENDIF
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        !
        DO iw = 1, nw_specfun
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k
            !
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - mu_t(itemp)
            !
            a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
               ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
            IF (opt_cond) THEN
              a_all_ibnd(iw, ik, ibnd, itemp) =  ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
               ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
            ENDIF
            !
          ENDDO
          !
        ENDDO
      ENDDO ! k pts
      !
      IF (opt_cond) CALL optical_conductivity()
      !
      DO ik = 1, nk_path
        CALL interpolate_path_scgd0(xkf_path(:, ik), xkf_all, a_all(:, :, itemp), a_path)
        !
        DO iw = 1, nw_specfun
          IF (ionode) THEN
            IF (ABS(a_path(iw) / ryd2mev) < 1.0d-99) a_path(iw) = 0.0d0
            WRITE(iospectral, '(2x, i7, 2x, f10.5, 2x, E12.5)') ik, ryd2ev * ww(iw), &
                                                                           a_path(iw) / ryd2mev
          ENDIF
        ENDDO  
        IF (ik ==1) THEN
          WRITE(stdout, '(/5x, a)') REPEAT('-', 67)
          WRITE(stdout, '(5x, "K-point        Energy (eV)         Spectral function (1/meV)")')
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          DO iw = 1, nw_specfun
            temp_a_all = a_path(iw) / ryd2mev
            IF (ABS(temp_a_all) < 1.0d-8) temp_a_all = 0.0d0
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), temp_a_all
          ENDDO
          WRITE(stdout, '(/5x, a)') 
          WRITE(stdout, '(5x, "Additional k-points are provided in the specfun.elself file")')
          WRITE(stdout, '(5x, a/)') 
        ENDIF
        IF (ionode) WRITE(iospectral, '(a)') ' '
      ENDDO
      !
      IF (ionode) CLOSE(iospectral)
      !
      DO ibnd = 1, nbndfst
        DO ik = 1, nk_path
          !
          CALL interpolate_path_scgd0_ekk(xkf_path(:, ik), xkf_all, etf_all(ibndmin -1 + ibnd, :), ekk)
          CALL interpolate_path_scgd0(xkf_path(:, ik), xkf_all, TRANSPOSE(esigmar_all(ibnd, : , :, itemp)), esigmarsc_interp)
          CALL interpolate_path_scgd0(xkf_path(:, ik), xkf_all, TRANSPOSE(esigmai_all(ibnd, : , :, itemp)), esigmaisc_interp)
          DO iw = 1, nw_specfun
            !
            IF (ionode) &
              WRITE(iospectral_sup, 102) ik, ibndmin - 1 + ibnd, ryd2ev * (ekk - mu_t(itemp)), ryd2ev * ww(iw), &
                  ryd2mev * esigmarsc_interp(iw), ryd2mev * esigmaisc_interp(iw)
            !
          ENDDO
          !
        ENDDO
          !
      ENDDO
      !
      IF (ionode) CLOSE(iospectral_sup)
      !
    ENDDO ! itemp
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error deallocating etf_all', 1)
    DEALLOCATE(xkf_path, STAT = ierr)
    IF (ierr /= 0) CALL errore('spectral_func_el_interpolate', 'Error deallocating xkf_path', 1)
    !
    101 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', ES13.5, ' meV^-1')
    102 FORMAT(2i9, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_el_interpolate
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_recompute_efermi(n_old, mu, wkf_, ekf_, nbnd, nkf_)
    !-----------------------------------------------------------------------
    ! This suborutine uses the bisection method to compute the Fermi level 
    ! from the spectral function. 
    !-----------------------------------------------------------------------
    !
    USE io_global,     ONLY : stdout
    USE input,            ONLY : wmax_specfun, wmin_specfun, nw_specfun, lwfpt, efermi_read,&
                                 nstemp, specfun_el, fsthick, ahc_win_max, ahc_win_min, fermi_energy
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : gtemp, esigmai_all, esigmar_all, ibndmin, &
                                 iter_scgd0, efnew, nbndfst, nktotf
    USE utilities,        ONLY : fermi_dirac
    USE supercond_common, ONLY : ef0, ixkff
    USE ep_constants,     ONLY : pi, ryd2mev
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(INOUT) :: mu(nstemp)
    ! given chemical potential from a previous run 
    REAL(KIND = DP), INTENT(IN)    :: n_old
    ! the number of electrons we want to achieve
    INTEGER, INTENT(IN)            :: nbnd
    !! nuumber of bands
    INTEGER, INTENT(IN)            :: nkf_
    !! number of k points
    REAL(KIND = DP), INTENT(IN)    :: wkf_(nkf_)
    !! k point weight array
    REAL(KIND = DP), INTENT(IN)    :: ekf_(nbnd, nkf_)
    !! energy array
    INTEGER :: itemp
    !! temperature index
    INTEGER :: iter
    !! bisection iteration index
    INTEGER :: ik
    !! k point counter
    INTEGER :: ikk
    !! k point counter
    INTEGER :: ib
    !! band index
    INTEGER :: ibnd
    !! band index
    INTEGER :: iw
    !! energy index 
    !
    REAL(KIND = DP) :: mu_lo
    !! lower boundary for bisection
    REAL(KIND = DP) :: mu_hi
    !! upper boundary for bisection
    REAL(KIND = DP) :: mu_mid
    !! middle value for bisection
    REAL(KIND = DP) :: n_lo
    !! lower boundary on the number of electrons
    REAL(KIND = DP) :: n_hi
    !! upper boundary on the number of electrons
    REAL(KIND = DP) :: n_mid
    !! middle value for the number
    REAL(KIND = DP) :: mu_old
    !! variable to save the mu value
    REAL(KIND = DP) :: a_dos(nw_specfun)
    !! DOS computed from the spectral function
    REAL(KIND = DP) :: alpha
    !! linear mixing factor
    REAL(KIND = DP)  :: dw
    !! energy increment
    REAL(KIND = DP)  :: ekk
    !! kohn-sham energy at n,k
    REAL(KIND = DP)  :: a_w
    !! spectral function
    REAL(KIND = DP)  :: mu_, n0
    !! chemical potential, temporary variable 
    REAL(KIND = DP)  :: ww(nw_specfun)
    !
    ! DOS-projected spectral function: A_dos(iw) = sum_{k,n} wk * A_{kn}(iw)
    ! Accumulated once per temperature; bisection then only integrates over iw.
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    alpha = 0.2
    !
    IF (specfun_el .OR. iter_scgd0 ==0) THEN
      alpha =1
      ! frmi level to center fsthick
     IF (efermi_read) THEN
        ef0 = fermi_energy
      ELSE
        ef0 = efnew
      ENDIF
    ENDIF
    DO itemp = 1, nstemp
      ! --- initial bracket
      a_dos(:) = 0.0
      !
      DO ik = 1, nktotf
        IF (specfun_el .OR. iter_scgd0 == 0) THEN
          ikk = 2 * ik - 1
        ELSE
          ikk = ik
        ENDIF
        !
        DO ib = 1, nbndfst
          IF (specfun_el .OR. iter_scgd0 == 0) THEN
            ibnd = ibndmin - 1 + ib
          ELSE
            ibnd = ib
          ENDIF
          !
          ! Only include states inside the fsthick energy window
          IF (ABS(ekf_(ibnd, ikk) - ef0) > fsthick) CYCLE
          !
          IF (lwfpt) THEN
            !
            ! Skip active states outside the ahc window
            IF (ekf_(ibnd, ikk) < ahc_win_min .OR. ekf_(ibnd, ikk) > ahc_win_max) CYCLE
            !
          ENDIF
          ekk = ekf_(ibnd, ikk) - mu(itemp)
          !
          DO iw = 1, nw_specfun
            a_w = ABS(esigmai_all(ib, ik, iw, itemp)) / (pi * ((ww(iw) - ekk - esigmar_all(ib, ik, iw, itemp))**2 &
                      + esigmai_all(ib, ik, iw, itemp)**2))
            a_dos(iw) = a_dos(iw) +  a_w  * wkf_(ikk)
          ENDDO
        ENDDO ! ibnd
      ENDDO   ! ik
      !
      mu_old = mu(itemp)
      mu_lo  =  -ABS(wmin_specfun) * 0.8
      mu_hi  =  wmax_specfun * 0.8
      !
      n_lo = compute_N(a_dos, ww, mu_lo, gtemp(itemp), nw_specfun, dw) - n_old
      n_hi = compute_N(a_dos, ww, mu_hi, gtemp(itemp), nw_specfun, dw) - n_old
      !
      n0 = compute_N(a_dos, ww, 0.d0, gtemp(itemp), nw_specfun, dw)
      !
      WRITE(stdout, '(5x,"Number of electrons computed:  N =",f15.8, " N_target=",f15.8)') &
           n0, n_old
      !
      IF (n_lo * n_hi > 0.0) THEN
        WRITE(stdout, '(5x,"WARNING: Failed to bracket Fermi level – ", &
                        "frequency window too small. Skipping update.")')
        CYCLE
      ENDIF
      !
      DO iter = 1, 500
        mu_mid = 0.5 * (mu_lo + mu_hi)
        n_mid  = compute_N(a_dos, ww, mu_mid, gtemp(itemp), nw_specfun, dw) - n_old
        IF (ABS(n_mid) < 1.0E-4) EXIT  
        IF (n_mid * n_lo < 0.0) THEN
          mu_hi = mu_mid
          n_hi  = n_mid
        ELSE
          mu_lo = mu_mid
          n_lo  = n_mid
        ENDIF
      ENDDO
      !
      mu_mid = mu_mid + mu_old
      mu_ = alpha *  mu_mid + (1 - alpha) * mu_old
      IF (ABS(mu_ - mu_old) > 1/ryd2mev) THEN
        mu(itemp) = mu_ 
      ELSE
        mu(itemp) = mu_old
      ENDIF 
    ENDDO
    !
    CONTAINS
    !
    PURE FUNCTION compute_N(a_dos, ww, mu, kT, nw, dw) RESULT(n)
      !! Integrate  N = sum_iw  A_dos(iw) * f(ww(iw) - mu, kT) * dw
      !! A_dos is the k/band-summed spectral DOS (precomputed outside bisection).
      !! This is an O(nw) operation with no k-point or band loops.
      USE kinds,     ONLY : DP
      USE utilities, ONLY : fermi_dirac
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nw
      !! number of energy points 
      INTEGER  :: iw
      !! energy counter
      REAL(KIND = DP), INTENT(IN) :: a_dos(nw)
      !! DOS computed from the spectral function
      REAL(KIND = DP), INTENT(IN) :: ww(nw)
      !! energy grid
      REAL(KIND = DP), INTENT(IN) :: mu
      !! chemical potential shift
      REAL(KIND = DP), INTENT(IN) :: kt
      !! temperature
      REAL(KIND = DP), INTENT(IN) :: dw
      !! energy increment
      REAL(KIND = DP) :: n
      !! number of electrons
      !
      n = 0.0
      DO iw = 1, nw
        n = n + a_dos(iw) * fermi_dirac(ww(iw) - mu, kT) * dw
      ENDDO
    !
    END FUNCTION compute_N
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_recompute_efermi
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE optical_conductivity()
    !-----------------------------------------------------------------------
    !!
    !! Compute the optical conductivity from the spectral function.
    !! can be used with specfun_el or specfun_el_scgd0.
    !!
    !-----------------------------------------------------------------------
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout, ionode
    USE io_var,            ONLY : iufilsigma
    USE input,             ONLY : nbndsub, wmin_specfun, wmax_specfun, nw_specfun, mp_mesh_k,  &
                                  efermi_read, fermi_energy, nstemp, lsda, lwfpt, ahc_win_min, &
                                  specfun_el_scgd0, lwfpt, ncarrier, system_2d, ahc_win_max,   &
                                  carrier, specfun_el, fsthick, nkf2, nkf1, nkf3
    USE noncollin_module,  ONLY : noncolin
    USE global_var,        ONLY : gtemp, etf, ibndmin, nkqf, nktotf, efnew, bztoibz, s_bztoibz,&
                                  xkf, nkqtotf, a_all_ibnd, nbndfst, iter_scgd0, esigmai_all,  &
                                  lower_bnd, upper_bnd, wkf, vmef, nkf, mu_t, a_all, esigmar_all
    USE ep_constants,      ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, hbar, Ang2m, &
                                  bohr2ang, ang2cm, hbarJ, czero
    USE constants,         ONLY : electron_si
    USE mp,                ONLY : mp_sum
    USE mp_global,         ONLY : inter_pool_comm, inter_image_comm
    USE parallelism,       ONLY : poolgather2, fkbounds
    USE supercond_common,  ONLY : xkfs, ekfs, wkfs, wkfs_all, nkfs, nbndfs, nkfs_all, ekfs_all, &
                                  ixkf, ef0
    USE cell_base,         ONLY : alat, at, omega, bg
    USE symm_base,         ONLY : s
    USE selfen,            ONLY : hilbert_transform
    USE utilities,         ONLY : fermi_dirac
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: fileoptcond
    !! File name for supporting information
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    CHARACTER(LEN = 10) :: itr
    !! string for the number of scGD0 iterations
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: iw0
    !! first positive frequency in the grid
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: ie
    !! Counter on the frequency
    INTEGER :: iw_plus
    !! Index of w + e on the ww grid
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ik_global
    !! global k  index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd      
    !! Counter on bands
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    INTEGER :: i
    !! Counter on x,y,z components
    INTEGER :: j
    !! Counter on x,y,z components
    INTEGER :: ik_vbm
    !! k index for vbm
    INTEGER :: ik_cbm
    !! k index for cbm
    INTEGER :: ibnd_vbm
    !! band index for vbm
    INTEGER :: ibnd_cbm
    !! band index for cbm
    !
    REAL(KIND = DP) :: evbm
    !! VBM energy
    REAL(KIND = DP) :: ecbm
    !! CBM energy
    REAL(KIND = DP) :: ecut
    !! truncation energy 
    REAL(KIND = DP) :: e_vbm
    !! renormalized  VBM energy
    REAL(KIND = DP) :: e_cbm
    !! renormalized  CBM energy
    REAL(KIND = DP) :: max_val
    !! value to search for band extrema
    REAL(KIND = DP) :: nelec_w
    !! number of bare electrons
    REAL(KIND = DP) :: ef
    !! variable t store mu_t(itemp)
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: conv_factor
    !! conversion of optical conductivity and optical mobility to correct units
    REAL(KIND = DP) :: conv_factor1
    !! conversion of optical conductivity and optical mobility to correct units
    REAL(KIND = DP) :: mob_factor
    !! dimension dependent factor for mobility 
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP) :: inv_cell
    !! inverse cell volume
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! Fermi-Dirac distribution function derivative (when -99)
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: ww_plus
    !! Value of w + e
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    COMPLEX(KIND = DP) :: vkk(3)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(KIND = DP) :: sigma_opt(3, 3, nw_specfun)
    !! optical conductivity tensor
    REAL(KIND = DP) :: sigma_opt_im(3, 3, nw_specfun)
    !! optical conductivity tensor
    REAL(KIND = DP) :: velocity_factor(3 ,3)
    !! prefactor for conductivity v_kk_i * v_kk_j * wkf
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! Collect k point weights from all pools in parallel case
    COMPLEX(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Rotation matrix
    REAL(KIND = DP) :: sb(3, 3)
    !! Rotation matrix
    REAL(KIND = DP) :: sr(3, 3)
    !! Rotation matrix
    REAL(KIND=DP) :: v_in_re(3)
    !! real part of the velocity vkk
    REAL(KIND=DP) :: v_in_im(3)
    !! imaginary part of the velocity vkk
    REAL(KIND=DP) ::  v_out_re(3)
    !! real part of the rotated velocity
    REAL(KIND=DP) :: v_out_im(3)    
    !! imaginary part of the rotated velocity
    !
    LOGICAL :: fermi_in_gap
    !! check to see whether the truncation method is applicable or not
    !
    IF (specfun_el_scgd0) CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
      conv_factor1 = electron_si**2/ (bohr2ang * Ang2m * hbarJ) * inv_cell
      mob_factor = (bohr2ang * ang2cm)**3
    ELSE
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
      conv_factor1 = electron_si**2/ hbarJ * inv_cell
      mob_factor = (bohr2ang * ang2cm)**2
    ENDIF
    conv_factor = conv_factor1
    !
    ! The k points are distributed among pools: here we collect them
    !
    ALLOCATE(wkf_all(nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('optical_conductivity', 'Error allocating wkf_all', 1)
    wkf_all(:) = zero
    sigma_opt(:, :, :) = zero
    sigma_opt_im(:, :, :) = zero
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('optical_conductivity', 'Error allocating etf_all', 1)
    etf_all(:, :) = zero
    !
    IF (specfun_el_scgd0 .AND. iter_scgd0 > 0) THEN
      DO ik = 1, nkfs_all
        !
	ikk = 2 * ik - 1
        wkf_all(ikk) = wkfs_all(ik)
        DO ibnd = 1, nbndfst
          etf_all(ibndmin - 1 + ibnd, ikk) = ekfs_all(ibnd, ik)
        ENDDO
      ENDDO
      !
    ELSE
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL poolgather2(1, nkqtotf, nkqf, wkf, wkf_all)
    ENDIF
    !
    WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
    WRITE(stdout, '(5x, "Computing optical conductivity in the bubble approximation.")')
    WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
    !
    DO itemp = 1, nstemp
      ef = mu_t(itemp)
      !
      ! in the one-shot calculation, we need to impose a truncation scheme for mobility 
      ! first, we find the cut off energy (see Supp. S2 B in  Phys. Rev. Lett. 134, 186401 (2025))
      ! for that, we need VBM, CBM
      !
      IF (specfun_el .OR. iter_scgd0 == 0) THEN
        IF (efermi_read) THEN
          ef0 = fermi_energy
        ELSE
          ef0 = efnew
        ENDIF
        !
      ENDIF
      IF (carrier .AND. (specfun_el .OR. iter_scgd0 == 0)) THEN
        !
        evbm = -1.0d+8
        ecbm = 1.0d+8
        !
	DO ik = 1, nktotf
          ikk = 2 * ik - 1
          DO ibnd = 1, nbndsub
            IF (etf_all(ibnd, ikk) <= ef) THEN
              IF (etf_all(ibnd, ikk) > evbm) THEN
                evbm = etf_all(ibnd, ikk)
                ik_vbm = ik
                ibnd_vbm = ibnd
              ENDIF
            ELSE
              IF (etf_all(ibnd, ikk) < ecbm) THEN
                ecbm = etf_all(ibnd, ikk)
                ik_cbm = ik
                ibnd_cbm = ibnd
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        !
        ! this variable checks if the fermi level is inside the gap.
        ! 10 mev is a tolerance in case the fermi level in fact crosses a band.
        !
        fermi_in_gap = (ef > evbm + 10/ ryd2mev) .AND. (ef + 10/ ryd2mev< ecbm)
        max_val = 1E-7
        IF (fermi_in_gap) THEN
          ! In case of light  p-doping:
          IF (ncarrier < 0.d0) THEN
            DO iw = 1, nw_specfun  
              IF (a_all_ibnd(iw, ik_vbm, ibnd_vbm, itemp) > max_val) THEN
                max_val = a_all_ibnd(iw, ik_vbm, ibnd_vbm, itemp)
                e_vbm = ww(iw)
              ENDIF
            ENDDO
            evbm = evbm - ef
            ecut = 0.5d0 * (evbm +  e_vbm) + 2 * gtemp(itemp) + SQRT( 4 * gtemp(itemp)**2 + (0.5d0 *evbm- 0.5d0 * e_vbm)**2)
            !
            WRITE(stdout, '(/5x, "Truncation scheme is applied with Ecut (p-doped) = ", f10.6, "[eV]")') ecut !* ryd2ev
            !
       	    ! In case of light  n-doping:
          ELSEIF (ncarrier > 0.d0) THEN
            DO iw = 1, nw_specfun 
              IF (a_all_ibnd(iw, ik_cbm, ibnd_cbm, itemp) > max_val) THEN
                max_val = a_all_ibnd(iw, ik_cbm, ibnd_cbm, itemp)
                e_cbm = ww(iw)
              ENDIF
            ENDDO
            ecbm = ecbm - ef
            ecut = 0.5d0 * (ecbm + e_cbm) - 2 * gtemp(itemp) - SQRT( 4 * gtemp(itemp)**2 + (0.5d0 *ecbm - 0.5d0 * e_cbm)**2)
            WRITE(stdout, '(/5x, "Truncation scheme is applied with Ecut (n-doped) = ", f10.6, "[eV]")') ecut * ryd2ev
          ENDIF
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              IF (ncarrier < 0.0d0) THEN
                DO iw =1, nw_specfun
                  ! p-doped → truncate ABOVE cutoff
                  ! truncate also esigmai_all so that 'spectral_recompute_efermi' doesn't count electrons above ecut
                  IF (ww(iw) > ecut) THEN
                    a_all_ibnd(iw, ik, ibnd, itemp) = 0.0
                    esigmai_all(ibnd, ik, iw, itemp) = 0.0
                  ENDIF
                ENDDO
              ELSEIF (ncarrier > 0.0d0) THEN
                ! n-doped → truncate BELOW cutoff
                ! truncate also esigmai_all so that 'spectral_recompute_efermi' doesn't count electrons below ecut
                DO iw = 1, nw_specfun
       	          IF (ww(iw) < ecut) THEN
                    a_all_ibnd(iw, ik, ibnd, itemp) = 0.0
                    esigmai_all(ibnd, ik, iw, itemp) = 0.0
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF (specfun_el .OR. iter_scgd0 == 0) THEN
        ! now, recompute the fermi level, in the case of one-shot G0D0 calculation
        ! in this case, we need to have ef0 - mu_t(itemp) shift in the fermi dirac factors. 
        ! first, we count the bare electrons inside the energy window
        nelec_w = 0.0
        !
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          DO ibnd = 1, nbndfst
            IF (ABS(etf(ibndmin - 1 + ibnd, ikk) - ef0) > fsthick) CYCLE
            IF (lwfpt) THEN
              !
              ! Skip active states outside the ahc window
             IF (etf(ibnd -1+ibndmin, ikk) < ahc_win_min .OR. etf(ibndmin-1+ibnd, ikk) > ahc_win_max) CYCLE
             !
             ENDIF
              ekk = etf(ibndmin - 1 + ibnd, ikk) - ef
              nelec_w = nelec_w + fermi_dirac(ekk, gtemp(itemp)) * wkf(ikk)
          ENDDO
        ENDDO
	CALL mp_sum(nelec_w, inter_pool_comm)
        CALL spectral_recompute_efermi(nelec_w, mu_t, wkf_all, etf_all, nbndsub, nkqtotf)
        WRITE(stdout, '(/5x, "Chemical potential is changed from ", f20.6, "[eV] to ",f20.6,"[eV] ")') &
          ef * ryd2ev, mu_t(itemp)*ryd2ev
      ENDIF
      IF (carrier .AND. ncarrier < -1E5) THEN
        carrier_density = 0.0
        !
        DO ik = 1, nktotf
          ikk = 2 * ik - 1
          DO ibnd = 1, nbndfst
            IF (ABS(etf_all(ibndmin - 1 + ibnd, ikk) - ef0) > fsthick) CYCLE
            IF (etf_all(ibndmin - 1 + ibnd, ikk) < mu_t(itemp)) THEN
              IF (lwfpt) THEN
                !
                ! Skip active states outside the ahc window
                IF (etf_all(ibnd - 1+ ibndmin, ikk) < ahc_win_min .OR. etf_all(ibndmin-1+ibnd, ikk) > ahc_win_max) CYCLE
                !
              ENDIF
              DO iw = 1, nw_specfun
                ! The wkf(ikk) already include a factor 2
                ! fermi_dirac has a shifted argumnet, because in G0D0 one shot case, ef is renormalized to mu_t(itemp)
                ! in the scGD0 case, this difference is 0.
                carrier_density = carrier_density + (1.0d0 - fermi_dirac(ww(iw) + ef - mu_t(itemp), gtemp(itemp))) * &
                a_all_ibnd(iw, ik, ibnd, itemp) * dw  * wkf_all(ikk)
              ENDDO 
            ENDIF
          ENDDO
        ENDDO
	conv_factor = conv_factor1 / (carrier_density * electron_si * inv_cell) * mob_factor
        IF (system_2d == 'no') THEN
          WRITE(stdout, '(/5x, "Carrier density recomputed from the spectral function is ", f30.4, "[*10^8 cm^-3]")') &
        carrier_density / mob_factor * inv_cell * 1.0d-8
        ELSE
          WRITE(stdout, '(/5x, "Carrier density recomputed from the spectral function is ", f30.4, "[*10^8 cm^-2]")')&
        carrier_density / mob_factor * inv_cell * 1.0d-8
        ENDIF
      ELSEIF (carrier .AND. ncarrier > 1E5) THEN 
        carrier_density = 0.0
        !
        DO ik = 1, nktotf
          ikk = 2 * ik - 1
          DO ibnd = 1, nbndfst
            IF (ABS(etf_all(ibndmin - 1 + ibnd, ikk) - ef0) > fsthick) CYCLE
            ! This selects only conduction bands for electron conduction
            IF (etf_all(ibndmin - 1 + ibnd, ikk) > mu_t(itemp)) THEN
              IF (lwfpt) THEN
                !
                ! Skip active states outside the ahc window
                IF (etf_all(ibnd -1+ibndmin, ikk) < ahc_win_min .OR. etf_all(ibndmin-1+ibnd, ikk) > ahc_win_max) CYCLE
                !
              ENDIF
              DO iw = 1, nw_specfun
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + fermi_dirac(ww(iw) + ef - mu_t(itemp), gtemp(itemp)) * &
                a_all_ibnd(iw, ik, ibnd, itemp) * dw * wkf_all(ikk)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        conv_factor = conv_factor1 / (carrier_density * electron_si * inv_cell) * mob_factor
        IF (system_2d == 'no') THEN
          WRITE(stdout, '(/5x, "Carrier density recomputed from the spectral function is ", f30.4, "[*10^8 cm^-3]")') &
        carrier_density / mob_factor * inv_cell * 1.0d-8
       	ELSE
          WRITE(stdout, '(/5x, "Carrier density recomputed from the spectral function is ", f30.4, "[*10^8 cm^-2]")') &
        carrier_density / mob_factor * inv_cell * 1.0d-8
       	ENDIF
      ENDIF
      IF (ionode) THEN
        fnm = ''
        IF (TRIM(lsda) == 'down') fnm = '.down'
        WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
        IF (specfun_el_scgd0 .AND. iter_scgd0 > 0) THEN
          IF (carrier) THEN
             fileoptcond = 'AC_mobility.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0'
          ELSE
             fileoptcond = 'AC_conductivity.' // trim(adjustl(tp)) // 'K' // TRIM(fnm) // '_scGD0'
          ENDIF
        ELSE
          IF (carrier) THEN
            fileoptcond = 'AC_mobility.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
       	  ELSE
            fileoptcond = 'AC_conductivity.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
       	  ENDIF
        ENDIF
        OPEN(UNIT = iufilsigma, FILE = TRIM(fileoptcond) )
        IF (carrier) THEN
          IF (system_2d == 'no') THEN
            WRITE(iufilsigma, '(a)') "# Optical mobility in cm^3 / Vs"
          ELSE
            WRITE(iufilsigma, '(a)') "# Optical mobility in cm^2 / Vs"
          ENDIF
            WRITE(iufilsigma, '(a)') "#      Omega (eV)     Re mobility_xx    Im mobility_xx    Re mobility_xy" // &
                         "    Im mobility_xy     Re mobility_xz    Im mobility_xz    Re mobility_yx    Im mobility_yx" // &
                         "    Re mobility_yy    Im mobility_yy    Re mobility_yz    Im mobility_yz     Re mobility_zx" // &
                         "    Im mobility_zx    Re mobility_zy    Im mobility_zy    Re mobility_zz    Im mobility_zz"
        ELSE
          IF (system_2d == 'no') THEN
            WRITE(iufilsigma, '(a)') "# Optical conductivity in 1/(Ohm * m)"
          ELSE
            WRITE(iufilsigma, '(a)') "# Optical conductivity in 1/(Ohm)"
       	  ENDIF
            WRITE(iufilsigma, '(a)') "#      Omega (eV)     Re Sigma_xx    Im Sigma_xx    Re Sigma_xy    Im Sigma_xy" // &
                         "    Re Sigma_xz    Im Sigma_xz    Re Sigma_yx    Im Sigma_yx" // &
                         "    Re Sigma_yy    Im Sigma_yy    Re Sigma_yz    Im Sigma_yz" // &
                         "    Re Sigma_zx    Im Sigma_zx    Re Sigma_zy    Im Sigma_zy  Re Sigma_zz    Im Sigma_zz"
        ENDIF
      ENDIF
      !
      iw0 = 0
      DO iw = 1, nw_specfun
        sigma_opt(:, :, iw) = zero
        sigma_opt_im(:, :, iw) = zero
        ! Compute real part of conductivity
        DO ik = 1, nkf
          !
          ik_global = ik + lower_bnd - 1
          ikk = 2 * ik - 1
          ! Intraband + Interband loops
          DO ibnd = 1, nbndfst
            !                
            IF (ABS(etf_all(ibndmin - 1 + ibnd, ik_global* 2 - 1) - ef0) > fsthick) CYCLE
            DO jbnd = 1, nbndfst
              IF (ABS(etf_all(ibndmin - 1 + jbnd, ik_global * 2 - 1) - ef0) > fsthick) CYCLE
              IF (lwfpt) THEN
                !
                ! Skip active states outside the ahc window
                IF (etf_all(ibnd -1 + ibndmin, ik_global * 2 - 1) < ahc_win_min .OR. &
                   etf_all(ibndmin - 1 + ibnd, ik_global * 2 - 1) > ahc_win_max) CYCLE
                IF (etf_all(jbnd -1 + ibndmin, ik_global * 2 - 1) < ahc_win_min .OR. &
                   etf_all(ibndmin - 1 + jbnd, ik_global * 2 - 1) > ahc_win_max) CYCLE
              ENDIF
              ! Velocity matrix elements
              IF (specfun_el_scgd0 .AND. iter_scgd0 > 0) THEN
                vkk(:) = vmef(:, ibnd, jbnd, ixkf(ik_global))
              ELSE
                vkk(:) = vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikk)
              ENDIF
              !here we precompute the volocity prefactor because it is the only difference between
              ! mp_mesh_k and no symmetries case
              velocity_factor(:, :) = zero
              vk_cart(:) = zero
              IF (mp_mesh_k) THEN
                DO ikbz = 1, nkf1 * nkf2 * nkf3
                  ! If the k-point from the full BZ is related by a symmetry operation
                  ! to the current k-point, then take it.
                  IF (bztoibz(ikbz) == ik_global) THEN
                    ! Transform the symmetry matrix from Crystal to cartesian
                    sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
                    sb    = MATMUL(bg, sa)
                    sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                    !Rotate the Real and Imaginary parts separately
                    v_in_re(:) = REAL(vkk(:))
                    v_in_im(:) = AIMAG(vkk(:))
                    v_out_re = 0.d0
                    v_out_im = 0.d0
                    CALL DGEMV('n', 3, 3, 1.d0, sr, 3, v_in_re, 1, 0.d0, v_out_re, 1) 
                    CALL DGEMV('n', 3, 3, 1.d0, sr, 3, v_in_im, 1, 0.d0, v_out_im, 1)
                    !Reconstruct the rotated complex vector
                    v_rot(:) = CMPLX(v_out_re(:), v_out_im(:), KIND=DP)
                    DO i = 1, 3
                      DO j = 1,3
                        IF (noncolin) THEN
                          velocity_factor(i, j) = velocity_factor(i, j) + & 
                           1.d0 / (nkf1 * nkf2 * nkf3) * REAL(CONJG(v_rot(i)) * v_rot(j))
                        ELSE
                          velocity_factor(i, j) = velocity_factor(i, j) + &
                           2.d0 / (nkf1 * nkf2 * nkf3) * REAL(CONJG(v_rot(i)) * v_rot(j))
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ELSE
                DO i = 1, 3
                  DO j = 1, 3
                    velocity_factor(i, j) = wkf_all(ik_global * 2 - 1) *REAL(CONJG(vkk(i)) * vkk(j))
                  ENDDO
                ENDDO
              ENDIF
              DO i = 1, 3
                DO j = 1, 3
                  DO ie = 1, nw_specfun
                    ww_plus = ww(iw) + ww(ie)
                    iw_plus = NINT((ww_plus - ww(1)) / dw) + 1
                    IF (iw_plus > 0 .AND. iw_plus < 1 + nw_specfun) THEN
                      IF ( ABS(ww(iw)) < 1.0D-5) THEN
                        sigma_opt(i, j, iw) = sigma_opt(i, j, iw) + velocity_factor(i, j) * dw  * pi *&
                        w0gauss( (ww(ie) + ef - mu_t(itemp))/ gtemp(itemp), -99) / gtemp(itemp) * &
                        a_all_ibnd(ie, ik_global, ibnd, itemp) * a_all_ibnd(iw_plus, ik_global, jbnd, itemp)
                      ELSE
                        sigma_opt(i, j, iw) = sigma_opt(i, j, iw) +  velocity_factor(i, j) *dw  * pi *&  
                        (fermi_dirac(ww(ie) + ef - mu_t(itemp), gtemp(itemp)) - &
                        fermi_dirac(ww_plus + ef - mu_t(itemp), gtemp(itemp))) / ww(iw) * &
                        a_all_ibnd(ie, ik_global, ibnd, itemp) * a_all_ibnd(iw_plus, ik_global, jbnd, itemp)
                      ENDIF
                    ENDIF
                  ENDDO ! energy integral
                ENDDO ! i
              ENDDO ! j
            ENDDO ! jbnd
          ENDDO !ibnd
        ENDDO ! k pts
        CALL mp_sum(sigma_opt(:,:,iw), inter_pool_comm)
        !
      ENDDO ! frequency loop
      DO i = 1, 3
        DO j = 1, 3
          CALL hilbert_transform(ww, nw_specfun, -sigma_opt(i, j, :), sigma_opt_im(i, j, :))
        ENDDO
      ENDDO
      DO iw = 1, nw_specfun
        IF (ww(iw) >= 0.0d0) THEN
          IF (iw0 == 0) iw0 = iw
          IF (ionode)  WRITE(iufilsigma, 103) ryd2ev * ww(iw), &
            sigma_opt(1,1,iw) * conv_factor, sigma_opt_im(1,1,iw) * conv_factor, &
            sigma_opt(1,2,iw) * conv_factor, sigma_opt_im(1,2,iw) * conv_factor, &
            sigma_opt(1,3,iw) * conv_factor, sigma_opt_im(1,3,iw) * conv_factor, &
            sigma_opt(2,1,iw) * conv_factor, sigma_opt_im(2,1,iw) * conv_factor, &
            sigma_opt(2,2,iw) * conv_factor, sigma_opt_im(2,2,iw) * conv_factor, &
            sigma_opt(2,3,iw) * conv_factor, sigma_opt_im(2,3,iw) * conv_factor, &
            sigma_opt(3,1,iw) * conv_factor, sigma_opt_im(3,1,iw) * conv_factor, &
            sigma_opt(3,2,iw) * conv_factor, sigma_opt_im(3,2,iw) * conv_factor, &
            sigma_opt(3,3,iw) * conv_factor, sigma_opt_im(3,3,iw) * conv_factor
        ENDIF ! positive freqs
      ENDDO
      ! print the static limit in stdout
      IF (ionode) THEN
        IF (carrier) THEN
          IF (system_2d == 'no') THEN
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_xx = ', sigma_opt(1,1,iw0) * conv_factor, '[cm^3 / Vs]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_yy = ', sigma_opt(2,2,iw0) * conv_factor, '[cm^3 / Vs]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_zz = ', sigma_opt(3,3,iw0) * conv_factor, '[cm^3 / Vs]'
          ELSE
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_xx = ', sigma_opt(1,1,iw0) * conv_factor, '[cm^2 / Vs]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_yy = ', sigma_opt(2,2,iw0) * conv_factor, '[cm^2 / Vs]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude mobility_zz = ', sigma_opt(3,3,iw0) * conv_factor, '[cm^2 / Vs]'
          ENDIF
        ELSE
          IF (system_2d == 'no') THEN
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_xx = ', sigma_opt(1,1,iw0) * conv_factor, '[1/(Ohm m)]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_yy = ', sigma_opt(2,2,iw0) * conv_factor, '[1/(Ohm m)]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_zz = ', sigma_opt(3,3,iw0) * conv_factor, '[1/(Ohm m)]'
          ELSE
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_xx = ', sigma_opt(1,1,iw0) * conv_factor, '[1/Ohm]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_yy = ', sigma_opt(2,2,iw0) * conv_factor, '[1/Ohm]'
            WRITE(stdout,'(/,5x,a,ES18.8,a)') 'Drude conductivity_zz = ', sigma_opt(3,3,iw0) * conv_factor, '[1/Ohm]'
          ENDIF
        ENDIF
      ENDIF
      IF (ionode) CLOSE(iufilsigma)
      !
      ! fermi level will be adjusted in scgd0, so return the old value
      IF (specfun_el_scgd0) mu_t(itemp) = ef
    ENDDO ! itemp
    DEALLOCATE(wkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('optical_conductivity', 'Error deallocating wkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('optical_conductivity', 'Error deallocating etf_all', 1)
    !
    103 FORMAT(ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, &
        ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5, 2x, ES13.5)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE optical_conductivity
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE interpolate_path_scgd0(kpt, q0, a_all, a_path)
    !-----------------------------------------------------------------------
    !!
    !! This subroutine is needed for scGD0 if the user want the spectral
    !! funciton to be computed on a k-path.
    !! The subroutine interpolates the spectral function or the self-energy
    !! computed on a full k grid, onto a desired k point.
    !! It does so by computing the distance of the desired point to all the 
    !! points in the full grid, and then it searches for 8 nearest neighbours.
    !! By assigning them a weight, the value at the desired k point is obtained.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode
    USE io_var,        ONLY : iufilsigma
    USE input,         ONLY : nbndsub, wmin_specfun, wmax_specfun, nw_specfun, &
                              efermi_read, fermi_energy, nstemp, lsda,         &
                              specfun_el_scgd0, lwfpt, ncarrier, system_2d,    &
                              opt_cond, carrier
    USE global_var,    ONLY : gtemp, etf, ibndmin, nkqf, nktotf, efnew, &
                              xkf, nkqtotf, nbndfst, lower_bnd, upper_bnd, wkf,&
                              vmef, nkf
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi,  &
                              bohr2ang, ang2cm, hbar, Ang2m, hbarJ, czero
    USE constants,        ONLY : electron_si
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE parallelism,   ONLY : poolgather2, poolgatherc4
    USE supercond_common, ONLY : xkfs, ekfs, ixkff, ekfs, wkfs, nkfs, nbndfs
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: kpt(3)
    !! k point to interpolate the quantity on
    REAL(KIND = DP), INTENT(IN) :: q0(3, nkqtotf)
    !! list of all computed points 
    REAL(KIND = DP), INTENT(IN) :: a_all(nw_specfun, nktotf)
    !! arrax computed on a grid of k points
    REAL(KIND = DP), INTENT(OUT) :: a_path(nw_specfun)
    !! interpolated array to kpt point
    !
    INTEGER :: ik
    !! counter on k points
    INTEGER :: ikk
    !! counter on k points
    INTEGER :: iw
    !! counter on frequencies
    INTEGER :: j
    !! counter on k points
    INTEGER :: i
    !! counter on k points
    INTEGER :: idx(nktotf)
    !! k point indices
    !
    REAL(KIND = DP) :: dist(nktotf)
    !! distance between kpt and all the k points
    REAL(KIND = DP) :: w(nktotf)
    !! weight; 1/distance
    REAL(KIND = DP) ::  wsum
    !! sum of weights
    REAL(KIND = DP) :: temp_dist
    !! temporary variable for k point index
    !
    a_path = 0.0
    ! Compute distance to all k-points
    DO ik = 1, nktotf
      !
      ikk = 2 * ik -1
      !
      dist(ik) = SQRT( (kpt(1)-q0(1, ikk))**2 + (kpt(2)-q0(2, ikk))**2 + (kpt(3)-q0(3, ikk))**2 )
      idx(ik) = ik
      ! 
    ENDDO
    ! Sort distances to find nearest neighbors. we search 8 neareset neighbors.  
    DO i = 1, 8
      DO j = i + 1, nktotf
        IF (dist(idx(j)) < dist(idx(i))) THEN
          temp_dist = idx(i)
          idx(i) = idx(j)
          idx(j) = temp_dist
        ENDIF
      ENDDO
    ENDDO
    ! Compute weighted sum
    wsum = 0.0
    DO i = 1, 8
      IF (dist(idx(i)) == 0.0) THEN
        ! Exactly on a grid point
        a_path = a_all(:, idx(i))
        RETURN
      ENDIF
      w(i) = 1.0D0 / dist(idx(i))
      wsum = wsum + w(i)
    ENDDO
    DO i = 1, 8
      DO iw = 1, nw_specfun
        a_path(iw) = a_path(iw) + w(i) * a_all(iw, idx(i))
      ENDDO
    ENDDO
    ! Normalize
    a_path = a_path / wsum
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE interpolate_path_scgd0
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE interpolate_path_scgd0_ekk(kpt, q0, ek, ekk)
    !-----------------------------------------------------------------------
    !!
    !! This subroutine is needed for scGD0 if the user want the spectral
    !! funciton to be computed on a k-path.
    !! The subroutine interpolates the  energy
    !! computed on a full k grid, onto a desired k point.
    !! It does so by computing the distance of the desired point to all the 
    !! points in the full grid, and then it searches for 8 nearest neighbours.
    !! By assigning them a weight, the value at the desired k point is obtained.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode
    USE global_var,    ONLY : gtemp, etf, ibndmin, nkqf, nktotf, efnew, &
                              xkf, nkqtotf, nbndfst, lower_bnd,         &
                              upper_bnd, wkf, vmef, nkf
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE parallelism,   ONLY : poolgather2, poolgatherc4
    USE supercond_common, ONLY : xkfs, ekfs, ixkff, wkfs, nkfs, nbndfs
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: kpt(3)
    !! k point to interpolate the quantity on
    REAL(KIND = DP), INTENT(IN) :: q0(3, nkqtotf)
    !! list of all computed points 
    REAL(KIND = DP), INTENT(IN) :: ek(nkqtotf)
    !! list of all computed energies
    REAL(KIND = DP), INTENT(OUT) :: ekk
    !! interpolated energy value
    !
    INTEGER :: ik
    !! counter on k points
    INTEGER :: ikk
    !! counter on k points
    INTEGER :: iw
    !! counter on frequencies
    INTEGER :: j
    !! counter on k points
    INTEGER :: i
    !! counter on k points
    INTEGER :: idx(nktotf)
    !! k point indices
    !
    REAL(KIND = DP) :: dist(nktotf)
    !! distance between kpt and all the k points
    REAL(KIND = DP) :: w(nktotf)
    !! weight; 1/distance
    REAL(KIND = DP) ::  wsum
    !! sum of weights
    REAL(KIND = DP) :: temp_dist
    !! temporary variable for k point index
    !
    ekk = 0.0
    ! Compute distance to all k-points
    DO ik = 1, nktotf
      !
      ikk = 2 * ik -1
      !
      dist(ik) = SQRT( (kpt(1)-q0(1, ikk))**2 + (kpt(2)-q0(2, ikk))**2 + (kpt(3)-q0(3, ikk))**2 )
      idx(ik) = ik
      ! 
    ENDDO
    ! Sort distances to find nearest neighbors (simple selection sort for small nneigh)
    DO i = 1, 8
      DO j = i + 1, nktotf
        IF (dist(idx(j)) < dist(idx(i))) THEN
          temp_dist = idx(i)
          idx(i) = idx(j)
          idx(j) = temp_dist
        ENDIF
      ENDDO
    ENDDO
    ! Compute weighted sum
    wsum = 0.0
    DO i = 1, 8
      IF (dist(idx(i)) == 0.0) THEN
        ! Exactly on a grid point
        ekk = ek(idx(i) * 2 - 1)
        RETURN
      ENDIF
      w(i) = 1.0D0 / dist(idx(i))
      wsum = wsum + w(i)
    ENDDO
    DO i = 1, 8
        ekk = ekk + w(i) * ek(idx(i) * 2 - 1)
    ENDDO
    ! Normalize
    ekk = ekk / wsum
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE interpolate_path_scgd0_ekk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_ph_q(iqq, iq, totq, nbndsub_tot, wkf_tot, etf_tot)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation. This corresponds to
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation.  This routine is similar to the one above
    !! but it is ONLY called from within ephwann_shuffle and calculates
    !! the selfenergy for one phonon at a time.  Much smaller footprint on the disk
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode
    USE io_var,    ONLY : iospectral_sup, iospectral
    USE modes,     ONLY : nmodes
    USE input,     ONLY : fsthick, shortrange, ngaussw, degaussw, &
                          nsmear, delta_smear, eps_acoustic, nstemp, &
                          wmin_specfun, wmax_specfun, nw_specfun, &
                          phonselfen
    USE pwcom,     ONLY : nelec, ef
    USE input,     ONLY : isk_dummy, lsda
    USE global_var,ONLY : gtemp, epf17, ibndmin, nbndfst, &
                          nkqf, wf, a_all_ph, pi_0, &
                          pir_all, gammai_all
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, cone, ci, eps8
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE selfen,        ONLY : selfen_phon_q
    !
    IMPLICIT NONE
    !
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: totq
    !! Total q-points in selecq window
    INTEGER, INTENT(in) :: nbndsub_tot
    !! Total number of bands 
    REAL(KIND = DP), INTENT(in) :: wkf_tot(nkqf)
    !! Integration weights
    REAL(KIND = DP), INTENT(in) :: etf_tot(nbndsub_tot,nkqf)
    !! Interpolated eigenvalues
    !
    ! Local variables
    INTEGER :: ismear
    !! Number of smearing values for the Gaussian function
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: iw
    !! Counter on frequency for the phonon spectra
    INTEGER :: itemp
    !! Counter on temperature
    INTEGER :: iqq_write
    !! Counter on q-point to write files
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: inv_pi
    !! Inverse pi
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: dwq(nmodes)
    !! $2\omega_{q\nu}$ defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: degaussw0
    !! degaussw0 = (ismear-1) * delta_smear + degaussw
    REAL(KIND = DP) :: eta
    !! artificial broadening of 0.1 meV for the visualization of the spectral function
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    fnm = ''
    IF (TRIM(lsda) == 'down') fnm = '.down'
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    eta = 0.1d0 / ryd2mev
    !
    !
    ! Now pre-treat phonon modes for efficiency
    ! Treat phonon frequency and Bose occupation
    wq(:)    = zero
    dwq(:)   = zero
    DO imode = 1, nmodes
      wq(imode) = wf(imode, iq)
      dwq(imode) = two * wq(imode)
      IF (wq(imode) > eps_acoustic) THEN
        inv_wq(imode) = one / (two * wq(imode))
      ELSE
        inv_wq(imode) = zero
      ENDIF
    ENDDO
    !
    !
    IF (.NOT. phonselfen) THEN
      CALL selfen_phon_q(iqq, iq, totq, nbndsub_tot, wkf_tot, etf_tot)
    ENDIF
    !
    ! collect contributions from all pools (sum over k-points) this finishes the integral over the BZ  (k)
    !
    CALL mp_sum(gammai_all, inter_pool_comm)
    CALL mp_sum(pir_all, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    !
    DO itemp = 1, nstemp
      ! Thomas-Fermi screening according to Resta PRB 1977
      ! Here specific case of Diamond
      !eps0   = 5.7
      !rtf    = 2.76
      !qtf    = 1.36
      !qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2) * tpiba2
      !epstf =  (qtf**2 + qsquared) / (qtf**2/eps0 * sin (sqrt(qsquared)*rtf)/(sqrt(qsquared)*rtf)+qsquared)
      !
      IF (iqq == 1) THEN
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Phonon Spectral Function Self-Energy in the Migdal Approximation (on the fly)")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
      ENDIF
      !
      DO ismear = 1, nsmear
        !
        degaussw0 = (ismear - 1) * delta_smear + degaussw
        !
        ! SP: Multiplication is faster than division ==> Important if called a lot
        !     in inner loops
        inv_pi       = one / pi
        !
        WRITE(stdout, '(5x, a)')
        !
        IF (iqq == 1 .and. itemp == 1) THEN
          IF (ionode) THEN
            WRITE(tp, "(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
            filespecsup = 'specfun_sup.phon' // TRIM(fnm) ! // trim(adjustl(tp)) // 'K'
            OPEN(UNIT = iospectral_sup, FILE = TRIM(filespecsup))
            WRITE(iospectral_sup, '(2x, a)') '#Phonon eigenenergies + real and im part of phonon self-energy (meV)'
            WRITE(iospectral_sup, '(2x, a)') '#Q-point    Mode      Temp.[K]       smearing[eV]       w_q[eV]    &    
                                          &w[eV]     Real Pi(w, T_low)[meV]   Real Pi(w=0,T_high)[meV]     Im Pi(w, T_low)[meV]'
          ENDIF
        ENDIF
        !
        ! Write to output file
        WRITE(stdout, '(/5x, a)') 'Real and Imaginary part of the phonon self-energy (omega=0) without gamma0.'
        DO imode = 1, nmodes
          ! Real and Im part of Phonon self-energy at 0 freq.
          WRITE(stdout, 105) imode, ryd2ev * wq(imode), ryd2mev * pir_all(1, imode, itemp, ismear), &
                             ryd2mev * gammai_all(1, imode, itemp, ismear)
        ENDDO
        !
        ! Write to support files
        DO iw = 1, nw_specfun
          !
          DO imode = 1, nmodes
            !
            !a_all(iw,iq) = a_all(iw,iq) + ABS(gammai_all(imode,iq,iw) ) / pi / &
            !      ( ( ww - wq - pir_all (imode,iq,iw) + pi_0 (imode))**two + (gammai_all(imode,iq,iw) )**two )
            ! SP: From Eq. 16 of PRB 9, 4733 (1974)
            !    Also in Eq.2 of PRL 119, 017001 (2017).
            ! in this spectral function, phonons are unscreened with the static part calculated at the high smearing values from the DFPT calculation
            a_all_ph(iw, iqq, itemp, ismear) = a_all_ph(iw, iqq, itemp, ismear) + inv_pi * dwq(imode) * &
                                    (two * ww(iw) * eta - dwq(imode) * gammai_all(iw, imode, itemp, ismear)) / &
                                    ((ww(iw)**two - eta**two - wq(imode)**two - dwq(imode) * &
                                    (pir_all(iw, imode, itemp, ismear) - pi_0(imode)))**two + (two*ww(iw) * eta - dwq(imode) * &
                                    gammai_all(iw, imode, itemp, ismear))**two)
            !
            IF (ionode) THEN
              WRITE(iospectral_sup, 102) iq, imode, gtemp(itemp) * ryd2ev / kelvin2eV, degaussw0 * ryd2ev, ryd2ev * wq(imode), &
                                       ryd2ev * ww(iw), ryd2mev * pir_all(iw, imode, itemp, ismear), ryd2mev * pi_0(imode), &
                                       ryd2mev * gammai_all(iw, imode, itemp, ismear)
            ENDIF
            !
          ENDDO
          !
!         IF (ionode) THEN
!           WRITE(iospectral, 103) iq, ryd2ev * ww(iw), a_all_ph(iw, iqq) / ryd2mev ! print to file
!         ENDIF
          !
        ENDDO !iw
      ENDDO !ismear
    ENDDO ! itemp
    !
    IF (iqq == totq) THEN
      IF (ionode) THEN
        DO itemp = 1, nstemp
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespec = 'specfun.phon.' // trim(adjustl(tp)) // 'K'//TRIM(fnm)
          OPEN(UNIT = iospectral, FILE = TRIM(filespec))                
          WRITE(iospectral, '(/2x, a)') '#Phonon spectral function (meV)'
          WRITE(iospectral, '(/2x, a)') '#Q-point    Energy[eV]      smearing[eV]     A(q,w)[meV^-1]'
          DO ismear = 1, nsmear
            degaussw0 = (ismear - 1) * delta_smear + degaussw
            DO iqq_write = 1, totq
              DO iw = 1, nw_specfun
                WRITE(iospectral, 103) iqq_write, ryd2ev * ww(iw), degaussw0 * ryd2ev, &
                                       a_all_ph(iw, iqq_write, itemp, ismear) / ryd2mev ! print to file
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iospectral)
        ENDDO
        CLOSE(iospectral_sup)
      ENDIF
    ENDIF
    WRITE(stdout, '(5x, a/)') REPEAT('-',67)
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6,' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    102 FORMAT(2i9, 2x, f8.3, 2x, f12.5, 2x, f12.5, 2x, f12.5, 2x, E22.14, 2x, E22.14, 2x, E22.14)
    103 FORMAT(2x, i7, 2x, f12.5, 2x, f12.5, 2x, E22.14)
    105 FORMAT(5x, 'Omega( ', i3, ' )=', f9.4,' eV   Re[Pi]=', f15.6, ' meV Im[Pi]=', f15.6, ' meV')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_ph_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_pl_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation.
    !!
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode
    USE io_var,        ONLY : iospectral_sup, iospectral
    USE input,         ONLY : nbndsub, fsthick, ngaussw, degaussw, nw_specfun, &
                              wmin_specfun, wmax_specfun, efermi_read, fermi_energy, &
                              nstemp, nel, meff, epsiheg, restart, restart_step, lsda
    USE pwcom,         ONLY : ef
    USE global_var,    ONLY : etf, ibndmin, nkqf, nbndfst, nkf, wqf, xkf, &
                              nkqtotf, xqf, vmef, esigmar_all, esigmai_all, a_all, &
                              gtemp, nktotf, lower_bnd, efnew
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, ci, eps6
    USE ep_constants,  ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega, alat, bg
    USE selfen,        ONLY : get_eps_mahan
    USE io_selfen,     ONLY : spectral_write
    USE parallelism,  ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-point index in selecq
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in fsthick window
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekk1
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq
    !! Plasmon frequency
    REAL(KIND = DP) :: wgq
    !! Bose occupation factor $n_{q wpl}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: weight
    !! SE factor
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature defined for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of degaussw defined for efficiency reasons
    REAL(KIND = DP) :: sq_degaussw
    !! Squared degaussw defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp
    !! Temporary variable defined for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Spectral frequency increment
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP) :: esigmar0
    !! static SE
    REAL(KIND = DP) :: tpiba_new
    !! 2 \pi / alat
    REAL(KIND = DP) :: kf
    !! Fermi wave-vector
    REAL(KIND = DP) :: vf
    !! Fermi velocity
    REAL(KIND = DP) :: fermiheg
    !! Fermi energy of a homageneous electron gas
    REAL(KIND = DP) :: qnorm
    !! |q|
    REAL(KIND = DP) :: qin
    !! (2 \pi / alat) |q|
    REAL(KIND = DP) :: sq_qin
    !! Squared qin defined for efficiency reasons
    REAL(KIND = DP) :: wpl0
    !! Plasmon frequency
    REAL(KIND = DP) :: eps0
    !! Dielectric function at zero frequency
    REAL(KIND = DP) :: deltaeps
    !!
    REAL(KIND = DP) :: qcut
    !! Cut-off of the maximum wave-vector of plasmon modes (qcut = wpl0 / vf)
    REAL(KIND = DP) :: qtf
    !! Thomas-Fermi screening wave-vector
    REAL(KIND = DP) :: dipole
    !! Dipole
    REAL(KIND = DP) :: rs
    !! Spherical radius used to describe the density of an electron gas
    REAL(KIND = DP) :: degen
    !! Degeneracy of the electron gas
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: q(3)
    !! The q-point in cartesian unit.
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! Spectral function
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    COMPLEX(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekq - wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekq + wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmpw1
    !! Temporary variable to store etmpw1 = ww - etmp1
    COMPLEX(KIND = DP) :: etmpw2
    !! Temporary variable to store etmpw1 = ww - etmp2
    COMPLEX(KIND = DP) :: fact
    !! SE factor
    !
    ! loop over temperatures can be introduced
    !
    fnm = ''
    IF (TRIM(lsda) == 'down') fnm = '.down'
    !
    inv_degaussw = one / degaussw
    sq_degaussw = degaussw * degaussw
    DO itemp = 1, nstemp
      inv_eptemp   = one / gtemp(itemp)
      ! energy range and spacing for spectral function
      !
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      !
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron Spectral Function in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
        !
      ENDIF
      !
      ! Fermi level
      !
      IF (efermi_read) THEN
        ef0 = fermi_energy
      ELSE
        ef0 = efnew
      ENDIF
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      ! SP: Sum rule added to conserve the number of electron.
      IF (iqq == 1) THEN
        WRITE(stdout, '(5x, a)') 'The sum rule to conserve the number of electron is enforced.'
        WRITE(stdout, '(5x, a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
        WRITE(stdout, '(5x, a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      !nel      =  0.01    ! this should be read from input - # of doping electrons
      !epsiheg  =  12.d0   ! this should be read from input - # dielectric constant at zero doping
      !meff     =  0.25    ! this should be read from input - effective mass
      !
      tpiba_new = two * pi / alat
      degen     = one
      !
      ! Based on Eqs. (5.3)-(5.6) and (5.127) of Mahan 2000.
      !
      ! omega is the unit cell volume in Bohr^3
      rs = (3.d0 / (4.d0 * pi * nel / omega / degen))**(1.d0 / 3.d0) * meff * degen
      kf = (3.d0 * (pi**2.d0) * nel / omega / degen)**(1.d0 / 3.d0)
      vf = (1.d0 / meff) * kf
      !
      ! fermiheg in [Ry] (multiplication by 2 converts from Ha to Ry)
      fermiheg = 2.d0 * (1.d0 / (2.d0 * meff)) * kf**2.d0
      ! qtf in ! [a.u.]
      qtf = DSQRT(6.d0 * pi * nel / omega / degen / (fermiheg / 2.d0))
      ! wpl0 in [Ry] (multiplication by 2 converts from Ha to Ry)
      wpl0 = two * DSQRT(4.d0 * pi * nel / omega / meff / epsiheg)
      wq = wpl0
      !
      q(:) = xqf(:, iq)
      CALL cryst_to_cart(1, q, bg, 1)
      qnorm = DSQRT(q(1)**two + q(2)**two + q(3)**two)
      qin = qnorm * tpiba_new
      sq_qin = qin * qin
      !
      ! qcut in [Ha] (1/2 converts from Ry to Ha)
      qcut = wpl0 / vf / tpiba_new / 2.d0
      !
      !IF (.TRUE.) qcut = qcut / 2.d0 ! renormalize to account for Landau damping
      !
      ! qin should be in atomic units for Mahan formula
      CALL get_eps_mahan(qin, rs, kf, eps0)
      deltaeps = -(1.d0 / (epsiheg + eps0 - 1.d0) - 1.d0 / epsiheg)
      !
      g2_tmp = 4.d0 * pi * (wq * deltaeps / 2.d0) / omega * 2.d0
      !
      IF (iqq == 1) THEN
        WRITE(stdout, '(12x, " nel       = ", E15.6)') nel
        WRITE(stdout, '(12x, " meff      = ", E15.6)') meff
        WRITE(stdout, '(12x, " rs        = ", E15.6)') rs
        WRITE(stdout, '(12x, " kf        = ", E15.6)') kf
        WRITE(stdout, '(12x, " vf        = ", E15.6)') vf
        WRITE(stdout, '(12x, " fermi_en  = ", E15.6)') fermiheg
        WRITE(stdout, '(12x, " qtf       = ", E15.6)') qtf
        WRITE(stdout, '(12x, " wpl       = ", E15.6)') wpl0
        WRITE(stdout, '(12x, " qcut      = ", E15.6)') qcut
        WRITE(stdout, '(12x, " eps0      = ", E15.6)') eps0
        WRITE(stdout, '(12x, " epsiheg   = ", E15.6)') epsiheg
        WRITE(stdout, '(12x, " deltaeps  = ", E15.6)') deltaeps
      ENDIF
      !
      IF (restart) THEN
        ! Make everythin 0 except the range of k-points we are working on
        esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmar_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmai_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        !
      ENDIF
      !
      ! In the case of a restart do not add the first step
      IF (first_cycle .and. itemp == nstemp) THEN
        first_cycle = .FALSE.
      ELSE
        IF (qnorm < qcut) THEN
          !
          ! wq is the plasmon frequency
          ! Bose occupation
          wgq = wgauss(-wq * inv_eptemp, -99)
          wgq = wgq / (one - two * wgq)
          !
          ! loop over all k points of the fine mesh
          !
          fermicount = 0
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle (but in this case they are the same)
            !
            IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
                (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
              !
              fermicount = fermicount + 1
              !
              DO ibnd = 1, nbndfst
                !
                !  the energy of the electron at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
                !
                DO jbnd = 1, nbndfst
                  !
                  ekk1 = etf(ibndmin - 1 + jbnd, ikk) - ef0
                  ! the energy of the electron at k+q (relative to Ef)
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  ! the Fermi occupation at k+q
                  wgkq = wgauss(-ekq * inv_eptemp, -99)
                  !
                  ! Computation of the dipole
                  IF (ibnd == jbnd) THEN
                    IF (qnorm > eps6) THEN
                      dipole = one / sq_qin
                    ELSE
                      dipole = zero
                    ENDIF
                  ELSE
                    IF (ABS(ekq - ekk1) > eps6) THEN
                      ! TODO: Check the expression to confirm that division by 2 is correct.
                      dipole = REAL(      vmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0 *  &
                                    CONJG(vmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0) / &
                                    ((ekk1 - ekk)**two + sq_degaussw))
                    ELSE
                      dipole = zero
                    ENDIF
                  ENDIF
                  !
                  ! The q^-2 is cancelled by the q->0 limit of the dipole.
                  ! See e.g., pg. 258 of Grosso Parravicini.
                  ! electron-plasmon scattering matrix elements squared
                  g2 = dipole * g2_tmp
                  !
                  fact1 =       wgkq + wgq
                  fact2 = one - wgkq + wgq
                  etmp1 = ekq - wq + ci * degaussw
                  etmp2 = ekq + wq + ci * degaussw
                  !
                  DO iw = 1, nw_specfun
                    !
                    etmpw1 = ww(iw) - etmp1
                    etmpw2 = ww(iw) - etmp2
                    !
                    fact = (fact1 / etmpw1) + (fact2 / etmpw2)
                    !
                    weight = wqf(iq) * REAL(fact)
                    !
                    ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
                    ! SP : Application of the sum rule
                    esigmar0 = - g2 *  wqf(iq) * REAL((fact1 / etmp1) + (fact2 / etmp2))
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) - &
                                                                       esigmar0
                    !
                    weight = wqf(iq) * AIMAG(fact)
                    !
                    ! \Im\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
!                    WRITE(stdout, '(5x, f8.5, f8.5)') g2_tmp, dipole
                  ENDDO
                ENDDO !jbnd
              ENDDO !ibnd
            ENDIF ! endif  fsthick
          ENDDO ! end loop on k
        ENDIF ! endif qnorm
        !
        ! Creation of a restart point
        IF (restart) THEN
          IF (MOD(iqq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iqq
            CALL mp_sum(esigmar_all, inter_pool_comm)
            CALL mp_sum(esigmai_all, inter_pool_comm)
            CALL mp_sum(fermicount, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      ! Collect pools and write the spectral function
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! Note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      DO itemp = 1, nstemp
        inv_eptemp = one / gtemp(itemp)
        !
        ! Output electron spectral function here after looping over all q-points (with their contributions summed in a etc.)
        !
        WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
        !
        ! construct the trace of the spectral function (assume diagonal selfenergy
        ! and constant matrix elements for dipole transitions)
        !
        IF (ionode) then
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespec = 'specfun.plself.' // trim(adjustl(tp)) // 'K'// TRIM(fnm)
          filespecsup = 'specfun_sup.plself.' // trim(adjustl(tp)) // 'K' // TRIM(fnm)
          OPEN(UNIT = iospectral, FILE = TRIM(filespec) )
          OPEN(UNIT = iospectral_sup, FILE = TRIM(filespecsup) )
          WRITE(iospectral, '(/2x, a/)') '#Electron-plasmon spectral function (meV)'
          WRITE(iospectral_sup, '(/2x, a/)') '#KS eigenenergies + real and im part of electron-plasmon self-energy (meV)'
          WRITE(iospectral, '(/2x, a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
          WRITE(iospectral_sup, '(/2x, a/)') '#K-point    Band   e_nk[eV]   w[eV]       Real Sigma[meV]  Im Sigma[meV]'
        ENDIF
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout, '(/5x, "ik = ", i5, " coord.: ", 3f12.7, " Temp.: ", f8.3 )') ik, xkf_all(:, ikk), &
                                                                                      gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO iw = 1, nw_specfun
            !
            DO ibnd = 1, nbndfst
              !
              !  the energy of the electron at k
              ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
              !
              a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp) ) / pi / &
                   ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
              !
            ENDDO
            !
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), a_all(iw, ik, itemp) / ryd2mev
            !
          ENDDO
          !
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          !
        ENDDO
        !
        DO ik = 1, nktotf
          !
          ! The spectral function should integrate to 1 for each k-point
          specfun_sum = zero
          !
          DO iw = 1, nw_specfun
            !
            fermi(iw) = wgauss(-ww(iw) * inv_eptemp, -99)
            specfun_sum = specfun_sum + a_all(iw, ik, itemp) * fermi(iw) * dw !/ ryd2mev
            !
           IF (ionode) WRITE(iospectral, '(2x, i7, 2x, f10.5, 2x, E12.5)') ik, ryd2ev * ww(iw), &
                                                                                 a_all(iw, ik, itemp) / ryd2mev
          ENDDO
          !
          IF (ionode) WRITE(iospectral, '(a)') ' '
          IF (ionode) WRITE(iospectral, '(2x, a, 2x, E12.5)') '# Integrated spectral function ', specfun_sum
        ENDDO
        !
        IF (ionode) CLOSE(iospectral)
        !
        DO ibnd = 1, nbndfst
          !
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            !  the energy of the electron at k
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            DO iw = 1, nw_specfun
              !
              WRITE(stdout, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
              IF (ionode) &
              WRITE(iospectral_sup, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
            ENDDO
            !
          ENDDO
          !
          WRITE(stdout, *) ' '
          !
        ENDDO
        !
        IF (ionode) CLOSE(iospectral_sup)
        !
      ENDDO ! itemp
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating etf_all', 1)
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', e12.5, ' meV^-1')
    102 FORMAT(2i9, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_pl_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE a2f_main()
    !-----------------------------------------------------------------------
    !!
    !! Compute the Eliasberg spectral function
    !! in the Migdal approximation.
    !!
    !! If the q-points are not on a uniform grid (i.e. a line)
    !! the function will not be correct
    !!
    !! 02/2009 works in serial on ionode at the moment.  can be parallelized
    !! 03/2009 added transport spectral function -- this involves a v_k dot v_kq term
    !!         in the quantities coming from selfen_phon.f90.  Not fully implemented
    !! 10/2009 the code is transitioning to 'on-the-fly' phonon selfenergies
    !!         and this routine is not currently functional
    !! 10/2015 RM: added calcution of Tc based on Allen-Dynes formula
    !! 09/2019 SP: Cleaning
    !!
    !
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE cell_base,     ONLY : omega
    USE input,         ONLY : degaussq, delta_qsmear, nqsmear, nqstep, nsmear, eps_acoustic, &
                          nstemp, delta_smear, degaussw, fsthick, nc, lsda
    USE global_var,    ONLY : gtemp, nqtotf, wf, wqf, lambda_all, lambda_v_all
    USE ep_constants,  ONLY : ryd2mev, ryd2ev, kelvin2eV, one, two, zero, kelvin2Ry
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE io_global,     ONLY : ionode, stdout
    USE io_var,        ONLY : iua2ffil, iudosfil, iua2ftrfil, iures
    USE io_files,      ONLY : prefix
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    !! Temperature
    CHARACTER(LEN = 256) :: fila2f
    !! File name for Eliashberg spectral function
    CHARACTER(LEN = 256) :: fila2ftr
    !! File name for transport Eliashberg spectral function
    CHARACTER(LEN = 256) :: fildos
    !! File name for phonon density of states
    CHARACTER(LEN = 256) :: filres
    !! File name for resistivity
    CHARACTER(LEN = 256) :: fnm
    !! Buffer file name
    !
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: iq
    !! Counter on the q-point index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ismear
    !! Counter on smearing values (phonons)
    INTEGER :: isig
    !! Counter on smearing values (electrons)
    INTEGER :: i
    !! Counter on mu
    INTEGER :: itemp
    !! Counter on temperature
    INTEGER :: itemprho
    !! Counter on temparture for rho
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Factor in a2f
    REAL(KIND = DP) :: temp
    !! Temperature
    REAL(KIND = DP) :: n
    !! Carrier density
    REAL(KIND = DP) :: be
    !! Bose-Einstein distribution
    REAL(KIND = DP) :: prefact
    !! Prefactor in resistivity
    REAL(KIND = DP) :: lambda_tot
    !! Total e-ph coupling strength (summation)
    REAL(KIND = DP) :: lambda_tr_tot
    !! Total transport e-ph coupling strength (summation)
    REAL(KIND = DP) :: degaussq0
    !! Phonon smearing
    REAL(KIND = DP) :: inv_degaussq0
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: a2f_tmp
    !! Temporary variable for Eliashberg spectral function
    REAL(KIND = DP) :: a2f_tr_tmp
    !! Temporary variable for transport Eliashberg spectral function
    REAL(KIND = DP) :: om_max
    !! max phonon frequency increased by 10%
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: w0
    !! Current frequency w(imode, iq)
    REAL(KIND = DP) :: l
    !! Temporary variable for e-ph coupling strength
    REAL(KIND = DP) :: l_tr
    !! Temporary variable for transport e-ph coupling strength
    REAL(KIND = DP) :: tc
    !! Critical temperature
    REAL(KIND = DP) :: mu
    !! Coulomb pseudopotential
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: ww(nqstep)
    !! Current frequency
    REAL(KIND = DP), ALLOCATABLE :: a2f_(:, :)
    !! Eliashberg spectral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: a2f_tr(:, :)
    !! Transport Eliashberg spectral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: l_a2f(:)
    !! total e-ph coupling strength (a2f_ integration) for different ismear
    REAL(KIND = DP), ALLOCATABLE :: l_a2f_tr(:)
    !! total transport e-ph coupling strength (a2f_tr integration) for different ismear
    REAL(KIND = DP), ALLOCATABLE :: dosph(:, :)
    !! Phonon density of states for different for different ismear
    REAL(KIND = DP), ALLOCATABLE :: logavg(:)
    !! logavg phonon frequency for different ismear
    REAL(KIND = DP), ALLOCATABLE :: rho(:, :)
    !! Resistivity for different for different ismear
    !
    CALL start_clock('a2F')
    !
    fnm = ''
    IF (TRIM(lsda) == 'down') fnm = '.down'
    IF (ionode) THEN
      !
      ALLOCATE(a2f_(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating a2f_', 1)
      ALLOCATE(a2f_tr(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating a2f_tr', 1)
      ALLOCATE(dosph(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating dosph', 1)
      ALLOCATE(l_a2f(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating l_a2f', 1)
      ALLOCATE(l_a2f_tr(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating l_a2f_tr', 1)
      ALLOCATE(logavg(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating logavg', 1)
      ! The resitivity is computed for temperature between 0K-1000K by step of 10
      ! This is hardcoded and needs to be changed here if one wants to modify it
      ALLOCATE(rho(100, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating rho', 1)
      !
      DO itemp = 1, nstemp
        DO isig = 1, nsmear
          !
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          IF (isig < 10) THEN
            WRITE(fila2f,   '(a, a6, i1, a, a)') TRIM(prefix) // TRIM(fnm), '.a2f.0', isig, '.', trim(adjustl(tp))
            WRITE(fila2ftr, '(a, a9, i1, a, a)') TRIM(prefix) // TRIM(fnm), '.a2f_tr.0', isig, '.', trim(adjustl(tp))
            WRITE(filres,   '(a, a6, i1, a, a)') TRIM(prefix) // TRIM(fnm), '.res.0', isig, '.', trim(adjustl(tp))
            WRITE(fildos,   '(a, a8, i1, a, a)') TRIM(prefix) // TRIM(fnm), '.phdos.0', isig, '.', trim(adjustl(tp))
          ELSE
            WRITE(fila2f,   '(a, a5, i2, a, a)') TRIM(prefix) // TRIM(fnm), '.a2f.', isig, '.', trim(adjustl(tp))
            WRITE(fila2ftr, '(a, a8, i2, a, a)') TRIM(prefix) // TRIM(fnm), '.a2f_tr.', isig, '.', trim(adjustl(tp))
            WRITE(filres,   '(a, a5, i2, a, a)') TRIM(prefix) // TRIM(fnm), '.res.', isig, '.', trim(adjustl(tp))
            WRITE(fildos,   '(a, a7, i2, a ,a)') TRIM(prefix) // TRIM(fnm), '.phdos.', isig, '.', trim(adjustl(tp))
          ENDIF
          OPEN(UNIT = iua2ffil, FILE = fila2f, FORM = 'formatted')
          OPEN(UNIT = iua2ftrfil, FILE = fila2ftr, FORM = 'formatted')
          OPEN(UNIT = iures, FILE = filres, FORM = 'formatted')
          OPEN(UNIT = iudosfil, FILE = fildos, FORM = 'formatted')
          !
          WRITE(stdout, '(/5x, a)') REPEAT('=',67)
          WRITE(stdout, '(5x, "Eliashberg Spectral Function in the Migdal Approximation")')
          WRITE(stdout, '(5x, a/)') REPEAT('=',67)
          !
          om_max = 1.1d0 * MAXVAL(wf(:, :)) ! increase by 10%
          dw = om_max / DBLE(nqstep)
          DO iw = 1, nqstep  !
            ww(iw) = DBLE(iw) * dw
          ENDDO
          !
          lambda_tot    = zero
          l_a2f(:)      = zero
          a2f_(:, :)    = zero
          lambda_tr_tot = zero
          l_a2f_tr(:)   = zero
          a2f_tr(:, :)  = zero
          dosph(:, :)   = zero
          logavg(:)     = zero
          !
          DO ismear = 1, nqsmear
            !
            degaussq0 = degaussq + (ismear - 1) * delta_qsmear
            inv_degaussq0 = one / degaussq0
            !
            DO iw = 1, nqstep  ! loop over points on the a2F(w) graph
              !
              DO iq = 1, nqtotf ! loop over q-points
                DO imode = 1, nmodes ! loop over modes
                  w0 = wf(imode, iq)
                  !
                  IF (w0 > eps_acoustic) THEN
                    !
                    l = lambda_all(imode, iq, isig, itemp)
                    IF (lambda_all(imode, iq, isig, itemp) < 0.d0) l = zero ! sanity check
                    !
                    a2f_tmp = wqf(iq) * w0 * l / two
                    !
                    weight = w0gauss((ww(iw) - w0) * inv_degaussq0, 0) * inv_degaussq0
                    a2f_(iw, ismear) = a2f_(iw, ismear) + a2f_tmp * weight
                    dosph(iw, ismear) = dosph(iw, ismear) + wqf(iq) * weight
                    !
                    l_tr = lambda_v_all(imode, iq, isig, itemp)
                    IF (lambda_v_all(imode, iq, isig, itemp) < 0.d0) l_tr = zero !sanity check
                    !
                    a2f_tr_tmp = wqf(iq) * w0 * l_tr / two
                    !
                    a2f_tr(iw, ismear) = a2f_tr(iw, ismear) + a2f_tr_tmp * weight
                    !
                  ENDIF
                ENDDO
              ENDDO
              !
              ! output a2f
              !
              IF (ismear == nqsmear) WRITE(iua2ffil,   '(f12.7, 15f12.7)') ww(iw) * ryd2mev, a2f_(iw, :)
              IF (ismear == nqsmear) WRITE(iua2ftrfil, '(f12.7, 15f12.7)') ww(iw) * ryd2mev, a2f_tr(iw, :)
              IF (ismear == nqsmear) WRITE(iudosfil,   '(f12.7, 15f12.7)') ww(iw) * ryd2mev, dosph(iw, :) / ryd2mev
              !
              ! do the integral 2 int (a2F(w)/w dw)
              !
              l_a2f(ismear) = l_a2f(ismear) + two * a2f_(iw, ismear) / ww(iw) * dw
              l_a2f_tr(ismear) = l_a2f_tr(ismear) + two * a2f_tr(iw, ismear) / ww(iw) * dw
              logavg(ismear) = logavg(ismear) + two *  a2f_(iw, ismear) * LOG(ww(iw)) / ww(iw) * dw
              !
            ENDDO
            !
            logavg(ismear) = EXP(logavg(ismear) / l_a2f(ismear))
            !
          ENDDO
          !
          DO iq = 1, nqtotf ! loop over q-points
            DO imode = 1, nmodes ! loop over modes
              IF (lambda_all(imode, iq, isig, itemp) > 0.d0 .AND. wf(imode, iq) > eps_acoustic ) &
                lambda_tot = lambda_tot + wqf(iq) * lambda_all(imode, iq, isig, itemp)
              IF (lambda_v_all(imode, iq, isig, itemp) > 0.d0 .AND. wf(imode, iq) > eps_acoustic) &
                lambda_tr_tot = lambda_tr_tot + wqf(iq) * lambda_v_all(imode, iq, isig, itemp)
            ENDDO
          ENDDO
          WRITE(stdout, '(5x, a, f12.7)') "lambda : ", lambda_tot
          WRITE(stdout, '(5x, a, f12.7)') "lambda_tr : ", lambda_tr_tot
          WRITE(stdout, '(a)') " "
          !
          !
          ! Allen-Dynes estimate of Tc for ismear = 1
          !
          WRITE(stdout, '(5x, a, f12.7, a)') "Estimated Allen-Dynes Tc"
          WRITE(stdout, '(a)') " "
          WRITE(stdout, '(5x, a, f12.7, a, f12.7)') "logavg = ", logavg(1), " l_a2f = ", l_a2f(1)
          DO i = 1, 6
            !
            mu = 0.1d0 + 0.02d0 * DBLE(i - 1)
            tc = logavg(1) / 1.2d0 * EXP(-1.04d0 * (1.d0 + l_a2f(1)) / (l_a2f(1) - mu * ( 1.d0 + 0.62d0 * l_a2f(1))))
            ! tc in K
            !
            tc = tc * ryd2ev / kelvin2eV
            !SP: IF Tc is too big, it is not physical
            IF (tc < 1000.0) THEN
              WRITE(stdout, '(5x, a, f6.2, a, f22.12, a)') "mu = ", mu, " Tc = ", tc, " K"
            ENDIF
            !
          ENDDO
          !
          rho(:, :) = zero
          ! Now compute the Resistivity of Metal using the Ziman formula
          ! rho(T,smearing) = 4 * pi * me/(n * e**2 * kb * T) int dw hbar w a2F_tr(w,smearing) n(w,T)(1+n(w,T))
          ! n is the number of electron per unit volume and n(w,T) is the Bose-Einstein distribution
          ! Usually this means "the number of electrons that contribute to the mobility" and so it is typically 8 (full shell)
          ! but not always. You might want to check this.
          !
          n = nc / omega
          WRITE(iures, '(a)') '# Temperature [K]                &
                              Resistivity [micro Ohm cm] for different Phonon smearing (meV)        '
          WRITE(iures, '("#     ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          DO ismear = 1, nqsmear
            DO itemprho = 1, 100 ! Per step of 10K
              temp = itemprho * 10.d0 * kelvin2Ry
              ! omega is the volume of the primitive cell in a.u.
              !
              prefact = 4.d0 * pi / (temp * n)
              DO iw = 1, nqstep  ! loop over points on the a2F(w)
                !
                be = one / (EXP(ww(iw) / temp) - one)
                ! Perform the integral with rectangle.
                rho(itemprho, ismear) = rho(itemprho, ismear) + prefact * ww(iw) * a2f_tr(iw, ismear) * be * (1.d0 + be) * dw
                !
              ENDDO
              ! From a.u. to micro Ohm cm
              ! Conductivity 1 a.u. = 2.2999241E6 S/m
              ! Now to go from Ohm*m to micro Ohm cm we need to multiply by 1E8
              rho(itemprho, ismear) = rho(itemprho, ismear) * 1E8 / 2.2999241E6
              IF (ismear == nqsmear) WRITE (iures, '(i8, 15f12.7)') itemprho * 10, rho(itemprho, :)
            ENDDO
          ENDDO
          CLOSE(iures)
          !
          WRITE(iua2ffil, *) "Integrated el-ph coupling"
          WRITE(iua2ffil, '("  #         ", 15f12.7)') l_a2f(:)
          WRITE(iua2ffil, *) "Phonon smearing (meV)"
          WRITE(iua2ffil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          WRITE(iua2ffil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
          WRITE(iua2ffil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
          WRITE(iua2ffil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
          CLOSE(iua2ffil)
          !
          WRITE(iua2ftrfil, *) "Integrated el-ph coupling"
          WRITE(iua2ftrfil, '("  #         ", 15f12.7)') l_a2f_tr(:)
          WRITE(iua2ftrfil, *) "Phonon smearing (meV)"
          WRITE(iua2ftrfil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          WRITE(iua2ftrfil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
          WRITE(iua2ftrfil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
          WRITE(iua2ftrfil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
          CLOSE(iua2ftrfil)
          !
          CLOSE(iudosfil)
          !
        ENDDO ! isig
        !
      ENDDO ! itemp
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2f', 1)
      DEALLOCATE(l_a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2f_tr', 1)
      DEALLOCATE(a2f_, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2f', 1)
      DEALLOCATE(a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2f_tr', 1)
      DEALLOCATE(rho, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating rho', 1)
      DEALLOCATE(dosph, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating dosph', 1)
      DEALLOCATE(logavg, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating logavg', 1)
      !
    ENDIF
    !
    CALL stop_clock('a2F')
    CALL print_clock('a2F')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE a2f_main
    !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE spectral
  !-----------------------------------------------------------------------
