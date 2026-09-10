  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE scgd0_driver
  !---------------------------------------------------------------------  
  !!
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------- 
    SUBROUTINE scgd0_run()
    !---------------------------------------------------------------------  
    !! This is the main driver to compute electron spectral functions self-consistently
    !! based on Phys. Rev. Lett. 134 186401 (2025).
    !! The method is triggered with the flag 'specfun_el_scgd0' and requires also ephwrite = .true.. 
    !! The calculation needs to be performed on full fine commensurate k and q grids, 
    !! but the user can in addition specify a k-path file with 'filkf' and the output will 
    !! contain and interpolated spectral function on this path.
    !! The driver first reads the files created by 'ephwrite' and the self-energy components from the previous iteration.
    !! Then it starts the iterations. 
    !---------------------------------------------------------------------
    !  
    USE io_global,         ONLY : stdout, ionode
    USE io_supercond,      ONLY : read_frequencies, read_eigenvalues, read_kqmap,   &
                                  read_ephmat
    USE supercond,         ONLY : deallocate_eliashberg_elphon
    USE supercond_common,  ONLY : nkfs, ekfs, ef0, g2, nbndfs_all, nbndfs, wkfs_all, nkfs_all, ekfs_all
    USE kinds,             ONLY : DP
    USE selfen,            ONLY : selfen_elec_scgd0, re_selfen_scgd0,               &
                                  check_convergence
    USE spectral,          ONLY : spectral_func_el_print, spectral_recompute_efermi,&
                                  spectral_func_el_interpolate
    USE io_selfen,         ONLY : selfen_el_read, spectral_read,                    &
                                  selfen_el_read_wfpt, spectral_read_scgd0,         &
                                  spectral_write_scgd0
    USE mp,                ONLY : mp_sum, mp_max, mp_min
    USE mp_global,         ONLY : inter_pool_comm
    USE ep_constants,      ONLY : ryd2ev, zero, pi, ryd2mev
    USE global_var,        ONLY : wf, etf, nqf, vmef, efnew, nbndfst, totq,         &
                                  esigmar_all, esigmai_all, a_all, gtemp,    &
                                  a_all_ibnd, esigmaisc_all, iter_scgd0, mu_t,      &
                                  sigma_ahc_hdw, sigma_ahc_uf, sigmar_dw_all,     &
                                  s_bztoibz
    USE input,             ONLY : fsthick, nw_specfun, nstemp, restart, opt_cond,   &
                                  restart_step, carrier, lwfpt, filkf, wmin_specfun,&
                                  wmax_specfun, mp_mesh_k
    USE pwcom,             ONLY : ef
    USE io_var,            ONLY : iospectral
    USE modes,             ONLY : nmodes
    USE Utilities,         ONLY : fermi_dirac
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: iter_restart
    !! current scgd0 iteration
    INTEGER :: iq_restart
    !! Counter on coarse q-point grid
    INTEGER :: conv_counter
    !! convergence criterion. once it is 2, iterations stop
    INTEGER :: ibnd
    !! counter on bands
    INTEGER :: ik
    !! counter on k points
    INTEGER :: iw, iw_qp
    !!  counter on frequencies
    LOGICAL :: first_cycle
    !! Check wheter this is the first cycle after a restart.
    LOGICAL :: exst
    !! If the file exist
    REAL (KIND = DP) :: nelec_w
    !! number of electrons
    REAL(KIND = DP) :: ww(nw_specfun)
    !! frequency aray
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: relerr
    !! relative error
    REAL(KIND = DP) :: im_max
    !! maximal value of the imaginary part of the electron self-energy
    REAL(KIND = DP) :: re_max
    !! maximal value of the real part of the electron self-energy
    REAL(KIND = DP) :: ekk
    !! electron energy
    CHARACTER (LEN = 256) :: filespec
    !! name of the binary file
    !
    CALL start_clock('scgd0')
    !
    CALL read_frequencies()
    CALL read_eigenvalues()
    CALL read_kqmap()
    CALL read_ephmat()
    !
    ALLOCATE(esigmar_all(nbndfst, nkfs_all, nw_specfun, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating esigmar_all', 1)
    ALLOCATE(esigmai_all(nbndfst, nkfs_all, nw_specfun, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating esigmai_all', 1)
    esigmar_all(:, :, :, :) = zero
    esigmai_all(:, :, :, :) = zero
    ALLOCATE(mu_t(nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating mu_t', 1)
    mu_t(:) = zero
    ALLOCATE(a_all(nw_specfun, nkfs_all, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating a_all', 1)
    a_all(:, :, :) = zero
    ALLOCATE(esigmaisc_all(nbndfst, nkfs_all, nw_specfun, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating esigmaisc_all', 1)
    esigmaisc_all(:, :, :, :) = zero
    IF (opt_cond) THEN
      ALLOCATE(a_all_ibnd(nw_specfun, nkfs_all, nbndfst, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating a_all_ibnd', 1)
      a_all_ibnd(:, :, :, :) = zero
    ENDIF
    IF (lwfpt) THEN
      ALLOCATE(sigmar_dw_all(nbndfst, nkfs_all, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating sigmar_dw_all', 1)
      ALLOCATE(sigma_ahc_hdw(nbndfst, nkfs_all, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating sigma_ahc_hdw', 1)
      ALLOCATE(sigma_ahc_uf(nbndfst, nkfs_all, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating sigma_ahc_uf', 1)
      sigmar_dw_all = zero
      sigma_ahc_hdw = zero
      sigma_ahc_uf = zero
    ENDIF
    !
    iter_restart = 1
    !
    IF (restart) THEN
      CALL spectral_read(iq_restart, totq, nkfs_all, esigmar_all, esigmai_all)
      IF (totq - iq_restart < restart_step) THEN
        iq_restart = 1
        esigmai_all(:, :, :, :) = zero
        esigmar_all(:, :, :, :) = zero
      ELSE
        CALL mp_sum(esigmai_all, inter_pool_comm)
        esigmaisc_all(:, :, :, :) = esigmai_all(:, :, :, :)
        esigmai_all(:, :, :, :) = zero
      ENDIF
      CALL spectral_read_scgd0(nkfs_all, esigmar_all, esigmai_all, iter_restart, nelec_w, mu_t)
    ELSE
      CALL spectral_read_scgd0(nkfs_all, esigmar_all, esigmai_all, iter_restart, nelec_w, mu_t)
    ENDIF
    !
    ! ephwrite option by default saves everything in eV, so here we go back to Rydberg!
    ekfs = ekfs / ryd2ev
    ekfs_all = ekfs_all / ryd2ev
    g2 = g2 / ryd2ev**2
    wf = wf / ryd2ev
    ef0 = ef0 / ryd2ev
    ef = ef / ryd2ev
    fsthick = fsthick / ryd2ev
    iter_scgd0 = iter_restart
    conv_counter = 0
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    WRITE(stdout, '(/5x, a)') REPEAT('-', 67)
    WRITE(stdout, '(5x, a)') 'Electron spectral function will be computed using the self-consistent GD0 method.'
    WRITE(stdout, '(5x, a)') REPEAT('-', 67)
    WRITE(stdout, '(/5x, a)') 'Maximum number of iterations is 100 and a convergence threshold is 0.01.'
    WRITE(stdout, '(5x, a)') 'Final result will be printed out in the specfun.elself.TTT.TTTK_scGD0 file'
    WRITE(stdout, '(5x, a)') 'specfun_iter.bin file contains all the iterations. '
    DO WHILE (iter_scgd0 <= 100 .AND. conv_counter < 2 ) 
      !
      ! recompute the fermi level using the spectral function
      WRITE(stdout, '(/5x, a)') 'Fermi level is recomputed from the spectral function using bisection. '
      CALL spectral_recompute_efermi(nelec_w, mu_t, wkfs_all, ekfs_all, nbndfs_all, nkfs_all)
      esigmaisc_all(:, :, :, :) = zero
      CALL selfen_elec_scgd0()
      !
      CALL mp_sum(esigmaisc_all, inter_pool_comm)
      ! compare esigmai_all and esigmaisc_all
      CALL check_convergence(relerr)
      IF (relerr < 0.01d0) THEN
        conv_counter = conv_counter + 1
      ELSE
        conv_counter = 0
      ENDIF
      ! redefine esigmai_all
      esigmai_all(:, :, :, :) = esigmaisc_all(:, :, :, :)
      ! ! compute esigmar_all from esigmai_all
      CALL re_selfen_scgd0()   
      IF (restart) THEN
        CALL spectral_write_scgd0(totq, nkfs_all, esigmar_all, esigmaisc_all, nelec_w) 
        ! write down the full real part (DW+FM) and the imaginary part after each iteration
      ENDIF
      !
      IF (conv_counter > 1 ) THEN
        WRITE(stdout,'(/5x, a, i6)') 'Convergence achieved after iteration number:', iter_scgd0 
        WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
        IF (filkf /= '') THEN
           CALL spectral_func_el_interpolate()
        ELSE
          CALL spectral_func_el_print()
        ENDIF
      ELSE
        DO itemp = 1, nstemp
          DO iw = 1, nw_specfun
            DO ik = 1, nkfs
              DO ibnd = 1, nbndfst
                 a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
                 ((ww(iw) - ekfs(ibnd, ik) + mu_t(itemp) - esigmar_all(ibnd, ik, iw, itemp))**2 + &
                 (esigmai_all(ibnd, ik, iw, itemp))**2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        ! write iteration to binary file 
        IF (ionode) THEN
          filespec = 'specfun_iter.bin'
          OPEN(unit = iospectral, FILE = TRIM(filespec), access = 'stream', form = 'unformatted', &
          STATUS = 'unknown', POSITION = 'append')
          WRITE(iospectral) iter_scgd0
          WRITE(iospectral) a_all
          CLOSE(iospectral)
        ENDIF
        ! initialize
        a_all(:, :, :) = zero
      ENDIF
      !
      iter_scgd0 = iter_scgd0 + 1
    ENDDO
    IF (conv_counter <= 1) THEN
      WRITE(stdout,'(/5x, a)') 'Maximum iterations number (100) reached without convergence'
    ENDIF
    !
    DEALLOCATE(esigmaisc_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating esigmaisc_all', 1)
    DEALLOCATE(a_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating a_all', 1)
    DEALLOCATE(esigmai_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating esigmai_all', 1)
    DEALLOCATE(esigmar_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating esigmar_all', 1)
    IF (opt_cond) THEN
      DEALLOCATE(a_all_ibnd, STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating a_all_ibnd', 1)
      DEALLOCATE(vmef, STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating vmef', 1)
      IF (mp_mesh_k) THEN 
        DEALLOCATE(s_bztoibz, STAT = ierr)
        IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating s_bztoibz', 1)
      ENDIF
    ENDIF
    IF (lwfpt) THEN
      DEALLOCATE(sigmar_dw_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating sigmar_dw_all', 1)
      DEALLOCATE(sigma_ahc_hdw, STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating sigmar_ahc_hdw', 1)
      DEALLOCATE(sigma_ahc_uf, STAT = ierr)
      IF (ierr /= 0) CALL errore('scgd0_run', 'Error deallocating sigma_ahc_uf', 1)
    ENDIF
    DEALLOCATE(mu_t, STAT = ierr)
    IF (ierr /= 0) CALL errore('scgd0_run', 'Error allocating mu_t', 1)
    CALL deallocate_eliashberg_elphon()
    CALL stop_clock('scgd0')
    !
    RETURN
    !---------------------------------------------------------------------  
    END SUBROUTINE scgd0_run
    !---------------------------------------------------------------------   

  !-----------------------------------------------------------------------
  END MODULE scgd0_driver
  !-----------------------------------------------------------------------


