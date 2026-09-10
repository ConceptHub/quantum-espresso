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
  MODULE polaron_scf_driver
  !--------------------------------------------------------------------------
  !!
  !! Drives the polaron SCF cycle: initial guess, build/diagonalize/update loop,
  !! convergence, formation energies, and band-edge helpers.
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : is_mirror_q, nbnd_plrn, band_pos, ik_edge, nRp, Rp, &
                             select_bands_plrn, etf_all, xkf_all, eigvec, Bmat, berry_phase
  USE polaron_grid,           ONLY : index_shift, ikqLocal2Global, isGVec
  USE polaron_hamiltonian,    ONLY : build_plrn_bmat, build_plrn_hamil, norm_plrn_wf
  USE polaron_diagonalization, ONLY : diag_serial, diag_parallel
#if defined(__ELPA)
  USE polaron_diagonalization, ONLY : diag_elpa
#endif
  USE io_polaron,     ONLY : read_plrn_wf, read_plrn_dtau, read_Rp_in_S, calc_den_of_state, &
                             write_plrn_wf, write_plrn_bmat, write_plrn_dtau_xsf,            &
                             scell_write_plrn_dtau_xsf
  USE polaron_interpolation, ONLY : check_time_rev_sym, plrn_eigvec_tran, scell_plrn_eigvec_tran, &
                             plrn_bmat_tran, scell_plrn_bmat_tran

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: polaron_scf, find_band_extreme, gather_band_eigenvalues

  CONTAINS

    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_random(eigvec_init)
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nktotf
    USE ep_constants,  ONLY : ci, cone
    USE input,         ONLY : nstate_plrn
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_init(:, :)
    !! Polaron wf coefficients, Ank, upon initialization
    !
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: rmat_tmp(:, :)
    !! Temporary variable for random number matrix
    !
    CALL RANDOM_SEED()
    ALLOCATE(rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error allocating rmat_tmp', 1)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigvec_init(1:nktotf * nbnd_plrn, 1:nstate_plrn) = cone * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigvec_init(1:nktotf * nbnd_plrn, 1:nstate_plrn) = eigvec_init(1:nktotf*nbnd_plrn, 1:nstate_plrn) + &
                                                ci * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    DEALLOCATE(rmat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error deallocating rmat_tmp', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_random
    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_gaussian(r0, k_all, k0, eigvec_init)
    !-----------------------------------------------------------------------
    !! Initialize Ank coefficients with a Gaussian lineshape
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : cone, ci, twopi, one, zero
    USE input,         ONLY : init_sigma_plrn
    USE global_var,    ONLY : nktotf
    USE cell_base,     ONLY : bg, alat
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: r0(3)
    !! Center of Gaussian in real space
    REAL(KIND = DP), INTENT(in) :: k0(3)
    !! Center of Gaussian in k-space
    REAL(KIND = DP), INTENT(in) :: k_all(:, :)
    !! List with all k-point coordinates
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_init(:, :)
    !! Polaron wf coefficients, Ank, upon initialization
    !
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: indexkn1
    !! Combined k-point and band counter
    INTEGER :: ishift
    !! Shift counter
    REAL(KIND = DP) :: qcart(3)
    !! q-point coordinates in cartesian
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinates
    REAL(KIND = DP) :: shift(3)
    !! shift coordinates
    REAL(KIND = DP) :: disK
    !! Value of Gaussian in neighbor BZs
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor from Fourier transform of Gaussian
    !
    DO ik = 1, nktotf
      xxq = k_all(1:3, ik) - (k0(:) - INT(k0(:))) ! shift k0 to 1BZ
      CALL dgemv('n', 3, 3, one, bg, 3, xxq, 1, zero, qcart, 1)
      ctemp = EXP(-ci * twopi * DOT_PRODUCT(qcart, r0))
      disK = -1
      ! Ensure periodicity of Ank checking distance to other equivalent BZ points
      DO ishift = 1, 27
        shift(1:3) = REAL(index_shift(ishift), KIND = DP)
        CALL dgemv('n', 3, 3, one, bg, 3, xxq + shift, 1, zero, qcart, 1)
        disK = MAX(disK, EXP(-init_sigma_plrn * NORM2(qcart) * twopi / alat)) ! for sigma to be in bohr
      ENDDO
      !
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        eigvec_init(indexkn1, :) = cone * disK * ctemp
      ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_gaussian
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE polaron_scf (nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
    !-----------------------------------------------------------------------
    !!
    !! Self consistency calculation of polaron wavefunction.
    !! Rewritten by Chao Lian based on the implementation by Denny Sio.
    !! SP: cleaning (Nov 2023)
    !!
    !
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : ryd2mev, one, ryd2ev, two, zero, twopi,           &
                              czero, cone, pi, ci, twopi, eps6, eps8, eps5
    USE input,         ONLY : type_plrn, full_diagon_plrn, debug_plrn,          &
                              init_sigma_plrn, init_k0_plrn, nstate_plrn,       &
                              conv_thr_plrn, init_plrn, niter_plrn,             &
                              nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, r0_plrn,      &
                              init_ntau_plrn, nbndsub, as, time_rev_A_plrn,     &
                              model_phfreq_plrn, omega_LO_plrn, scell_mat_plrn, &
                              acoustic_plrn, cal_acous_plrn, dtau_max_plrn,     &
                              eigen_solver_plrn, istate_relax_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id, ionode_id
    USE io_var,        ONLY : iufileigplrn
    USE global_var,    ONLY : nqtotf, nktotf, wf
    USE cell_base,     ONLY : alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_global,     ONLY : my_pool_id
    USE parallelism,   ONLY : poolgather2
    USE mp_world,      ONLY : world_comm
    USE input,         ONLY : ethrdg_plrn
    USE control_flags, ONLY : use_gpu
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) ::  dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    !
    ! local variables
    CHARACTER(LEN = 256) :: filename
    !! Output file name
    CHARACTER(LEN = 256) :: tmpch
    !! Temporary character to assign name to different displacement files
    LOGICAL :: debug
    !! .true. if extra output is to be printed for debugging
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! band counter
    INTEGER :: ierr
    !! Error status
    INTEGER :: itau
    !! Atom index
    INTEGER :: iplrn
    !! Polaron index
    INTEGER :: iter
    !! Iteration counter
    INTEGER :: indexkn1
    !! Combined k-point and band index
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nbndsub_p
    !! Number of bands in polaron calculation
    INTEGER :: nktotf_p
    !! Number of k-points in the fine grid
    !REAL(KIND = DP) :: estmteRt(nstate_plrn)
    !! Polaron eigenvalue in diagonalization
    !REAL(KIND = DP) :: eigval(nstate_plrn)
    REAL(KIND = DP), ALLOCATABLE :: eigval(:)
    !! Polaron eigenvalue
    REAL(KIND = DP) :: esterr
    !! Difference of displacements between scf loops
    REAL(KIND = DP) :: eplrnelec
    !! Electron part of polaron formation energy
    REAL(KIND = DP) :: eplrnphon
    !! Phonon part of polaron formation energy
    REAL(KIND = DP) :: r_cry(1:3)
    !! Polaron center
    REAL(KIND = DP) :: dtau_diff
    !! Max difference between displacements between scf loops
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_acoustic(:, :)
    !! Polaron displacements by acoustic phonon modes
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_save(:, :)
    !! Polaron displacements from previous iteration to compare
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_list(:, :, :)
    !! List of displacements from which scf is to be initiated
    !
    ALLOCATE(dtau(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau', 1)
    ALLOCATE(dtau_save(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_save', 1)
    dtau = czero
    dtau_save = czero
    debug = debug_plrn
    IF (cal_acous_plrn) THEN
      ALLOCATE(dtau_acoustic(nktotf, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_acoustic', 1)
      dtau_acoustic = czero
    ENDIF
    !
    CALL start_clock('main_prln')
    ! Gather all the eigenvalues to determine the EBM/VBM,
    CALL start_clock('re_omega')
    ! Recalculate the frequency, when restart from save g
    CALL cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wf)
    !
    IF(model_phfreq_plrn) THEN
      wf = zero
      wf(nmodes, :) = omega_LO_plrn
    ENDIF
    !
    CALL stop_clock('re_omega')
    !
    !! Initialize Ac(k) based on profile
    !! TODO: ik_bm should be user-adjustable
    CALL start_clock('init_Ank')
    eigvec = czero
    SELECT CASE (init_plrn)
      CASE (1)
        ! If k0 has not been set on input, center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        !
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
           &packet with a width of", ES14.6)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f14.6)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigvec)
      CASE (3)
        ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
        WRITE(stdout, '(5x, a)') "Initializing the polaron wavefunction with previously saved Amp.plrn file"
        CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nstate_plrn, 'Amp.plrn')
        CALL plrn_eigvec_tran('Wan2Bloch', time_rev_A_plrn, eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
           nrr_k, ndegen_k, irvec_r, dims, eigvec)
        DEALLOCATE(eigvec_wan, STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
      CASE (6)
        WRITE(stdout, '(5x, a, I6)') "Starting from displacements read from file; number of displacements:", init_ntau_plrn
        ALLOCATE(dtau_list(init_ntau_plrn, nktotf, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_list', 1)
        dtau_list = CMPLX(0.d0, 0.d0)
        !
        IF (init_ntau_plrn == 1) THEN
          filename = 'dtau_disp.plrn'
          CALL read_plrn_dtau(dtau, nqtotf, nmodes, filename, scell_mat_plrn)
          dtau_list(1, :, :) = dtau(:, :)
        ELSE
          DO itau = 1, init_ntau_plrn
            WRITE(tmpch,'(I4)') itau
            filename = TRIM('dtau_disp.plrn_'//ADJUSTL(tmpch))
            CALL read_plrn_dtau(dtau, nqtotf, nmodes, filename, scell_mat_plrn)
            dtau_list(itau, :, :) = dtau(:, :)
          ENDDO
        ENDIF
        !
        CALL mp_bcast(dtau_list, meta_ionode_id, world_comm)
        ! Initialize Ank wavefunction for iterative diagonalization Gaussian
        ! If k0 has not been set on input,     center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
        &packet with the width of", f15.7)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f15.7)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigvec)
        CALL norm_plrn_wf(eigvec, REAL(nktotf, DP))
      CASE DEFAULT
        CALL errore('polaron_scf','init_plrn not implemented!', 1)
    END SELECT
    !
    ! Only keep the coefficients in lowest/highest band,
    ! since the electron/hole localized at this band will be more stable.
    IF (init_plrn <= 2) THEN
      DO ik = 1, nktotf
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          IF (select_bands_plrn(ibnd) /= band_pos) eigvec(indexkn1, 1:nstate_plrn) = czero
        ENDDO
      ENDDO
      CALL norm_plrn_wf(eigvec, REAL(nktotf, DP))
    ENDIF
    CALL stop_clock('init_Ank')
    !
    IF ( istate_relax_plrn .NE. 1 ) THEN
      WRITE(stdout, '(5x, "Warning: the relaxed polaron is an excited state since istate_relax_plrn != 1")')
      WRITE(stdout, '(5x, a, i8)') 'istate_relax_plrn = ', istate_relax_plrn
    ENDIF
    WRITE(stdout, '(5x, "Starting the SCF cycles")')
    IF (full_diagon_plrn) THEN
      WRITE(stdout, '(5x, a)') "Using serial direct diagonalization"
    ELSE
      WRITE(stdout, '(5x, a)') "Using parallel iterative diagonalization"
      WRITE(stdout, '(5x, "Diagonalizing polaron Hamiltonian with a threshold of ", ES18.6)') ethrdg_plrn
      WRITE(stdout, '(5x, "Please check the results are convergent with this value")')
    ENDIF
    !
#if defined(__ELPA)
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
    WRITE(stdout,'(5x,a)')  'Full diagonalization can be accelerated with ELPA Library                '
    WRITE(stdout,'(5x,a)')  'Refer: https://elpa.mpcdf.mpg.de/ELPA_PUBLICATIONS.html for citation     '
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
#else
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
    WRITE(stdout,'(5x,a)')  'ELPA Library not found                                                   '
    WRITE(stdout,'(5x,a)')  'Please install EPW with ELPA for faster and denser calculations          '
    WRITE(stdout,'(5x,a)')  'Visit: https://elpa.mpcdf.mpg.de for installing ELPA                     '
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
#endif
    !
    WRITE(stdout, '(/5x, a)') "Starting the self-consistent process"
    WRITE(stdout, '( 5x, a)') REPEAT('-',80)
    WRITE(stdout, '(5x, " iter", 60a15)') "  Eigval/eV", "Phonon/eV", "Electron/eV", &
                                          "Formation/eV", "Error/eV"
    ALLOCATE(Bmat(nqtotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating Bmat', 1)
    !
    IF (scell_mat_plrn) THEN
      CALL read_Rp_in_S()
    ENDIF
    !
    ALLOCATE(eigval(nktotf*nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigval', 1)
    !JLB: possibility of multiple displacements read from file, to calculate polaron energy landscape.
    !     Calculate and print the energies; .plrn files written to disk for last calculation only.
    DO itau = 1, init_ntau_plrn ! ntau_plrn=1 by default
      !
      IF (init_plrn == 6) dtau(:, :) = dtau_list(itau, :, :)
      !
      !estmteRt = 1E3
      eigval(:) = 1E3
      esterr = 1E5
      DO iter = 1, niter_plrn
        ! Enforce the relation A_k = A*_{G-k} and normalize |A| = 1
        ! Calculating $$ B_{\bq\nu} = \frac{1}{\omega_{\bq,\nu} N_p} \sum_\bk A^\dagger_{\bk+\bq} g_\nu(\bk,\bq) A_\bk $$
        CALL start_clock('cal_bqu')
        IF (init_plrn == 6 .AND. iter==1) THEN
          Bmat = czero
          IF (scell_mat_plrn) THEN
            CALL scell_plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
          ELSE
            CALL plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
          ENDIF
        ELSE
          CALL build_plrn_bmat(Bmat)
          !
          IF (scell_mat_plrn) THEN
            CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
          ELSE
            CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
          ENDIF
          !
        ENDIF
        !
        dtau_diff = MAXVAL(ABS(REAL(dtau - dtau_save)))
        esterr = dtau_diff
        IF (dtau_diff < conv_thr_plrn .AND. iter > 1) THEN
          !! IF(MAXVAL(ABS(REAL(dtau))) > alat / 2.d0) THEN
          !! KL: dtau_max_plrn is a criteria of maximual polaron displacements
          !! of a physical solution, which is 0.5d0 by default.
          !! When acoustic phonon modes are found to be important, one can try to increase.
          IF(MAXVAL(ABS(REAL(dtau))) > alat * dtau_max_plrn) THEN
            CALL errore("polaron_scf","Non-physical solution, check initial guess and convergence.", 1)
          ENDIF
          ! converged, write the final value of eigenvalue
          WRITE(stdout,'(5x,a)') REPEAT('-',80)
          WRITE(stdout, '(5x,a,f10.6,a)' )  'End of self-consistent cycle'
          EXIT
        ELSE
          dtau_save = dtau
        ENDIF
        !
        CALL stop_clock('cal_bqu')
        !
        CALL start_clock('Setup_H')
        !
        CALL build_plrn_hamil(Bmat)
        CALL stop_clock('Setup_H')
        CALL start_clock('DiagonH')
        ! For hole polaron (type_plrn = 1),
        ! we need the highest eigenvalues instead of the lowest eigenvalues
        ! To use KS_solver, which only gives the lowest eigenvalues,
        ! we multiply -1 to the Hamiltonian to get the lowest eigenvalues
        IF (full_diagon_plrn .OR. use_gpu) THEN
          IF (eigen_solver_plrn == 'lapack') THEN
            ! Diagonalize Hamiltonian with Serial LAPACK subroutine
            ! Used for testing or robust benchmark
            !CALL diag_serial(estmteRt, eigvec)
            CALL diag_serial(eigval, eigvec)
#if defined(__ELPA)
          ELSEIF (eigen_solver_plrn == 'elpa') THEN
            !CALL diag_elpa(estmteRt, eigvec)
            CALL diag_elpa(eigval, eigvec)
#endif
          ENDIF
          !
          !CALL mp_bcast(estmteRt, meta_ionode_id, world_comm)
          CALL mp_bcast(eigval, meta_ionode_id, world_comm)
          CALL mp_bcast(eigvec, meta_ionode_id, world_comm)
        ELSE
          ! Diagonalize Hamiltonian with Davidson Solver
          !CALL diag_parallel(estmteRt, eigvec)
          CALL diag_parallel(eigval, eigvec)
        ENDIF
        !
        CALL stop_clock('DiagonH')
        !
        ! Reverse the eigenvalues if it is the hole polaron
        !estmteRt(1:nstate_plrn) = (-type_plrn) * estmteRt(1:nstate_plrn)
        eigval(:) = (-type_plrn) * eigval(:)
        !
        ! impose the time-reversal symmetry: A^T_k = A_k + A^*_{-k}
        IF (time_rev_A_plrn) CALL check_time_rev_sym(eigvec)
        CALL norm_plrn_wf(eigvec, REAL(nktotf, dp))
        !
        eigvec = (- type_plrn) * eigvec
        !
        CALL start_clock('cal_E_Form')
        CALL calc_form_energy(eplrnphon, eplrnelec)
        CALL stop_clock('cal_E_Form')
        !
        ! TODO : use exact number instead of 20 in 20e15.7
        r_cry(1:3) = IMAG(LOG(berry_phase(1:3) * EXP(- twopi * ci * r0_plrn(1:3)))) / twopi
        r_cry(1:3) = r_cry(1:3) - NINT(r_cry(1:3))
        !WRITE(stdout, '(5x, i5, 60e15.4)') iter, estmteRt(1:nstate_plrn) * ryd2ev, eplrnphon * ryd2ev, &
        !   - eplrnelec * ryd2ev, (eplrnelec + eplrnphon) * ryd2ev, esterr
        !! We only print the first eigenvalue during the iteration
        !WRITE(stdout, '(5x, i5, 60e15.4)') iter, eigval(1) * ryd2ev, eplrnphon * ryd2ev, &
        !   - eplrnelec * ryd2ev, (eplrnelec + eplrnphon) * ryd2ev, esterr
        !! the selected excited-state polron eigenvalue is printed
        WRITE(stdout, '(5x, i5, 60e15.4)') iter, eigval(istate_relax_plrn) * ryd2ev, eplrnphon * ryd2ev, &
           - eplrnelec * ryd2ev, (eplrnelec + eplrnphon) * ryd2ev, esterr
        !eigval = estmteRt(1:nstate_plrn)
      ENDDO
      !
      ! Calculate and write the energies
      IF (nstate_plrn .EQ. 1) THEN
          !WRITE(stdout, '(5x, a, 50f16.7)') '      Eigenvalue (eV): ', eigval(1) * ryd2ev
          WRITE(stdout, '(5x, a, 50f16.7)') '      Eigenvalue (eV): ', eigval(istate_relax_plrn) * ryd2ev
      ELSE
          WRITE(stdout, '(5x, a)') '      Eigenvalue (eV): '
          DO iplrn = 1, nstate_plrn, 5
              WRITE(stdout, '(5x, a, 5f16.7)') '      ', eigval(iplrn: MIN(iplrn+4, nstate_plrn)) * ryd2ev
          ENDDO
      ENDIF
      WRITE(stdout, '(5x, a, f16.7)')   '     Phonon part (eV): ', eplrnphon * ryd2ev
      WRITE(stdout, '(5x, a, f16.7)')   '   Electron part (eV): ', eplrnelec * ryd2ev
      IF (init_plrn == 6) THEN
        !WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', ((-type_plrn) * eigval - eplrnphon) * ryd2ev
        !WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', ((-type_plrn) * eigval(1) - eplrnphon) * ryd2ev
        WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', &
                ((-type_plrn) * eigval(istate_relax_plrn) - eplrnphon) * ryd2ev
      ELSE
        WRITE(stdout, '(5x, a, f16.7)')   'Formation Energy (eV): ', (eplrnelec + eplrnphon) * ryd2ev
      ENDIF
    ENDDO ! init_ntau_plrn
    !
    ! Save all eigenvalues
    IF ((full_diagon_plrn) .AND. (eigen_solver_plrn == 'elpa') .AND. (my_pool_id == ionode_id)) THEN
      OPEN(UNIT = iufileigplrn, FILE = 'eigenvalues.elpa.plrn')
      WRITE(iufileigplrn, '(5x, a, 2i12)') 'polaron eigenvalues (eV). ', nktotf, nbnd_plrn
      DO iplrn = 1, nktotf*nbnd_plrn, 5
        WRITE(iufileigplrn, '(5x, 5f16.7)') eigval(iplrn:MIN(iplrn+4, nktotf*nbnd_plrn)) * ryd2ev
      ENDDO
    ENDIF
    !
    ! Calculate and write Density of State of Bqnu and Ank
    WRITE(stdout, '(5x, a)') "Calculating density of states to save in dos.plrn"
    CALL calc_den_of_state(eigvec, Bmat)
    !
    ! Do Bloch to Wannier transform, with U matrix
    CALL start_clock('Ank_trans')
    WRITE(stdout, '(5x, a)') "Generating the polaron wavefunction in Wannier basis to save in Amp.plrn"
    !
    ALLOCATE(eigvec_wan(nbndsub * nktotf, nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
    eigvec_wan = czero
    IF (scell_mat_plrn) THEN
      ! JLB: t_rev set to .true. before
      CALL scell_plrn_eigvec_tran('Bloch2Wan', time_rev_A_plrn, eigvec, nktotf, nRp, Rp, nbndsub, nrr_k, &
              ndegen_k, irvec_r, dims, eigvec_wan)
    ELSE
      ! JLB: t_rev set to .true. before
      CALL plrn_eigvec_tran('Bloch2Wan', time_rev_A_plrn, eigvec, nkf1, nkf2, nkf3, nbndsub, nrr_k, &
         ndegen_k, irvec_r, dims, eigvec_wan)
    ENDIF
    CALL stop_clock('Ank_trans')
    !
    ! Calculate displacements of ions dtau, which is B matrix in Wannier basis
    CALL start_clock('Bqu_tran')
    dtau = czero
    WRITE(stdout, '(5x, a)') "Generating the ionic displacements to save in dtau.plrn and dtau.plrn.xsf"
    IF (scell_mat_plrn) THEN
      CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
    ELSE
      CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
      IF (cal_acous_plrn) THEN
        CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, &
                            dtau_acoustic, acoustic_plrn=acoustic_plrn)
      ENDIF
    ENDIF
    CALL stop_clock('Bqu_tran')
    CALL start_clock('write_files')
    !
    IF (ionode) THEN
      ! Write Amp in Wannier basis
      CALL write_plrn_wf(eigvec_wan, 'Amp.plrn')
      ! Write Ank in Bloch basis
      CALL write_plrn_wf(eigvec, 'Ank.plrn',  etf_all)
      ! Write Bqnu
      CALL write_plrn_bmat(Bmat, 'Bmat.plrn', wf)
      ! Write dtau
      CALL write_plrn_bmat(dtau, 'dtau.plrn')
      IF (cal_acous_plrn) THEN
        CALL write_plrn_bmat(dtau_acoustic, 'dtau.acoustic.plrn')
      ENDIF
      ! Write dtau in a user-friendly format for visulization
      IF (scell_mat_plrn) THEN
        CALL scell_write_plrn_dtau_xsf(dtau, nqtotf, nRp, Rp, as, 'dtau.plrn.xsf')
      ELSE
        CALL write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, 'dtau.plrn.xsf')
      IF (cal_acous_plrn) THEN
        CALL write_plrn_dtau_xsf(dtau_acoustic, nqf1, nqf2, nqf3, 'dtau.acoustic.plrn.xsf')
      ENDIF
      ENDIF
    ENDIF
    CALL stop_clock('write_files')
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau', 1)
    DEALLOCATE(dtau_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau', 1)
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
    DEALLOCATE(Bmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating Bmat', 1)
    IF (cal_acous_plrn) THEN
      DEALLOCATE(dtau_acoustic, STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau_acoustic', 1)
    ENDIF
    DEALLOCATE(eigval, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigval', 1)
    !
    CALL stop_clock('main_prln')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE polaron_scf
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE calc_form_energy(eplrnphon, eplrnelec)
    !-----------------------------------------------------------------------
    !!
    !! Computes the polaron formation energy
    !! Note: Require etf_all to be properly initialized
    !!
    !-----------------------------------------------------------------------
    USE ep_constants,    ONLY : zero, czero, twopi, ci, two
    USE modes,           ONLY : nmodes
    USE global_var,      ONLY : xqf, wf, nqtotf, nktotf, nkf
    USE input,           ONLY : type_plrn, nstate_plrn, istate_relax_plrn
    USE mp,              ONLY : mp_sum
    USE mp_global,       ONLY : inter_pool_comm
    USE control_flags,   ONLY : iverbosity
    USE io_global,       ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: eplrnphon
    !! Phonon part of polaron formation energy
    REAL(KIND = DP), INTENT(out) :: eplrnelec
    !! Electron part of polaron formation energy
    !
    ! Local variables
    INTEGER :: iq
    !! q point index
    INTEGER :: ik
    !! k point index
    INTEGER :: ik_global
    !! k point global index
    INTEGER :: iq_global
    !! q point global index
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! Phonon mode index
    INTEGER :: ibnd
    !! Electron band index
    !INTEGER :: iplrn
    !! Polaron state index
    INTEGER :: indexkn1
    !! Combined k-point and band index
    !
    ! Based on Eq. 41 of Ref. 2:
    ! E_{f,ph} = 1/N_p \sum_{q\nu}|B_{q\nu}|^2\hbar\omega_{q\nu}
    ! iq -> q, nqtotf -> N_p
    ! Bmat(iq, inu) -> B_{q\nu}
    ! wf(inu, iq) -> \hbar\omega_{q\nu}
    eplrnphon = zero
    DO iq = 1, nkf
      iq_global = ikqLocal2Global(iq, nqtotf)
      ! JLB - Swapped indices!
      ! I think it would be better to discard modes by looking at the frequencies, i.e. discard negative or zero frequency modes.
      IF(isGVec(xqf(1:3, iq_global))) THEN
        start_mode = 4
      ELSE
        start_mode = 1
      ENDIF
      DO inu = start_mode, nmodes
        eplrnphon = eplrnphon - ABS(Bmat(iq_global, inu))**2 * (wf(inu, iq_global) / nqtotf)
      ENDDO
    ENDDO
    CALL mp_sum(eplrnphon, inter_pool_comm)
    !
    ! E_{f,el} = 1/N_p \sum_{nk}|A_{nk}|^2(\epsilon_{nk}-\epsilon_{F})
    ! indexkn1 -> nk, nktotf -> N_p
    ! eigvec(indexkn1, iplrn) -> A_{nk}
    ! etf_all(select_bands_plrn(ibnd), ik) - ef -> \epsilon_{nk}-\epsilon_{F}
    IF (iverbosity == 5) WRITE(stdout, '(5x, a, i8)') 'istate_relax_plrn = ', istate_relax_plrn
    eplrnelec = zero
    ! TODO: what should we do in iplrn
    !DO iplrn = 1, nstate_plrn
    !iplrn = 1 !! KL: only the first (ground state) polaron is counted
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nqtotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik_global - 1) * nbnd_plrn + ibnd
          !eplrnelec = eplrnelec - type_plrn * ABS(eigvec(indexkn1, iplrn))**2 / nktotf *&
          !   etf_all(select_bands_plrn(ibnd), ik_global)
          eplrnelec = eplrnelec - type_plrn * ABS(eigvec(indexkn1, istate_relax_plrn))**2 &
                 / nktotf * etf_all(select_bands_plrn(ibnd), ik_global)
        ENDDO
      ENDDO
    !ENDDO
    CALL mp_sum(eplrnelec, inter_pool_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calc_form_energy
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE find_band_extreme(type_plrn, enk_all, ik_bm, band_loc, efermi)
    !-----------------------------------------------------------------------
    !!
    !! Determine the Fermi energy, read from the input or calculated from band structure
    !!
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, ryd2ev
    USE input,         ONLY : efermi_read, fermi_energy
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : nkf, nktotf
    USE mp,            ONLY : mp_max, mp_min, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: type_plrn
    !! Whether electron polaron (-1) or hole polaron (+1) is to be calculated
    INTEGER, INTENT(out) :: ik_bm
    !! k-point index for CBM or VBM
    INTEGER, INTENT(out) :: band_loc
    !! band index for CBM or VBM
    REAL(KIND = DP), INTENT(in) :: enk_all(:,:)
    !! KS eigenvalues
    REAL(KIND = DP), INTENT(out) :: efermi
    !! Fermi energy
    !
    ! Local variable
    INTEGER :: ik
    !! k-point index
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: k_extreme_local(npool)
    !! Auxiliary index of CBM or VBM k-point at different pools
    INTEGER :: ipool(1)
    !! Index to locate CBM or VBM
    REAL(KIND = DP) :: band_edge
    !! Temporary KS eigenvalue
    REAL(KIND = DP) :: extreme_local(npool)
    !! Auxiliary CBM or VBM energy at different pools
    !
    IF (type_plrn == -1 ) THEN
      band_loc = select_bands_plrn(1)
    ELSE IF ( type_plrn == 1 ) THEN
      band_loc = select_bands_plrn(nbnd_plrn)
    ENDIF
    !
    WRITE(stdout, '(5x, "The band extremes are at band ",  i0)') band_loc
    !
    ! Determine the Fermi energy, read from the input or calculated from band structure
    ! = 1E4*(-type_plrn)
    ik_bm = 0
    k_extreme_local = 0
    extreme_local = zero
    IF (efermi_read) THEN
      efermi = fermi_energy
      WRITE(stdout, '(5x, "Polaron Reference energy (VBM or CBM) is read from the input file: ",&
         &f16.6, " eV.")') efermi * ryd2ev
    ELSE
      IF(type_plrn == 1) THEN
        efermi = -1E5
      ELSE IF (type_plrn == -1) THEN
        efermi =  1E5
      ELSE
        CALL errore('','Wrong type_plrn, should be 1 or -1', 1)
      ENDIF
      !
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        band_edge = enk_all(band_loc, ik_global)
        !
        IF (type_plrn == 1) THEN
          ! For hole polaron (type_plrn = 1), find the highest eigenvalue
          IF (band_edge > efermi) THEN
            efermi = band_edge
            ik_bm = ik_global
          ENDIF
        ELSE IF (type_plrn == -1) THEN
          ! For electron polaron (type_plrn = -1), find the lowest eigenvalue
          IF (band_edge < efermi) THEN
            efermi = band_edge
            ik_bm = ik_global
          ENDIF
        ENDIF
      ENDDO
      !
      k_extreme_local(my_pool_id + 1) = ik_bm
      extreme_local(my_pool_id + 1) = efermi
      CALL mp_sum(k_extreme_local, inter_pool_comm)
      CALL mp_sum(extreme_local, inter_pool_comm)
      !
      IF (type_plrn == 1) THEN
        ipool = MAXLOC(extreme_local)
        ik_bm = k_extreme_local(ipool(1))
        efermi = MAXVAL(extreme_local)
      ELSE IF (type_plrn == -1) THEN
        ipool = MINLOC(extreme_local)
        ik_bm = k_extreme_local(ipool(1))
        efermi = MINVAL(extreme_local)
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE find_band_extreme
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE gather_band_eigenvalues(etf, enk_all)
    !-----------------------------------------------------------------------
    !! Gather all the eigenvalues to determine the EBM/VBM,
    !! and calculate the density state of Ank and Bqnu
    !-----------------------------------------------------------------------
    USE input,         ONLY : nbndsub
    USE global_var,    ONLY : nkqf, nktotf
    USE ep_constants,  ONLY : zero
    USE parallelism,  ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etf(:, :)
    !! Eigenvalues per cpu
    REAL(KIND = DP), INTENT(out) :: enk_all(:, :)
    !! Eigenvalues (total)
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:, :)
    !! Temporary variable to gather eigenvalues across pools
    !
    ALLOCATE(rtmp2(nbndsub, nktotf*2), STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error allocating rtmp2', 1)
    rtmp2 = zero
    !
    CALL poolgather2 ( nbndsub, nktotf*2, nkqf, etf, rtmp2  )
    enk_all(1:nbndsub, 1:nktotf) = rtmp2(1:nbndsub, 1:nktotf*2:2)
    !
    DEALLOCATE(rtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error deallocating rtmp2', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gather_band_eigenvalues
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wfreq)
    !-----------------------------------------------------------------------
    !! Calculate the phonon eigen frequencies. This is needed when restarting the polaron
    !! calculation with recalculating el-ph vertex
    !-----------------------------------------------------------------------
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : xqf, nqtotf
    USE ep_constants,  ONLY : zero, eps8, czero
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc
    USE io_global,     ONLY : ionode, stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    REAL(KIND = DP), INTENT(out) :: wfreq(:, :)
    !! Phonon frequencies
    !
    ! Local variables
    INTEGER  :: inu
    !! Phonon mode counter
    INTEGER  :: iq
    !! q-point counter
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinates
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    uf = czero
    w2 = zero
    !
    !TODO: make this part parallel over q
    wfreq = zero
    DO iq = 1, nqtotf
      ! iq -> q
      xxq = xqf(1:3, iq)
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2, is_mirror_q(iq))
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2, is_mirror_q(iq))
      ENDIF
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wfreq(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          IF (ionode) THEN
            WRITE(stdout, '(5x, "WARNING: Imaginary frequency mode ",&
            &I6, " at iq=", I6)') inu, iq
          ENDIF
          wfreq(inu, iq) = 0.d0
        ENDIF
      ENDDO ! inu
    ENDDO ! iq
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_phonon_eigenfreq

  END MODULE polaron_scf_driver
