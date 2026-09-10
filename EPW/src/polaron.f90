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
  !-----------------------------------------------------------------------
  MODULE polaron
  !-----------------------------------------------------------------------
  !! This module contains the main drivers and entry points for polaron.
  !! It imports and re-exports subroutines from modular polaron sub-files.
  !!
  !! The refactored modules fit into the EPW flow as follows:
  !!
  !! epw.x
  !! └── use_wannier
  !!     └── polaron                              facade: prepare, dispatch, close
  !!         ├── q-point loop → io_polaron         save and collect g
  !!         └── plrn_flow_select                  SCF or post-processing, never both
  !!             ├── polaron_scf_driver            run the SCF calculation
  !!             │   ├── polaron_hamiltonian       build and apply the Hamiltonian
  !!             │   ├── polaron_diagonalization   solve the Hamiltonian
  !!             │   └── io_polaron                write results
  !!             ├── polaron_interpolation         Ank/Bqu post-processing
  !!             │   └── io_polaron                read stored, write interpolated
  !!             └── io_polaron                    real-space output
  !!
  !! Supporting modules:
  !!   polaron_common   shared calculation state
  !!   polaron_grid     grids, mappings, and indices
  !!
  !! Authored by Chao Lian, Weng Hong (Denny) Sio, and Jon Lafuente-Bartolome
  !! Partial cleaning by SP (Nov 2023)
  !! Partial cleaning by STiwari (Nov 2023)
  !! Partial cleaning by JLB (Aug 2024)
  !! Adding variables by KL (Oct 2024)
  !! Optimization by DK, TYK, JLB (May 2025)
  !! Implemention of ELPA by STiwari and KL (Aug 2025)
  !! Supporting q-parallelism by KL and STiwari (Sep 2025)
  !! Supporting excited states by KL (Oct 2025)
  !! Modularization refactoring by DK (Jul 2026)
  !!
  USE kinds,          ONLY : DP
  USE polaron_common, ONLY : is_mirror_k, is_mirror_q, is_tri_k, is_tri_q, nbnd_plrn, &
                              nbnd_g_plrn, lword_h, lword_g, lword_m, io_level_g_plrn, &
                              io_level_h_plrn, hblocksize, band_pos, ik_edge, Rp,       &
                              select_bands_plrn, kpg_map, etf_all, xkf_all, Hamil,      &
                              eigvec, gq_model, epf, epfall
  USE polaron_grid,   ONLY : ikq_all, isGVec, ikqLocal2Global, indexGamma, index_shift
  USE io_polaron,     ONLY : plrn_save_g_to_file, plrn_collect_image,                  &
                              write_real_space_wavefunction,                           &
                              scell_write_real_space_wavefunction
  USE polaron_interpolation, ONLY : interp_plrn_wf, interp_plrn_bq
  USE polaron_scf_driver, ONLY : polaron_scf, find_band_extreme, gather_band_eigenvalues
  USE buffers,        ONLY : open_buffer, close_buffer
  USE io_var,         ONLY : iepfall, ihamil
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: plrn_prepare, plrn_flow_select, plrn_save_g_to_file
  PUBLIC :: is_mirror_q, is_mirror_k, kpg_map, ikq_all
  PUBLIC :: plrn_collect_image
  
  CONTAINS
  
    !
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_prepare(totq, iq_restart)
    !-----------------------------------------------------------------------
    !!
    !! Routine to prepare quantities for polaron calculation
    !!
    !-----------------------------------------------------------------------
    !
    USE input,         ONLY : start_band_plrn, end_band_plrn, nbndsub, nstate_plrn, debug_plrn, &
                              cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn,   &
                              model_vertex_plrn, nhblock_plrn, g_start_band_plrn,               &
                              g_end_band_plrn, lrot, lphase, type_plrn, g_start_energy_plrn,    &
                              g_end_energy_plrn, model_enband_plrn, model_vertex_plrn,          &
                              g_power_order_plrn, io_lvl_plrn, nkf1, nkf2, nkf3, m_eff_plrn,    &
                              kappa_plrn, omega_LO_plrn, lfast_kmesh, scell_mat_plrn
    USE global_var,    ONLY : nkqf, nkf, nqf, nqtotf, nktotf, etf, xkf, xqf
    USE modes,         ONLY : nmodes
    USE cell_base,     ONLY : bg, omega, alat
    USE ep_constants,  ONLY : czero, cone, pi, ci, twopi, fpi, eps6, eps8, eps5, zero, ryd2ev
    USE parallelism,   ONLY : poolgather2
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE io_files,      ONLY : check_tempdir
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out)  :: iq_restart
    !! Restart q-point indexx
    INTEGER, INTENT(in)   :: totq
    !! Total q-point in the fsthick
    !
    ! Local variables
    LOGICAL :: debug
    !! Debug flag
    LOGICAL :: plrn_scf
    !! .true. if self-consistent polaron calculation is to be performed
    LOGICAL :: exst
    !! Checking on the file presence
    LOGICAL :: pfs
    !! FIXME
    INTEGER :: iq
    !! q-point index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: ikq
    !! k+q point index
    INTEGER :: ik_global
    !! k-point index in global list
    INTEGER :: ierr
    !! Error status
    INTEGER :: ikGamma
    !! Index of \Gamma point in k-point list
    INTEGER :: iqGamma
    !! Index of \Gamma point in q-point list
    INTEGER :: minNBlock
    !! FIXME
    INTEGER :: ishift
    !! Index of neighbor G-vectors
    INTEGER :: lword_h_tmp
    !! Temporary Hamiltonian record size for I/O
    INTEGER :: lword_g_tmp
    !! Temporary el-ph matrix element record size for I/O
    INTEGER, PARAMETER :: maxword = HUGE(1)
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: efermi
    !! Fermi level (VBM or CBM)
    REAL(KIND = DP) :: klen
    !! Distance to the nearest \Gamma point
    REAL(KIND = DP) :: shift(3)
    !! Shift G-vector to find nearest \Gamma point
    REAL(KIND = DP) :: rfac
    !! Numerator of matrix element in Frohlich model
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:,:)
    !! Temporary array for gathering eigenvalues across pools
    !
    CALL start_clock('plrn_prepare')
    !
    IF(lfast_kmesh) THEN
      CALL errore('polaron_prepare', 'Polaron module not working with lfast_kmesh', 1)
    ENDIF
    !
    WRITE(stdout, '(5x,"fsthick not working in polaron module, selecting all the k/q points.")')
    !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
    !! Legalize the type_plrn input,     in case that the user use an arbitrary number
    IF(type_plrn < 0) THEN
      type_plrn = -1
      WRITE(stdout, '(5x, "The electron polaron is calculated.")')
    ELSE
      type_plrn = 1
      WRITE(stdout, '(5x, "The hole polaron is calculated.")')
    ENDIF
    !
    lrot = .TRUE.
    lphase = .TRUE.
    !
    debug = debug_plrn
    IF (debug_plrn) CALL check_tempdir('test_out', exst, pfs)
    !
    WRITE(stdout,'(5x,a)') REPEAT('=',67)
    !
    IF (g_start_band_plrn == 0) g_start_band_plrn = 1
    IF (g_end_band_plrn == 0) g_end_band_plrn = nbndsub
    nbnd_g_plrn = g_end_band_plrn - g_start_band_plrn + 1
    !
    IF (start_band_plrn == 0) start_band_plrn = g_start_band_plrn
    IF (end_band_plrn == 0) end_band_plrn = g_end_band_plrn
    nbnd_plrn = end_band_plrn - start_band_plrn + 1
    !
    IF(g_start_band_plrn > start_band_plrn .OR. g_end_band_plrn < end_band_plrn) THEN
      CALL errore('polaron_prepare', 'Selecting more bands in polaron than saving g matrix', 1)
    ENDIF
    !
    ALLOCATE(select_bands_plrn(nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating select_bands_plrn', 1)
    !
    select_bands_plrn = 0
    DO ibnd = 1, nbnd_plrn
      select_bands_plrn(ibnd) = start_band_plrn + ibnd - 1
    ENDDO
    !
    !! copy q(x,y,z) to xkf_all, save the copy of all kpoints
    !! Note that poolgather2 has the dimension of nktotf*2,
    !! which has k at ik and k+q at ik+1
    !! This is because that xkf has the dimension of 3, nkf*2
    !! where the ik is k and ik+1 is k+q
    ALLOCATE(xkf_all(3, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating xkf_all', 1)
    ALLOCATE(rtmp2(3, nktotf * 2), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating rtmp2', 1)
    ALLOCATE(epf(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating epf', 1)
    epf = czero
    !
    xkf_all = zero
    rtmp2 = zero
    CALL poolgather2(3, nktotf * 2, nkqf, xkf, rtmp2)
    xkf_all(1:3, 1:nktotf) = rtmp2(1:3, 1:nktotf * 2:2)
    !
    DEALLOCATE(rtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error deallocating rtmp2', 1)
    !
    WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0)") start_band_plrn, end_band_plrn, nbnd_plrn
    WRITE(stdout, "(5x, 'Including bands: ', 10i3)") select_bands_plrn
    WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0, ' in saving g')") &
          g_start_band_plrn, g_end_band_plrn, nbnd_g_plrn
    !
    WRITE(stdout, "(5x, 'Gathering eigenvalues of ', i0, ' bands and ', i0, ' k points')") nbndsub, nktotf
    !
    ALLOCATE(etf_all(nbndsub, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating etf_all', 1)
    !
    IF(model_enband_plrn) THEN
      etf = zero
      ! Find the distance to the nearest Gamma point in crystal coordinates
      DO ik = 1, 2 * nkf
        klen = 1E3
        xxk = xkf(:, ik)
        CALL cryst_to_cart(1, xxk, bg, 1)
        DO ishift = 1, 27
          shift(1:3) = REAL(index_shift(ishift), KIND = DP)
          xxk = xkf(:, ik) + shift
          CALL cryst_to_cart(1, xxk, bg, 1)
          klen = MIN(klen, NORM2(xxk))
        ENDDO
        etf(1, ik) = 0.5 / m_eff_plrn * (klen * twopi / alat)**2
      ENDDO
    ENDIF
    !
    etf_all = zero
    CALL gather_band_eigenvalues(etf, etf_all)
    !
    IF(model_vertex_plrn) THEN
      WRITE(stdout, '(5x, a, f8.3)') "Using model g vertex, with order ", g_power_order_plrn
      ALLOCATE(gq_model(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating gq_model', 1)
      gq_model = zero
      !
      rfac = SQRT(fpi / omega * omega_LO_plrn / kappa_plrn)
      DO iq = 1, nqf
        klen = 1E3
        DO ishift = 1, 27
          shift(1:3) = REAL(index_shift(ishift), KIND = DP)
          xxk = xqf(:, iq) + shift
          CALL cryst_to_cart(1, xxk, bg, 1)
          klen = MIN(klen, NORM2(xxk))
        ENDDO
        IF(klen > eps8) THEN
          gq_model(iq) = rfac / ((klen * twopi / alat) ** g_power_order_plrn)
        ENDIF
      ENDDO
    ENDIF

    ! change unit from eV to Rydberg
    g_start_energy_plrn = g_start_energy_plrn / ryd2ev
    g_end_energy_plrn = g_end_energy_plrn / ryd2ev
    !
    CALL start_clock('find_EVBM')
    CALL find_band_extreme(type_plrn, etf_all, ik_edge, band_pos, efermi)
    !
    ! Determine the Fermi energy, read from the input or calculated from band structure
    WRITE(stdout, '(5x, "Fermi Energy is", f16.7, &
    &" (eV) located at kpoint ", i6, 3f8.3, " band ", i3)') efermi * ryd2ev, ik_edge, xkf_all(1:3, ik_edge), band_pos
    ! Shift the eigenvalues to make VBM/CBM zero
    etf_all(1:nbndsub, 1:nktotf) = etf_all(1:nbndsub, 1:nktotf) - efermi
    !
    CALL stop_clock('find_EVBM')
    !
    WRITE(stdout, "(5x, 'Allocating arrays and open files.')")
    !
    IF(interp_Ank_plrn .OR. interp_Bqu_plrn .OR. cal_psir_plrn) THEN
      plrn_scf = .FALSE.
      restart_plrn = .TRUE.
    ELSE
      plrn_scf = .TRUE.
    ENDIF
    !
    IF(restart_plrn) THEN
      iq_restart = totq + 1
    ELSE
      iq_restart = 1
    ENDIF
    !
    IF(interp_Ank_plrn) THEN
      ALLOCATE(eigvec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigvec', 1)
      eigvec = czero
    ELSE IF(plrn_scf) THEN
       CALL check_tempdir('plrn_tmp', exst, pfs)
       !
       io_level_g_plrn = 1
       lword_g_tmp = nbnd_g_plrn * nbnd_g_plrn * nmodes * nkf
       IF(lword_g_tmp > maxword) THEN
         CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores!', 1)
       ELSE
         lword_g = INT(lword_g_tmp)
       ENDIF
       !
       IF (io_lvl_plrn == 0) THEN
         ALLOCATE(epfall(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf, nqtotf), STAT = ierr)
         IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating epfall', 1)
         epfall = czero
       ELSE IF (io_lvl_plrn == 1) THEN
         CALL open_buffer( iepfall , 'ephf' , lword_g, io_level_g_plrn, exst, direc = 'plrn_tmp/')
       ENDIF
       !
       io_level_h_plrn = 1
       IF(nhblock_plrn < 1 .OR. nhblock_plrn > nkf * nbnd_plrn) THEN
         CALL errore('plrn_prepare','Illegal nhblock_plrn, should between 1 and nkf * nbnd_plrn', 1)
       ENDIF
       !
       minNBlock = CEILING(REAL(nkf * nbnd_plrn * nktotf * nbnd_plrn, dp) / maxword)
       !
       IF(minNBlock >  nhblock_plrn .AND. nhblock_plrn /= 1) THEN
         CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores!', 1)
       ENDIF
       !
       hblocksize = CEILING(REAL(nkf * nbnd_plrn, dp) / nhblock_plrn)
       !
       lword_h_tmp = nktotf * nbnd_plrn * hblocksize
       !
       IF(nhblock_plrn /= 1) THEN
         IF(lword_h_tmp > maxword) THEN
           CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores or larger nhblock_plrn!', 1)
         ELSE
           lword_h = lword_h_tmp
         ENDIF
       ENDIF
       !
       lword_m = nbnd_plrn * nbnd_plrn * nktotf * 3
       !
       ALLOCATE(Hamil(nktotf * nbnd_plrn, hblocksize), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating Hamil', 1)
       !
       ! Allocate and initialize the variables
       ALLOCATE(eigvec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigvec', 1)
       eigvec = czero
       !
       ! Check whether the input is legal, otherwise print warning and stop
       ! Check the input now, because if inputs are illegal, we can stop the calculation
       ! before the heavy el-ph interpolation begins.
       IF(nktotf /= nqtotf .OR. nktotf < 1) CALL errore('plrn_prepare','Not identical k and q grid. Do use same nkf and nqf!', 1)
       !
       IF( (.NOT. scell_mat_plrn) .AND. (nkf1 == 0 .or. nkf2 == 0 .or. nkf3 == 0) ) THEN
          CALL errore('plrn_prepare','Try to use nkf and nqf to generate k and q grid, &
             &IF you are using a manual grid, also provide this information.', 1)
       ENDIF
       !
       IF(nkf < 1) CALL errore('plrn_prepare','Some node has no k points!', 1)
       IF(nqtotf /= nqf) CALL errore('plrn_prepare','Parallel over q is not available for polaron calculations.', 1)
       !
       WRITE(stdout, '(5x, "Polaron wavefunction calculation starts with k points ",&
          &i0, ", q points ", i0, " and KS band ", i0)') nktotf,  nqtotf,  nbnd_plrn
       !
       ! check whether the k and q mesh are identical
       ! This may not be theoretically necessary,
       ! but necessary in this implementation
       DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         IF (ANY(ABS(xkf_all(1:3, ik_global) - xqf(1:3, ik_global)) > eps6)) THEN
           CALL errore('plrn_prepare', 'The k and q meshes must be exactly the same!', 1)
         ENDIF
       ENDDO
       !
       ! map iq to G-iq, ik to G-ik
       ! find the position of Gamma point in k and q grid
       ! Note that xqf is not a MPI-local variable, xqf = xqtotf otherwise
       ! the program gives wrong results
       ikGamma = indexGamma(xkf_all)
       iqGamma = indexGamma(xqf)
       IF (ikGamma == 0) CALL errore('plrn_prepare','k = 0 not included in k grid!', 1)
       IF (iqGamma == 0) CALL errore('plrn_prepare','q = 0 not included in q grid!', 1)
       WRITE(stdout, '(5x, "The index of Gamma point in k grid is ", i0, " &
          &and q grid IS ", i0)') ikGamma, iqGamma
       !
       ! Given k, find the index of -k+G, to impose the relation
       ! A_{n,-k+G} = A^*_{n,k} and B_{-q+G, \nu} = B^*_{q, \nu}
       ! the relation k1 + k2 = G is unique given k1
       ! Since both A and B have the dimension of nktotf/nqtotf,
       ! kpg_map should map all the k from 1 to nktotf (global)
       ALLOCATE(kpg_map(nqtotf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating kpg_map', 1)
       kpg_map = 0
       !
       ALLOCATE(is_mirror_k(nkf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_k', 1)
       is_mirror_k = .FALSE.
       ALLOCATE(is_mirror_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_mirror_q = .FALSE.
       ALLOCATE(is_tri_k(nkf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_k = .FALSE.
       ALLOCATE(is_tri_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_q = .FALSE.
       !
       WRITE(stdout, '(5x, a)') "Finding the index of -k for each k point."
       ! For two k points k1 and k2, if k1 = G - k2, G is any reciprocal vector
       ! then k2 is the mirror point of k1 if the index of k2 is larger than k1
       ! Same rule for q, while is_mirror_q(nqf) is global but is_mirror_k(nkf) is local
       DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         DO ikq = 1, nktotf
           ! -k+G = k', i.e. k' + k = G, G may be (0, 0, 0)
           xxk = xkf_all(1:3, ik_global) + xkf_all(1:3, ikq)
           IF (isGVec(xxk)) kpg_map(ik_global) = ikq
         ENDDO
         IF (kpg_map(ik_global) == 0) CALL errore('plrn_prepare', 'Not legal k/q grid!', 1)
         !
         ikq = kpg_map(ik_global)
         IF (ik_global > ikq) THEN
           is_mirror_k(ik) = .TRUE.
         ELSE IF (ik == ikq) THEN
           is_tri_k(ik) = .TRUE.
         ENDIF
       ENDDO
       !
       CALL mp_sum(kpg_map, inter_pool_comm)
       !
       DO iq = 1, nqf
         ikq = kpg_map(iq)
         IF (iq > ikq) THEN
           is_mirror_q(iq) = .TRUE.
         ELSE IF (iq == ikq) THEN
           is_tri_q(iq) = .TRUE.
         ENDIF
       ENDDO
       !
       WRITE(stdout, '(5x, a)') "Checking the k + q is included in the mesh grid for each k and q."
       ! find the global index of ik_global, ikq with vector k and k+q.
       ! Different from kpg_map, it is used in constructing Hamiltonian or hpsi
       ! ik goes over all the local k points, to parallel the program
       DO ik = 1, nkf
         DO iq = 1, nqtotf
           ik_global = ikqLocal2Global(ik, nktotf)
           xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)
           IF (ikq_all(ik, iq) == 0) THEN
             CALL errore('plrn_prepare','Not commensurate k and q grid!', 1)
           ENDIF
         ENDDO
       ENDDO
    ENDIF
    !
    WRITE(stdout, "(5x, 'End of plrn_prepare')")
    !
    CALL stop_clock('plrn_prepare')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_prepare
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_flow_select(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
    !-----------------------------------------------------------------------
    !!
    !! Driver which selects whether a self-consistent polaron
    !! or a post-processing calculation is to be performed.
    !!
    !-----------------------------------------------------------------------
    USE input,         ONLY : cal_psir_plrn,  interp_Ank_plrn, interp_Bqu_plrn, &
                              io_lvl_plrn, scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT (in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT (in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT (in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT (in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT (in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT (in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT (in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), INTENT (in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    !
    ! Local variable
    !
    ! Bqu Ank interpolation is not compatible with self-consistency process
    ! Added by Chao Lian for polaron calculations flow select
    ! If postprocess is ON, i.e. Bqu interpolation with saved dtau,
    ! Ank interpolation with saved Amp, and polaron visualization with saved Wannier function cube files,
    ! then self-consistent process is skipped.
    IF (.NOT. (interp_Bqu_plrn .OR. interp_Ank_plrn .OR. cal_psir_plrn)) THEN
      CALL polaron_scf(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
      IF (io_lvl_plrn == 1) CALL close_buffer(iepfall, 'KEEP')
      CALL close_buffer(ihamil, 'delete')
    ENDIF
    !
    IF (interp_Ank_plrn) THEN
      IF(ionode) WRITE(stdout, "(5x, 'Interpolating the Ank at given k-point set....')")
      CALL interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
    ENDIF
    !
    IF (interp_Bqu_plrn) THEN
      IF(ionode) THEN
        WRITE(stdout, "(5x, 'Interpolating the Bqu at given q-point set....')")
      ENDIF
      CALL interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
    ENDIF
    !
    IF (cal_psir_plrn) THEN
      IF(ionode) WRITE(stdout, "(5x, 'Calculating the real-space distribution of polaron wavefunction....')")
      IF (scell_mat_plrn) THEN
        CALL scell_write_real_space_wavefunction()
      ELSE
        CALL write_real_space_wavefunction()
      ENDIF
    ENDIF
    !
    ! Deallocate allocated arrays and close open files
    CALL plrn_close()
    !
    WRITE(stdout, '(/5x, "======================== Polaron Timers ===========================")')
    CALL print_clock('main_prln')
    CALL print_clock('plrn_prepare')
    CALL print_clock('write_files')
    CALL print_clock('Bqu_tran')
    CALL print_clock('Ank_trans')
    CALL print_clock('cal_E_Form')
    CALL print_clock('DiagonH')
    CALL print_clock('Setup_H')
    CALL print_clock('H_alloc')
    CALL print_clock('read_gmat')
    CALL print_clock('read_Hmat')
    CALL print_clock('Write_Hmat')
    CALL print_clock( 'cegterg' )
    CALL print_clock( 'cegterg:init' )
    CALL print_clock( 'cegterg:diag' )
    CALL print_clock( 'cegterg:update' )
    CALL print_clock( 'cegterg:overlap' )
    CALL print_clock( 'cegterg:last' )
    CALL print_clock('cal_bqu')
    CALL print_clock('init_Ank')
    CALL print_clock('find_EVBM')
    CALL print_clock('re_omega')
    CALL print_clock('cal_hpsi')
    WRITE(stdout, '(5x, "===================================================================")')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_flow_select
    !-----------------------------------------------------------------------------------
    !
    SUBROUTINE plrn_close()
    !-----------------------------------------------------------------------------------
    !! Deallocate arrays and close polaron calculations
    !-----------------------------------------------------------------------------------
    USE input,         ONLY : model_vertex_plrn, interp_Ank_plrn, &
                              interp_Bqu_plrn, cal_psir_plrn, scell_mat_plrn
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (.NOT. (interp_Ank_plrn .OR. interp_Bqu_plrn .OR. cal_psir_plrn)) THEN
      DEALLOCATE(is_mirror_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_mirror_k', 1)
      DEALLOCATE(is_mirror_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_mirror_q', 1)
      DEALLOCATE(is_tri_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_tri_k', 1)
      DEALLOCATE(is_tri_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_tri_q', 1)
      DEALLOCATE(kpg_map, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating kpg_map', 1)
      DEALLOCATE(Hamil, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Hamil', 1)
    ENDIF
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Hamil', 1)
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating xkf_all', 1)
    DEALLOCATE(select_bands_plrn, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating select_bands_plrn', 1)
    IF (interp_Ank_plrn .OR. (.NOT. (interp_Bqu_plrn .OR. cal_psir_plrn))) THEN
      DEALLOCATE(eigvec, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating eigvec', 1)
    END IF
    IF (scell_mat_plrn) THEN
      DEALLOCATE(Rp, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Rp', 1)
    END IF
    IF (model_vertex_plrn) THEN
      DEALLOCATE(gq_model, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating gq_model', 1)
    END IF
    !
    !-----------------------------------------------------------------------------------    
    END SUBROUTINE plrn_close
    !-----------------------------------------------------------------------------------
  !  
  !-----------------------------------------------------------------------------------
  END MODULE polaron
  !-----------------------------------------------------------------------------------
