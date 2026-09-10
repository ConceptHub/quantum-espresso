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
  MODULE io_polaron
  !--------------------------------------------------------------------------
  !!
  !! All disk I/O for polaron calculations: buffered g, Ank and dtau, .xsf and
  !! cube files, density of states. Formats only; transforms are in
  !! polaron_interpolation.
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : nbnd_plrn, nbnd_g_plrn, nRp, Rp, select_bands_plrn, &
                             etf_all, eigvec, epf, epfall
  USE polaron_grid,   ONLY : ikq_all, index_Rp, index_shift
  USE buffers,        ONLY : save_buffer
  USE io_var,         ONLY : iepfall

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: plrn_save_g_to_file, plrn_collect_image
  PUBLIC :: write_plrn_dtau_xsf, scell_write_plrn_dtau_xsf, write_plrn_wf
  PUBLIC :: read_plrn_wf_grid, read_plrn_wf, write_plrn_bmat
  PUBLIC :: read_plrn_dtau_grid, read_plrn_dtau, calc_den_of_state
  PUBLIC :: write_real_space_wavefunction, scell_write_real_space_wavefunction
  PUBLIC :: read_Rp_in_S

  CONTAINS

    !-----------------------------------------------------------------------
    SUBROUTINE cal_f_delta(energy, sigma, f_delta)
    !-----------------------------------------------------------------------
    !! Return a normalized Gaussian distribution.
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : twopi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: energy(:)
    !! Energy variable
    REAL(KIND = DP), INTENT(in) :: sigma
    !! Width of Gaussian
    REAL(KIND = DP), INTENT(out) :: f_delta(:)
    !! Gaussian function
    !
    f_delta = EXP(-energy**2 / (2*sigma**2))/SQRT(twopi*sigma**2)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_f_delta
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    SUBROUTINE plrn_save_g_to_file(iq, epfg, wfreq)
    !-----------------------------------------------------------------------
    !! Save el-ph matrix element to file
    !!
    USE modes,         ONLY : nmodes
    USE input,         ONLY : g_start_band_plrn, g_end_band_plrn, &
                              io_lvl_plrn, eps_acoustic, g_scale_plrn
    USE global_var,    ONLY : nkf
    USE ep_constants,  ONLY : two
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q-point index
    REAL(KIND = DP), INTENT(in) :: wfreq(:, :)
    !! Phonon frequency
    COMPLEX(KIND = DP), INTENT(inout) :: epfg(:, :, :, :)
    !! el-ph matrix element
    !
    INTEGER :: ik
    !! k-point index
    INTEGER :: ikq
    !! k+q point index
    INTEGER :: imode
    !! mode index
    !
    ! In polaron equations, g is not epf but epf/omega
    ! To ensure a Hermitian Hamiltonian, g_{mnu}(k, -q) is calculated as g*_{nmu}(k-q, q)
    DO ik = 1, nkf
      ikq = ikq_all(ik, iq)
      DO imode = 1, nmodes
        IF (wfreq(imode, iq) > eps_acoustic) THEN
          !epf(:, :, imode, ik) = &
          epf(:, :, imode, ik) = g_scale_plrn * &
          epfg(g_start_band_plrn:g_end_band_plrn, g_start_band_plrn:g_end_band_plrn, imode, ik) / &
          DSQRT(two * wfreq(imode, iq))
        ENDIF
      ENDDO ! imode
    ENDDO ! ik
    !
    IF(io_lvl_plrn == 0) THEN
      epfall(:, :, : ,: ,iq) = epf(:, :, :, :)
    ELSE IF (io_lvl_plrn == 1) THEN
      CALL save_buffer(epf(:, :, :, :), nbnd_g_plrn * nbnd_g_plrn * nmodes * nkf, iepfall, iq)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_save_g_to_file
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_collect_image()
    !-----------------------------------------------------------------------
    !! collect data from multiple images
    !-----------------------------------------------------------------------
    USE mp,            ONLY : mp_sum
    USE mp_images,     ONLY : inter_image_comm
    USE input,         ONLY : restart_plrn, io_lvl_plrn
    !
    IF ((.NOT. restart_plrn) .AND. (io_lvl_plrn == 0)) THEN
      CALL mp_sum(epfall, inter_image_comm)
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_collect_image

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    SUBROUTINE write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, filename, species)
    !------------------------------------------------------------------------
    !! Write ionic positions and displacements in XSF format
    !------------------------------------------------------------------------
    USE ep_constants,   ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE mp,             ONLY : mp_sum, mp_bcast
    USE ions_base,      ONLY : nat, ityp, tau, ntypx
    USE cell_base,      ONLY : at, alat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqf1
    !! Fine q-point grid along b1
    INTEGER, INTENT(in) :: nqf2
    !! Fine q-point grid along b2
    INTEGER, INTENT(in) :: nqf3
    !! Fine q-point grid along b3
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! Polaron displacements in real space
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! Atomic species in unit cell
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    INTEGER :: nat_all
    !! Number of atoms in supercell
    INTEGER :: nptotf
    !! Number of unit cells in supercell
    INTEGER :: nqf_s(1:3)
    !! q-point grid
    INTEGER :: iRp
    !! Counter for unit cell vectors in supercell
    INTEGER :: iatm
    !! Counter for atoms in unit cell
    INTEGER :: iatm_all
    !! Counter for atoms in supercell
    INTEGER :: ika
    !! Combined atom and cartesian direction index
    INTEGER :: isp
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER, ALLOCATABLE :: elements(:)
    !! Atomic species array
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell coordinates
    REAL(KIND = DP) :: shift(1:3)
    !! Shift vector coordinates
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! Atomic position coordinates in supercell
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! Displacement coordinates in supercell
    !
    ! total number of atoms is
    ! (number of atoms in unit cell) x (number of cells in the supercell)
    nptotf =  nqf1 * nqf2 * nqf3
    nqf_s  = (/nqf1, nqf2, nqf3/)
    nat_all = nat * nptotf
    !
    ALLOCATE(atoms(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating atoms', 1)
    ALLOCATE(elements(nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating elements', 1)
    ALLOCATE(displacements(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating displacements', 1)
    !
    atoms = zero
    elements = 0
    displacements = zero
    !
    cell(1:3, 1) = at(1:3, 1) * nqf1
    cell(1:3, 2) = at(1:3, 2) * nqf2
    cell(1:3, 3) = at(1:3, 3) * nqf3
    !
    iatm_all = 0
    DO isp = 1, ntypx
      DO iRp = 1, nptotf
        Rp_vec(1:3) = index_Rp(iRp, nqf_s)
        DO iatm = 1, nat
          IF(ityp(iatm) == isp) THEN
            iatm_all = iatm_all + 1
            ika = (iatm - 1) * 3 + 1
            !Rp(1:3) = (ix - nqf1/2) * at(1:3, 1) + (iy - nqf2/2) * at(1:3, 2) + (iz - nqf3/2) * at(1:3, 3)
            shift(1:3) = Rp_vec(1) * at(1:3, 1) + Rp_vec(2) * at(1:3, 2) + Rp_vec(3) * at(1:3, 3)
            IF(PRESENT(species)) THEN
              elements(iatm_all) = species(ityp(iatm))
            ELSE
              elements(iatm_all) = ityp(iatm)
            ENDIF
            atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
            displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika + 2))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    cell = cell * alat
    atoms = atoms * alat
    !
    CALL write_xsf_file(filename, cell * bohr2ang, elements, atoms * bohr2ang, displacements * bohr2ang)
    !
    DEALLOCATE(atoms, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating atoms', 1)
    DEALLOCATE(elements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating elements', 1)
    DEALLOCATE(displacements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating displacements', 1)
    !----------------------------------------------------------------------------------------
    END SUBROUTINE write_plrn_dtau_xsf
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    SUBROUTINE scell_write_plrn_dtau_xsf(dtau, nqtotf_p, nRp_p, Rp_p, as_p, filename, species)
    !----------------------------------------------------------------------------------------
    !! JLB: Write ionic positions and displacements for transformed supercell
    !----------------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE mp,            ONLY : mp_sum, mp_bcast
    USE ions_base,     ONLY : nat, ityp, tau, ntypx
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !! Output file name
    INTEGER, INTENT(in) :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Coordinates of lattice vectors in supercell
    REAL(KIND = DP), INTENT(in) :: as_p(3,3)
    !! Supercell lattice vectors
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! Polaron displacement coordinates
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! Atomic species in unit cell
    !
    ! Local variable
    INTEGER :: ierr
    !! Error index
    INTEGER :: nat_all
    !! Total number of atoms in supercell
    INTEGER :: iRp
    !! Lattice vector counter in supercell
    INTEGER :: iatm
    !! Atom counter in unit cell
    INTEGER :: iatm_all
    !! Atom counter in supercell
    INTEGER :: ika
    !! Combined atom and cartesian direction index
    INTEGER :: isp
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER, ALLOCATABLE :: elements(:)
    !! Atomic species
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP) :: shift(1:3)
    !! Shift vector coordinates
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! Atomic position coordinates in supercell
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! Atomic displacement coordinates in supercell
    !
    ! total number of atoms is
    ! (number of atoms in unit cell) x (number of cells in the supercell)
    nat_all = nat * nRp_p
    !
    ALLOCATE(atoms(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating atoms', 1)
    ALLOCATE(elements(nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating elements', 1)
    ALLOCATE(displacements(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating displacements', 1)
    !
    atoms = zero
    elements = 0
    displacements = zero
    !
    cell(1:3, 1) = as_p(1, 1:3)
    cell(1:3, 2) = as_p(2, 1:3)
    cell(1:3, 3) = as_p(3, 1:3)
    !
    iatm_all = 0
    DO isp = 1, ntypx
      DO iRp = 1, nRp_p
        Rp_vec(1:3) = Rp_p(1:3, iRp)
        DO iatm = 1, nat
          IF(ityp(iatm) == isp) THEN
            iatm_all = iatm_all + 1
            ika = (iatm - 1) * 3 + 1
            !Rp(1:3) = (ix - nqf1/2) * at(1:3, 1) + (iy - nqf2/2) * at(1:3, 2) + (iz - nqf3/2) * at(1:3, 3)
            shift(1:3) = Rp_vec(1) * at(1:3, 1) + Rp_vec(2) * at(1:3, 2) + Rp_vec(3) * at(1:3, 3)
            IF(PRESENT(species)) THEN
              elements(iatm_all) = species(ityp(iatm))
            ELSE
              elements(iatm_all) = ityp(iatm)
            ENDIF
            atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
            displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika + 2))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    cell = cell * alat
    atoms = atoms * alat
    !
    CALL write_xsf_file(filename, cell * bohr2ang, elements, atoms * bohr2ang, displacements * bohr2ang)
    !
    DEALLOCATE(atoms, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating atoms', 1)
    DEALLOCATE(elements, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating elements', 1)
    DEALLOCATE(displacements, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating displacements', 1)
    !---------------------------------------------------------------------------
    END SUBROUTINE scell_write_plrn_dtau_xsf
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    SUBROUTINE write_xsf_file(filename, cell, elements, atoms, forces, data_cube)
    !---------------------------------------------------------------------------
    !! Write xsf to file
    !---------------------------------------------------------------------------
    USE io_var,            ONLY : ixsfplrn
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in):: filename
    !! Output file name
    INTEGER, INTENT(in) :: elements(:)
    !! Atomic species
    REAL(KIND = DP), INTENT(in) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP), INTENT(in) :: atoms(:, :)
    !! Atomic position coordinates
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: forces(:, :)
    !! Atomic displacement coordinates
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: data_cube(:, :, :)
    !! Data for isosurface
    !
    ! Local variables
    INTEGER :: ix
    !! x-coordinate counter
    INTEGER :: iy
    !! y-coordinate counter
    INTEGER :: iz
    !! z-coordinate counter
    INTEGER :: iatm
    !! Atom counter
    INTEGER :: natm
    !! Number of atoms in supercell
    INTEGER :: shapeTemp(3)
    !! Shape of isosurface data
    !
    natm = UBOUND(elements, DIM = 1)
    !
    OPEN(UNIT = ixsfplrn, FILE = TRIM(filename), FORM = 'formatted', STATUS = 'unknown')
    !
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') '# Generated by the EPW polaron code'
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') 'CRYSTAL'
    WRITE(ixsfplrn, '(a)') 'PRIMVEC'
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 1)
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 2)
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 3)
    WRITE(ixsfplrn, '(a)') 'PRIMCOORD'
    ! The second number is always 1 for PRIMCOORD coordinates,
    ! according to http://www.xcrysden.org/doc/XSF.html
    WRITE(ixsfplrn, '(2i6)')  natm, 1
    !
    DO iatm = 1, natm
      IF (PRESENT(forces)) THEN
        WRITE(ixsfplrn,'(I3, 3x, 3f15.9, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm), forces(1:3, iatm)
      ELSE
        WRITE(ixsfplrn,'(I3, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm)
      ENDIF
    ENDDO
    !
    IF(PRESENT(data_cube)) THEN
      shapeTemp = SHAPE(data_cube)
      WRITE(ixsfplrn, '(/)')
      WRITE(ixsfplrn, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(ixsfplrn, '(3i6)') SHAPE(data_cube)
      WRITE(ixsfplrn, '(3f12.6)') 0.0, 0.0, 0.0
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 1)
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 2)
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 3)
      ! TODO: data cube is probably to large to take in the same way of lattice information
      ! May be usefull and implemented in the furture
      WRITE(ixsfplrn, *) (((data_cube(ix, iy, iz), ix = 1, shapeTemp(1)), &
        iy = 1, shapeTemp(2)), iz = 1, shapeTemp(3))
      WRITE(ixsfplrn, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
    ENDIF
    CLOSE(ixsfplrn)
    !----------------------------------------------------------------------------------------
    END SUBROUTINE write_xsf_file
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    SUBROUTINE write_plrn_wf(eigvec_coef, filename, enk_all)
    !----------------------------------------------------------------------------------------
    !! Write polaron wavefunction coefficients
    !----------------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev
    USE global_var,    ONLY : nktotf
    USE io_var,        ONLY : iwfplrn
    USE input,         ONLY : nstate_plrn, nkf1, nkf2, nkf3, nbndsub, scell_mat_plrn
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues on global fine grid
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch (Ank) or Wannier (Amp) basis
    COMPLEX(KIND = DP) :: phase
    !! phase of the first element of the first polaron eigenvector
    !
    ! Local variables
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: nbnd_out
    !! Number of bands in polaron wave function expansion
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: iplrn
    !! Polaron state counter
    !
    IF(PRESENT(enk_all)) THEN
      nbnd_out = nbnd_plrn
    ELSE
      nbnd_out = nbndsub
    ENDIF
    ! now rotate all the eigenvectors by a global phase
    ! to ensure the first element of the first eigenvector
    ! is pure real - KL
    phase = EXP(-CMPLX(0.0_DP, ATAN2(AIMAG(eigvec_coef(1, 1)), REAL(eigvec_coef(1,1)))))
    phase = phase / ABS(phase)
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    !WRITE(iwfplrn, '(a, 3f22.12)') 'The global phase', phase, ABS(phase)
    !
    IF (scell_mat_plrn) THEN
      WRITE(iwfplrn, '(a, 3I10)') 'Scell', nktotf, nbndsub, nstate_plrn
    ELSE
      WRITE(iwfplrn, '(6I10)') nkf1, nkf2, nkf3, nktotf, nbndsub, nstate_plrn
    ENDIF
    !
    DO ik = 1, nktotf
      DO ibnd = 1, nbnd_out
        DO iplrn = 1, nstate_plrn
          indexkn1 = (ik - 1) * nbnd_out + ibnd
          IF (PRESENT(enk_all)) THEN
            !WRITE(iwfplrn, '(2I5, 4f15.7)') ik, ibnd, enk_all(select_bands_plrn(ibnd), ik) * ryd2ev, &
            WRITE(iwfplrn, '(I12, I5, 4f15.7)') ik, ibnd, enk_all(select_bands_plrn(ibnd), ik) * ryd2ev, &
              ! eigvec_coef(indexkn1, iplrn), ABS(eigvec_coef(indexkn1, iplrn))
               eigvec_coef(indexkn1, iplrn) * phase, ABS(eigvec_coef(indexkn1, iplrn))
          ELSE
            !WRITE(iwfplrn, '(2f15.7)') eigvec_coef(indexkn1, iplrn)
            WRITE(iwfplrn, '(2f15.7)') eigvec_coef(indexkn1, iplrn) * phase
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    CLOSE(iwfplrn)
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE write_plrn_wf
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, &
               nbndsub_p, nplrn_p, filename, scell)
    !----------------------------------------------------------------------------
    !! Read k-point grid in which polaron wave function has been written to file.
    !! Needed for correct allocation when interp_plrn_wf.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    INTEGER, INTENT(out) :: nkf1_p
    !! Number of k-points in fine grid along b1
    INTEGER, INTENT(out) :: nkf2_p
    !! Number of k-points in fine grid along b2
    INTEGER, INTENT(out) :: nkf3_p
    !! Number of k-points in fine grid along b3
    INTEGER, INTENT(out) :: nktotf_p
    !! Total number of k-points in fine grid
    INTEGER, INTENT(out) :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER, INTENT(out) :: nplrn_p
    !! Number of polaron states
    !
    ! Local variables
    LOGICAL :: scell_
    !! Dummy variable to set a default value (.FALSE.) for scell
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variable to read from scell wf file
    !
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    scell_ = .FALSE.
    IF (PRESENT(scell)) scell_ = scell
    !
    IF (scell_) THEN
      READ(iwfplrn, '(a, 3I10)') dmmy, nktotf_p, nbndsub_p, nplrn_p
      !READ(iwfplrn, '(a, 3I10)') dmmy, nktotf_p, nbndsub_p, nPlrn_p
      ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
      ! Just assigning an arbitrary value
      nkf1_p = 0
      nkf2_p = 0
      nkf3_p = 0
    ELSE
      READ(iwfplrn, '(6I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p
      !READ(iwfplrn, '(6I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nPlrn_p
      IF(nkf1_p * nkf2_p * nkf3_p /= nktotf_p) THEN
        CALL errore("read_plrn_wf_grid", 'Amp.plrn'//'Not generated from the uniform grid!', 1)
      ENDIF
    ENDIF
    !
    CLOSE(iwfplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_wf_grid
    !-----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_wf(eigvec_coef, nkf1_p, nkf2_p, nkf3_p, nktotf_p, &
               nbndsub_p, nplrn_p, filename, enk_all)
    !----------------------------------------------------------------------------
    !! Read polaron wavefunction coefficients.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    INTEGER, INTENT(in) :: nkf1_p
    !! Number of k-points in fine grid along b1
    INTEGER, INTENT(in) :: nkf2_p
    !! Number of k-points in fine grid along b2
    INTEGER, INTENT(in) :: nkf3_p
    !! Number of k-points in fine grid along b3
    INTEGER, INTENT(in) :: nktotf_p
    !! Total number of k-points in fine grid
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER, INTENT(in) :: nplrn_p
    !! Number of polaron states
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues in global fine k-point grid
    !
    ! Local variables
    INTEGER:: ierr
    !! Error status
    INTEGER :: indexkn1
    !! Combined band and k-poin index
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: i1
    !! Dummy integer to read from file
    INTEGER :: i2
    !! Dummy integer to read from file
    REAL(KIND = DP) :: r1
    !! Dummy real to read from file
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    ! First line should have been read already,
    ! so that eigvec has been properly allocated
    READ(iwfplrn, *)
    !
    eigvec_coef = czero
    DO ik = 1, nktotf_p
      DO ibnd = 1, nbndsub_p
        DO iplrn = 1, nplrn_p
          indexkn1 = (ik - 1) * nbndsub_p + ibnd
          IF(PRESENT(enk_all)) THEN ! Ank.plrn is read
            !READ(iwfplrn, '(2I5, 3f15.7)') i1, i2, r1, eigvec_coef(indexkn1, iplrn)
            READ(iwfplrn, '(I12, I5, 3f15.7)') i1, i2, r1, eigvec_coef(indexkn1, iplrn)
          ELSE ! Amp.plrn is read
            READ(iwfplrn, '(2f15.7)') eigvec_coef(indexkn1, iplrn)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CLOSE(iwfplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE write_plrn_bmat(bqv_coef, filename, enk_all)
    !-----------------------------------------------------------------------------------
    !!
    !! Write Bqv coeffients (or dtau displacements) and phonon frequencies to filename
    !!
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev
    USE global_var,    ONLY : nqtotf
    USE io_var,        ONLY : idtauplrn
    USE input,         ONLY : nqf1, nqf2, nqf3, scell_mat_plrn
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(in) :: bqv_coef(:, :)
    !! Polaron displacement coefficients
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues in fine grid
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: imode
    !! Phonon mode counter
    !
    OPEN(UNIT = idtauplrn, FILE = TRIM(filename))
    IF (scell_mat_plrn) THEN
      WRITE(idtauplrn, '(a, 2I10)') 'Scell', nqtotf, nmodes
    ELSE
      WRITE(idtauplrn, '(5I10)') nqf1, nqf2, nqf3, nqtotf, nmodes
    ENDIF
    !
    DO iq = 1, nqtotf
      DO imode = 1, nmodes ! p
        IF (PRESENT(enk_all)) THEN ! write Bqv
          !JLB: Changed format for improved accuracy
          WRITE(idtauplrn, '(2I5, 4ES18.10)') iq, imode, enk_all(imode, iq) * ryd2mev, &
                                                  bqv_coef(iq, imode), ABS(bqv_coef(iq, imode))
        ELSE ! write \dtau
          !JLB: Changed format for improved accuracy
          WRITE(idtauplrn, '(2ES18.10)') bqv_coef(iq, imode)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(idtauplrn)
    !---------------------------------------------------------------------------------
    END SUBROUTINE write_plrn_bmat
    !---------------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_dtau_grid(nqf1_p, nqf2_p, nqf3_p, nqtotf_p, &
               nmodes_p, filename, scell)
    !----------------------------------------------------------------------------
    !! Read q-point grid in which polaron displacements have been written to file.
    !! Needed for correct allocation when interp_plrn_bq.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    INTEGER, INTENT(out) :: nqf1_p
    !! Number of q-points in fine grid along b1
    INTEGER, INTENT(out) :: nqf2_p
    !! Number of q-points in fine grid along b2
    INTEGER, INTENT(out) :: nqf3_p
    !! Number of q-points in fine grid along b3
    INTEGER, INTENT(out) :: nqtotf_p
    !! Total number of q-points in fine grid
    INTEGER, INTENT(out) :: nmodes_p
    !! Number of phonon modes
    !
    ! Local variables
    LOGICAL :: scell_
    !! Dummy variable to set a default value (.FALSE.) for scell
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variable to read from scell wf file
    !
    !
    OPEN(UNIT = idtauplrn, FILE = 'dtau.plrn')
    !
    scell_ = .FALSE.
    IF (PRESENT(scell)) scell_ = scell
    !
    IF (scell_) THEN
      READ(idtauplrn, '(a, 3I10)') dmmy, nqtotf_p, nmodes_p
      ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
      ! Just assigning an arbitrary value
      nqf1_p = 0
      nqf2_p = 0
      nqf3_p = 0
    ELSE
      READ(idtauplrn, '(6I10)') nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p
      IF(nqf1_p * nqf2_p * nqf3_p /= nqtotf_p) THEN
        CALL errore("read_plrn_dtau_grid", 'dtau.plrn'//'Not generated from the uniform grid!', 1)
      ENDIF
    ENDIF
    !
    CLOSE(idtauplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_dtau_grid
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_plrn_dtau(dtau_read, nqtotf_p, nmodes_p,&
                  filename, scell, wfreq)
    !-----------------------------------------------------------------------------------
    !! Read displacement coefficients in phonon (Bqv) or real-space (dtau) basis from file
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL :: scell
    !! .true. if non-diagonal supercell has been used in polaron calculation
    INTEGER, INTENT(in) :: nqtotf_p
    !! Total number of q-points in fine grid
    INTEGER, INTENT(in) :: nmodes_p
    !! Number of phonon modes
    REAL(KIND = DP), INTENT(in), OPTIONAL :: wfreq(:, :)
    !! Phonon frequencies
    COMPLEX(KIND = DP), INTENT(out) :: dtau_read(:, :)
    !! Polaron displacement coefficients
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    INTEGER :: iatm
    !! Atom counter
    INTEGER :: iq
    !! q-point counter
    INTEGER :: i1
    !! Dummy integer to read from file
    INTEGER :: i2
    !! Dummy integer to read from file
    REAL(KIND = DP) :: r1
    !! Dummy real to read from file
    !
    OPEN(UNIT = idtauplrn, FILE = TRIM(filename))
    !
    ! JLB:
    ! First line should have been read already,
    ! so that dtau has been properly allocated
    READ(idtauplrn, *)
    !
    dtau_read = czero
    DO iq = 1, nqtotf_p
      DO iatm = 1, nmodes_p
        IF(PRESENT(wfreq)) THEN ! read Bqv
          !JLB: Changed format for improved accuracy
          READ(idtauplrn, '(2I5, 3ES18.10)') i1, i2, r1, dtau_read(iq, iatm)
        ELSE ! read \dtau
          !JLB: Changed format for improved accuracy
          READ(idtauplrn, '(2ES18.10)') dtau_read(iq, iatm)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(idtauplrn)
    !--------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_dtau
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE calc_den_of_state(eigvec_coef, bqv_coef)
    !-----------------------------------------------------------------------------------
    !! Compute the DOS for Ank and Bqv coefficients
    !-----------------------------------------------------------------------------------
    USE input,         ONLY : nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn,   &
                              pdos_max_plrn, pdos_min_plrn, pdos_sigma_plrn,              &
                              istate_relax_plrn
    USE global_var,    ONLY : nqtotf, nktotf, wf
    USE io_var,        ONLY : idosplrn
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : ryd2mev, czero, one, ryd2ev, two, zero, cone, pi, ci, twopi,&
                              eps6, eps8, eps5
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in the Bloch basis, Ank
    COMPLEX(KIND = DP), INTENT(in) :: bqv_coef(:, :)
    !! Polaron displacement coefficients in the phonon basis, Bqv
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    INTEGER :: idos
    !! Counter for DOS
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: ik
    !! k-point counter
    !INTEGER :: iplrn
    !!! Polaron state counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    REAL(KIND = DP) :: temp
    !! Temporary max or min eigenvalue
    REAL(KIND = DP), ALLOCATABLE :: f_tmp(:)
    !! Temporary Gaussian function
    REAL(KIND = DP), ALLOCATABLE :: edos(:)
    !! Ank DOS
    REAL(KIND = DP), ALLOCATABLE :: pdos(:)
    !! Bqv DOS
    REAL(KIND = DP), ALLOCATABLE :: sdos(:)
    !! Sqv=wqv|Bqv|^2/2 DOS: Huang-Rhys factor
    REAL(KIND = DP), ALLOCATABLE :: edos_all(:)
    !! Electron DOS
    REAL(KIND = DP), ALLOCATABLE :: pdos_all(:)
    !! Phonon DOS
    REAL(KIND = DP), ALLOCATABLE :: e_grid(:)
    !! Grid of energy points to compute Ank DOS
    REAL(KIND = DP), ALLOCATABLE :: p_grid(:)
    !! Grid of energy points to compute Bqv DOS
    !
    !Calculating DOS
    ALLOCATE(f_tmp(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating f_tmp', 1)
    ALLOCATE(e_grid(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating e_grid', 1)
    ALLOCATE(edos(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos', 1)
    ALLOCATE(edos_all(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos_all', 1)
    temp = MAXVAL(etf_all) * ryd2ev + one
    IF (edos_max_plrn < temp) edos_max_plrn = temp
    temp = MINVAL(etf_all) * ryd2ev - one
    IF (edos_min_plrn > temp) edos_min_plrn = temp
    !
    e_grid = zero
    DO idos = 1, nDOS_plrn
      e_grid(idos) = edos_min_plrn + idos * (edos_max_plrn - edos_min_plrn) / (nDOS_plrn)
    ENDDO
    !
    edos = zero
    edos_all = zero
    DO ik = 1, nktotf
      DO ibnd = 1, nbnd_plrn
        ! TODO : iplrn
        !  KL: we should not run the loop over iplrn
        !DO iplrn = 1, 1
          CALL cal_f_delta(e_grid - (etf_all(select_bands_plrn(ibnd), ik) * ryd2ev), &
             edos_sigma_plrn, f_tmp)
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          !edos = edos + (ABS(eigvec_coef(indexkn1, iplrn))**2) * f_tmp
          edos = edos + (ABS(eigvec_coef(indexkn1, istate_relax_plrn))**2) * f_tmp
          edos_all = edos_all + f_tmp
        !ENDDO
      ENDDO
    ENDDO
    !
    ALLOCATE(p_grid(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating p_grid', 1)
    ALLOCATE(pdos(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos', 1)
    ALLOCATE(pdos_all(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos_all', 1)
    ALLOCATE(sdos(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating sdos', 1)
    !
    temp = MAXVAL(wf) * ryd2mev + 10.0_dp
    IF (pdos_max_plrn < temp) pdos_max_plrn = temp
    temp = MINVAL(wf) * ryd2mev - 10.0_dp
    IF (pdos_min_plrn > temp) pdos_min_plrn = temp
    !
    p_grid = zero
    DO idos = 1, nDOS_plrn
      p_grid(idos) = pdos_min_plrn + idos * (pdos_max_plrn - pdos_min_plrn) / (nDOS_plrn)
    ENDDO
    !
    pdos = zero
    pdos_all = zero
    sdos = zero
    DO iq = 1, nqtotf
      DO inu = 1, nmodes
        CALL cal_f_delta(p_grid - wf(inu, iq) * ryd2mev, pdos_sigma_plrn, f_tmp)
        pdos = pdos + (ABS(bqv_coef(iq, inu))**2) * f_tmp
        sdos = sdos + (wf(inu,iq)*ABS(bqv_coef(iq, inu))**2) * f_tmp
        pdos_all = pdos_all + f_tmp
      ENDDO
    ENDDO
    !
    OPEN(UNIT = idosplrn, FILE = 'dos.plrn')
    !WRITE(idosplrn, '(/2x, a/)') '#energy(ev)  A^2   edos  energy(mev)  B^2  pdos'
    WRITE(idosplrn, '(/2x, a/)') '#energy(ev)  A^2   edos  energy(mev)  B^2  pdos  sdos'
    DO idos = 1, nDOS_plrn
      WRITE(idosplrn, '(7f15.7)') e_grid(idos), edos(idos), &
         edos_all(idos), p_grid(idos), pdos(idos), pdos_all(idos), sdos(idos)
    ENDDO
    CLOSE(idosplrn)
    !
    DEALLOCATE(pdos_all)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos_all', 1)
    DEALLOCATE(f_tmp)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating f_tmp', 1)
    DEALLOCATE(e_grid)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating e_grid', 1)
    DEALLOCATE(edos)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos', 1)
    DEALLOCATE(edos_all)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos_all', 1)
    DEALLOCATE(p_grid)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating p_grid', 1)
    DEALLOCATE(pdos)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos', 1)
    DEALLOCATE(sdos)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating sdos', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE calc_den_of_state
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE write_real_space_wavefunction()
    !-----------------------------------------------------------------------------------
    !! Write polaron wave function in real space to file.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone, twopi, ci, bohr2ang, eps8
    USE input,         ONLY : nbndsub, step_wf_grid_plrn, &
                              plot_psir_plrn, lsign_psir_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE io_var,        ONLY : ipsirplrn
    USE mp_world,      ONLY : world_comm
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 60) :: plrn_file
    !! Name of output file
    INTEGER :: ierr
    !! Error status
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    INTEGER :: nqf_p(3)
    !! Fine q-point grid
    INTEGER :: nqtotf_p
    !! Number of k-points in fine grid
    INTEGER :: nmodes_p
    !! Number of phonon modes
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: idir
    !! Cartesian direction counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: nxx, nyy, nzz
    !! Number of grid points in real space supercell along cartesian directions
    INTEGER :: ip_min
    !! Initial lattice vector within this pool
    INTEGER :: ip_max
    !! Final lattice vector within this pool
    INTEGER :: ig_vec(1:3)
    !! Supercell latice vector coordinates
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: n_grid(3)
    !! Number of grid points in cell where Wannier functions are written
    INTEGER :: grid_start(3)
    !! Initial grid point within this pool
    INTEGER :: grid_end(3)
    !! Final grid point within this pool
    INTEGER :: n_grid_super(3)
    !! Number of grid points in supercellcell where polaron wave function is written
    INTEGER :: r_grid_vec(3)
    !! Grid point coordinates
    INTEGER :: rpc(1:3)
    !! Grid point coordinates for lattice vectors
    INTEGER :: species(50)
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER :: shift(1:3)
    !! Shift vector coordinates
    INTEGER :: ishift
    !! Shift counter
    REAL(KIND = DP) :: orig(3)
    !! Supercell origin coordinates
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP) :: r_cry(3)
    !! Polaron center in crystal coords
    REAL(KIND = DP) :: r_cart(3)
    !! Polaron center in cart coord
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! Wannier function in real space
    COMPLEX(KIND = DP) :: ctemp(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP) :: b_vec(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements in real space
    COMPLEX(KIND = DP), ALLOCATABLE :: cvec(:)
    !! Polaron wave function magnitude squared at each grid point
    REAL(KIND = DP), ALLOCATABLE :: sign_psir(:)
    !! Sign of the polaron wave function at each grid point
    REAL(DP)                        :: progress(11)
    !! Array to store the progress of the calculation
    REAL(DP)                        :: counter, counter_tot
    !! Counters for computing current progress
    INTEGER                         :: iprogress
    !! Counters for computing current progress
    CHARACTER(LEN = 256) :: tmpch
    !! Temporary character to assign name to psir files
    !
    ! read Amp.plrn, save eigvec_wan for the latter use
    IF (ionode) THEN
      CALL read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nplrn_p,  meta_ionode_id, world_comm)
    !
    ALLOCATE(eigvec_wan(nktotf_p * nbndsub_p, nplrn_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_wf', 'Error allocating eigvec_wan', 1)
    !
    IF (ionode) THEN
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !
    ! read dtau.plrn, get the displacement.
    IF (ionode) THEN
      CALL read_plrn_dtau_grid(nqf_p(1), nqf_p(2), nqf_p(3), nqtotf_p, nmodes_p, 'dtau.plrn')
    END IF
    CALL mp_bcast(nqf_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nqtotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nmodes_p, meta_ionode_id, world_comm)
    !
    ALLOCATE(dtau(nqtotf_p, nmodes_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating dtau', 1)
    !
    IF (ionode) THEN
      CALL read_plrn_dtau(dtau, nqtotf_p, nmodes_p, 'dtau.plrn')
    END IF
    CALL mp_bcast(dtau, meta_ionode_id, world_comm)
    !
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    cell(1:3, 1) = at(1:3, 1) * nqf_p(1) * alat
    cell(1:3, 2) = at(1:3, 2) * nqf_p(2) * alat
    cell(1:3, 3) = at(1:3, 3) * nqf_p(3) * alat
    !
    IF (plot_psir_plrn == 1) THEN
      plrn_file = 'psir_plrn.xsf'
    ELSE
      WRITE(tmpch,'(I6)') plot_psir_plrn
      plrn_file = 'psir_plrn_' // TRIM(ADJUSTL(tmpch)) // '.xsf'
    ENDIF
    ! Write the file head including information of structures
    IF (ionode) THEN
      IF (plot_psir_plrn == 1) THEN
        CALL write_plrn_dtau_xsf(dtau, nqf_p(1), nqf_p(2), nqf_p(3), plrn_file, species)
      ENDIF
    END IF
    !
    orig(1:3) = zero
    n_grid_super(1:3) = nqf_p(1:3) * n_grid(1:3)
    !
    IF (ionode) THEN
      OPEN(UNIT = ipsirplrn, FILE = TRIM(plrn_file), POSITION='APPEND')
      WRITE(ipsirplrn, '(/)')
      WRITE(ipsirplrn, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(ipsirplrn, '(3i6)')  n_grid_super / step_wf_grid_plrn
      WRITE(ipsirplrn, '(3f12.6)') zero, zero, zero
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 1) * bohr2ang
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 2) * bohr2ang
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 3) * bohr2ang
    ENDIF
    !
    b_vec(1:3) = twopi * ci / REAL(n_grid_super(1:3))
    !
    ALLOCATE(cvec(1:n_grid_super(1)), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating cvec', 1)
    ALLOCATE(sign_psir(1:n_grid_super(1)), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating sign_psir', 1)
    !
    CALL fkbounds(nqtotf_p, ip_min, ip_max)
    !
    ctemp = czero
    progress = (/0.0_DP, 0.1_DP, 0.2_DP, 0.3_DP, 0.4_DP, 0.5_DP, 0.6_DP, 0.7_DP, 0.8_DP, 0.9_DP, 1.0_DP/)
    iprogress = 1
    counter_tot = n_grid_super(2) * n_grid_super(3)
    IF(lsign_psir_plrn) THEN
        WRITE(stdout, "(5x, 'lsign_psir_plrn is turned on, ')")
        WRITE(stdout, "(5x, 'the sign of polaron wavefunction is retained.')")
    ENDIF
    DO nzz = 1, n_grid_super(3), step_wf_grid_plrn
      DO nyy = 1, n_grid_super(2), step_wf_grid_plrn
        cvec = czero
        sign_psir = zero
        !
        ! report progress
        counter =  REAL(nzz, KIND=DP) * nyy
        !
        IF ( counter / counter_tot > progress(iprogress)) THEN
          WRITE(stdout, '(5x, a, 1F5.2, a)') 'Current progress: ', &
                counter / counter_tot * 100, '%'
          iprogress = iprogress + 1
        ENDIF
        DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
          DO iRp = ip_min, ip_max !-nqf_p(3)/2, (nqf_p(3)+1)/2
            Rp_vec(1:3) = index_Rp(iRp, nqf_p)
            rpc(1:3) = Rp_vec(1:3) * n_grid(1:3) !- (n_grid_super(1:3))/2
            ! To make sure that all the nonzero points are included,
            ! we need to try from -1 to 1 neighbor supercells
            DO ishift = 1, 27
              shift(1:3) = index_shift(ishift)
              ig_vec(1:3) = (/nxx, nyy, nzz/) - rpc(1:3) + shift(1:3) * n_grid_super(1:3)
              IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                DO ibnd = 1, nbndsub !TODO change to nbndsub
                  indexkn1 = (iRp - 1) * nbndsub + ibnd
                  !TODO eigvec_wan(indexkn1, 1) should be eigvec_wan(indexkn1, iplrn)
                  !cvec(nxx) = cvec(nxx)  +  eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                  cvec(nxx) = cvec(nxx)  +  eigvec_wan(indexkn1, plot_psir_plrn) * &
                              wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                ENDDO !ibnd
              ENDIF
            ENDDO ! iscx
          ENDDO ! ipx
        ENDDO ! nxx
        CALL mp_sum(cvec, inter_pool_comm)
        IF (ionode) THEN
          !KL: retain the sign of psir
          IF(lsign_psir_plrn) THEN
            sign_psir(:) = DBLE(cvec(:))
            sign_psir(:) = SIGN(1.d0, sign_psir(:))
            WRITE (ipsirplrn, '(5e13.5)', ADVANCE='yes') &
              ABS(cvec(::step_wf_grid_plrn))**2 * sign_psir(::step_wf_grid_plrn)
          ELSE
            !JLB: Changed to |\Psi(r)|^{2}, I think it's physically more meaningful
            WRITE (ipsirplrn, '(5e13.5)', ADVANCE='yes') ABS(cvec(::step_wf_grid_plrn))**2
          ENDIF
        ENDIF
        ! Calculate the center of polaron
        ! TODO: not parallel, all the processors are doing the same calculations
        DO nxx = 1, n_grid_super(1)
          r_grid_vec(1:3) = (/nxx - 1, nyy - 1, nzz - 1/)
          ! |cvec|^2, not cvec**2: the center must not depend on the global phase.
          ctemp(1:3) = ctemp(1:3) + EXP(b_vec(1:3) * r_grid_vec(1:3)) * ABS(cvec(nxx))**2
        ENDDO
        ! End calculating the center of polaron
      ENDDO
    ENDDO
    !
    DO idir = 1, 3
      ! A vanishing sum has no phase, and ATAN2(0, 0) is processor dependent.
      IF (ABS(ctemp(idir)) <= eps8) THEN
        CALL errore('write_real_space_wavefunction', 'Real-space density has undefined center', idir)
      ENDIF
      r_cry(idir) = ATAN2(AIMAG(ctemp(idir)), REAL(ctemp(idir), KIND = DP)) / twopi
    ENDDO
    ! make crystal coordinates with 0 to 1
    r_cry(1:3) = r_cry - FLOOR(r_cry)
    ! cell is in bohr, the position is printed in Angstrom.
    r_cart(1:3) = zero
    DO idir = 1, 3
      r_cart(1:3) = r_cart(1:3) + r_cry(idir) * cell(1:3, idir) * bohr2ang
    ENDDO
    !
    IF (ionode) THEN
      WRITE(stdout, "(5x, 'The position of polaron:')")
      WRITE(stdout, "(5x, 3f9.4, ' in crystal coordinates')") r_cry(1:3)
      WRITE(stdout, "(5x, 3f9.4, ' in Cartesian coordinates (Angstrom)')") r_cart(1:3)
      WRITE(ipsirplrn, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
      CLOSE(ipsirplrn)
      WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
    ENDIF
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating dtau', 1)
    DEALLOCATE(eigvec_wan , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(cvec , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE write_real_space_wavefunction
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_write_real_space_wavefunction()
    !-----------------------------------------------------------------------------------
    !!JLB: write psir in transformed supercell
    !!     xsf format no longer compatible,
    !!     as .cube files are written in primitive coords.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone, twopi, ci, bohr2ang, eps8
    USE input,         ONLY : nbndsub, step_wf_grid_plrn, as
    USE io_var,        ONLY : iunpsirscell
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    USE low_lvl,       ONLY : matinv3
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 60) :: plrn_file
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    INTEGER :: species(50)
    !! Atomic species in unit cell
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: ig_vec(1:3)
    !! Supercell latice vector coordinates
    INTEGER :: ir1, ir2, ir3
    !! Real space grid counters
    INTEGER :: iRp1
    !! Lattice vector counter
    INTEGER :: iRp2
    !! Lattice vector counter
    INTEGER :: n_grid_total
    !! Number of points in real space grid within supercell
    INTEGER :: ip_min
    !! First lattice vector within this pool
    INTEGER :: ip_max
    !! Last lattice vector within this pool
    INTEGER :: n_grid(3)
    !! Number of grid points in cell where Wannier functions are written
    INTEGER :: grid_start(3)
    !! Initial grid point within this pool
    INTEGER :: grid_end(3)
    !! Final grid point within this pool
    INTEGER :: r_in_crys_p_sup(3)
    !! Real space grid point indices in supercell
    INTEGER :: ishift
    !! Shift counter
    INTEGER :: shift(3)
    !! Shift vector coordinates
    INTEGER :: idir
    !! Cartesian direction counter
    REAL(KIND = DP) :: r_cry(3)
    !! Polaron center in supercell crystal coordinates
    REAL(KIND = DP) :: r_cart(3)
    !! Polaron center in Cartesian coordinates
    COMPLEX(KIND = DP) :: ctemp(3)
    !! Phase sum used to locate the polaron center
    REAL(KIND = DP) :: r_in_crys_p(3)
    !! Unit cell grid points in crystal coordinates
    REAL(KIND = DP) :: r_in_crys_s(3)
    !! Supercell grid points in crystal coordinates
    REAL(KIND = DP) :: r_in_cart(3)
    !! Supercell grid points in cartesian coordinates
    REAL(KIND = DP) :: p2s(3, 3)
    !! Primitive to supercell coordinate transformation
    REAL(KIND = DP) :: s2p(3,3)
    !! Supercell to primitive coordinate transformation
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! Wannier functions in real space grid
    COMPLEX(KIND = DP) :: cvec
    !! Polaron wave function modulus squared
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficient in Wannier basis, Amp
    !
    ! Broadcast supercell lattice vectors
    CALL mp_bcast(as, meta_ionode_id, world_comm)
    !
    ! read Amp.plrn, save eigvec_wan
    IF (ionode) THEN
      ! Amp.plrn has the 'Scell' header in a non-diagonal supercell run.
      CALL read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn', scell = .TRUE.)
    END IF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nplrn_p,  meta_ionode_id, world_comm)
    !
    ALLOCATE(eigvec_wan(nktotf_p * nbndsub_p, nplrn_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_wf', 'Error allocating eigvec_wan', 1)
    !
    IF (ionode) THEN
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    ! Read list of Rp-s within supercell
    CALL read_Rp_in_S()
    WRITE(stdout, '(a, i12)') "     Number of unit cells within supercell:", nRp
    !
    ! Open file
    plrn_file = 'psir_plrn.scell.csv'
    IF (ionode) THEN
      OPEN(UNIT = iunpsirscell, FILE = TRIM(plrn_file), FORM = 'formatted', STATUS = 'unknown')
      WRITE(iunpsirscell, '(a)') "x , y , z, |\psi(r)|^2"
    ENDIF
    !
    ! Total number of grid points
    n_grid_total = nRp * n_grid(1) * n_grid(2) * n_grid(3)
    WRITE(stdout, '(a, i12)') "     Total grid points:", n_grid_total
    WRITE(stdout, '(a, i12)') "     Step:", step_wf_grid_plrn
    !
    ! Parallelize iRp
    CALL fkbounds(nRp, ip_min, ip_max)
    !
    ! Matrix to transform from primitive to supercell crystal coordinates
    p2s = matinv3(TRANSPOSE(as))
    p2s = MATMUL(p2s, at)
    ! Supercell to primitive coordinates
    s2p = matinv3(p2s)
    !
    ctemp(1:3) = czero
    !
    ! Loop over all the grid points
    DO iRp1 = 1, nRp
      !
      DO ir1 = 1, n_grid(1), step_wf_grid_plrn
        DO ir2 = 1, n_grid(2), step_wf_grid_plrn
          DO ir3 = 1, n_grid(3), step_wf_grid_plrn
            !
            r_in_crys_p(1:3) = (/REAL(ir1 - 1, DP) / n_grid(1), REAL(ir2 - 1, DP) / n_grid(2), REAL(ir3 - 1, DP) / n_grid(3)/) &
                               + REAL(Rp(1:3, iRp1), DP)
            !
            ! Wannier functions stored in (1:ngrid*iRp) list
            r_in_crys_p_sup(1:3) = (/ir1, ir2, ir3/) +  Rp(1:3, iRp1) * n_grid(1:3)
            !
            ! Move the r-point to the first supercell and store in cartesian coordinates for plotting
            r_in_crys_s = MATMUL(p2s, r_in_crys_p)
            r_in_crys_s = MODULO(r_in_crys_s, (/1.d0, 1.d0, 1.d0/))
            r_in_cart   = MATMUL(TRANSPOSE(as), r_in_crys_s) * alat * bohr2ang
            !
            ! Sum over p, PRB 99, 235139 Eq.(47)
            cvec = czero
            DO iRp2 = ip_min, ip_max !1, nRp
              !
              DO ishift = 1, 27
                !
                shift(1:3) = index_shift(ishift)
                ig_vec(1:3) = r_in_crys_p_sup(1:3) - Rp(1:3, iRp2) * n_grid(1:3) + MATMUL(s2p, shift(1:3)) * n_grid(1:3)
                !
                IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                  ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                  !
                  DO ibnd = 1, nbndsub ! sum over m
                    !
                    indexkn1 = (iRp2 - 1) * nbndsub + ibnd
                    cvec = cvec + eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                    !
                  ENDDO !ibnd
                  !
                ENDIF
                !
              ENDDO ! ishift
              !
            ENDDO ! iRp2
            CALL mp_sum(cvec, inter_pool_comm)
            !
            ! Write |\psi(r)|^2 data point to file
            IF (ionode) THEN
              WRITE(iunpsirscell, '(f12.6,", ", f12.6,", ", f12.6,", ", E13.5)') r_in_cart(1:3), ABS(cvec)**2
            ENDIF
            !
            ! Polaron center. r_in_crys_s is already reduced in the supercell,
            ! which is the periodic cell here.
            ctemp(1:3) = ctemp(1:3) + EXP(twopi * ci * r_in_crys_s(1:3)) * ABS(cvec)**2
            !
          ENDDO ! ir3
        ENDDO ! ir2
      ENDDO ! ir1
    ENDDO ! iRp1
    !
    DO idir = 1, 3
      ! A vanishing sum has no phase, and ATAN2(0, 0) is processor dependent.
      IF (ABS(ctemp(idir)) <= eps8) THEN
        CALL errore('scell_write_real_space_wavefunction', 'Real-space density has undefined center', idir)
      ENDIF
      r_cry(idir) = ATAN2(AIMAG(ctemp(idir)), REAL(ctemp(idir), KIND = DP)) / twopi
    ENDDO
    ! make crystal coordinates with 0 to 1
    r_cry(1:3) = r_cry - FLOOR(r_cry)
    r_cart(1:3) = MATMUL(TRANSPOSE(as), r_cry) * alat * bohr2ang
    !
    IF (ionode) THEN
      CLOSE(iunpsirscell)
      WRITE(stdout, "(5x, 'The position of polaron:')")
      WRITE(stdout, "(5x, 3f9.4, ' in crystal coordinates')") r_cry(1:3)
      WRITE(stdout, "(5x, 3f9.4, ' in Cartesian coordinates (Angstrom)')") r_cart(1:3)
      WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
    ENDIF
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_real_space_wavefunction', 'Error allocating eigvec_wan', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_write_real_space_wavefunction
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_wannier_cube(select_bands, wann_func, species, n_grid, &
               grid_start_min, grid_end_max)
    !-----------------------------------------------------------------------------------
    !! Read the nth Wannier function from prefix_0000n.cube file
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE input,         ONLY : nbndsub
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: select_bands(:)
    !! Wannier functions in which polaron wave function has been expanded
    INTEGER, INTENT(out) :: species(50)
    !! Atomic species
    INTEGER, INTENT(out) :: n_grid(3)
    !! Number of points in real space grid where Wannier functions are written
    INTEGER, INTENT(out) :: grid_start_min(3)
    !! Initial grid point within this pool
    INTEGER, INTENT(out) :: grid_end_max(3)
    !! Final grid point within this pool
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wann_func(:, :, :, :)
    !! Wannier function in real space
    !
    ! Local variables
    CHARACTER(LEN = 60) :: wancube
    !! Name of file containing Wannier function
    CHARACTER(LEN = 60) :: temp_str
    !! Temporary string
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: ie
    !! Atomic species counter
    INTEGER :: idir
    !! Cartesian direction counter
    INTEGER :: i_species
    !! Atomic species index
    INTEGER :: nbnd
    !! Number of Wannier funcions in which polaron wave function is expanded
    INTEGER :: iline
    !! Counter along line
    INTEGER :: nAtoms
    !! Total number of atoms
    INTEGER :: nxx, nyy, nzz
    !! Number of grid points in Cartesian directions
    INTEGER :: n_len_z
    !! Number of grid points within this pool
    INTEGER :: grid_start(3)
    !! Initial grid point within this loop
    INTEGER :: grid_end(3)
    !! Final grid point within this loop
    REAL(KIND = DP) :: rtempvec(4)
    !! Temporary vector to be read from .cube file
    REAL(KIND = DP) :: norm
    !! Wannier function normalization
    !
    nbnd = SIZE(select_bands)
    ! find the max and min of real space grid of Wannier functions of all Wannier orbitals
    IF (ionode) THEN
      grid_start_min(:) = 100000
      grid_end_max(:) = -100000
      DO ibnd = 1, nbndsub
        WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str !, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        DO idir = 1, 3
          IF (grid_start_min(idir) >= grid_start(idir)) grid_start_min(idir) = grid_start(idir)
          IF (grid_end_max(idir) <= grid_end(idir))   grid_end_max(idir) = grid_end(idir)
        ENDDO
        CLOSE(iun_plot)
      ENDDO
    ENDIF
    !
    CALL mp_bcast(n_grid,         meta_ionode_id, world_comm)
    CALL mp_bcast(grid_start_min, meta_ionode_id, world_comm)
    CALL mp_bcast(grid_end_max,   meta_ionode_id, world_comm)
    !
    ! Read the xth Wannier functions from prefix_0000x.cube in ionode
    ! and broadcast to all nodes
    ALLOCATE(wann_func(grid_start_min(1):grid_end_max(1), &
       grid_start_min(2):grid_end_max(2), &
       grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    wann_func = zero
    species = 0
    IF (ionode) THEN
      DO ibnd = 1, nbndsub
        WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        READ(iun_plot, *) nAtoms, rtempvec(1:3)
        !
        DO iline = 1, 3
          READ(iun_plot, '(8A)') temp_str
        ENDDO
        ie = 1
        DO iline = 1, nAtoms
          READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
          IF (iline == 1 ) THEN
            species(ie) = i_species
            ie = ie + 1
          ELSE IF (species(ie - 1) /= i_species) THEN
            species(ie) = i_species
            ie = ie + 1
          ENDIF
        ENDDO
        n_len_z = grid_end(3) - grid_start(3) + 1
        !
        DO nxx = grid_start(1), grid_end(1)
          DO nyy = grid_start(2), grid_end(2)
            DO nzz = grid_start(3), grid_end(3), 6
              IF (grid_end(3) - nzz < 6) THEN
                READ(iun_plot, *) wann_func(nxx, nyy, nzz:grid_end(3) - 1, ibnd)
              ELSE
                READ(iun_plot, '(6E13.5)') wann_func(nxx, nyy, nzz:nzz + 5, ibnd)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iun_plot)
        ! Wannier function is not well normalized
        ! Normalize here will make the calculations with Wannier functions easier
        norm = SUM(wann_func(:, :, :, ibnd) * wann_func(:, :, :, ibnd))
        wann_func(:, :, :, ibnd) = wann_func(:, :, :, ibnd) / SQRT(norm)
      ENDDO
    ENDIF
    CALL mp_bcast(wann_func, meta_ionode_id, world_comm)
    CALL mp_bcast(species, meta_ionode_id, world_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_wannier_cube
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_Rp_in_S()
    !-----------------------------------------------------------------------------------
    ! JLB
    !! Allocate and read list of Rp unit cell vectors contained on transformed supercell
    !-----------------------------------------------------------------------------------
    USE io_var,    ONLY : iunRpscell
    USE io_global, ONLY : ionode, meta_ionode_id
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    USE global_var,ONLY : nqtotf
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ierr
    !! Error status
    INTEGER :: nRp2
    !! Number of lattice vectors within supercell
    !
    IF (ionode) THEN
      OPEN(UNIT = iunRpscell, FILE = 'Rp.scell.plrn', FORM = 'formatted', STATUS = 'unknown')
      READ(iunRpscell, *) nRp
      IF (nRp /= nqtotf) CALL errore('read_Rp_in_S', 'nRp and nqtotf are not the same!',1)
      CLOSE(UNIT = iunRpscell)
    ENDIF
    CALL mp_bcast(nRp, meta_ionode_id, world_comm)
    ALLOCATE(Rp(3, nRp), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_Rp_in_S', 'Error allocating Rp', 1)
    Rp = 0
    IF (ionode) THEN
      OPEN(UNIT = iunRpscell, FILE = 'Rp.scell.plrn', FORM = 'formatted', STATUS = 'unknown')
      READ(iunRpscell, *) nRp2
      DO iRp = 1, nRp
        READ(iunRpscell, *) Rp(1:3, iRp)
      ENDDO
      CLOSE(UNIT = iunRpscell)
    ENDIF
    CALL mp_bcast(Rp, meta_ionode_id, world_comm)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_Rp_in_S

  END MODULE io_polaron
