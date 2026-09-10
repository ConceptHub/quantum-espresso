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
  MODULE polaron_interpolation
  !--------------------------------------------------------------------------
  !!
  !! Transforms Ank and Bqu between bases and grids, including supercell variants,
  !! plus the interp_plrn_wf / interp_plrn_bq post-processing drivers.
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : is_mirror_q, nbnd_plrn, nRp, select_bands_plrn, &
                             kpg_map, etf_all, eigvec, Bmat
  USE polaron_grid,       ONLY : ikqLocal2Global, index_Rp
  USE polaron_hamiltonian, ONLY : get_cfac
  USE io_polaron,     ONLY : read_plrn_wf_grid, read_plrn_wf, write_plrn_wf, &
                             read_plrn_dtau_grid, read_plrn_dtau, write_plrn_bmat

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: check_time_rev_sym, plrn_eigvec_tran, scell_plrn_eigvec_tran
  PUBLIC :: interp_plrn_wf, interp_plrn_bq, plrn_bmat_tran, scell_plrn_bmat_tran

  CONTAINS

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE check_time_rev_sym(eigvec_coef)
    !-----------------------------------------------------------------------
    !! Enforces TR symmetry on polaron wave function coefficients.
    !! Not used by default, only for testing purposes.
    !-----------------------------------------------------------------------
    USE global_var,     ONLY : nkf, nktotf
    USE ep_constants,   ONLY : czero, one, two, cone
    USE mp,             ONLY : mp_sum
    USE mp_global,      ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(inout) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    !! FIXME: nPlrn_l seems replicated with nstate_plrn
    INTEGER :: nPlrn_l
    !! Number of polaron states
    REAL(KIND = DP) :: norm
    !! Norm of polaron wave function
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_save(:, :)
    !! Auxiliary array for polaron wave function coefficients
    !
    ! nstate_plrn
    nPlrn_l = 1
    !
    ALLOCATE(eigvec_save(nktotf * nbnd_plrn, nPlrn_l), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error allocating eigvec_save', 1)
    eigvec_save = czero
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      ikpg = kpg_map(ik_global)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ikpg - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigvec_save(indexkn1, 1:nPlrn_l)  = CONJG(eigvec(indexkn2, 1:nPlrn_l))
      ENDDO
    ENDDO
    CALL mp_sum(eigvec_save, inter_pool_comm)
    eigvec_coef(:, 1:nPlrn_l) = (eigvec_coef(:, 1:nPlrn_l) + eigvec_save(:, 1:nPlrn_l))
    !
    DO iplrn = 1, nPlrn_l
      norm = REAL(DOT_PRODUCT(eigvec_coef(1:nbnd_plrn*nktotf, iplrn), eigvec_coef(1:nbnd_plrn*nktotf, iplrn)))!nktotf*nbnd_plrn*
      eigvec_coef(:, iplrn) = eigvec_coef(:, iplrn)/DSQRT(norm)
    ENDDO
    !
    DEALLOCATE(eigvec_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error deallocating eigvec_save', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE check_time_rev_sym
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    SUBROUTINE plrn_eigvec_tran(ttype, t_rev, eigvecin, nkf1_p, nkf2_p, nkf3_p, &
               nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigvecout, ip_center)
    !--------------------------------------------------------------------------------
    !! Fourier transform from eigvecin to eigvecout
    !! ttype is 'Bloch2Wan' or 'Wan2Bloch'
    !! Parallel version, each pool calculates its own k point set (nkf),
    !! then the mp_sum is used to sum over different pools.
    !! require the correct initialization of Rp_array
    !--------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, cone, two
    USE global_var,    ONLY : nkf, xkf, chw, nktotf
    USE input,         ONLY : nstate_plrn, nbndsub
    USE wannier2bloch, ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time reversal symmetry is to be imposed in the Ank coefficients
    INTEGER, INTENT(in) :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER, INTENT(in) :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER, INTENT(in) :: nkf3_p
    !! Fine k-point grid along b2
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! Center of polaron wave function, to shift supercell accordingly
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: eigvecout(:, :)
    !! Output wave function coefficients
    COMPLEX(KIND = DP), INTENT(in) :: eigvecin(:, :)
    !! Input wave function coefficients
    !
    ! Local variables
    LOGICAL :: is_mirror
    !! .true. if k-point is a time-reversal mirror point
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: ikglob
    !! Global inner k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: i_vec(3)
    !! Shifted lattice vector coordinates
    INTEGER :: center_shift(1:3)
    !! Shift coordinates
    INTEGER :: nkf_p(3)
    !! k-point coordinates for shift
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! Eigenvalues after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! U_{mn} matrices after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: cfac(nrr_k)
    !!! Exponential factor for Hamiltonian transformation
    !
    nkf_p(1:3) = (/nkf1_p, nkf2_p, nkf3_p/)
    IF (nbndsub_p /= nbndsub) CALL errore('plrnwfwan2bloch','Different bands included in last calculation!', 1)
    IF (ttype == 'Bloch2Wan') THEN
      itype =  1
    ELSE IF (ttype == 'Wan2Bloch') THEN
      itype = -1
    ELSE
      CALL errore('plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
    ENDIF
    !
    IF(PRESENT(ip_center)) THEN
      center_shift(1:3) = nkf_p / 2 - ip_center
    ELSE
      center_shift(1:3) = 0
    ENDIF
    !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
    !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
    !! ibnd -> m, jbnd -> n
    !! R_p from 1 to nkf1/2/3_p, note that loop in the sequence of ix, iy, and iz,
    !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
    eigvecout = czero
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      ik_global = ikqLocal2Global(ik, nktotf)
      !
      CALL get_cfac(xxk, nrr_k, irvec_r, cfac)
      !
      ! Only pass is_mirror where it has a meaning. hamwan2bloch guards every
      ! use with PRESENT, so omitting it is equivalent to passing .FALSE., and
      ! it leaves no local that could be read before it is assigned.
      IF (t_rev) THEN
        ikpg = kpg_map(ik_global)
        is_mirror = (ik_global > ikpg)
        CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
           etf_tmp, chw, cfac, is_mirror)
      ELSE
        CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
           etf_tmp, chw, cfac)
      END IF
      !
      IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
      DO iplrn = 1, nstate_plrn
        !ikglob = 0
        !! loop over all Wannier position p
        IF (nkf1_p == 0 .OR. nkf2_p == 0 .OR. nkf3_p == 0) THEN
          CALL errore('plrn_eigvec_tran','Wrong k grid, use nkf1/2/3 to give k grid!', 1)
        ENDIF
        DO ikglob = 1, nkf1_p * nkf2_p * nkf3_p
          i_vec(1:3) = MODULO(index_Rp(ikglob, nkf_p) + center_shift, nkf_p)
          ctemp = EXP(CMPLX(0.0_DP, twopi * DOT_PRODUCT(xxk, i_vec), KIND = DP))
          DO ibnd = 1, nbndsub_p ! loop over all Wannier state m
            DO jbnd = 1, nbnd_plrn ! loop over all Bloch state n
              indexkn1 = (ikglob - 1) * nbndsub + ibnd !mp
              indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd !nk
              SELECT CASE(itype)
                CASE(1)  ! Bloch2Wan !
                   eigvecout(indexkn1, iplrn) = eigvecout(indexkn1, iplrn) + &
                      eigvecin(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd)) !JLB: Conjugate transpose taken above!
                CASE(-1) ! Wan2Bloch !
                   eigvecout(indexkn2, iplrn) = eigvecout(indexkn2, iplrn) + &
                      eigvecin(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd) !JLB
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO ! ikglob
      ENDDO ! iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigvecout, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_plrn_eigvec_tran(ttype, t_rev, eigvecin, nktotf_p, nRp_p, Rp_p, &
                                nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigvecout)
    !-----------------------------------------------------------------------------------
    !! JLB: Fourier transform for non-diagonal supercells
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, cone, two
    USE global_var,    ONLY : nkf, xkf, chw, nktotf
    USE input,         ONLY : nstate_plrn, nbndsub
    USE wannier2bloch, ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed
    INTEGER, INTENT(in) :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Lattice vector coefficients in supercell
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands in polaron expansion
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: eigvecout(:, :)
    !! Output wave function coefficients
    COMPLEX(KIND = DP), INTENT(in) :: eigvecin(:, :)
    !! Input wave function coefficients
    !
    ! Local Variables
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! Eigenvalues after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! U_{mn} matrices after interpolated KS Hamiltonian diagonalization
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    LOGICAL :: is_mirror
    !! .true. if k-point is time-reversal mirror point
    INTEGER :: ierr
    !! Error status
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:, :, :)
    !! Exponential prefactor
    !
    IF (nbndsub_p /= nbndsub) CALL errore('scell_plrn_eigvec_tran','Different bands included in last calculation!',1)
    IF (ttype == 'Bloch2Wan') THEN
      itype =  1
    ELSEIF (ttype == 'Wan2Bloch') THEN
      itype = -1
    ELSE
      CALL errore('scell_plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
    ENDIF
    !
    ALLOCATE(cfac(nrr_k, dims, dims), STAT = ierr)
    IF(ierr /= 0) CALL errore('scell_plrn_eigvec_tran', 'Error allocating cfac', 1)
    !
    !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
    !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
    !! ibnd -> m, jbnd -> n
    !! R_p from 1 to nktotf_p
    !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
    eigvecout = czero
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      ik_global = ikqLocal2Global(ik, nktotf)
      CALL get_cfac(xxk, nrr_k, irvec_r, cfac)
      ! See plrn_eigvec_tran: pass is_mirror only where it is defined.
      IF (t_rev) THEN
        ikpg = kpg_map(ik_global)
        is_mirror = (ik_global > ikpg)
        CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
           etf_tmp, chw, cfac, is_mirror)
      ELSE
        CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
           etf_tmp, chw, cfac)
      ENDIF
      IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
      !
      DO iplrn = 1, nstate_plrn
        !icount = 0
        !! loop over all Wannier position p
        DO iRp = 1, nRp_p
          ctemp = EXP(twopi * ci * DOT_PRODUCT(xxk, Rp_p(1:3, iRp)))
          DO ibnd = 1, nbndsub_p ! loop over all Wannier state m
            DO jbnd = 1, nbnd_plrn ! loop over all Bloch state n
              indexkn1 = (iRp - 1) * nbndsub + ibnd !mp
              indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd !nk
              SELECT CASE(itype)
                CASE(1)  ! Bloch2Wan !
                   eigvecout(indexkn1, iplrn) = eigvecout(indexkn1, iplrn) + &
                      eigvecin(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd))
                CASE(-1) ! Wan2Bloch !
                   eigvecout(indexkn2, iplrn) = eigvecout(indexkn2, iplrn) + &
                      eigvecin(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd)
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO
      ENDDO !iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigvecout, inter_pool_comm)
    !
    DEALLOCATE(cfac, STAT = ierr)
    IF(ierr /= 0) CALL errore('scell_plrn_eigvec_tran', 'Error deallocating cfac', 1)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
    !-----------------------------------------------------------------------------------
    !! Interpolate polaron wave function coeffcients (Ank) and write to Ank.band.plrn.
    !! Mostly used to visualize contributions from different bands and k-points.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, ryd2ev, czero
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_bcast
    USE input,         ONLY : nstate_plrn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    !
    ! Local variables
    INTEGER :: ierr
    !! Error code when reading file
    INTEGER :: i_center(2)
    !! Index of polaron center lattice vector
    INTEGER :: ip_center(3)
    !! Coordinates of polaron center
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nbndsub_p
    !! Number of bands in polaron expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    !
    IF (ionode) WRITE(stdout, "(5x, a)") "Start of interpolation of electronic band structure."
    !
    ! read Amp.plrn, save eigvec_wan for the latter use
    IF(ionode) THEN
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
    i_center = MAXLOC(ABS(eigvec_wan))
    !
    ! i_center(1) is the flattened (ik, ibnd) index, recover ik from it.
    ip_center = index_Rp((i_center(1) - 1) / nbndsub_p + 1, (/nkf1_p, nkf2_p, nkf3_p/))
    WRITE(stdout, '(5x, a, i8, 3i5)') "The largest Amp ", i_center(1), ip_center
    WRITE(stdout, '(5x, a, i8)') "The number of polaron states in Amp.plrn: ", nplrn_p
    WRITE(stdout, '(5x, a, i8)') "The number of polaron states to be interpolated: ", nstate_plrn
    !
    ! JLB: kpg_map cannot be generally defined in interpolation k-paths,
    !      thus t_rev set to .false.
    CALL plrn_eigvec_tran('Wan2Bloch', .FALSE., eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
       nrr_k, ndegen_k, irvec_r, dims, eigvec, ip_center)
    !
    WRITE(stdout, '(5x, a)') "Polaron states have been interpolated to Bloch basis."
    !
    CALL write_plrn_wf(eigvec, 'Ank.band.plrn', etf_all)
    !
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF(ierr /= 0) CALL errore('interp_plrn_wf', 'Error deallocating eigvec_wan', 1)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
    !-----------------------------------------------------------------------------------
    !! Interpolate polaron displacements coefficients (Bqv) and write to Bmat.band.plrn.
    !! Mostly used to visualize contributions from each phonon mode and q-point.
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : wf, nqtotf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER,  INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    !
    ! Local variables
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variables read from file
    INTEGER :: nqf1_p
    !! Fine q-point grid along b1
    INTEGER :: nqf2_p
    !! Fine q-point grid along b2
    INTEGER :: nqf3_p
    !! Fine q-point grid along b3
    INTEGER :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER :: nmodes_p
    !! Number of phonon modes
    INTEGER :: ierr
    !! Error status
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ina
    !! Atom counter
    INTEGER :: i_center(2)
    !! Index of polaron center
    INTEGER :: ip_center(3)
    !! Coordinates of polaron center
    COMPLEX(KIND = DP), ALLOCATABLE :: bqv_coef(:,:)
    !! Polaron displacement coefficients in phonon basis, Bqv
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements in real space
    REAL(KIND = DP),    ALLOCATABLE :: dtau_r(:, :)
    !! Auxiliary polaron displacements to find polaron center
    !
    IF(ionode) THEN
      CALL read_plrn_dtau_grid(nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p, 'dtau.plrn')
    ENDIF
    CALL mp_bcast(nqf1_p,   meta_ionode_id, world_comm)
    CALL mp_bcast(nqf2_p,   meta_ionode_id, world_comm)
    CALL mp_bcast(nqf3_p,   meta_ionode_id, world_comm)
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
    ! Locate max displacement to center supercell
    ALLOCATE(dtau_r(nqtotf_p, nmodes/3), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating dtau_e', 1)
    dtau_r = czero
    DO iRp = 1, nqtotf_p
      DO ina = 1, nmodes / 3 ! ika -> kappa alpha
        dtau_r(iRp, ina) = NORM2(REAL(dtau(iRp, (ina - 1) * 3 + 1:ina * 3)))
      ENDDO
    ENDDO
    i_center = MAXLOC(ABS(dtau_r))
    ip_center = index_Rp(i_center(1), (/nqf1_p, nqf2_p, nqf3_p/))
    DEALLOCATE(dtau_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau_r', 1)
    !
    ALLOCATE(bqv_coef(nqtotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
    bqv_coef = czero
    !
    CALL plrn_bmat_tran('Dtau2Bmat', .false., dtau, nqf1_p, nqf2_p, nqf3_p, &
       nrr_q, ndegen_q, irvec_q, rws, nrws, bqv_coef, ip_center)
    !
    IF (ionode) CALL write_plrn_bmat(bqv_coef, 'Bmat.band.plrn', wf)
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau', 1)
    DEALLOCATE(bqv_coef, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating Bmat', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_bq
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE plrn_bmat_tran(ttype, t_rev, mat_in, nqf1_p, nqf2_p, nqf3_p, &
          nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out, ip_center, acoustic_plrn)
    !-----------------------------------------------------------------------------------
    !! Fourier transform between Bmat and dtau,
    !! Eq.(39) of PRB 99, 235139 (2019).
    !! Dtau2Bmat : B_{q\nu} = -1/N_p\sum_{\kappa\alpha p}C_{q\kappa \nu}\Delta\tau_{\kappa\alpha p}  e_{\kappa\alpha\nu}(q)\exp(iqR_p)
    !! Bmat2Dtau : \Delta \tau_{\kappa\alpha p} = -\sum_{q\nu} 1/(C_{q\kappa \nu}) B^*_{q\nu} e_{\kappa\alpha,\nu}(q) \exp(iqR_p)
    !! C_{q\kappa \nu} = N_p\left(\frac{M_k\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}} = N_p(M_k)^{\frac{1}{2}}D_{q\nu}
    !! D_{q \nu} = \left(\frac{\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}}
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : xqf, nqtotf, wf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc, type_plrn, eps_acoustic
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed in the Bqv coefficients
    INTEGER, INTENT(in) :: nqf1_p
    !! Fine q-point grid along b1
    INTEGER, INTENT(in) :: nqf2_p
    !! Fine q-point grid along b2
    INTEGER, INTENT(in) :: nqf3_p
    !! Fine q-point grid along b3
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! Coordinates of polaron center, for shifting supercell
    REAL(KIND = DP), INTENT(in), OPTIONAL :: acoustic_plrn
    !! the cutoff frequency of acoustic phonon modes in dtau.acoustic.plrn.xsf
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! Input matrix with Bqv/dtau coefficients
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! Output matrix with dtau/Bqv coefficients
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! .true. if q1 is TR mirror of another q2 point
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ika
    !! Combined atom and cartesian direction counter
    INTEGER :: ip_start
    !! Initial lattice vector in this pool
    INTEGER :: ip_end
    !! Final lattice vector in this pool
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: nqf_p(1:3)
    !! Fine q-point grid
    INTEGER :: ina
    !! Atom counter
    INTEGER :: nptotf
    !! Lattice vector counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER :: center_shift(1:3)
    !! Coordinates of shift to center supercell around polaron
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinate
    REAL(KIND = DP) :: xxq_r(3)
    !! q-point coordinate in case time-reversal has to be taken
    REAL(KIND = DP) :: ctemp
    !! Prefactor
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    COMPLEX(KIND = DP) :: dtemp
    !! Prefactor
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    nptotf = nqf1_p * nqf2_p * nqf3_p
    nqf_p(1:3) = (/nqf1_p, nqf2_p, nqf3_p/)
    !
    IF (nptotf <= 0) CALL errore('plrn_eigvec_tran', 'Use correct .plrn file with nqf1_p \= 0!', 1)
    IF (ttype == 'Bmat2Dtau') THEN
      itype =  1
    ELSE IF (ttype == 'Dtau2Bmat') THEN
      itype = -1
    ELSE
      CALL errore('plrn_eigvec_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
    ENDIF
    !
    uf = czero
    w2 = zero
    wf = zero
    !
    mat_out = czero
    !
    CALL fkbounds(nptotf, ip_start, ip_end)
    !
    DO iq = 1, nqtotf ! iq -> q
      xxq = xqf(1:3, iq)
      xxq_r = xxq(1:3)
      mirror_q = .false.
      ! if we need to force the time-rev symmetry, we have to ensure that the phase of uf is fixed
      ! i.e. uf = uf*(-q)
      IF (t_rev) THEN
        IF (is_mirror_q (iq)) THEN
          xxq_r = xqf(1:3, kpg_map(iq))
          mirror_q = .true.
        ENDIF
      ENDIF
      !
      ! Get phonon eigenmode and eigenfrequencies
      IF (.NOT. lifc) THEN
        ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
        ! parallel q is not working any more due to mp_sum in rgd_blk
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
      ENDIF
      !
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wf(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          wf(inu, iq) = 0.d0
        ENDIF
      ENDDO
      !
      IF (PRESENT(ip_center)) THEN
        center_shift(1:3) = nqf_p / 2 - ip_center
      ELSE
        center_shift(1:3) = 0
      ENDIF
      ! For mirror q, calculate the time-symmetric q' and get uf from q'
      ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
      !!IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
      DO inu = 1, nmodes ! inu -> nu
        IF (wf(inu, iq) < eps_acoustic) CYCLE !JLB - cycle zero and imaginary frequency modes
        IF (PRESENT(acoustic_plrn)) THEN
          IF (wf(inu, iq) > acoustic_plrn) CYCLE
          ! KL - include only modes with frequency lower than acoustic_plrn
        ENDIF
        DO ika = 1, nmodes ! ika -> kappa alpha
          ina = (ika - 1) / 3 + 1
          ctemp = DSQRT(two / (wf(inu, iq) * amass(ityp(ina))))
          ! Parallel run, only calculate the local cell ip
          ! Note that, ip_end obtained from fkbounds should be included
          ! If you have 19 kpts and 2 pool,
          ! lower_bnd= 1 and upper_bnd=10 for the first pool
          ! lower_bnd= 1 and upper_bnd=9 for the second pool
          DO iRp = ip_start, ip_end !, (nqf1_p + 1)/2
            Rp_vec(1:3) = MODULO(index_Rp(iRp, nqf_p) + center_shift, nqf_p)
            ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
            dtemp = uf(ika, inu) * EXP(CMPLX(0.0_DP, twopi * DOT_PRODUCT(xxq, Rp_vec), KIND = DP))
            IF (itype == 1) THEN ! Bqv -> dtau
              ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B^*_{q\nu}
              ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
              mat_out(iRp, ika) = mat_out(iRp, ika) -  cone / REAL(nptotf, dp) * dtemp * ctemp &
                 * (-type_plrn) * CONJG(mat_in(iq, inu))
            ELSE IF(itype == -1) THEN
              !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
              mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp / ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
              !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! sum all the cell index ip
    CALL mp_sum(mat_out, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE plrn_bmat_tran
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_plrn_bmat_tran(ttype, t_rev, mat_in, nqtotf_p, nRp_p, Rp_p, &
          nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out)
    !-----------------------------------------------------------------------------------
    !! JLB: Fourier transform between Bmat and dtau for non-diagonal supercells
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : xqf, wf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc, type_plrn, eps_acoustic
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed in Bqv coefficients
    INTEGER, INTENT(in) :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Lattice vector coordinates within supercell
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
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! Input matrix with Bqv/dtau coefficients
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! Output matrix with dtau/Bqv coefficients
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! .true. if q1 is a TR mirror point of another q2 point
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: itype
    !! Atom type counter
    INTEGER :: ika
    !! Combined atom and cartesian direction counter
    INTEGER :: ip_start
    !! Initial lattice vector within this pool
    INTEGER :: ip_end
    !! Final lattice vector within this pool
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ina
    !! Atom counter
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinate
    REAL(KIND = DP) :: xxq_r(3)
    !! auxiliary q-point coordinate in case TR is to be imposed
    REAL(KIND = DP) :: ctemp
    !! Prefactor
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    COMPLEX(KIND = DP) :: dtemp
    !! Prefactor
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    IF (ttype == 'Bmat2Dtau') THEN
      itype =  1
    ELSE IF (ttype == 'Dtau2Bmat') THEN
      itype = -1
    ELSE
      CALL errore('scell_plrn_bmat_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
    ENDIF
    !
    uf = czero
    w2 = zero
    wf = zero
    !
    mat_out = czero
    !
    CALL fkbounds(nRp_p, ip_start, ip_end)
    !
    DO iq = 1, nqtotf_p ! iq -> q
      xxq = xqf(1:3, iq)
      xxq_r = xxq(1:3)
      mirror_q = .false.
      ! if we need to force the time-rev symmetry, we have to ensure that the phase of uf is fixed
      ! i.e. uf = uf*(-q)
      IF (t_rev) THEN
        IF (is_mirror_q (iq)) THEN
          xxq_r = xqf(1:3, kpg_map(iq))
          mirror_q = .TRUE.
        ENDIF
      ENDIF
      !
      ! Get phonon eigenmode and eigenfrequencies
      IF (.NOT. lifc) THEN
        ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
        ! parallel q is not working any more due to mp_sum in rgd_blk
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
      ENDIF
      !
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wf(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          wf(inu, iq) = 0.d0
        ENDIF
      ENDDO
      ! For mirror q, calculate the time-symmetric q' and get uf from q'
      ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
      !IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
      DO inu = 1, nmodes ! inu -> nu
        IF (wf(inu, iq) < eps_acoustic) CYCLE !JLB - cycle zero and imaginary frequency modes
        DO ika = 1, nmodes ! ika -> kappa alpha
          ina = (ika - 1) / 3 + 1
          ctemp = DSQRT(two / (wf(inu, iq) * amass(ityp(ina))))
          !
          DO iRp = ip_start, ip_end
            ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
            dtemp = uf(ika, inu) * EXP( twopi * ci * DOT_PRODUCT(xxq(1:3), Rp_p(1:3, iRp)) )
            IF (itype == 1) THEN ! Bqv -> dtau
              ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B^*_{q\nu}
              ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
              mat_out(iRp, ika) = mat_out(iRp, ika) -  cone / REAL(nRp, dp) * dtemp * ctemp &
                 * (-type_plrn) * CONJG(mat_in(iq, inu))
            ELSE IF(itype == -1) THEN
              !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
              mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp / ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
              !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! sum all the cell index ip
    CALL mp_sum(mat_out, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_plrn_bmat_tran

  END MODULE polaron_interpolation
