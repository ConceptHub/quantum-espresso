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
  MODULE polaron_grid
  !--------------------------------------------------------------------------
  !!
  !! Index arithmetic for polaron grids: k+q lookup, local/global maps, G-vector
  !! and Gamma tests, Rp indexing. Pure functions, no state.
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : nbnd_plrn, xkf_all

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ikq_all, isGVec, ikqLocal2Global, indexGamma
  PUBLIC :: find_ik, ikGlobal2Local, index_Rp, index_shift

  CONTAINS

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikq_all(ik, iq)
    !-----------------------------------------------------------------------
    !!
    !! find the global index of k+q for the local ik and global iq
    !!
    USE global_var, ONLY : nktotf
    USE input,      ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! k-point index
    INTEGER, INTENT(in) :: iq
    !! q-point index
    !
    ! Local variables
    INTEGER :: ikq
    ! k+q index
    INTEGER :: ik_global
    !! k-point index in global list
    INTEGER :: ikq_all
    !! k+q point index in global list
    INTEGER :: index_target(1:3)
    !! Auxiliary index
    INTEGER :: index_kq
    !! Auxiliary index of k+q point
    INTEGER :: ikq_loop
    !! Loop counter
    REAL(KIND = DP) :: xxk(1:3)
    !! k+q point coordinates
    REAL(KIND = DP) :: xxk_target(1:3)
    !! k+q point coordinates in 1BZ
    !
    ikq_all = 0
    !
    ik_global = ikqLocal2Global(ik, nktotf)
    xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)
    !
    xxk_target(1:3) = xxk(1:3) - INT(xxk(1:3))
    index_target(1:3) = NINT(xxk_target(1:3) * (/nkf1, nkf2, nkf3/))
    !
    index_kq = index_target(1) * nkf1 * nkf2 + index_target(2) * nkf2 + index_target(3) + 1
    !
    DO ikq_loop = index_kq - 1, nktotf + index_kq
      ! ik (local) + iq (global) = ikq (global)
      ! get ikq to locate the column of the Hamiltonian
      ikq = MOD(ikq_loop, nktotf) + 1
      IF (isGVec(xxk - xkf_all(1:3, ikq))) THEN
        ikq_all = ikq
        EXIT
      ENDIF
    ENDDO
    !
    IF (ikq_all == 0) CALL errore('ikq_all','k + q not found', 1)
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikq_all
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION find_ik(xxk, xkf_global)
    !-----------------------------------------------------------------------
    !!
    !! Find k-point index
    !!
    USE global_var, ONLY : nktotf
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(1:3)
    !! K-point position per cpu
    REAL(KIND = DP), INTENT(in) :: xkf_global(1:3, 1:nktotf)
    !! global k-points
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    INTEGER :: find_ik
    !! k-point index of found k-point
    REAL(KIND = DP) :: xkq(1:3)
    !! k-point coordinates
    !
    CALL start_clock('find_k')
    !
    find_ik = 0
    DO ik = 1, nktotf
      xkq(1:3) = xkf_global(1:3, ik) - xxk(1:3)
      IF(isGVec(xkq)) THEN
        find_ik = ik
        EXIT
      ENDIF
    ENDDO
    !
    IF (find_ik == 0) CALL errore('find_ik','k not found', 1)
    !
    CALL stop_clock('find_k')
    !
    !-----------------------------------------------------------------------
    END FUNCTION find_ik
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION isGVec(xxk)
    !-----------------------------------------------------------------------
    !! Return true if xxk integer times of the reciprocal vector
    !! if xxk is the difference of two vectors, then return true if these
    !! two vector are the same
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : eps6
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(3)
    !! k-point coordinate
    !
    ! Local variable
    LOGICAL :: isGVec
    !! .true. if k-point is a G-vector
    !
    isGVec = &
      ABS(xxk(1) - NINT(xxk(1))) < eps6 .AND. &
      ABS(xxk(2) - NINT(xxk(2))) < eps6 .AND. &
      ABS(xxk(3) - NINT(xxk(3))) < eps6
    !-----------------------------------------------------------------------
    END FUNCTION isGVec
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikqLocal2Global(ikq, nkqtotf)
    !-----------------------------------------------------------------------
    !! Return the global index of the local k point ik
    !-----------------------------------------------------------------------
    USE parallelism, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ikq
    !! k+q point counter
    INTEGER, INTENT(in) :: nkqtotf
    !! number of k-points in the fine grid
    !
    ! Local variable
    INTEGER :: ikqLocal2Global
    !! Index of k+q point in global list
    INTEGER :: startn
    !! Lower bound for k-points in pools
    INTEGER :: lastn
    !! Upper bound for k-points in pools
    !
    CALL fkbounds(nkqtotf, startn, lastn)
    !
    ikqLocal2Global = startn + ikq - 1
    IF (ikqLocal2Global > lastn) THEN
      CALL errore('ikqLocal2Global', 'Index of k/q is beyond this pool.', 1)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikqLocal2Global
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikGlobal2Local(ik_g, nktotf)
    !-----------------------------------------------------------------------
    !! Return the global index of the local k point ik
    !-----------------------------------------------------------------------
    USE parallelism, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik_g
    !! Global k-point index
    INTEGER, INTENT(in) :: nktotf
    !! Number of k-points in global fine grid
    !
    ! Local variable
    INTEGER :: ikGlobal2Local
    !! Index of k-point in local pool list
    INTEGER :: startn
    !! Lower bound for k-points in pools
    INTEGER :: lastn
    !! Upper bound for k-points in pools
    !
    CALL fkbounds(nktotf, startn, lastn)
    !
    ikGlobal2Local = ik_g - startn + 1
    !
    IF(ikGlobal2Local <= 0) THEN
      ikGlobal2Local = 0
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION ikGlobal2Local
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION indexGamma(k_all)
    !-----------------------------------------------------------------------
    !! Find the index of Gamma point i.e. (0, 0, 0) in xkf_all
    !! which contains all the crystal coordinates of the k/q points
    !! if Gamma point is not included, return 0
    !
    !-----------------------------------------------------------------------
    USE global_var,   ONLY : nkf, nktotf
    USE mp,           ONLY : mp_sum
    USE mp_global,    ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: k_all(:, :)
    !! crystal coordinates of k/q points.
    !! Renamed from xkf_all to avoid variable shadowing.
    !
    ! Local variable
    INTEGER :: indexGamma
    !! Index of \Gamma point in global k-point list
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    !
    indexGamma = 0
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      IF(isGVec(k_all(1:3, ik_global))) THEN
        indexGamma = ik_global
      ENDIF
    ENDDO
    CALL mp_sum(indexGamma, inter_pool_comm)
    !
    IF (.NOT. isGVec(k_all(1:3, indexGamma))) THEN
      CALL errore('indexGamma','The index of Gamma point is wrong!', 1)
    ENDIF
    !----------------------------------------------------------------------
    END FUNCTION indexGamma
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    FUNCTION index_Rp(iRp, nqfs)
    !-----------------------------------------------------------------------------------
    !! Index
    !-----------------------------------------------------------------------------------
    USE input,      ONLY: nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iRp
    !! Lattice vector index
    INTEGER, INTENT(in), OPTIONAL  :: nqfs(1:3)
    !! k/q points along each direcion in grid
    !
    ! Local variable
    INTEGER  :: index_Rp(1:3)
    !! Lattice vector in crystal coords
    INTEGER  :: nqf_c(1:3)
    !! k/q points along each direcion in grid
    !
    IF (PRESENT(nqfs)) THEN
      nqf_c(1:3) = nqfs(1:3)
    ELSE
      nqf_c(1:3) = (/nqf1, nqf2, nqf3/)
    ENDIF
    !
    index_Rp(1) = (iRp - 1)/(nqf_c(2) * nqf_c(3))
    index_Rp(2) = MOD(iRp - 1, nqf_c(2) * nqf_c(3))/nqf_c(3)
    index_Rp(3) = MOD(iRp - 1, nqf_c(3))
    !
    IF (ANY(index_Rp < 0) .OR. ANY(index_Rp >= nqf_c)) THEN
      CALL errore('index_Rp','index_Rp not correct!',1)
    ENDIF
    !-----------------------------------------------------------------------------------
    END FUNCTION index_Rp
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    FUNCTION index_shift(ishift)
    !-----------------------------------------------------------------------------------
    !! Find supercell lattice vector in crystal coords for shift loop around neighbors
    !-----------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ishift
    !! Shift vector index
    !
    ! Local variable
    INTEGER  :: index_shift(1:3)
    !! Shift vector crystal coords.
    !
    index_shift(1) = (ishift - 1) / 9 - 1
    index_shift(2) = MOD(ishift - 1, 9) / 3 - 1
    index_shift(3) = MOD(ishift - 1, 3) - 1
    !
    IF (ANY(index_shift < -1) .OR. ANY(index_shift > 1)) THEN
      CALL errore('index_shift', 'index_shift not correct!', 1)
    ENDIF
    !-----------------------------------------------------------------------------------
    END FUNCTION index_shift

  END MODULE polaron_grid
