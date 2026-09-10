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
  MODULE polaron_diagonalization
  !--------------------------------------------------------------------------
  !!
  !! Solves H A = E A. Three backends: diag_serial (LAPACK), diag_parallel
  !! (QE Davidson), diag_elpa (ELPA, built only under -D__ELPA).
  !!
  USE kinds,     ONLY : DP
  USE polaron_common, ONLY : nbnd_plrn, lword_h, hblocksize, Hamil
  USE polaron_grid,        ONLY : ikqLocal2Global
  USE polaron_hamiltonian, ONLY : h_psi_plrn, s_psi_plrn, g_psi_plrn
  USE buffers,        ONLY : get_buffer
  USE io_var,         ONLY : ihamil

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: diag_serial, diag_parallel
#if defined(__ELPA)
  PUBLIC :: diag_elpa
#endif

  CONTAINS

#if defined(__ELPA)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_elpa(estmteRt, eigvec_coef)
    !-----------------------------------------------------------------------
    !! This subroutine diagonalizes the full polaron Hamiltonian
    !! using ELPA library
    !! Implemented first version on 08/11/2025: S. Tiwari and K.F. Luo
    !-----------------------------------------------------------------------
    USE ep_constants,        ONLY : czero, twopi, ci, cone, zero, ryd2ev
    USE global_var,          ONLY : nkf, nktotf, Hamil_save, index_save
    USE input,               ONLY : nstate_plrn, &
                                    type_plrn, nhblock_plrn
    USE io_global,           ONLY : ionode, meta_ionode_id, stdout, ionode_id
    USE mp_world,            ONLY : world_comm
    USE mp_global,           ONLY : inter_pool_comm, my_pool_id, npool
    USE mp,                  ONLY : mp_sum, mp_bcast, mp_barrier
    USE control_flags,       ONLY : iverbosity
    USE elpa
    !
    IMPLICIT NONE
    !
    class(elpa_t), pointer :: elp
    !! Elpa API
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! Polaron eigenvalue
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron eigenvector coefficients
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: info
    !! FIXME
    INTEGER :: mm
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    INTEGER :: i, j, k, l, ipool
    !! Counters
    INTEGER :: tot, totx, toty, totn
    !! Total number of matrix elements
    INTEGER, EXTERNAL                    :: numroc, INDXG2L, INDXG2P, INDXL2G
    !! SCALAPACK variables
    REAL(KIND = c_double), ALLOCATABLE   :: ev(:), RWORK(:)
    !! Eigenvaues and RWORK for ZGEEV
    COMPLEX(KIND = c_double) :: inter
    !! Interaction strength e-p
    COMPLEX(KIND = c_double),ALLOCATABLE :: al(:,:), zl(:,:), WORK(:), &
                                            buff(:,:), H_f(:,:)
    !! Hamiltonian distributed over processors, eigenvectors, ZGEEV Work,
    !! Buffer array, full Hamiltonian
    !
    ! --------------------------------------------------------------------------------!!
    !! All variables are for diagonalization
    CHARACTER(len=8)                   :: task_suffix
    !! ELPA related task
    INTEGER                            :: success
    !! Success for ELPA diagonalization or not
    INTEGER                            :: nblk
    !! Block length for cyclic distribution
    INTEGER                            :: np_rows, np_cols, na_rows, na_cols
    !! Processor distribution in rows and columns
    INTEGER                            :: my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
    !! Present processor row and column, row/column communicators
    INTEGER                            :: my_blacs_ctxt, sc_desc(9), nprow,npcol,mpierr
    !! BLACS contexts, number of processor rows
    INTEGER                            :: il, jl, kl, kp, jp, jg, kg
    !! indices, l:local, g:global, p:processor
    INTEGER, ALLOCATABLE               :: buff_ind(:,:)
    !! BUffer array for redistribution
    INTEGER, ALLOCATABLE               :: na_rows_max(:), na_cols_max(:)
    !! maximum number of elements(row/column) in each processor
    INTEGER, ALLOCATABLE               :: loc(:)
    !! location in an array
    INTEGER                            :: STATUS
    !! Status if diagonalization
    INTEGER, PARAMETER                 :: error_units = 0
    !! Error unit
    INTEGER                            :: size_vec
    !! Size of the leading dimension of the largest eigenvector
    !!---------------------------------------------------------------------------------!!
    nblk=8
    tot = nktotf * nbnd_plrn
    !
    eigvec_coef = czero
    estmteRt = zero
    !
    IF (iverbosity == 5) THEN
       WRITE(stdout, '(/5x,a)'), '------------------------------'
       WRITE(stdout, '(5x,a)'), 'iverbosity = 5, printing more info ...'
       WRITE(stdout, '(5x,a)'), 'Using ELPA for diagonalization'
       WRITE(stdout, '(5x, a, 3I10)'), 'nkf, nktotf, nbnd_plrn', nkf, nktotf, nbnd_plrn
       WRITE(stdout, '(5x, a)'), 'WARNING: please ensure nkf (=nktotf/ncore) is integer'
    ENDIF
    !
    ALLOCATE(Hamil_save(nkf * nbnd_plrn, nktotf * nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating Hamil_save', 1)
    ALLOCATE(index_save(nkf * nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating index_save', 1)
    !
    ALLOCATE(na_rows_max(npool), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating na_rows_max(npool)', 1)
    ALLOCATE(na_cols_max(npool), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating na_cols_max(npool)', 1)
    !
    DO np_cols = NINT(SQRT(REAL(npool))), 2, -1
      IF (MOD(npool, np_cols) == 0 ) exit
    ENDDO
    ! at the end of the above loop, npools is always divisible by np_cols
    np_rows = npool / np_cols
    totn = tot
    ! initialise BLACS
    my_blacs_ctxt = inter_pool_comm
    !
    CALL BLACS_Gridinit(my_blacs_ctxt, 'C', np_rows, np_cols)
    CALL BLACS_Gridinfo(my_blacs_ctxt, nprow, npcol, my_prow, my_pcol)
    !
    !
    IF ((my_pool_id == ionode_id) .AND. (iverbosity == 5) ) THEN
      WRITE(stdout,'(/5x,a)'),'| Past BLACS_Gridinfo.'
    ENDIF
    ! determine the neccessary size of the distributed matrices,
    ! we use the scalapack tools routine NUMROC
    !
    na_rows = numroc(totn, nblk, my_prow, 0, np_rows)
    na_cols = numroc(totn, nblk, my_pcol, 0, np_cols)
    na_rows_max = 0
    na_cols_max = 0
    na_rows_max(my_pool_id + 1) = na_rows
    IF (na_cols > 0) THEN
      !
      na_cols_max(my_pool_id + 1) = na_cols
      !
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    CALL mp_sum(na_rows_max, inter_pool_comm)
    CALL mp_sum(na_cols_max, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    totx = INT(maxval(na_rows_max))
    toty = INT(maxval(na_cols_max))
    !
    IF ((my_pool_id == ionode_id) .AND. (iverbosity == 5))  THEN
      WRITE(stdout, '(/5x,a,3I10)'), 'na_rows, na_cols, my_prow', na_rows, na_cols, my_prow
      WRITE(stdout,'(5x,a,3I10)'), 'nprow, npcol, my_pcol', nprow, npcol, my_pcol
    ENDIF
    !
    ! set up the scalapack descriptor for the checks below
    ! For ELPA the following restrictions hold:
    ! - block sizes in both directions must be identical (args 4 a. 5)
    ! - first row and column of the distributed matrix must be on
    !   row/col 0/0 (arg 6 and 7)
    !
    CALL descinit(sc_desc, totn, totn, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info)
    !
    IF (info .NE. 0) THEN
      WRITE(error_units,*) 'Error in BLACS descinit! info=',info
      WRITE(error_units,*) 'Most likely this happend since you want to use'
      WRITE(error_units,*) 'more MPI tasks than are possible for your'
      WRITE(error_units,*) 'problem size (matrix size and blocksize)!'
      WRITE(error_units,*) 'The blacsgrid can not be set up properly'
      WRITE(error_units,*) 'Try reducing the number of MPI tasks...'
     ! call MPI_ABORT(inter_pool_comm, 1, mpierr)
    ENDIF
    !
    ! Error is not assigned for these because there are some issus
    ALLOCATE(al (na_rows, na_cols))
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating al,', 1)
    ALLOCATE(zl (na_rows, na_cols))
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating zl,', 1)
    ALLOCATE(ev (totn))
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error allocating ev,', 1)
    al=(0.d0, 0.d0)
    zl=(0.d0, 0.d0)
    !
    CALL mp_barrier(inter_pool_comm)
    IF (iverbosity == 5) WRITE(stdout, '(/5x,a)') 'Passed initial array allocation'
    !
    ! Distribute the Hamiltonian to different processors
    !
    DO ipool = 1,npool
      Hamil_save(:, :) = czero
      index_save(:) = 0
      l = 1
      !WRITE(stdout,'(/5x,a,I90)') 'Finished core: ', ipool
      IF ((MOD(ipool, 20) == 0) .AND. (iverbosity == 5)) THEN
        WRITE(stdout, '(5x, a, i10, a, i10)' ) 'Distributing Hamiltonian: ', ipool, '/', npool
      ENDIF
      !
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          !
          index_loc = MOD(indexkn1 - 1, hblocksize) + 1
          index_blk = INT((indexkn1 - 1) / hblocksize) + 1
          IF (index_loc == 1 .AND. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
          !
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          IF (my_pool_id + 1 == ipool) THEN
            index_save(l) = indexkn2
            Hamil_save(l, 1:nktotf * nbnd_plrn) = - type_plrn * Hamil(1:nktotf * nbnd_plrn, index_loc)
            l = l+1
          ENDIF
        ENDDO
      ENDDO
      CALL mp_sum(Hamil_save,inter_pool_comm)
      CALL mp_sum(index_save,inter_pool_comm)
      DO i = 1, nkf * nbnd_plrn
        DO j = 1, tot
          !
          jg = index_save(i)
          kg = j
          !
          inter = Hamil_save(i, j)
          jl = INDXG2L(jg, nblk, 0, 0, nprow) !FLOOR((j-1)/REAL(nprow))!
          jp = INDXG2P(jg, nblk, 0, 0, nprow)   !mod(j-1,nprow)
          !
          kl = INDXG2L(kg, nblk, 0, 0, npcol)!FLOOR((k-1)/REAL(npcol))
          kp = INDXG2P(kg, nblk, 0, 0, npcol) !mod(k-1,npcol)
          IF ((jp == my_prow) .AND. (kp == my_pcol)) THEN
            !
            IF ((jl <= na_rows) .AND. (kl <= na_cols)) THEN
              !
              al(jl, kl) = inter
              !
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_barrier(inter_pool_comm)
    !
    IF (iverbosity == 5) WRITE(stdout, '(5x,a)'),'Finished Hamiltonian setup for ELPA'
    !
    ! Start ELPA diagonalization
    !
    IF (elpa_init(20180501) /= elpa_ok) THEN
      WRITE(stdout, *), "ELPA API version not supported"
      STOP
    ENDIF
    elp => elpa_allocate()
    !
    ! set parameters decribing the matrix and it's MPI distribution
    CALL elp%set("na", totn, success)
    CALL elp%set("nev", totn, success)
    CALL elp%set("local_nrows", na_rows, success)
    CALL elp%set("local_ncols", na_cols, success)
    CALL elp%set("nblk", nblk, success)
    CALL elp%set("mpi_comm_parent", inter_pool_comm, success)
    CALL elp%set("process_row", my_prow, success)
    CALL elp%set("process_col", my_pcol, success)
    success = elp%setup()
    CALL elp%set("solver", elpa_solver_2stage, success)
    CALL elp%set("complex_kernel", elpa_solver_2stage, success)
    ! Calculate eigenvalues/eigenvectors
    !
    IF (iverbosity == 5) WRITE(stdout, '(/5x,a)')'| Entering one-step ELPA solver ... '
    !
    CALL mp_barrier(inter_pool_comm) ! for correct timings only
    CALL elp%eigenvectors(al, ev, zl, success)
    CALL mp_barrier(inter_pool_comm) ! for correct timings only
    !
    IF (iverbosity == 5) THEN
      WRITE(stdout, '(5x,a)') '| One-step ELPA solver complete.'
      WRITE(stdout, '(5x,a)'), 'The lowest 10 polaron eigenvalues (eV):'
      DO i = 1, 10
        WRITE(stdout,'(5x, i8, f12.6)') i, ev(i) * ryd2ev
      ENDDO
      WRITE(stdout, '(5x,a)'), '------------------------------'
    ENDIF
    !
    estmteRt(:) = ev(:)
    !
    DO jl = 1, totx
      IF ((MOD(jl, 500) == 0) .AND. (iverbosity == 5)) THEN
        WRITE(stdout,'(5x, a, I10, a, I10)'),'Redistributing eigenvectors: ', jl, '/' , totx
      ENDIF
      DO kl = 1, toty
        IF ((jl <= na_rows) .AND. (kl <= na_cols)) THEN
          !
          kg = INDXL2G(kl, nblk, my_pcol, 0, npcol)
          jg = INDXL2G(jl, nblk, my_prow, 0, nprow)
          IF ((kg > 0) .AND. (kg <= nstate_plrn)) THEN
            eigvec_coef(jg, kg) = zl(jl, kl)
          ENDIF
          !
        ENDIF
      ENDDO
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    CALL mp_sum(eigvec_coef, inter_pool_comm)
    !
    DEALLOCATE(al, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error deallocating al,', 1)
    DEALLOCATE(ev, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error deallocating ev,', 1)
    DEALLOCATE(zl, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error deallocating zl,', 1)
    DEALLOCATE(Hamil_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error deallocating Hamil_save,', 1)
    DEALLOCATE(index_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_elpa', 'Error deallocating index_save,', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE diag_elpa
#endif
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_serial(estmteRt, eigvec_coef)
    !-----------------------------------------------------------------------
    !! Serial diagonalization using LAPACK library
    !-----------------------------------------------------------------------
    USE ep_constants,        ONLY : czero, twopi, ci, cone, zero
    USE global_var,          ONLY : nkf, nktotf
    USE input,               ONLY : nstate_plrn, &
                                    type_plrn, nhblock_plrn
    USE io_global,           ONLY : ionode, meta_ionode_id
    USE mp_world,            ONLY : world_comm
    USE mp_global,           ONLY : inter_pool_comm
    USE mp,                  ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! Polaron eigenvalue
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron eigenvector coefficients
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: lwork
    !! FIXME
    INTEGER :: info
    !! FIXME
    INTEGER :: mm
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    INTEGER, ALLOCATABLE :: iwork(:)
    !! FIXME
    INTEGER, ALLOCATABLE :: ifail(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: work(:)
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: Hamil_save(:,:)
    !! Auxiliary array for storing Hamiltonian for all k-points
    COMPLEX(KIND = DP),  ALLOCATABLE :: Identity(:,:)
    !! FIXME
    !
    ALLOCATE(Hamil_save(nktotf * nbnd_plrn, nktotf * nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating Hamil_save', 1)
    !
    Hamil_save = czero
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        !
        index_loc = MOD(indexkn1 - 1, hblocksize) + 1
        index_blk = INT((indexkn1 - 1) / hblocksize) + 1
        IF (index_loc == 1 .AND. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
        !
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        Hamil_save(indexkn2, 1:nktotf * nbnd_plrn) = - type_plrn * Hamil(1:nktotf * nbnd_plrn, index_loc)
      ENDDO
    ENDDO
    CALL mp_sum(Hamil_save, inter_pool_comm)
    IF (ionode) THEN
      ALLOCATE(Identity(nktotf * nbnd_plrn, nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating Identity', 1)
      lwork = 5 * nktotf * nbnd_plrn
      ALLOCATE(rwork(7 * nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating rwork', 1)
      ALLOCATE(iwork(5 * nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating iwork', 1)
      ALLOCATE(ifail(nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating ifail', 1)
      ALLOCATE(work(lwork), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating work', 1)
      Identity = czero
      !
      DO ibnd = 1, nbnd_plrn * nktotf
        Identity(ibnd, ibnd) = cone
      ENDDO
      !
      eigvec_coef = czero
      estmteRt = zero
      !
      ! TODO: check out what is mm
      CALL ZHEGVX( 1, 'V', 'I', 'U', nktotf * nbnd_plrn, Hamil_save, nktotf * nbnd_plrn, Identity,&
         nktotf * nbnd_plrn, zero, zero, 1, nstate_plrn, zero, mm, estmteRt(1:nstate_plrn), &
         eigvec_coef, nktotf * nbnd_plrn, work, lwork, rwork, iwork, ifail, info)
      !
      IF (info /= 0) CALL errore('diag_serial','Polaron: diagonal error.', 1)
      DEALLOCATE(rwork, iwork, ifail, work, Identity, STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating rwork,', 1)
      !
    ENDIF
    !
    DEALLOCATE(Hamil_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating Hamil_save,', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE diag_serial
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_parallel(estmteRt, eigvec_coef)
    !-----------------------------------------------------------------------
    !! Parallel diagonalization using Davidson library from QE
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, eps5, eps6, eps4, eps2, eps8, eps10
    USE global_var,    ONLY : nkf, nktotf
    USE input,         ONLY : nstate_plrn, ethrdg_plrn, &
                              adapt_ethrdg_plrn, init_ethrdg_plrn, nethrdg_plrn, &
                              david_ndim_plrn
    USE io_global,     ONLY : stdout, ionode
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum, mp_bcast, mp_size, mp_max
    USE mp_bands,      ONLY : inter_bgrp_comm, mp_start_bands
    USE mp_bands_util, ONLY : intra_bgrp_comm_ => intra_bgrp_comm, &
                              inter_bgrp_comm_ => inter_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! Polaron eigenvalue
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron eigenvector coefficients
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER  :: ik
    !! k-point counter
    INTEGER  :: ik_global
    !! Global k-point index
    INTEGER  :: ibnd
    !! Electron band counter
    INTEGER  :: itemp
    !! FIXME
    INTEGER  :: jtemp
    !! FIXME
    INTEGER  :: indexkn1
    !! Combined band and k-point index
    INTEGER  :: indexkn2
    !! Combined band and k-point index
    INTEGER  :: ithr
    !! Counter for incremental threshold
    INTEGER  :: nthr
    !! Number of incremental threshold steps
    INTEGER  :: npw
    !! FIXME
    INTEGER  :: npwx
    !! FIXME
    INTEGER  :: dav_iter
    !! FIXME
    INTEGER  :: notcnv
    !! Number of non-converged eigenvalues
    INTEGER  :: btype(nstate_plrn)
    !! FIXME
    INTEGER  :: nhpsi
    !! FIXME
    REAL(KIND = DP) :: ethrdg_init
    !! Initial incremental threshold
    REAL(KIND = DP) :: ethrdg
    !! Threshold for convergence in diagonalization
    COMPLEX(KIND = DP), ALLOCATABLE :: psi(:, :)
    !! Eigenvector
    !
    npw = nkf * nbnd_plrn
    npwx = npw
    CALL mp_max(npwx, inter_pool_comm)
    !
    ALLOCATE(psi(1:npwx, 1:nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_parallel', 'Error allocating psi', 1)
    !
    ! JLB: Option for adaptive threshold
    IF (adapt_ethrdg_plrn) THEN
      ethrdg_init = init_ethrdg_plrn
      nthr = nethrdg_plrn
      IF(ionode) THEN
        WRITE(stdout, "(a)") "     Adaptive threshold on iterative diagonalization activated:"
        WRITE(stdout, "(a)") "     threshold, # iterations, eigenvalue(Ry)"
      ENDIF
    ELSE
      nthr = 1
    ENDIF
    !
    DO ithr = 1, nthr
      psi = czero
      btype(1:nstate_plrn) = 1
      !
      IF (adapt_ethrdg_plrn) THEN
        ethrdg = 10**(LOG10(ethrdg_init) + (ithr - 1) * (LOG10(ethrdg_plrn) - LOG10(ethrdg_init)) /(nthr - 1))
      ELSE
        ethrdg = ethrdg_plrn
      ENDIF
      !
      ! split eigvector (nqtotf) into parallel pieces psi (nkf), contains corresponding part with Hpsi
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          psi(indexkn1, 1:nstate_plrn) = eigvec_coef(indexkn2, 1:nstate_plrn)
        ENDDO
      ENDDO
      ! inter_bgrp_comm should be some non-existing number,
      ! to make the nodes in bgrp equal to 1
      ! intra_bgrp_comm is parallel PW in pwscf
      ! but here it should be parallel K.
      ! Save them before change them
      itemp = intra_bgrp_comm_
      jtemp = inter_bgrp_comm_
      !
      intra_bgrp_comm_ = inter_pool_comm
      inter_bgrp_comm_ = inter_bgrp_comm
      !
      CALL start_clock('cegterg_prln')
      CALL cegterg( h_psi_plrn, s_psi_plrn, .FALSE., g_psi_plrn, &
        !npw, npwx, nstate_plrn, nstate_plrn * david_ndim_plrn, 1, psi, ethrdg, &
        npw, npwx, nstate_plrn, david_ndim_plrn, 1, psi, ethrdg, &
        estmteRt, btype, notcnv, .FALSE., dav_iter, nhpsi)
      CALL start_clock('cegterg_prln')
      IF(adapt_ethrdg_plrn .AND. ionode) WRITE(stdout, "(a, E14.6, I6, E14.6)") "   ", ethrdg, dav_iter, estmteRt(1)
      IF(notcnv > 0 .AND. ionode) WRITE(stdout, "(a)") "   WARNING: Some eigenvalues not converged, &
      &check initialization, ethrdg_plrn or try adapt_ethrdg_plrn"
      !
      intra_bgrp_comm_ = itemp
      inter_bgrp_comm_ = jtemp
      !
      eigvec_coef = czero
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          eigvec_coef(indexkn2, 1:nstate_plrn) = psi(indexkn1, 1:nstate_plrn)
        ENDDO
      ENDDO
      CALL mp_sum(eigvec_coef, inter_pool_comm)
      !
    ENDDO
    !
    DEALLOCATE(psi, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_parallel', 'Error deallocating psi', 1)
    !------------------------------------------------------------------------
    END SUBROUTINE diag_parallel

  END MODULE polaron_diagonalization
