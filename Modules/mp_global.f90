!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
  ! ... Wrapper module, for compatibility. Contains a few "leftover" variables
  ! ... used for checks (all the *_file variables, read from data file),
  ! ... plus the routine mp_startup initializing MPI and the command line, 
  ! ... plus the routine mp_global_end stopping MPI.
  ! ... Do not use this module to reference variables (e.g. communicators)
  ! ... belonging to each of the various parallelization levels:
  ! ... use the specific modules instead. 
  ! ... PLEASE DO NOT ADD NEW STUFF TO THIS MODULE. Removing stuff is ok.
  !
  USE mp_world, ONLY: world_comm, mp_world_start, mp_world_end
  USE mp_images
  USE mp_pools
  USE mp_bands
  USE mp_bands_TDDFPT
  USE mp_exx
  USE mp_diag
  USE mp_orthopools
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... number of processors for the various groups: values read from file
  !
  INTEGER :: nproc_file = 1
  INTEGER :: nproc_image_file = 1
  INTEGER :: nproc_pool_file  = 1
  INTEGER :: nproc_ortho_file = 1
  INTEGER :: nproc_bgrp_file  = 1
  INTEGER :: ntask_groups_file= 1
  INTEGER :: nyfft_file= 1
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_startup ( my_world_comm, start_images, diag_in_band_group, what_band_group )
    !-----------------------------------------------------------------------
    ! ... This wrapper subroutine initializes all parallelization levels.
    ! ... If option with_images=.true., processes are organized into images,
    ! ... each performing a quasi-indipendent calculation, such as a point
    ! ..  in configuration space (NEB) or a phonon irrep (PHonon)
    ! ... Within each image processes are further subdivided into various
    ! ... groups and parallelization levels.
    ! ... IMPORTANT NOTICE 1: since the command line is read here, it may be
    ! ...                     convenient to call it in serial execution as well
    ! ... IMPORTANT NOTICE 2: most parallelization levels are initialized here 
    ! ...                     but they should be moved to a later stage
    !
    USE command_line_options, ONLY : get_command_line, &
        nimage_, npool_, ndiag_, nband_, ntg_, nyfft_
    USE parallel_include
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN), OPTIONAL :: my_world_comm
    LOGICAL, INTENT(IN), OPTIONAL :: start_images
    LOGICAL, INTENT(IN), OPTIONAL :: diag_in_band_group
    INTEGER, INTENT(IN), OPTIONAL :: what_band_group
    LOGICAL :: do_images
    LOGICAL :: do_diag_in_band
    INTEGER :: my_comm
    INTEGER :: what_band_group_
    LOGICAL :: do_distr_diag_inside_bgrp
    !
    my_comm = MPI_COMM_WORLD
    IF ( PRESENT(my_world_comm) ) my_comm = my_world_comm
    !
    what_band_group_ = 1
    IF( PRESENT( what_band_group ) ) THEN
       what_band_group_ = what_band_group
    END IF
    !
    CALL mp_world_start( my_comm )
    CALL get_command_line ( )
    !
    do_images = .FALSE.
    IF ( PRESENT(start_images) ) do_images = start_images
    IF ( do_images ) THEN
       CALL mp_start_images ( nimage_, world_comm )
    ELSE
       CALL mp_init_image ( world_comm  )
    END IF
    !
    CALL mp_start_pools ( npool_, intra_image_comm )
    ! Init orthopools is done during EXX bootstrap but,
    ! if they become more used, do it here:
    ! CALL mp_start_orthopools ( intra_image_comm )
    CALL mp_start_bands ( nband_, ntg_, nyfft_, intra_pool_comm )
    CALL mp_start_exx ( nband_, ntg_, intra_pool_comm )
    !
    do_diag_in_band = .FALSE.
    IF ( PRESENT(diag_in_band_group) ) do_diag_in_band = diag_in_band_group
    !
    IF( negrp.gt.1 .or. do_diag_in_band ) THEN
       ! used to be the default : one diag group per bgrp
       ! with strict hierarchy: POOL > BAND > DIAG
       ! if using exx groups from mp_exx still use this diag method
       my_comm = intra_bgrp_comm
    ELSE
       ! new default: one diag group per pool ( individual k-point level )
       ! with band group and diag group both being children of POOL comm
       my_comm = intra_pool_comm
    END IF
    do_distr_diag_inside_bgrp = (negrp.gt.1) .or. do_diag_in_band
    CALL mp_start_diag ( ndiag_, world_comm, my_comm, do_distr_diag_inside_bgrp )
    !
    call set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE mp_startup
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_end ( )
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY : mp_comm_free
    !
    CALL unset_mpi_comm_4_solvers()
    CALL mp_comm_free ( intra_bgrp_comm )
    CALL mp_comm_free ( inter_bgrp_comm )
    CALL mp_comm_free ( intra_pool_comm )
    CALL mp_comm_free ( inter_pool_comm )
    CALL mp_stop_orthopools( ) ! cleans orthopools if used in exx
    CALL mp_world_end( )
    !
    RETURN
    !
  END SUBROUTINE mp_global_end
  !
END MODULE mp_global
