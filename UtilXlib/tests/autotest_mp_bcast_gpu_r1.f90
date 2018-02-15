#if defined(__CUDA)
PROGRAM test_mp_bcast_r1_gpu
!
! Simple program to check the functionalities of test_mp_bcast_i1.
!
    USE cudafor
#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    REAL(8), DEVICE :: r1_d
    REAL(8) :: r1_h
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    r1_h = mpime
    r1_d = r1_h
    CALL mp_bcast(r1_d, root, world_comm)
    r1_h = r1_d
    !
    CALL test%assert_equal((r1_h .eq. 0) , .true. , fail=.true.)
    !
    r1_h = mpime
    r1_d = r1_h
    CALL mp_bcast(r1_d, nproc-1, world_comm)
    r1_h = r1_d
    !
    CALL test%assert_equal((r1_h .eq. nproc-1) , .true. , fail=.true.)
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_r1_gpu
#else
PROGRAM test_mp_bcast_r1_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_r1_gpu
#endif
