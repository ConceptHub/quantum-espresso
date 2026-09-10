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
  MODULE polaron_common
  !--------------------------------------------------------------------------
  !!
  !! Shared state for the polaron modules: flags, dimensions, k/q maps, and the
  !! working arrays Hamil, eigvec, Bmat, epf. Data only; polaron.f90 allocates.
  !!
  USE kinds,     ONLY : DP
  
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: test_tags_plrn, mem_save_h, is_mirror_k, is_mirror_q, is_tri_k, is_tri_q
  PUBLIC :: nbnd_plrn, nbnd_g_plrn, lword_h, lword_g, lword_m
  PUBLIC :: io_level_g_plrn, io_level_h_plrn, hblocksize, band_pos, ik_edge, nRp
  PUBLIC :: Rp, select_bands_plrn, kpg_map, wq_model, etf_model, etf_all, xkf_all
  PUBLIC :: Hamil, eigvec, Bmat, gq_model, epf, epfall, berry_phase
  
  !
  LOGICAL :: test_tags_plrn(20) = .FALSE.
  !! The B matrix Bqu
  LOGICAL :: mem_save_h = .FALSE.
  !! The B matrix Bqu
  LOGICAL, ALLOCATABLE :: is_mirror_k(:)
  !! .true. if k is a mirror point, used for time-reversal symmetry
  LOGICAL, ALLOCATABLE :: is_mirror_q(:)
  !! .true. if q is a mirror point, used for time-reversal symmertry
  LOGICAL, ALLOCATABLE :: is_tri_k(:)
  !! .true. if k is a time-reversal invariant point, used for time-reversal symmetry
  LOGICAL, ALLOCATABLE :: is_tri_q(:)
  !! .true. if q is a time-reversal invariant point, used for time-reversal symmetry
  INTEGER :: nbnd_plrn
  !! Number of bands used in polaron calculations
  INTEGER :: nbnd_g_plrn
  !! Number of bands in which g is to be interpolated
  INTEGER :: lword_h
  !! Hamiltonian record size for I/O
  INTEGER :: lword_g
  !! el-ph matrix element record size for I/O
  INTEGER :: lword_m
  !! FIXME
  INTEGER :: io_level_g_plrn
  !! Write el-ph matrix elements to disk or store in memory
  INTEGER :: io_level_h_plrn
  !! Write Hamiltonian to disk or store in memory
  INTEGER :: hblocksize
  !! FIXME
  INTEGER :: band_pos
  !! Band in which CBM or VBM is located
  INTEGER :: ik_edge
  !! k-point in which CBM or VBM is located
  INTEGER :: nRp
  !! Number of unit cells on supercell for non-diagonal supercells
  INTEGER, ALLOCATABLE :: Rp(:,:)
  !! List of unit cell vectors within supercell
  INTEGER, ALLOCATABLE :: select_bands_plrn(:)
  !! Map from {start_band_plrn, end_band_plrn} to {1, end_band_plrn - start_band_plrn}
  INTEGER, ALLOCATABLE :: kpg_map(:)
  !! Map from a given k1 to its mirror point k2 = -k1 + G
  REAL(KIND = DP) :: wq_model
  !! Phonon freq in Frohlich model
  REAL(KIND = DP), ALLOCATABLE :: etf_model(:)
  !! Band structure in Frohlich model
  REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
  !! Gathered KS eigenvalues over the pools
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
  !! Gathered k-point coordinates over the pools
  COMPLEX(KIND = DP), ALLOCATABLE :: Hamil(:, :)
  !! Effective polaron Hamiltonian
  COMPLEX(KIND = DP), ALLOCATABLE :: eigvec(:, :)
  !! Polaron wave function coefficients in Bloch basis, Ank
  COMPLEX(KIND = DP), ALLOCATABLE :: Bmat(:,:)
  !! Polaron displacement coefficients in phono basis, Bqv
  COMPLEX(KIND = DP), ALLOCATABLE :: gq_model(:)
  !! el-ph matrix element in a simplified Fr\"ohlich model
  COMPLEX(KIND = DP), ALLOCATABLE :: epf(:, :, :, :)
  !! el-ph matrix element in for a given q
  COMPLEX(KIND = DP), ALLOCATABLE :: epfall(:, :, :, :, :)
  !! el-ph matrix element for all q
  COMPLEX(KIND = DP) :: berry_phase(1:3)
  !! FIXME
  !! Subroutine preparing variables for polaron calculations
  !! Subroutine selecting polaron scf or post-processing
  !

  END MODULE polaron_common
