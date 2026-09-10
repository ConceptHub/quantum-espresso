! Copyright (C) 2016-2023 Marios Zacharias, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
! This module is identical to Module ifconstants in PHonon/PH/matdyn.f90.
! It is named ifconstants_zg (rather than ifconstants) to avoid a module-name
! collision with matdyn.f90: that file sits in EPW/ZG/src's dependency search
! path, so "make depend" would otherwise resolve "USE ifconstants" to matdyn.o
! instead of this file. It is USEd by both ZG.f90 and disca.f90, which are two
! separate programs in this directory; keeping the module in its own file
! (rather than duplicated inside each program) prevents them from concurrently
! writing the same .mod file during a parallel build.
!
! One should make sure this file is updated when PHonon/PH/matdyn.f90 is updated.
!
Module ifconstants_zg
  !
  !! All variables read from file that need dynamical allocation.
  !
  USE kinds, ONLY: DP
  !
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:)
  !! interatomic force constants in real space
  REAL(DP), ALLOCATABLE :: frc_lr(:,:,:,:,:,:,:)
  !! long-range part of interatomic force constants in real space
  REAL(DP), ALLOCATABLE :: tau_blk(:,:)
  !! atomic positions for the original cell
  REAL(DP), ALLOCATABLE :: zeu(:,:,:)
  !! effective charges for the original cell
  REAL(DP), ALLOCATABLE :: m_loc(:,:)
  !! the magnetic moments of each atom
  INTEGER, ALLOCATABLE  :: ityp_blk(:)
  !! atomic types for each atom of the original cell
  !
  CHARACTER(LEN=6), ALLOCATABLE :: atm(:)
  !
end Module ifconstants_zg
