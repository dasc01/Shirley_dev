!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM diag_epsilon
  !----------------------------------------------------------------------- 
  ! 
  ! David Prendergast
  ! UCB, Feb 2007
  !
  ! Diagonalizes the dielectric matrix (static) using Shirley Brillouin zone 
  ! interpolation of electronic wave functions 
  !
#include "f_defs.h" 
  USE parameters, ONLY : ntypx, npk, lmaxx, nchix, ndmx, nqfx, nbrx
  USE io_global,  ONLY : stdout, ionode
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir 
  USE mp_global,  ONLY : npool
  USE control_flags,        ONLY : lscf
  USE basis,                ONLY : startingpot, startingwfc
  !
  USE ions_base, ONLY : ntyp=>nsp
  use gvect
  use uspp, only : deeq
  use uspp_param, only : nh
  !
  use shirley_epsilon_input, only : get_input
  use hamq_shirley, only : init_stdout, read_hamq
  use epsilon_shirley
  !
  !
  IMPLICIT NONE 
  character(255) :: hamqfile
  integer :: ios
  integer :: iunhq
  integer,external :: freeunit
  ! 
  CALL start_shirley (nd_nmbr) 
  !
  if( npool /= 1 ) then
    call errore('diag_epsilon','number of pools should be 1',abs(npool))
  endif
  !
  ! read input and distribute
  !call get_input

  call scala_zheev( 1000, 1000, 1000 )

  call stop_shirley
  ! 

END PROGRAM diag_epsilon
!
