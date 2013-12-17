!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM shirley_basis 
  !----------------------------------------------------------------------- 
  ! 
  ! David Prendergast
  ! UCB, Dec 2006
  !
  ! Generates the optimal basis set for Shirley Brillouin zone interpolation
  ! of electronic wave functions for use in the parametrized Shirley Hamiltonian
  !
#include "f_defs.h" 
  USE parameters,         ONLY : ntypx, npk, lmaxx
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev 
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : degauss, ngauss, lgauss, nkstot
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir 
  USE noncollin_module, ONLY : noncolin
  USE mp,         ONLY : mp_bcast, mp_barrier      
  USE mp_global,  ONLY : me_pool, root_pool, my_pool_id, npool, mpime
  USE control_flags,        ONLY : lscf, twfcollect
  USE basis,                ONLY : starting_pot, starting_wfc
  use gvect, only: g, ngm, gcutm
  use wvfct, only: ecutwfc
  USE fixed_occ,          ONLY : tfixed_occ
  use dfunct, only : newd
  use ldaU,   only : lda_plus_u, U_projection
  !
  use shirley_basis_input, only : get_input, &
                                  prefix_input=>prefix, &
                                  tmp_dir_input=>outdir, &
                                  debug
  !
  !
  !
  ! 
  CALL start_shirley (nd_nmbr) 
  !
  if( npool /= 1 ) then
    call errore('shirley_basis','number of pools should be 1',abs(npool))
  endif
  !
  IF ( ionode ) THEN
     !     
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
     !     
     WRITE( unit = stdout, FMT = 9010 ) & 
         ntypx, npk, lmaxx
     !     
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !

  END IF  
  !
  ! read input and distribute
  call get_input
  !
  if( debug .and. .not. ionode ) stdout = 500+mpime
  !
  prefix = prefix_input
  tmp_dir = tmp_dir_input
  !
  ! this is a nonselfconsistent run
  lscf = .false.
  starting_pot = 'file'
  starting_wfc = 'file'
  !
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  call read_file_shirley( 0 )
  !
  call openfil
  !
  CALL hinit0()
  CALL potinit()
  !
  CALL newd()
  !
  CALL wfcinit()
  !
  ! initialize basis set by GS orthogonalizing periodic functions over k
  !
  call build_optimal_basis()
  !
  ! stop
  call stop_shirley


END PROGRAM shirley_basis 
!
