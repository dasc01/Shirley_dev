!----------------------------------------------------------------------- 
PROGRAM change_rhogrid
  !----------------------------------------------------------------------- 
  ! 
  ! David Prendergast
  !
  ! Changes the FFT grid used to represent the charge density .rho
  !
#include "f_defs.h" 
  USE parameters,         ONLY : ntypx, npk, lmaxx, nchix, ndmx, nqfx, nbrx
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev 
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : degauss, ngauss, lgauss, nkstot
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir 
  USE noncollin_module, ONLY : noncolin
  USE mp,         ONLY : mp_bcast, mp_barrier      
  USE mp_global,  ONLY : me_pool, root_pool, my_pool_id, npool, mpime
  USE control_flags,        ONLY : lscf, twfcollect
  USE basis,                ONLY : startingpot, startingwfc
  use gvect, only: g, ngm, gcutm, ecutwfc
  USE fixed_occ,          ONLY : tfixed_occ
  !
  !use shirley_basis_input, only : get_input, &
  !                                prefix_input=>prefix, &
  !                                tmp_dir_input=>outdir, &
  !                                debug
  use pot_rho_G
  use wfc_shirley
  !
  real(dp) :: ecutwfc_new
  character(255) :: cecutwfc_new
  !
  !
  ! 
  CALL start_shirley (nd_nmbr) 
  !
  if( npool /= 1 ) then
    call errore('shirley_basis','number of pools should be 1',abs(npool))
  endif
  !
  ! read input and distribute
  !call get_input
  if( iargc() .ne. 3 ) &
    call errore('change_rhogrid','usage: change_rhogrid.x prefix tmp_dir cut_off',1)
  call getarg( 1, prefix )
  call getarg( 2, tmp_dir )
  call getarg( 3, cecutwfc_new )
  read(cecutwfc_new,*) ecutwfc_new
  write(stdout,*) ' generating rho using new cut-off ', ecutwfc_new
  !
  !prefix = prefix_input
  !tmp_dir = tmp_dir_input
  !
  ! this is a nonselfconsistent run
  lscf = .false.
  startingpot = 'file'
  !
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  call read_file_shirley
  !
  ! read the potential and density from file and transform to G-space
  call read_pot_rho_G

  ! redefine cut-off
  ecutwfc = ecutwfc_new

  call update_shirley_fft

  ! transform the potential and density from G-space and write to file
  write(stdout,*) ' dump charge and potential to file using updated FFT grid'
  call write_pot_rho_G
  !
  ! stop
  call stop_shirley

END PROGRAM change_rhogrid 
!
