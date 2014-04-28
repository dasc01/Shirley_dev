  module shirley_epsilon_input

  use kinds, only : dp
  use constants, only : rytoev, convert_E_to_temp

  implicit none

  character(255) :: outdir
  ! band subset to read and process
  integer :: band_subset(2)
  ! k-point grid and offset
  integer :: nkgrid(3), ikgrid(3)
  ! q-point grid and offset
  integer :: nqgrid(3), iqgrid(3)
  ! epsilon cut-off
  real(dp) :: ecuteps
  ! nstride
  integer :: nstride(3)
  ! block size
  integer :: block_size
  ! fermi energy (eV) and electron temperature (K)
  real(dp) :: efermi, etemp
  ! debug flag
  logical :: debug
  !
  !
  contains

  !----------------------------------------------------------------------- 
    subroutine get_input
  !----------------------------------------------------------------------- 
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir 
  USE mp,         ONLY : mp_bcast, mp_barrier      
  !
  INTEGER :: ios
  integer :: i
  ! 
  NAMELIST / input / outdir, prefix, band_subset, &
                     ecuteps, efermi, etemp, debug, &
                     nstride, block_size
  !
  !   set default values for variables in namelist 
  ! 
  outdir = './' 
  prefix = 'pwscf' 
  band_subset = 0
  ecuteps = 0.d0
  nstride = 1  ! default - no reduction in resolution
  block_size = 1  ! default - each proc only deals with a single r,r' pair
  efermi = 0.d0
  etemp = -1.d0
  debug = .false.
  !
  IF ( ionode )  THEN  
     !
     CALL input_from_file ( )
     !
     READ (5, input, err = 200, iostat = ios) 
200  CALL errore ('shirley_ham', 'reading input namelist', ABS (ios) ) 
     ! 
     tmp_dir = TRIM(outdir) 
     efermi = efermi / rytoev
     etemp  = etemp  / convert_E_to_temp ! K to Ryd
     !
     write(stdout,*)
     !
     ! read k-point list and shift integers
     read(5,*) ! heading
     read(5,*) nkgrid(1:3), ikgrid(1:3)
     write(stdout,*) 'k-point grid:'
     write(stdout,*) nkgrid
     write(stdout,*) ikgrid
     !
     ! read q-point list and shift integers
     read(5,*) ! heading
     read(5,*) nqgrid(1:3), iqgrid(1:3)
     write(stdout,*) 'q-point grid:'
     write(stdout,*) nqgrid
     write(stdout,*) iqgrid
     !
     !
  END IF 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir,  ionode_id ) 
  CALL mp_bcast( prefix,   ionode_id ) 
  CALL mp_bcast( ecuteps,   ionode_id ) 
  CALL mp_bcast( nstride,   ionode_id ) 
  CALL mp_bcast( block_size,   ionode_id ) 
  CALL mp_bcast( efermi,   ionode_id ) 
  CALL mp_bcast( etemp,   ionode_id ) 
  CALL mp_bcast( debug, ionode_id )
  CALL mp_bcast( band_subset, ionode_id )
  ! 
  CALL mp_bcast( nkgrid,   ionode_id )
  CALL mp_bcast( ikgrid,   ionode_id )
  !
  CALL mp_bcast( nqgrid,   ionode_id )
  CALL mp_bcast( iqgrid,   ionode_id )
  !
  return
  end subroutine get_input
! ---------------------------------------------------------------------- 

  end module shirley_epsilon_input

