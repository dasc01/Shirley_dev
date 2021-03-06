!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE ioneb()
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : autoev, eV_to_kelvin, pi, rytoev, &
                            uakbar, amconv, bohr_radius_angs, eps8
  USE io_global,     ONLY : stdout
  !
  !
  USE io_files,      ONLY : tmp_dir 
  !
  USE path_variables, ONLY : lsteep_des, lquick_min, &
                             lbroyden, lbroyden2, &
                             llangevin, &
                             nstep_path_      => nstep_path, &    
                             ds_              => ds, &
                             use_masses_      => use_masses, &
                             CI_scheme_       => CI_scheme, &
                             fixed_tan_       => fixed_tan, &
                             use_freezing_    => use_freezing, &
                             k_max_           => k_max, &
                             k_min_           => k_min, &
                             num_of_images_   => num_of_images, &
                             first_last_opt_  => first_last_opt, &
                             temp_req_        => temp_req, &
                             path_thr_        => path_thr
  !
  !
  USE path_variables, ONLY : lneb, lsmd
  !
  USE path_input_parameters_module, ONLY : restart_mode
  !
  USE path_variables, ONLY : restart
  !                             
  USE path_input_parameters_module, ONLY : nstep_path, string_method, &
                               num_of_images, path_thr, CI_scheme, opt_scheme, &
                               use_masses, first_last_opt, temp_req, k_max,    &
                               k_min, ds, use_freezing, fixed_tan
  !
  !
  ! ... "path" specific
  !
  !
  IMPLICIT NONE
  !
  INTEGER  :: ia, image, nt
  REAL(DP) :: theta, phi
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  !
  !
  SELECT CASE(trim( string_method ))
     !
  CASE( 'neb' )
     !
     lneb  = .true.
     !
  CASE( 'smd' )
     !
     lsmd  = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'ioneb', 'string_method ' // &
                & trim( string_method ) // ' not implemented', 1 )
     !
  END SELECT
  !
  SELECT CASE( trim( restart_mode ) )
  CASE( 'from_scratch' )
     !
     restart        = .false.
     !
  CASE( 'restart' )
     !
     IF ( lneb .or. lsmd ) THEN
        !
        ! ... "path" specific
        !
        restart = .true.
        !
     ENDIF
     !
  CASE DEFAULT
     !
     CALL errore( 'ioneb', &
                & 'unknown restart_mode ' // trim( restart_mode ), 1 )
     !
  END SELECT
  !
     !
!
! check da mettere dopo iosys del pw
!
!     IF( io_level < 0) CALL errore ( 'ioneb', &
!                       'NEB, SMD do not work with "disk_io" set to "none"', 1)
     !
     !
     IF ( num_of_images < 2 ) &
        CALL errore( 'ioneb', 'string_method=' // trim( string_method ) // &
                   & ': num_of_images must be at least 2', 1 )
     !
     IF ( ( CI_scheme /= "no-CI"  ) .and. &
          ( CI_scheme /= "auto"   ) .and. &
          ( CI_scheme /= "manual" ) ) THEN
        !
        CALL errore( 'ioneb', 'string_method=' // trim( string_method ) // &
                   & ': unknown CI_scheme', 1 )
        !
     ENDIF
     !
     ! ... initialization of logical variables
     !
     lsteep_des  = .false.
     lquick_min  = .false.
     lbroyden    = .false.
     lbroyden2   = .false.
     !
     SELECT CASE( opt_scheme )
     CASE( "sd" )
        !
        lsteep_des = .true.
        !
     CASE( "quick-min" )
        !
        lquick_min = .true.
        !
     CASE( "broyden" )
        !
        lbroyden     = .true.
        !
     CASE( "broyden2" )
        !
        lbroyden2    = .true.
        !
     CASE( "langevin" )
        !
        llangevin = .true.
        !
        IF ( lneb ) &
           CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                      & ': langevin dynamics not implemented', 1 )
        !
        temp_req = temp_req / ( eV_to_kelvin * autoev )
        !
        IF ( temp_req <= 0.D0 ) &
           CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                      & ': tepm_req has not been set', 1 )
        !
        IF ( use_freezing ) &
           WRITE( UNIT = stdout, &
                  FMT = '(5X,"warning: freezing cannot be used in langevin")' )
        !
        use_freezing = .false.
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                   & ': unknown opt_scheme', 1 )
        !
     END SELECT
     !
  !
  ! ... "path"-optimization variables
  !
  nstep_path_     = nstep_path
  ds_             = ds
  num_of_images_  = num_of_images
  first_last_opt_ = first_last_opt
  use_masses_     = use_masses
  use_freezing_   = use_freezing
  temp_req_       = temp_req
  path_thr_       = path_thr
  CI_scheme_      = CI_scheme
  k_max_          = k_max
  k_min_          = k_min
  fixed_tan_      = fixed_tan
  !
  !
  !
  CALL verify_neb_tmpdir( tmp_dir )
  !
  !
  RETURN
  !
END SUBROUTINE ioneb
!
!-----------------------------------------------------------------------
SUBROUTINE verify_neb_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE wrappers,         ONLY : f_mkdir
  USE path_input_parameters_module, ONLY : restart_mode
  USE io_files,         ONLY : prefix, xmlpun, &
                               delete_if_present
  USE path_variables,   ONLY : num_of_images
  USE mp_image_global_module,        ONLY : mpime, nproc, nimage
  USE io_global,        ONLY : meta_ionode
  USE mp,               ONLY : mp_barrier
  USE xml_io_base,      ONLY : copy_file
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(inout) :: tmp_dir
  !
  INTEGER             :: ios, image, proc, nofi
  LOGICAL             :: exst
  CHARACTER (len=256) :: file_path, filename
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !
  file_path = trim( tmp_dir ) // trim( prefix )
  !
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     ! ... let us try to create the scratch directory
     !
     CALL parallel_mkdir ( tmp_dir )
     !
  ENDIF
  !
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )
  !
  IF ( meta_ionode ) THEN
     !
     ! ... files needed by parallelization among images are removed
     !
     CALL delete_if_present( trim( file_path ) // '.newimage' )
     !
     ! ... file containing the broyden's history
     !
     IF ( restart_mode == 'from_scratch' ) THEN
        !
        CALL delete_if_present( trim( file_path ) // '.broyden' )
        !
     ENDIF
     !
  ENDIF ! end if ionode
  !
  nofi = num_of_images
  !
  DO image = 1, nofi
     !
     file_path = trim( tmp_dir ) // trim( prefix ) //"_" // &
                 trim( int_to_char( image ) ) // '/'
     !
     CALL parallel_mkdir ( file_path )
     !
     ! ... if starting from scratch all temporary files are removed
     ! ... from tmp_dir ( by all the cpus in sequence )
     !
     IF ( restart_mode == 'from_scratch' ) THEN
        !
        DO proc = 0, nproc - 1
           !
           IF ( proc == mpime ) THEN
              !
              ! ... extrapolation file is removed
              !
              CALL delete_if_present( trim( file_path ) // &
                                    & trim( prefix ) // '.update' )
              !
              ! ... standard output of the self-consistency is removed
              !
              CALL delete_if_present( trim( file_path ) // 'PW.out' )
              !
           ENDIF
           !
           CALL mp_barrier()
           !
        ENDDO
        !
     ENDIF ! end restart_mode
     !
  ENDDO ! end do image
  !
  RETURN
  !
END SUBROUTINE verify_neb_tmpdir
