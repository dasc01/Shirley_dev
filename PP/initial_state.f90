!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM initial_state
  !-----------------------------------------------------------------------
  !
  !  compute forces on atoms as a post-process
  !
  ! input: namelist "&inputpp", with variables
  !   prefix      prefix of input files saved by program pwscf
  !   outdir      temporary directory where files resides
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE kinds,      ONLY : DP
  USE io_files,   ONLY : prefix, tmp_dir, iunwfc, nwordwfc, trimcheck
  USE ions_base,  ONLY : nat
  USE klist,      ONLY : nks, xk
  USE wvfct,      ONLY : npw, igk
  USE uspp,       ONLY : nkb, vkb
  USE wavefunctions_module, ONLY : evc
  USE parameters, ONLY : ntypx
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start
  !
  IMPLICIT NONE
  CHARACTER(len=256) :: outdir
  INTEGER :: ios, ik, excite(ntypx)
  NAMELIST / inputpp / outdir, prefix, excite
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'initstate' )
  !
  !   set default values for variables in namelist
  !
  excite(:) = 0
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id )
  !
  IF ( ios /= 0) &
     CALL errore ('postforces', 'reading inputpp namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( excite, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  CALL hinit0
  CALL hinit1
  IF ( nks == 1 ) THEN
     ik = 1
     CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
     IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
  ENDIF

  !CALL sum_band
  !
  CALL do_initial_state (excite)
  !
  CALL stop_pp
  !

END PROGRAM initial_state
