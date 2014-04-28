!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine start_shirley (nodenumber)
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N]
  !
  !  Wrapper routine for shirley initialization
  !
  USE global_version, ONLY: version_number
  use environment, only : environment_start
  use mp_global, only : mp_startup, mpime
  USE mp_global,        ONLY : nimage, me_image, root_image
  USE image_io_routines,  ONLY : io_image_start
  implicit none
  character(len=3) :: nodenumber
  character(len=9) :: code = 'SHIRLEY'
  !
#ifdef __PARA
  CALL mp_startup ( )
  write(nodenumber,'(i3)') mpime
  ! reset IO nodes
  ! (do this to make each "image head node" an ionode)
  ! Has to be used ONLY to run nimage copies of pwscf
  !
  IF ( nimage > 1 ) CALL io_image_start( )
#endif
  CALL environment_start ( code )
  !
  return
end subroutine start_shirley
