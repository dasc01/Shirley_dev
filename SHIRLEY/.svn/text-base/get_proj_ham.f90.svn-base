!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE get_nproj_ham( local_channel, nkb_ham )
  !----------------------------------------------------------------------------
  !
  ! With knowledge of the local channel for each species,
  ! this routine generates the set of indices corresponding to
  ! non-local projectors only, for use in constructing the
  ! Hamiltonian
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : current_spin
  USE uspp,       ONLY : nkb, nhtol
  USE uspp_param, ONLY : nh
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  integer,intent(in) :: local_channel(ntyp)
  integer,intent(out) :: nkb_ham  
  !
  ! ... here the local variables
  !
  INTEGER :: ikb, ih, na, nt
  !
  write(stdout,*) ' nhtol = '
  write(stdout,'(4x,3i6)') ((ih, nt, nhtol(ih,nt), ih=1,nh(nt)), nt=1,ntyp)
  write(stdout,*) ' local_channel = '
  write(stdout,'(4x,2i6)') (nt, local_channel(nt), nt=1,ntyp)
  !
  nkb_ham=0
  !
  DO nt = 1, ntyp
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           DO ih = 1, nh(nt)
              if( nhtol(ih,nt) /= local_channel(nt) .or. local_channel(nt) < 0 ) then
                nkb_ham = nkb_ham + 1
              endif
           END DO
        END IF
     END DO
  END DO
  !
  RETURN
  !
END SUBROUTINE get_nproj_ham
!
!----------------------------------------------------------------------------
SUBROUTINE get_proj_ham( local_channel, nkb_ham, ikb_ham )
  !----------------------------------------------------------------------------
  !
  ! With knowledge of the local channel for each species,
  ! this routine generates the set of indices corresponding to
  ! non-local projectors only, for use in constructing the
  ! Hamiltonian
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : current_spin
  USE uspp,       ONLY : nkb, nhtol
  USE uspp_param, ONLY : nh
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  integer,intent(in) :: local_channel(ntyp)
  integer,intent(in) :: nkb_ham 
  integer,intent(out) :: ikb_ham(nkb_ham)
  !
  ! ... here the local variables
  !
  INTEGER :: ikb, ih, na, nt, nkb_ham_
  !
  write(stdout,*) ' nhtol = '
  write(stdout,'(4x,3i6)') ((ih, nt, nhtol(ih,nt), ih=1,nh(nt)), nt=1,ntyp)
  write(stdout,*) ' local_channel = '
  write(stdout,'(4x,2i6)') (nt, local_channel(nt), nt=1,ntyp)
  !
  nkb_ham_=0
  ikb=0
  !
  DO nt = 1, ntyp
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           DO ih = 1, nh(nt)
              ikb=ikb+1
              if( nhtol(ih,nt) /= local_channel(nt) .or. local_channel(nt) < 0 ) then
                nkb_ham_ = nkb_ham_ + 1
                ikb_ham(nkb_ham_) = ikb
              endif
           END DO
        END IF
     END DO
  END DO
  !
  if( ikb /= nkb ) call errore('get_proj_ham','wrong counting of number of projectors',1)
  !
  RETURN
  !
END SUBROUTINE get_proj_ham

