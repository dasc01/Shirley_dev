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
SUBROUTINE get_dnl( dnl )
  !----------------------------------------------------------------------------
  !
  !    This routine applies the sqrt of the atomic term D
  !    to the projections on each band, where
  !      Vnl = sum_ikb,jkb | beta(ikb) > D(ikb,jkb) < beta(jkb) |
  !    If we want to apply this 
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector psi and puts the result in hpsi.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated by ccalbec)
  ! input:
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     hpsi  V_US|psi> is added to hpsi
  !
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : current_spin
  USE uspp,       ONLY : nkb, deeq
  USE uspp_param, ONLY : nh
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  COMPLEX(DP), INTENT(OUT) :: dnl(nkb,nkb)  
  !
  ! ... here the local variables
  !
  INTEGER :: jkb, ikb, ih, jh, na, nt, ijkb0
  !
  write(stdout,*) ' get_dnl : current_spin = ', current_spin
  write(stdout,'(i,e)') (ih, deeq(ih,ih,1,current_spin), ih=1,nh(1))
  IF ( nkb == 0 ) RETURN
  !
  ijkb0 = 0
  iblock=0
  !
  DO nt = 1, ntyp
     !
     DO na = 1, nat
        !
        IF ( ityp(na) == nt ) THEN
           !
           iblock=iblock+1
           dnl(iblock)%block_dspl = ijkb0
           dnl(iblock)%block_size = nh(nt)
           if( allocated( dnl(iblock)%block ) ) deallocate( dnl(iblock)%block )
           allocate( dnl(iblock)%block(nh(nt),nh(nt)) )
           dnl(iblock)%block = deeq(1:nh(nt),1:nh(nt),na,current_spin)
           !
           ijkb0 = ijkb0 + nh(nt)
           !
        END IF
        !
     END DO
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE get_dnl
