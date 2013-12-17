!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rgrid( r )
  !----------------------------------------------------------------------------
  !
  USE kinds,        ONLY : DP
  use fft_base,     only : dfftp, grid_scatter
  USE mp_global,    ONLY : inter_pool_comm, me_pool, root_pool
  USE mp,           ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  REAL(DP)     :: r(dfftp%nnr,3)
#if defined (__PARA)
  REAL(DP), ALLOCATABLE :: allr(:,:)
  integer :: ixyz
#endif
  integer :: irx1, irx2, irx3, n
  !
  ! 
#if defined (__PARA)
  !
  ! ... parallel case
  !
  ALLOCATE( allr( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3 ) )
  !
  n=0
  do irx3=1,dfftp%nr3x
  do irx2=1,dfftp%nr2x
  do irx1=1,dfftp%nr1x
    n=n+1
    allr(n,1) = dble(irx1-1)/dble(dfftp%nr1x)
    allr(n,2) = dble(irx2-1)/dble(dfftp%nr2x)
    allr(n,3) = dble(irx3-1)/dble(dfftp%nr3x)
  enddo
  enddo
  enddo
  !
  IF( me_pool == root_pool ) &
    CALL mp_bcast( allr, root_pool, inter_pool_comm )
  !
  DO ixyz = 1, 3
    !        
    CALL grid_scatter( allr(:,ixyz), r(:,ixyz) )
    !
  END DO
  !
  DEALLOCATE( allr )
  !
#else
  !
  ! ... serial case
  !
  n=0
  do irx3=1,dfftp%nr3x
  do irx2=1,dfftp%nr2x
  do irx1=1,dfftp%nr1x
    n=n+1
    r(n,1) = dble(irx1-1)/dble(dfftp%nr1x)
    r(n,2) = dble(irx2-1)/dble(dfftp%nr2x)
    r(n,3) = dble(irx3-1)/dble(dfftp%nr3x)
  enddo
  enddo
  enddo
  if( n /= dfftp%nnr ) call errore('rgrid','the size of r is not the product of dimensions',1)
  !
#endif
  !
  RETURN
  !
END SUBROUTINE rgrid
!----------------------------------------------------------------------------
SUBROUTINE scatter( f_in, f_out, dfft )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... REAL*8  f_in  = gathered variable (nrx1*nrx2*nrx3)
  ! ... REAL*8  f_out = distributed variable (nxx)
  !
  USE fft_types, ONLY: fft_dlay_descriptor
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL (DP) :: f_in(:), f_out(:)
  type(fft_dlay_descriptor),intent(in) :: dfft
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( dfft%nproc - 1 )
     !
     sendcount(proc) = dfft%nnp * dfft%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( dfft%comm )
  !  
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_REAL8,   &
                     f_out, sendcount(dfft%mype), MPI_REAL8, &
                     dfft%root, dfft%comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(dfft%mype) /= dfft%nnr ) f_out(sendcount(dfft%mype)+1:dfft%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE scatter
!
