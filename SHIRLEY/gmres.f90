! *************************************************************************
! **                 DRIVER EXAMPLE FOR THE GMRes CODE
! *************************************************************************
      subroutine gmres( n, a, b, x, restart ) 

      implicit none

      integer,parameter :: dp=kind(1.d0)

      integer,intent(in) :: n
      complex(dp),intent(in)  :: a(n,n)
      complex(dp),intent(in)  :: b(n)
      complex(dp),intent(out) :: x(n)
      integer :: restart

      complex(dp),parameter :: ZERO = (0.0d0, 0.0d0), &
                                ONE = (1.0d0, 0.0d0)

      complex(dp),allocatable :: work(:)

      integer,parameter :: ldstrt = 60
      integer :: lwork

      integer irc(5), icntl(8), info(3)

      integer matvec, precondLeft, precondRight, dotProd
      parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)

      real(dp) :: cntl(5), rinfo(2)

      integer revcom, colx, coly, colz, nbscal

!**************************************************************
!* Generate the test matrix a and set the right-hand side
!* in positions (n+1) to 2n of the array work.
!**************************************************************

      write(*,*) '***********************************************'
      write(*,*) 'This code is an example of use of GMRES'
      write(*,*) '***********************************************'
      write(*,*)

      if( restart > ldstrt ) write(*,*) ' Warning: restart = ', restart

!
!******************************************************
!* Initialize the control parameters to default value
!******************************************************
!
      call init_zgmres(icntl,cntl)
!
!************************
!c Tune some parameters
!************************
!
! Save the convergence history on standard output
      icntl(3) = 6
! Maximum number of iterations
      icntl(7) = 10000 
!
! preconditioner location
      icntl(4) = 0
! orthogonalization scheme
      icntl(5)=0
! initial guess
      icntl(6) = 0
! residual calculation strategy at restart
      icntl(8) = 1

!* Initialise the right hand side
      lwork = restart**2 + restart*(n+5) + 5*n + 1
      if( icntl(5)==2 .or. icntl(5)==3 ) then
        lwork = lwork + restart
      else
        lwork = lwork + 1
      endif
      if( icntl(8)/=1 ) lwork = lwork + n

      allocate( work(lwork) )

      work(n+1:2*n) = b
      if( any(x /= ZERO) ) then
        work(1:n) = x
      else
        work(1:n) = cmplx(1.0d-4)
      endif


!****************************************
!* Reverse communication implementation
!****************************************
!
10     call drive_zgmres(n,n,restart,lwork,work, &
               irc,icntl,cntl,info,rinfo)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)
!
       if (revcom.eq.matvec) then
! perform the matrix vector product
!        work(colz) <-- A * work(colx)
         call zgemv('N',n,n,ONE,a,n,work(colx),1, &
                  ZERO,work(colz),1)
         goto 10
!
       else if (revcom.eq.precondLeft) then
! perform the left preconditioning
!         work(colz) <-- M^{-1} * work(colx)
         call zcopy(n,work(colx),1,work(colz),1)
         call ztrsm('L','L','N','N',n,1,ONE,A,n,work(colz),n)
         goto 10
!
       else if (revcom.eq.precondRight) then
! perform the right preconditioning
         call zcopy(n,work(colx),1,work(colz),1)
         call ztrsm('L','U','N','N',n,1,ONE,A,n,work(colz),n)
         goto 10
!
       else if (revcom.eq.dotProd) then
!      perform the scalar product
!      work(colz) <-- work(colx) work(coly)
!
         call zgemv('C',n,nbscal,ONE,work(colx),n, &
                     work(coly),1,ZERO,work(colz),1)
         goto 10
       endif
!
!******************************
! dump the solution on a file
!******************************
!
      x = work(1:n)


      if (icntl(5).eq.0) then
        write(*,*) 'Orthogonalisation : MGS'
      elseif (icntl(5).eq.1) then
        write(*,*) 'Orthogonalisation : IMGS'
      elseif (icntl(5).eq.2) then
        write(*,*) 'Orthogonalisation : CGS'
      elseif (icntl(5).eq.3) then
        write(*,*) 'Orthogonalisation : ICGS'
      endif
      write(*,*) 'Restart : ', restart
      write(*,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
      write(*,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
      write(*,*) 'Optimal workspace = ', info(3)
      write(*,*) ' ************************************************ '
!
!
100    continue
!
      if( allocated(work) ) deallocate( work )

      return
      end subroutine gmres
