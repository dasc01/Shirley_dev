C
C     ***********************
C     SIMPLE DRIVER FOR LBFGS
C     ***********************
C
C     Example of driver for LBFGS routine, using a
C     simple test problem. The solution point is at 
C     X=(1,...,1) and the optimal function value of 0.
C
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
      PROGRAM SDRIVE

      use lbfgs_module

      DOUBLE PRECISION F,EPS,XTOL,T1,T2
      DOUBLE PRECISION, ALLOCATABLE :: X(:), G(:), DIAG(:), W(:)
      INTEGER IPRINT(2),IFLAG,ICALL,N,M,J
      LOGICAL DIAGCO
C
      write(0,*) ' number of parameters?'
      read(*,*) N
      if( N < 1 ) then
        write(*,*) ' error: N must be positive'
        write(*,*) ' N = ', N
        stop
      endif

c     number of steps
      M=3

C     allocate space
      allocate( X(N), G(N), DIAG(N), W( N*(2*M+1)+2*M ) )

      IPRINT(1)= 0
      IPRINT(2)= 0
      LP=0
      MP=0
C
C     We do not wish to provide the diagonal matrices Hk0, and 
C     therefore set DIAGCO to FALSE.
C
      DIAGCO= .FALSE.

      write(0,*) ' requested tolerance in renormalized gradient?'
      read(*,*) EPS

C     estimate of machine precision
      XTOL=1.d-16

      ICALL=0
      IFLAG=0

      write(0,*) ' initial estimate of x?'
      read(*,*) X
C
      write(0,*) ' entering lbfgs ... '

      STPMAX=1.d-1

 20   CONTINUE

      write(0,*) ' want to calculate the function at x = '
      write(*,'(e)') X

      write(0,*) ' please provide the function value f = '
      read(*,*) F
      
      write(0,*) ' please provide the ', N, ' gradient values g = '
      read(*,*) G

      CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
C     We allow at most 2000 evaluations of F and G
      IF(ICALL.GT.2000) GO TO 50
      STPMAX=1.d+20
      GO TO 20
  50  CONTINUE
      write(*,*) ' iflag = ', iflag
      END
C
C     ** LAST LINE OF SIMPLE DRIVER (SDRIVE) **

