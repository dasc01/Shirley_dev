  program cellrelax

  implicit none
  integer,parameter :: dp=kind(1.d0)

  integer :: n, lwa, info
  real(dp),allocatable :: x(:), fvec(:), wa(:)
  real(dp) :: tol
  
  external hybrd1

  ! provide the initial estimate of x
  write(0,*) ' number of parameters?'
  read(*,*) n
  if( n < 1 ) then
    write(*,*) ' error: n must be positive'
    write(*,*) ' n = ', n
    stop
  endif

  lwa = (n*(3*n+13))/2
  allocate( x(n), fvec(n), wa(lwa) )

  write(0,*) ' requested tolerance in solution x?'
  read(*,*) tol

  write(0,*) ' initial estimate of x?'
  read(*,*) x

  write(0,*) ' entering hybrd1 ... '

  call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)

  write(0,*) ' hybrd1 finished with info = ', info
  write(0,*) ' fvec = '
  write(0,'(e)') fvec
  write(0,*) ' x = '
  write(0,'(e)') x

  stop

  contains

    subroutine fcn(n,x,fvec,iflag)
    integer :: n,iflag
    real(dp) :: x(n),fvec(n)

    character :: cont
!   ----------
!   calculate the functions at x and
!   return this vector in fvec.
!   ---------
    write(0,*) ' want to calculate the function at x = '
    write(*,'(e)') x

    write(0,*) ' please provide the ', n, ' function values fvec = '
    read(*,*) fvec

!    write(*,*) ' Continue (y/n)?'
!    read(*,*) cont
!    if( cont /= 'y' ) then
!      iflag = -1
!    endif

    return
    end subroutine fcn

  end program cellrelax
