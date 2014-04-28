  program addrho

  implicit none

  integer,parameter :: maxchar=255
  integer,parameter :: dp=kind(1.d0)
  integer,parameter :: iun1=20
  integer,parameter :: iun2=21
  integer,parameter :: iun12=22

  character(maxchar) :: frho1, frho2, frho12
  character(maxchar) :: cr1, cr2, cr3
  integer :: narg
  integer,external :: iargc
  integer :: nr1, nr2, nr3
  real(dp),allocatable,dimension(:) :: rho1, rho2

  integer :: ierr, i

  narg=iargc()
  if( narg /= 6 ) then
    write(0,*) ' usage: addrho nr1 nr2 nr3 rho1 rho2 rho12'
    stop
  endif

  call getarg( 1, cr1 )
  call getarg( 2, cr2 )
  call getarg( 3, cr3 )
  call getarg( 4, frho1 )
  call getarg( 5, frho2 )
  call getarg( 6, frho12 )

  read(cr1,*) nr1
  read(cr2,*) nr2
  read(cr3,*) nr3

  allocate( rho1(nr1*nr2*nr3), rho2(nr1*nr2*nr3), stat=ierr )
  if( ierr/=0 ) then
    write(0,*) ' problem allocating memory'
    stop
  endif
  
  open(iun1,file=trim(frho1),form='unformatted',status='old',iostat=ierr)
  if( ierr/=0 ) then
    write(0,*) ' problem opening rho1: ', trim(frho1)
    stop
  endif
  open(iun2,file=trim(frho2),form='unformatted',status='old',iostat=ierr)
  if( ierr/=0 ) then
    write(0,*) ' problem opening rho2: ', trim(frho2)
    stop
  endif
  open(iun12,file=trim(frho12),form='unformatted',status='unknown',iostat=ierr)
  if( ierr/=0 ) then
    write(0,*) ' problem opening rho12: ', trim(frho12)
    stop
  endif
  
  read(iun1) rho1
  read(iun2) rho2
  forall( i=1:nr1*nr2*nr3 ) rho2(i)=rho1(i)+rho2(i)
  write(iun12) rho2

  close(iun12)
  close(iun1)
  close(iun2)

  stop
  end program addrho
