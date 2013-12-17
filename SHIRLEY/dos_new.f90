  program dos

  implicit none

  integer,parameter :: dp=kind(1.d0)
  real(dp),parameter :: evtory=7.3498649d-2

  integer :: narg

  integer :: nener
  real(dp) :: e1, e2, sigma
  real(dp),allocatable :: ener(:), spec(:)

  integer :: i, n, m, nk
  real(dp) :: de

  real(dp) :: kvec(3), wk

  integer :: nblk
  type block_type
    complex(dp),pointer :: block(:,:)
  end type
  type(block_type),allocatable :: datablk(:)

  integer :: iunf
  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma

  integer,external :: freeunit
  integer,external :: iargc

  narg = iargc()
  if( narg /= 5 ) then
    write(*,*) ' usage: dos_new e1 e2 nener sigma filename'
    stop
  endif

  call getarg( 1, ce1 )
  call getarg( 2, ce2 )
  call getarg( 3, cnener )
  call getarg( 4, csigma )
  call getarg( 5, filename )

  read(ce1,*) e1
  read(ce2,*) e2
  read(cnener,*) nener
  read(csigma,*) sigma

  e1 = e1*evtory
  e2 = e2*evtory
  sigma = sigma*evtory

  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec=0.d0

  iunf=freeunit()
  open(iunf,file=trim(filename),form='unformatted')

  nk=0
  do while (.true.)
    read(iunf,end=901,err=911) kvec, wk
    read(iunf,end=901,err=911) nblk
    if( .not. allocated( datablk ) ) &
      allocate( datablk(nblk) )
    do i=1,nblk
      read(iunf,end=901,err=911) n, m
      if( .not. associated(datablk(i)%block) ) &
        allocate( datablk(i)%block(n,m) )
      read(iunf,end=901,err=911) datablk(i)%block(1:n,1:m)
    enddo

    do i=1,n
      call add_gauss( nener, ener, real(datablk(1)%block(i,1)), &
                      sigma, wk, spec )
    enddo

    nk=nk+1
  enddo

  901 continue

  close(iunf)
  
  write(*,*) '# ', nk, ' k-points'
  do i=1,nener
    write(*,*) ener(i)/evtory,spec(i)
  enddo

  911 stop

  contains

    subroutine add_gauss( nener, ener, e, sigma, weight, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, weight, spec(nener)

    real(dp) :: its2, pre
    real(dp) :: arg(nener)

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2) * weight
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    spec = spec + pre * exp( arg )
     
    end subroutine add_gauss

  end program dos

