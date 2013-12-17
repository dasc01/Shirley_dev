  program dos

  implicit none

  integer,parameter :: dp=kind(1.d0)

  integer :: nk, nbnd, nener, icol
  real(dp) :: e1, e2, sigma
  real(dp),allocatable :: eigval(:,:)
  real(dp),allocatable :: ener(:), spec(:)

  integer :: i, ik, ibnd
  real(dp) :: de

  integer :: idum
  real(dp),allocatable :: rdum(:)
  integer :: nfiles
  integer,allocatable :: iunf(:)
  character(255),allocatable :: filename(:)

  integer,external :: freeunit

  read(*,*) nk, nbnd, icol
  read(*,*) e1, e2, nener, sigma
  read(*,*) nfiles
  allocate( filename(nfiles), iunf(nfiles) )
  do i=1,nfiles
    read(*,*) filename(i)
  enddo
  !write(*,*) '#', nk, nbnd, icol, e1, e2, nener, sigma, nfiles, filename
  allocate( eigval(nbnd,nk) )
  allocate( rdum(icol-1) )
  do i=1,nfiles
    iunf(i)=freeunit()
    open(iunf(i),file=trim(filename(i)),form='formatted')
  enddo
  i=1
  do ik=1,nk
    read(iunf(i),*,end=911,err=911) idum, rdum, eigval(1:nbnd,ik)
    cycle
    911 i=i+1
    if( ik==nk ) exit
    if( i > nfiles ) stop
    read(iunf(i),*,end=911) idum, rdum, eigval(1:nbnd,ik)
  enddo
  do i=1,nfiles
    close(iunf(i))
  enddo
  !write(*,*) 'done reading'
  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec=0.d0
  do ik=1,nk
    do ibnd=1,nbnd
      call add_gauss( nener, ener, eigval(ibnd,ik), sigma, spec )
    enddo
  enddo
  spec = spec *2.d0 / dble(nk)
  
  do i=1,nener
    write(*,*) ener(i),spec(i)
  enddo

  stop

  contains

    subroutine add_gauss( nener, ener, e, sigma, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, spec(nener)

    real(dp) :: its2, pre
    real(dp) :: arg(nener)

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2)
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    spec = spec + pre * exp( arg )
     
    end subroutine add_gauss

  end program dos

