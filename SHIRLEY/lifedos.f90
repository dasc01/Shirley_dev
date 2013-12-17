  program lifedos

  implicit none

  integer,parameter :: dp=kind(1.d0)

  integer :: ndata, nk, npedv, nener, icol
  real(dp) :: e1, e2, sigma
  real(dp),allocatable :: eigval(:)
  real(dp),allocatable :: f(:,:)
  real(dp),allocatable :: ener(:), spec(:,:), dos(:)

  integer :: i, ik, ibnd
  real(dp) :: de

  integer :: ipedv
  integer :: nfiles
  integer,allocatable :: iunf(:)
  character(255),allocatable :: filename(:)
  character(255) :: fmtstr

  integer,external :: freeunit

  ! ndata is total number of data points (may include a factor of nk)
  read(*,*) ndata, nk, npedv
  read(*,*) e1, e2, nener, sigma
  read(*,*) nfiles
  allocate( filename(nfiles), iunf(nfiles) )
  do i=1,nfiles
    read(*,*) filename(i)
  enddo
  ! write(*,*) '#', ndata, nbnd, e1, e2, nener, sigma, nfiles, (trim(filename(i)), i=1,nfiles)
  allocate( eigval(ndata) )
  allocate( f(ndata,npedv) )
  do i=1,nfiles
    iunf(i)=freeunit()
    open(iunf(i),file=trim(filename(i)),form='formatted')
  enddo
  i=1
  do ik=1,ndata
    !!write(*,*) '# file '//trim(filename(i)), ik
    read(iunf(i),*,end=911,err=911) eigval(ik), f(ik,1:npedv)
    cycle
    911 i=i+1
    if( ik==ndata ) exit
    if( i > nfiles ) then
      write(*,*) ' only found ', ik, ' data points'
      stop
    endif
    read(iunf(i),*,end=911,err=911) eigval(ik), f(ik,1:npedv)
  enddo
  do i=1,nfiles
    close(iunf(i))
  enddo
  !write(*,*) 'done reading'
  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener,npedv), dos(nener) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec=0.d0
  dos=0.d0
  do ipedv=1,npedv
  do ik=1,ndata
    call add_gaussw( nener, ener, f(ik,ipedv), eigval(ik), sigma, spec(1,ipedv) )
    call add_gaussw( nener, ener, 1.d0, eigval(ik), sigma, dos )
  enddo
  enddo
  spec = spec *2.d0 / dble(nk)
  dos = dos *2.d0 / dble(nk)
  
  write(fmtstr,*) '(4e,', npedv, 'e)'
  do i=1,nener
    write(*,trim(fmtstr)) ener(i), dos(i), &        !DOS
                          sum(spec(i,1:npedv)), &   !Total scattering
                          sum(spec(i,4:npedv)), &   !Non-acoustic scattering
                          spec(i,1:npedv)           !Scattering per mode
  enddo

  stop

  contains

    subroutine add_gaussw( nener, ener, w, e, sigma, spec )

    integer :: nener
    real(dp) :: ener(nener), w, e, sigma, spec(nener)

    real(dp) :: its2, pre
    real(dp) :: arg(nener)

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2)
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    where( arg >= -30.d0 ) arg = exp(arg)
    where( arg < -30.d0 ) arg=0.d0
    spec = spec + w * pre * arg
     
    end subroutine add_gaussw

  end program lifedos

