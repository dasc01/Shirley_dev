  program dos

#include "f_defs.h"
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum

  implicit none

  integer,parameter :: dp=kind(1.d0)
  real(dp),parameter :: evtory=7.3498649d-2

  character(len=3) :: nodenumber

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

  integer :: iunout
  character(255) :: fout

  integer,external :: freeunit
  integer,external :: iargc

  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

  narg = iargc()
  if( narg /= 5 ) then
    write(stdout,*) ' usage: dos_para e1 e2 nener sigma filename'
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

  endif

  call mp_bcast( e1, ionode_id )
  call mp_bcast( e2, ionode_id )
  call mp_bcast( nener, ionode_id )
  call mp_bcast( sigma, ionode_id )
  call mp_bcast( filename, ionode_id )

  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec=0.d0

  if( ionode ) then
    iunout=freeunit()
    fout=trim(filename)//'.dos'
    open(iunout,file=trim(fout),form='formatted')
    write(stdout,*) '    output in '//trim(fout)
  endif

  iunf=freeunit()
  if( trim(nodenumber) /= '' ) then
    filename = trim(filename) // '.' // trim(nodenumber)
  endif
  open(iunf,file=trim(filename),form='unformatted')

  write(stdout,*) '    running...'

  nk=0
  do while (.true.)
    read(iunf,end=901,err=911) kvec, wk
    read(iunf,end=901,err=911) nblk

    if( .not. allocated( datablk ) ) then
      allocate( datablk(nblk) )
      do i=1,nblk
        nullify( datablk(i)%block )
      enddo
    endif

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
  
  if( ionode ) write(stdout,*) '    waiting for other processors...'
  call mp_barrier

  call mp_sum( nk )
  if( ionode ) write(stdout,*) '    in total: ', nk, ' k-points'

  call mp_sum( spec )
  if( ionode ) then
    do i=1,nener
      write(iunout,*) ener(i)/evtory,spec(i)
    enddo
    close(iunout)
  endif

  911 continue
  call mp_end
  stop

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

