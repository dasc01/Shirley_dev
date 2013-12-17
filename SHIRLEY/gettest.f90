  program gettest

  use kinds, only : dp
  use mp, only : mp_get, mp_scatter, mp_barrier, mp_scatter_size
  use mp_global, only : mpime, nproc
  use io_global, only : ionode, ionode_id

  implicit none

  character(3) :: nodenumber

  real(dp),allocatable :: global(:)
  real(dp),allocatable :: local(:)

  integer :: n, i
  integer,allocatable :: n_local(:), i_local(:)
  integer :: n_l

  call start_shirley( nodenumber )

  n=19

  ! construct global data set
  if( ionode ) then
    allocate( global(n) )
    do i=1,n
      global(i) = i
    enddo
  endif

  call mp_scatter_size( n, n_l, ionode_id )

  allocate( local(n_l) )

  call mp_scatter( global, local, ionode_id )

!  if( ionode ) then
!    write(*,*) ' global = ', global
!  endif
  call mp_barrier
  do i=0,nproc-1
    if( i==mpime ) then
      write(*,*) ' proc = ', mpime, ' local = ', local
      call flush(6)
    endif
    call mp_barrier
  enddo

  end program gettest
