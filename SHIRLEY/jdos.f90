  program main

  ! computes for a given DOS(E) the integral:
  !  JDOS(E) = int_{Evbm-E}^{Evbm} DOS(Ep) * DOS(Ep+E) dEp

  implicit none

  integer,parameter :: dp=kind(1.d0)

  character(255) :: cne, cevbm
  integer :: ne
  real(dp) :: evbm, de
  real(dp),allocatable :: ener(:), dos(:), jdos(:)
  integer :: i, nvbm, nej, j

  call getarg( 1, cne )
  call getarg( 2, cevbm )
  read(cne,*) ne 
  read(cevbm,*) evbm 
  allocate( ener(ne), dos(ne) )
  do i=1,ne
    read(*,*) ener(i), dos(i)
  enddo

  ! assuming a uniform grid
  de = ener(2)-ener(1)

  ! find valence band max position
  do i=1,ne
    if( ener(i) > evbm ) exit 
  enddo
  nvbm = i
  write(*,*) '# Evbm = ', ener(nvbm), evbm
  if( nvbm == ne ) stop

  nej=min(ne-nvbm,nvbm-1)
  allocate( jdos(nej) )
  do i=1,nej
    jdos(i)=0.d0
    do j=nvbm-i,nvbm
      jdos(i) = jdos(i) + dos(j)*dos(j+i)
    enddo
    write(*,*) ener(i+1)-ener(1), jdos(i)*de
  enddo
  
  deallocate( jdos, dos, ener )

  end program main
