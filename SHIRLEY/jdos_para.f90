  program jdos

#include "f_defs.h"
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum

  implicit none

  integer,parameter :: dp=kind(1.d0)
  real(dp),parameter :: evtory=7.3498649d-2
  real(dp),parameter :: ktoev=8.6173425d-5

  character(len=3) :: nodenumber

  integer :: narg

  integer :: nener
  real(dp) :: e1, e2, sigma, efermi, temp
  real(dp),allocatable :: ener(:), spec(:,:)

  integer :: n, m, nk
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
  character(255) :: cefermi
  character(255) :: ctemp

  integer :: iunout
  character(255) :: fout

  integer :: ierr
  real(dp),allocatable :: eig(:), occ(:), deig(:,:), weig(:,:)
  real(dp),allocatable :: fosc(:,:,:)
  real(dp) :: w
  real(dp),parameter :: eps=1.d-12
  integer :: i, j, iblk

  integer,external :: freeunit
  integer,external :: iargc

  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

  narg = iargc()
  if( narg /= 7 ) then
    write(stdout,*) ' usage: jdos_para e1 e2 nener sigma efermi temp filename'
    stop
  endif

  call getarg( 1, ce1 )
  call getarg( 2, ce2 )
  call getarg( 3, cnener )
  call getarg( 4, csigma )
  call getarg( 5, cefermi )
  call getarg( 6, ctemp )
  call getarg( 7, filename )

  read(ce1,*) e1
  read(ce2,*) e2
  read(cnener,*) nener
  read(csigma,*) sigma
  read(cefermi,*) efermi
  read(ctemp,*) temp

  e1 = e1*evtory
  e2 = e2*evtory
  sigma = sigma*evtory
  efermi = efermi*evtory
  temp = temp*ktoev*evtory

  endif

  call mp_bcast( e1, ionode_id )
  call mp_bcast( e2, ionode_id )
  call mp_bcast( nener, ionode_id )
  call mp_bcast( sigma, ionode_id )
  call mp_bcast( efermi, ionode_id )
  call mp_bcast( temp, ionode_id )
  call mp_bcast( filename, ionode_id )

  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener,7) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec=0.d0

  if( ionode ) then
    iunout=freeunit()
    fout=trim(filename)//'.jdos'
    open(iunout,file=trim(fout),form='formatted')
    write(stdout,*) '    output in '//trim(fout)
  endif

  iunf=freeunit()
  if( trim(nodenumber) /= '' ) then
    filename = trim(filename) // '.' // trim(nodenumber)
  endif
  open(iunf,file=trim(filename),form='unformatted',status='old',iostat=ierr)
  if( ierr /= 0 ) then
    call errore('jdos_para','problem opening file '//trim(filename),1)
  endif

  write(stdout,*) ' using file ', trim(filename)
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


    if( .not. allocated(eig) ) allocate(eig(n))
    if( .not. allocated(occ) ) allocate(occ(n))
    eig = real(datablk(1)%block(:,1))
    do i=1,n
      occ(i) = fermi_occ(eig(i), efermi, temp )
    enddo

    if( .not. allocated(deig) ) allocate(deig(n,n))
    if( .not. allocated(weig) ) allocate(weig(n,n))
    do j=1,n
      do i=1,n
        deig(i,j) = eig(j)-eig(i)
        weig(i,j) = occ(i)-occ(j) ! reversed order
      enddo
    enddo

    if( .not. allocated(fosc) ) allocate(fosc(n,n,7))
    fosc(:,:,1) = 1.d0
    do j=1,n
      do i=1,n
        do iblk=2,7
          fosc(i,j,iblk) = real( conjg(datablk(iblk)%block(i,j)) &
                                     * datablk(iblk)%block(i,j) )
        enddo
      enddo
    enddo

    do j=1,n
    do i=1,n
      if( abs(weig(i,j)) > eps ) then
        do iblk=1,7
          w = fosc(i,j,iblk)*weig(i,j)*wk
          call add_gauss( nener, ener, deig(i,j), &
                          sigma, w, spec(:,iblk) )
        enddo
      endif
    enddo
    enddo

    nk=nk+1
    write(stdout,*) ' nk = ', nk
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
      write(iunout,'(e,7e)') ener(i)/evtory,(spec(i,iblk), iblk=1,7)
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


    function fermi_occ( ener, mu, temp )

    real(dp) :: fermi_occ
    real(dp),intent(in) :: ener, mu, temp
    real(dp) :: f

    f = exp((ener-mu)/temp)+1.d0
    fermi_occ = 1.d0 / f
    
    end function fermi_occ

  end program jdos

