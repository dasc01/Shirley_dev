  program cond_para

#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, world_comm
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum, mp_min
  use fermi, only : fermifunc, fermideriv
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2
  real(dp),parameter :: rytoev=1.d0/evtory
  real(dp),parameter :: kelvin2rydberg=6.333630d-6

  character(len=3) :: nodenumber

  integer :: narg

  integer :: nener, ib1, ib2
  real(dp) :: e1, e2, sigma
  real(dp) :: kT, ef, tol
  real(dp),allocatable :: ener(:), spec(:,:)

  integer :: ik, i, j, k, ij, nij, ideriv, niter

  real(dp) :: de
  real(dp),allocatable :: eigval(:,:), wg(:,:), wk(:), kvec(:,:)
  logical :: cartesian
  complex(dp),allocatable :: cond(:)
  real(dp),allocatable :: df(:), deig(:), dsig(:)
  real(dp) :: efermi
  integer :: nktot
  real(dp) :: wktot, wkfac
  integer,allocatable :: index_wg(:,:)
  integer,allocatable :: index_fs(:)
  integer :: nfs
  complex(dp),allocatable :: dhdk(:,:,:,:)
  logical :: nonuniform
  real(dp) :: prefac, deltak

  integer :: iunf
  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma
  character(255) :: ckT
  character(255) :: cef
  character(255) :: ctol
  character(255) :: cib1,cib2

  integer :: iunout
  character(255) :: fout
  integer :: iuninf, ierr
  integer :: fheigval, fhdhdk
  character(255) :: eigval_file, dhdk_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval, reqdhdk
  integer :: nk_loc
  real(dp),allocatable :: wk_loc(:)
  real(dp) :: demin, dh

  integer,external :: freeunit
  integer,external :: iargc

  integer :: nk, nbnd
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


  integer :: nspin
  integer,allocatable :: isk(:)
  integer :: ntetra
  integer,allocatable :: tetra(:,:)
  real(dp),external :: efermit
  logical :: lda_plus_u
  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

  narg = iargc()
  if( narg /= 7 .and. narg /= 8 .and. narg /= 9 .and. narg /= 10  ) then
    write(stdout,*) ' usage: cond_para e1 e2 nener sigma temp ef [tol] filename [ib1 ib2]'
    stop
  endif
  write(*,*) "narg=",narg
  call getarg( 1, ce1 )
  call getarg( 2, ce2 )
  call getarg( 3, cnener )
  call getarg( 4, csigma )
  call getarg( 5, ckT )
  call getarg( 6, cef )
  if( narg==7 ) then
     call getarg( 7, filename )
  else if( narg==8 ) then
     call getarg( 7, ctol )
     call getarg( 8, filename )
  else if (narg==9) then
     call getarg( 7, filename )
     call getarg( 8, cib1)
     call getarg( 9, cib2)
  else if (narg==10) then
     call getarg( 7, ctol )
     call getarg( 8, filename )
     call getarg( 9, cib1)
     call getarg(10, cib2)
  endif


  read(ce1,*) e1
  read(ce2,*) e2
  read(cnener,*) nener
  read(csigma,*) sigma
  read(ckT,*) kT
  read(cef,*) ef
  tol=0.d0
  if( narg==8 ) then
    read(ctol,*) tol
  endif
  if(narg==9) then
     read(cib1,*) ib1
     read(cib2,*) ib2     
  else if(narg==10) then
     read(ctol,*) tol
     read(cib1,*) ib1
     read(cib2,*) ib2
  endif
  e1 = e1*evtory
  e2 = e2*evtory
  sigma = sigma*evtory
  kT = kT*kelvin2rydberg
  ef = ef*evtory
  tol = tol*evtory

  endif

  call mp_bcast( e1, ionode_id )
  call mp_bcast( e2, ionode_id )
  call mp_bcast( nener, ionode_id )
  call mp_bcast( sigma, ionode_id )
  call mp_bcast( kT, ionode_id )
  call mp_bcast( ef, ionode_id )
  call mp_bcast( tol, ionode_id )
  call mp_bcast( filename, ionode_id )

  nonuniform = ( sigma < 0.d0 )
  sigma = abs(sigma)

  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener,8) )
  do i=1,nener
    ener(i) = e1 + dble(i)*de
  enddo
  spec=0.d0


  if( ionode ) then
    iunout=freeunit()
    fout=trim(filename)//'.cond'
    open(iunout,file=trim(fout),form='formatted')
    write(stdout,*) '    output in '//trim(fout)

    iuninf=freeunit()
    open(iuninf,file=trim(filename)//'.info',form='formatted')
    read(iuninf,nml=info)
    allocate( wk(nk), kvec(1:3,nk) )
    read(iuninf,*) wk
    read(iuninf,*) cartesian
    read(iuninf,*) kvec
    close(iuninf)
    write(stdout,*) ' reading info from ', trim(filename)
    write(stdout,*) ' fermi_energy (stored) = ', fermi_energy
  endif

  call mp_bcast( nk, ionode_id )
  call mp_bcast( nbnd, ionode_id )
  call mp_bcast( nelec, ionode_id )
  call mp_bcast( alat, ionode_id )
  call mp_bcast( volume, ionode_id )
  call mp_bcast( at, ionode_id )
  call mp_bcast( bg, ionode_id )
  call mp_bcast( tpiba, ionode_id )
  call mp_bcast( fermi_energy, ionode_id )

  if( .not. ionode ) allocate( wk(nk), kvec(3,nk) )
  call mp_bcast( wk, ionode_id )
  call mp_bcast( cartesian, ionode_id )
  call mp_bcast( kvec, ionode_id )

  if( ionode ) then
     if(narg<9) then
        ib1=1
        ib2=nbnd
     endif
  endif
  call mp_bcast( ib1, ionode_id )
  call mp_bcast( ib2, ionode_id )

  wktot = sum(wk)
  wkfac=2.d0/wktot

  ! copy the input ef
  efermi = ef 

! compute deltak - linear spacing between k-vectors in grid
  deltak = 2.d0*acos(-1.d0)/(volume*dble(nk))**(1.d0/3.d0)
  write(stdout,*) ' linear spacing between k-vectors in grid'
  write(stdout,*) ' deltak = ', deltak

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'
  dhdk_file=trim(filename)//'.dhdk'
  
  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )
  call mp_file_open_dp( dhdk_file, fhdhdk, ionode_id, world_comm )


  ! Read Ambrosch-Draxl and Sofo, Comput Phys Comm 175, 1 (2006)
  ! This is the prefactor from Eq. 20 but converted back to a discrete sum
  ! over k-points, i.e. multiplying by (2 pi)^3/volume
  prefac = 8.d0 * acos(-1.d0)**2 
  prefac = prefac / volume

  write(stdout,*) '    running ...'

  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
  enddo

!  write(mpime+100,*) ' nk_loc = ', nk_loc
!  write(mpime+100,*) nk_loc*(nbnd*(2+nbnd*3))

  allocate( eigval(nbnd,nk_loc), &
            wg(nbnd,nk_loc), wk_loc(nk_loc), &
            dhdk(nbnd,nbnd,3,nk_loc), stat=ierr )
  if( ierr/=0 ) call errore('allocation error',ierr)

  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle

    write(stdout,*) ' reading ik= ', ik, ' of ', nk

    nk_loc = nk_loc + 1

    ! read eigenvalues
    offset = (ik-1)*nbnd
    call mpi_file_iread_at( fheigval, offset, &
                           eigval(1,nk_loc), nbnd, &
                           MPI_DOUBLE_PRECISION, reqeigval, ierr )
    ! read dhdk
    offset = (ik-1)*nbnd*nbnd*3*2
    call mpi_file_iread_at( fhdhdk, offset, &
                           dhdk(1,1,1,nk_loc), nbnd*nbnd*3*2, &
                           MPI_DOUBLE_PRECISION, reqdhdk, ierr )
  enddo

  ! close - otherwise I don't know if I have the data or not
  do ik=1,min(nk,nproc)
    if( mod(ik,nproc)/=mpime ) cycle
    call mp_wait( reqeigval, status, ierr )
    call mp_wait( reqdhdk, status, ierr )
  enddo
  call mpi_file_close( fheigval, ierr )
  call mpi_file_close( fhdhdk, ierr )

  write(stdout,*) '    input Fermi energy ...'
  write(stdout,*) ' Fermi Energy = ', efermi*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT/kelvin2rydberg, ' K'
  write(stdout,*) ' nelec = ', nelec

  ! distribute weights
  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
    wk_loc(nk_loc)=wk(ik)
  enddo

  write(stdout,*) ' Fermi Energy = ', efermi*rytoev, ' eV'

  forall( i=1:nbnd, ik=1:nk_loc ) &
    wg(i,ik) = fermifunc( eigval(i,ik), efermi, kT )*wk_loc(ik)*wkfac

  ! now wg contain weights
  ! find pairs of energies which have useful weight differences
  allocate( index_wg(2,nbnd*nbnd) )

  do ik=1,nk_loc
    ! establish non-zero contributions
    ij=0
    do j=1,nbnd
    do i=ib1,ib2
      de=eigval(j,ik)-eigval(i,ik)
      ! use the input tolerance to determine when we include transitions
      ! ideally, one should treat transitions of energy less than tol using
      ! an approximation like the Drude term
      if( abs(de) > tol ) then
        if( de > ener(1)-4.d0*sigma .and. de < ener(nener)+4.d0*sigma ) then
          ij=ij+1
          index_wg(:,ij) = (/ i, j /)
        endif
      endif
    enddo
    enddo
    nij=ij
    if( nij > 0 ) then
      write(stdout,'(i,a,i)') ik, ' nij= ', nij
    else
      cycle
    endif

    allocate( cond(nij), deig(nij), df(nij), dsig(nij), stat=ierr )
    if( ierr/=0 ) call errore('failed to allocate',ierr)

    dsig = sigma

    do ideriv=1,3
      do ij=1,nij
        cond(ij)=dhdk(index_wg(1,ij),index_wg(2,ij),ideriv,ik)
      enddo
      do ij=1,nij
        i = index_wg(1,ij)
        j = index_wg(2,ij)
        deig(ij) = eigval(j,ik)-eigval(i,ik)
        if( deig(ij) > 1.d-6 ) then
          df(ij)=(wg(i,ik)-wg(j,ik))/deig(ij)
        else
          df(ij)=-wk_loc(ik) * wkfac * 0.5d0 * ( &
            fermideriv( eigval(i,ik), efermi, kT ) &
          + fermideriv( eigval(j,ik), efermi, kT ) )
        endif

        ! reweighted broadening
        demin=1.d20
        if( nonuniform ) then
          dh=deltak * sqrt( &
               dhdk(i,i,ideriv,ik)*conjg(dhdk(i,i,ideriv,ik)) + &
               dhdk(j,j,ideriv,ik)*conjg(dhdk(j,j,ideriv,ik)) )
          demin=min(demin,dh)
          dsig(ij) = max( sigma, dh )
        endif
      enddo

      cond = (conjg(cond) * cond) * df
      cond = cond * prefac

      do ij=1,nij
        ! Gaussian broadening
        call add_gauss( nener, ener, deig(ij), &
                        dsig(ij), real(cond(ij)), spec(1:nener,ideriv+1) )

        ! Lorentzian broadening
        call add_lorentz( nener, ener, deig(ij), &
                        dsig(ij), real(cond(ij)), spec(1:nener,ideriv+5) )
      enddo
    enddo
    deallocate( cond, df, deig, dsig )

  enddo
  call mp_min( demin )
  write(stdout,*) ' minimum energy resolution = ', demin

  if( nonuniform ) then
    if( ionode ) write(stdout,*) '    spectrum broadened non-uniformly using'
    if( ionode ) write(stdout,*) '      local derivate dH/dk'
  endif
  if( ionode ) write(stdout,*) '    waiting for other processors...'
  call mp_barrier

  call mp_sum( spec )
  if( ionode ) then
    ! construct spatial average
    forall( i=1:nener) spec(i,1) = sum(spec(i,2:4))/3.d0
    forall( i=1:nener) spec(i,5) = sum(spec(i,6:8))/3.d0
    do i=1,nener
      write(iunout,'(9e)') ener(i)/evtory,spec(i,1:8)
    enddo
    close(iunout)
  endif

  call mp_barrier
  stop

  contains

    subroutine add_gauss( nener, ener, e, sigma, weight, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, weight, spec(:)

    real(dp) :: its2, pre
    real(dp) :: arg(nener)
    integer :: i

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2) * weight
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    forall( i=1:nener ) &
      spec(i) = spec(i) + pre * exp( arg(i) )
     
    end subroutine add_gauss

    subroutine add_lorentz( nener, ener, e, sigma, weight, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, weight, spec(:)

    real(dp) :: its2, pre
    real(dp) :: arg(nener)
    integer :: i

    its2=sigma*sigma
    pre=1.d0/acos(-1.d0) * weight
    arg = ener - e
    arg = arg**2.d0 + its2
    arg = sigma / arg
    forall( i=1:nener ) &
      spec(i) = spec(i) + pre * arg(i)
     
    end subroutine add_lorentz

  end program cond_para

