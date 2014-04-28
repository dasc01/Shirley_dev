  program plasma_freq

#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, world_comm
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use fermi, only : fermifunc, fermideriv

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2
  real(dp),parameter :: rytoev=1.d0/evtory
  real(dp),parameter :: kelvin2rydberg=6.333630d-6

  character(len=3) :: nodenumber

  integer :: narg

  real(dp) :: kT

  integer :: ik, i, j, k, ij, nij, ideriv, niter

  real(dp) :: wq, dw, de, wt, wkfac
  real(dp),allocatable :: eigval(:,:), wg(:,:), wk(:), wk_loc(:), kvec(:,:)
  logical :: cartesian
  real(dp) :: wpl2(1,4)
  complex(dp) :: condpl
  real(dp) :: efermi
  integer :: nktot, nbndf, wktot, wkzero
  integer :: nfs
  complex(dp),allocatable :: dhdk(:,:,:,:)
  logical :: nonuniform
  real(dp) :: prefac

  integer :: iunf
  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma
  character(255) :: ckT, cefermi

  integer :: iunout
  character(255) :: fout
  integer :: iuninf, ierr
  integer :: fheigval, fhdhdk
  character(255) :: eigval_file, dhdk_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval, reqdhdk
  integer :: nk_loc

  integer,external :: freeunit
  integer,external :: iargc

  integer :: nk, nbnd
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy

  integer :: nspin
  integer,allocatable :: isk(:)
  integer :: ntetra
  integer,allocatable :: tetra(:,:)
  real(dp),external :: efermit

  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

  narg = iargc()
  if( narg /= 3 ) then
    write(stdout,*) ' usage: plasma_freq.x temp efermi filename'
    stop
  endif

  call getarg( 1, ckT )
  call getarg( 2, cefermi )
  call getarg( 3, filename )

  read(ckT,*) kT
  read(cefermi,*) efermi

  kT = kT*kelvin2rydberg
  efermi = efermi*evtory

  endif

  call mp_bcast( kT, ionode_id )
  call mp_bcast( efermi, ionode_id )
  call mp_bcast( filename, ionode_id )

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
!    read(iuninf,*) ntetra
!    allocate( tetra(4,ntetra) )
!    read(iuninf,*) tetra
    close(iuninf)
    write(stdout,*) ' reading info from ', trim(filename)
    write(stdout,*) ' fermi_energy = ', fermi_energy
    write(stdout,*) '      sum(wk) = ', sum(wk)
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

!  call mp_bcast( ntetra, ionode_id )
!  if( .not. ionode ) allocate( tetra(4,ntetra) )
!  call mp_bcast( tetra, ionode_id )

  wktot = sum(wk)
  wkfac=2.d0/wktot

  ! how many are local?
  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
  enddo

  allocate( wk_loc(nk_loc) )
  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
    wk_loc(nk_loc) = wk(nk)
  enddo
  deallocate( wk )

!  call mp_bcast( ntetra, ionode_id )
!  if( .not. ionode ) allocate( tetra(4,ntetra) )
!  call mp_bcast( tetra, ionode_id )

! get this value from the command line instead
!  efermi = fermi_energy / rytoev

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


  allocate( eigval(nbnd,nk_loc), &
            wg(nbnd,nk_loc), &
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
    ! read weights
    offset = (nk+ik-1)*nbnd
    call mpi_file_iread_at( fheigval, offset, &
                           wg(1,nk_loc), nbnd, &
                           MPI_DOUBLE_PRECISION, reqeigval, ierr )
    ! read dhdk
    offset = (ik-1)*nbnd*nbnd*3*2
    call mpi_file_iread_at( fhdhdk, offset, &
                           dhdk(1,1,1,nk_loc), nbnd*nbnd*3*2, &
                           MPI_DOUBLE_PRECISION, reqdhdk, ierr )
  enddo

  ! close - otherwise I don't know if I have the data or not
  call mp_wait( reqeigval, status, ierr )
  call mpi_file_close( fheigval, ierr )
  call mp_wait( reqdhdk, status, ierr )
  call mpi_file_close( fhdhdk, ierr )

  write(stdout,*) '    input Fermi energy ...'
  write(stdout,*) ' Fermi Energy = ', efermi*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT/kelvin2rydberg, ' K'
  write(stdout,*) ' nelec = ', nelec

!  ! redetermine fermi energy
!  nspin=1
!  allocate( isk(nk) )
!  isk=1
!  efermi = efermit (eigvalk, nbnd, nk, nelec, nspin, ntetra, tetra, 0, isk)

  write(stdout,*) ' Fermi Energy = ', efermi*rytoev, ' eV'
  write(stdout,*) '       volume = ', volume

!  do ik=1,nk_loc
!  do i=1,nbnd
!    write(stdout,*) ik, i, wg(i,ik)
!  enddo
!  enddo
!  write(stdout,*)

! weights read in from file should be good - don't overwrite them
!  forall( i=1:nbnd, ik=1:nk_loc ) &
!    wg(i,ik) = fermifunc( eigval(i,ik), efermi, kT )*2.d0/dble(nk)

!  do ik=1,nk_loc
!  do i=1,nbnd
!    write(stdout,*) ik, i, wg(i,ik)
!  enddo
!  enddo
!  write(stdout,*)


  wpl2=0.d0
  do ik=1,nk_loc
    do ideriv=1,3
      ! plasma frequency squared
      do i=1,nbnd
        condpl = dhdk(i,i,ideriv,ik)
        condpl = condpl * conjg(condpl)
        condpl = condpl * prefac

        wt = -1.d0 * fermideriv( eigval(i,ik), efermi, kT )

        !write(stdout,'(2i,2f)') i, ideriv, wt, real(condpl)
        wpl2(1,ideriv+1)=wpl2(1,ideriv+1)+ wt*condpl*wk_loc(ik)*wkfac
      enddo
    enddo
  enddo

  if( ionode ) write(stdout,*) '    waiting for other processors...'
  call mp_barrier

  call mp_sum( wpl2 )
  wpl2 = wpl2
  if( ionode ) then
    wpl2(1,1) = sum(wpl2(1,2:4))/3.d0
    write(stdout,*) ' plasma frequencies'
    write(stdout,*) sqrt(wpl2(1,2:4))/evtory
    write(stdout,*) ' ave: ', sqrt(wpl2(1,1))/evtory
  endif

  call mp_barrier
  stop

  end program plasma_freq

