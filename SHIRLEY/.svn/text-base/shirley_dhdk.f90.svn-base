  program shirley_dhdk

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and generate dH/dk, i.e. derivatives
  ! of the Hamiltonian in k-space - useful matrix elements for
  ! computing optical conductivity

  ! David Prendergast, UCB, Feb 2008

#include "f_defs.h"
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, world_comm
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_stop,  &
                 mp_scatter
  use mpio
  use kpt_module
  use corerepair_module
  use shirley_input_module, nelec_in=>nelec
  use diag_module

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)
  integer(dp),parameter :: bytes_per_realdp=8
  integer(dp),parameter :: bytes_per_cplxdp=16

  character(255) :: fmtstr, filename

  character(255) :: eigval_file, info_file, dhdk_file

  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: fheigval, fheigvec, fhdhdk
  integer :: reqeigval, reqreigvec, reqdhdk
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  integer :: group_ranks(1), group_size, ionode_group, ionode_comm
  integer :: wcomm, wgroup
  

  complex(dp),allocatable :: zeigval(:), dhamk(:,:,:), dhamkx(:,:,:)
  real(dp) :: fermi_energy

  integer :: ik, i, j, ispin
  integer :: ibasis, jbasis, ideriv
  integer :: nband_subset, nbnd

  real(dp) :: nelec, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin
  integer :: nk

  integer :: iuninf

  integer,external :: freeunit

  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


  call shirley_input

  call dump_system( nelec, alat, volume, at, bg, tpiba, nspin, lda_plus_u  )
  write(stdout,*)
  write(stdout,*) ' System details: '
  write(stdout,*) '       nelec = ', nelec
  write(stdout,*) '        alat = ', alat
  write(stdout,*) ' cell volume = ', volume
  write(stdout,*)

  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis

  if( band_subset(1) /= 1 ) then
    write(stdout,*) ' to determine the fermi energy I need all occupied bands'
    band_subset(1)=1
  endif
  write(stdout,*) ' band_subset = ', band_subset
  nband_subset = band_subset(2)-band_subset(1)+1

  ! close iunout as formatted and reopen as unformatted
!  write(6,*) 'filename=', filename
!  inquire(unit=iunout,name=filename)
!  close(iunout,status='delete')

  call diagx_init( band_subset(1), band_subset(2) )

  write(stdout,*) ' shirley_dhdk'
  write(stdout,*)

  ! store total number of k-points
  nk=kpt%list%nk

  ! MPI-IO
  if( ionode ) then
    eigval_file=trim(outfile)//'.eigval'
    dhdk_file=trim(outfile)//'.dhdk'
    info_file=trim(outfile)//'.info'
  endif
  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )
  call mp_file_open_dp( dhdk_file, fhdhdk, ionode_id, world_comm )


  allocate( dhamk(nbasis,nbasis,3), dhamkx(nband_subset,nband_subset,3) )

! ======================================================================
  do ispin=1,nspin
  do ik=1,kpt%list%nk
! ======================================================================

    if( mod(ik,nproc)/=mpime ) cycle

    write(stdout,'(a,3f,i,a,i,a,a)') ' k-point ', &
      kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                  ' on node ', trim(nodenumber)

    ! build the Hamiltonian for this q-point
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )
    
    ! diagonalize matrix to find eigenvalues and vectors
    call diagx_ham
    
    ! write eigenvalues
    offset = ((ispin-1)*kpt%list%nk + ik-1)*nband_subset
    call mpi_file_write_at( fheigval, offset, &
                             eigval(band_subset(1)), nband_subset, &
                             MPI_DOUBLE_PRECISION, status, ierr )

    do ideriv=1,3
      ! evaluate derivative along ideriv-th direction (crystal coordinates)
      call diag_build_dhamk( ideriv, kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, dhamk(:,:,ideriv) )
    enddo

    ! transform to Cartesian coordinates
    forall(i=band_subset(1):band_subset(2),j=band_subset(1):band_subset(2)) &
      dhamkx(i-band_subset(1)+1,j-band_subset(1)+1,:) &
        = matmul( dhamk(i,j,:), trkin2nlp )

    ! write dhdk
    offset = ((ispin-1)*kpt%list%nk + ik-1)*nband_subset*nband_subset*3*2
    call mpi_file_write_at( fhdhdk, offset, &
                             dhamkx, nband_subset*nband_subset*3*2, &
                             MPI_DOUBLE_PRECISION, status, ierr )

    if( debug ) then
      if( kpt%param%cartesian ) then
        write(mpime+100,'(i,3f)') ik, matmul( transpose(at), kpt%list%kvec(:,ik) )/tpiba
        write(mpime+100,'(i,3f)') ik, kpt%list%kvec(1:3,ik)/tpiba
      else
        write(mpime+100,'(i,3f)') ik, kpt%list%kvec(1:3,ik)
        write(mpime+100,'(i,3f)') ik, matmul( bg, kpt%list%kvec(:,ik) )
      endif
      write(mpime+100,*) eigval(band_subset(1):band_subset(2))*rytoev
    endif

! ======================================================================
  enddo ! loop over k-points ik
  enddo ! loop over spin ispin
! ======================================================================

!  if( ionode) then
     ! close binary dump
     call mpi_file_close( fheigval, ierr )
     
     ! close dhdk file
     call mpi_file_close( fhdhdk, ierr )
!  endif
  
  if( ionode ) then
    nbnd=nband_subset

    fermi_energy = efermi

    write(stdout,*) ' Fermi energy = ', fermi_energy

    iuninf=freeunit()
    open(iuninf,file=info_file,form='formatted')
    write(iuninf,nml=info)
    write(iuninf,*) kpt%list%wk
    write(iuninf,*) kpt%param%cartesian
    write(iuninf,*) kpt%list%kvec
    if( trim(kpt%param%grid_type) == 'tetrahedra' ) then
      write(iuninf,*) kpt%tetra%ntetra
      write(iuninf,*) kpt%tetra%tetra
    endif
    close(iuninf)
  endif
  
  ! deallocate
  deallocate( dhamk, dhamkx )

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  write(stdout,*) ' end shirley_dhdk'
  
  end program shirley_dhdk
