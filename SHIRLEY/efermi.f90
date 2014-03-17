  program efermi

#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, world_comm, mp_global_end
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum, mp_max, mp_min
  use fermi, only : fermifunc, fermideriv, fermi_energy_sub=>fermi_energy
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2
  real(dp),parameter :: rytoev=1.d0/evtory
  real(dp),parameter :: kelvin2rydberg=6.333630d-6

  character(len=3) :: nodenumber

  integer :: narg

  real(dp) :: kT

  integer :: ik, i, j, k, ij, nij, ideriv, niter, ibnd

  real(dp),allocatable :: eigval(:,:), wg(:,:),tmpeig(:)
  real(dp),allocatable :: wk(:), kvec(:,:)
  logical :: cartesian
  real(dp) :: ef, cbm, vbm
  integer :: nk_loc
  real(dp),allocatable :: wk_loc(:)

  character(255) :: filename
  character(255) :: ckT
  character(255) :: cnelec,cmag

  integer :: iuninf, ierr
  integer :: fheigval, fhdhdk
  character(255) :: eigval_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval

  integer,external :: freeunit
#ifdef __PGI
  integer,external :: iargc
#endif

  character(255) :: smearing, grid_type

  integer :: nk, nbnd, ncp, nspin,ispin,iband
  logical :: lda_plus_u
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
!DASb
  logical::gwrun
!DASe  
  namelist /info/ nk, nbnd, ncp, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u, gwrun
  real(dp) :: nelec_in, mag_in, nele_sp

! call sleep (30)
  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

     narg = iargc()
     if( narg /= 2 .and. narg /= 3 .and. narg /= 4) then
        write(stdout,*) ' usage: efermi temp filename [nelec] [mag]'
        stop
     endif

     call getarg( 1, ckT )
     call getarg( 2, filename )
     if( narg>=3 ) then
        call getarg( 3, cnelec )
        read(cnelec,*) nelec_in
     endif
     mag_in=0
     if( narg==4 ) then
        call getarg( 4, cmag )
        read(cmag,*) mag_in
     endif
     
     read(ckT,*) kT
     
     kT = kT*kelvin2rydberg
     
  endif

 
  call mp_bcast( kT, ionode_id )
  call mp_bcast( filename, ionode_id )

  if( ionode ) then
    iuninf=freeunit()
    open(iuninf,file=trim(filename)//'.info',form='formatted')
    read(iuninf,nml=info)
    allocate( wk(nk), kvec(1:3,nk) )
    read(iuninf,*) cartesian
    read(iuninf,*) grid_type
    read(iuninf,*) wk
    read(iuninf,*) kvec
    close(iuninf)
    write(stdout,*) ' reading info from ', trim(filename)
    write(stdout,*) ' fermi_energy (stored) = ', fermi_energy

    if( narg>=3 ) then
      write(stdout,*) ' number of electrons modified based on input:'
      write(stdout,*) ' nelec from Hamiltonian = ', nelec
      nelec=nelec_in
      write(stdout,*) ' nelec from input       = ', nelec
    endif
    if( narg==4 ) then
       write(stdout,*) ' mag from input       = ', mag_in
    endif
  endif

  call mp_bcast( nk, ionode_id )
  call mp_bcast( nbnd, ionode_id )
  call mp_bcast( nelec, ionode_id )
  call mp_bcast( mag_in, ionode_id )
  call mp_bcast( alat, ionode_id )
  call mp_bcast( volume, ionode_id )
  call mp_bcast( at, ionode_id )
  call mp_bcast( bg, ionode_id )
  call mp_bcast( tpiba, ionode_id )
  call mp_bcast( fermi_energy, ionode_id )
  call mp_bcast( nspin, ionode_id )
  call mp_bcast( lda_plus_u, ionode_id )

  if( .not. ionode ) allocate( wk(nk), kvec(3,nk) )
  call mp_bcast( cartesian, ionode_id )
  call mp_bcast( wk, ionode_id )
  call mp_bcast( kvec, ionode_id )

  ef = fermi_energy / rytoev

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'
  
  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )


  write(stdout,*) '    running ...'

  nk_loc=0
  do ik=1,nk
    if( mod(ik-1,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
  enddo

! debug
!  write(mpime+100,*) nk_loc
! debug

  allocate( eigval(nbnd,nk_loc), tmpeig(nbnd), &
            wg(nbnd,nk_loc), wk_loc(nk_loc), stat=ierr )
  if( ierr/=0 ) call errore('allocation error',ierr)

  ! distribute weights
  nk_loc=0
  do ik=1,nk
    if( mod(ik-1,nproc)/=mpime ) cycle

    nk_loc = nk_loc + 1
    wk_loc(nk_loc) = wk(ik)
  enddo

  do ispin=1,nspin

     nele_sp=(nelec+((-1)**(ispin-1))*mag_in)/nspin
     write(stdout,*)  ' nele_sp=', nele_sp
     
     nk_loc=0
     do ik=1,nk
        if( mod(ik-1,nproc)/=mpime ) cycle
        
        nk_loc = nk_loc + 1
        
        write(stdout,*) ' reading ik= ', ik, ' of ', nk, ' wt= ', wk_loc(nk_loc)
        
        ! read eigenvalues
        offset = ((ispin-1)*nk + ik-1)*nbnd
        call mpi_file_iread_at( fheigval, offset, &
             eigval(1,nk_loc), nbnd, &
             MPI_DOUBLE_PRECISION, reqeigval, ierr )
     enddo

  ! close - otherwise I don't know if I have the data or not
     do ik=1,min(nk,nproc)
        if( mod(ik-1,nproc)/=mpime ) cycle
        call mp_wait( reqeigval, status, ierr )
     enddo

     ! debug
     !  nk_loc=0
     !  do ik=1,nk
     !    if( mod(ik,nproc)/=mpime ) cycle
     !    nk_loc = nk_loc + 1
     !    if( cartesian ) then
     !      write(mpime+100,'(i,3f)') ik, matmul( transpose(at), kvec(:,ik) )/tpiba
     !      write(mpime+100,'(i,3f)') ik, kvec(1:3,ik)/tpiba
     !    else
     !      write(mpime+100,'(i,3f)') ik, kvec(1:3,ik)
     !      write(mpime+100,'(i,3f)') ik, matmul( bg, kvec(:,ik) )
     !    endif
     !    write(mpime+100,*) eigval(1:nbnd,nk_loc)*rytoev
     !  enddo
     ! debug
     
     write(stdout,*) '    input Fermi energy ...'
     write(stdout,*) ' Fermi Energy = ', ef*rytoev, ' eV'
     write(stdout,*) ' temperature = ', kT*rytoev, ' eV'
     write(stdout,*) ' temperature = ', kT/kelvin2rydberg, ' K'
     write(stdout,*) ' nelec = ', nele_sp
     
     ! determine fermi energy
     smearing='fermi-dirac'
     call fermi_energy_sub( nele_sp, nk_loc, nbnd, wk_loc, eigval, kT, smearing, ef , nspin)
     
     ! find VBM and CBM if present
     cbm= 9.d+20
     vbm=-9.d+20
     do ik=1,nk_loc
        do ibnd=1,nbnd
           if( eigval(ibnd,ik) >= ef ) exit
        enddo
        cbm = min( cbm, eigval(min(ibnd,nbnd),ik) )
        vbm = max( vbm, eigval(ibnd-1,ik) )
     enddo
     call mp_max( vbm )
     call mp_min( cbm )
     
     write(stdout,*) '          VBM = ', vbm*rytoev, ' eV'
     write(stdout,*) '          CBM = ', cbm*rytoev, ' eV'
     write(stdout,*) ' Fermi Energy = ', ef*rytoev, ' eV'
     
  enddo
  call mpi_file_close( fheigval, ierr )
  call mp_barrier
  call mp_global_end()  
  stop

  end program efermi

