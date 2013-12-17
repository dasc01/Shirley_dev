  program xas_para

#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, world_comm
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2

  character(len=3) :: nodenumber

  integer :: narg

  integer :: nener
  real(dp) :: e1, e2, sigma, delta, efermi,estretch
  logical :: have_efermi, have_stretch
  real(dp),allocatable :: ener(:), spec(:,:,:),tmp_spec(:,:)
  complex(dp),allocatable :: xas(:,:,:)
  real(dp) :: xas_xyz(3)


  integer :: i, j, k, n, m, ik, ispin
  real(dp) :: de
  real(dp),allocatable :: eigval(:)

  logical :: cartesian
  character(255) :: grid_type
  real(dp),allocatable :: kvec(:,:), wk(:)
  real(dp) :: wktot, wkfac

  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma
  character(255) :: cdelta
  character(255) :: cefermi
  character(255) :: cstretch

  integer :: iunout
  character(255) :: fout
  integer :: iunstick
  character(255) :: stick_file, fmtstr
  integer :: iuninf, ierr
  integer :: fheigval, fhxmat
  character(255) :: eigval_file, xmat_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval, reqxmat

  integer,external :: freeunit
#ifdef __PGI
  integer,external :: iargc
#else
  integer,intrinsic :: iargc
#endif

  integer :: nk, nbnd, ncp
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  integer :: nspin
  logical :: lda_plus_u
!DAS-->
  logical::gwrun
  integer::becunit
  real(dp),allocatable :: xas_xyz_sp(:,:)

  namelist /info/ nk, nbnd, ncp, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u, gwrun

  data gwrun/.false./
  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then

  narg = iargc()
  if( narg /= 6 .and. narg /= 7 .and. narg /= 8) then
    write(stdout,*) ' usage: xas_para e1 e2 nener sigma delta [efermi] [stretch] filename'
    stop
  endif
  have_efermi=.false.
  have_stretch=.false.
  call getarg( 1, ce1 )
  call getarg( 2, ce2 )
  call getarg( 3, cnener )
  call getarg( 4, csigma )
  call getarg( 5, cdelta )
  if( narg==8 ) then
     have_efermi=.true.
     call getarg( 6, cefermi )
     have_stretch=.true.
     call getarg( 7, cstretch )
     call getarg( 8, filename )     
  else if( narg==7 ) then
     have_efermi=.true.
     call getarg( 6, cefermi )
     call getarg( 7, filename )
  else
     have_efermi=.false.
     call getarg( 6, filename )
  endif

  read(ce1,*) e1
  read(ce2,*) e2
  read(cnener,*) nener
  read(csigma,*) sigma
  read(cdelta,*) delta
  if( have_efermi ) then
    read(cefermi,*) efermi
  endif

  if( have_efermi ) then
    write(stdout,*) ' Excluding states with eigenenergies below efermi = ', efermi, ' eV'
  endif

  estretch=1.0d0
  if(have_stretch) then
    read(cstretch,*) estretch
  endif

  e1 = e1*evtory
  e2 = e2*evtory
  sigma = sigma*evtory
  delta = delta*evtory
  efermi = efermi*evtory

  endif

  call mp_bcast( e1, ionode_id )
  call mp_bcast( e2, ionode_id )
  call mp_bcast( nener, ionode_id )
  call mp_bcast( sigma, ionode_id )
  call mp_bcast( delta, ionode_id )
  call mp_bcast( have_efermi, ionode_id )
  call mp_bcast( efermi, ionode_id )
  call mp_bcast( have_stretch, ionode_id )
  call mp_bcast( estretch, ionode_id )
  call mp_bcast( filename, ionode_id )
!  write(*,*) efermi, estretch
  if( ionode ) then
    iunout=freeunit()
    fout=trim(filename)//'.xas'
    open(iunout,file=trim(fout),form='formatted')
    write(stdout,*) '    output in '//trim(fout)

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


    if(gwrun) then
       becunit=freeunit()
       open(becunit,file='bands.dat',form='formatted',status='old')
       read(becunit,'(i8)') nbnd
       close(becunit)
    endif
  endif
  call mp_bcast( nk, ionode_id )
  call mp_bcast( nbnd, ionode_id )
  call mp_bcast( ncp, ionode_id )
  call mp_bcast( nelec, ionode_id )
  call mp_bcast( alat, ionode_id )
  call mp_bcast( volume, ionode_id )
  call mp_bcast( at, ionode_id )
  call mp_bcast( bg, ionode_id )
  call mp_bcast( tpiba, ionode_id )
  call mp_bcast( fermi_energy, ionode_id )
  call mp_bcast( nspin, ionode_id )

  if( .not. ionode ) allocate( wk(nk), kvec(3,nk) )
  call mp_bcast( wk, ionode_id )
  call mp_bcast( cartesian, ionode_id )
  call mp_bcast( kvec, ionode_id )


  wktot = sum(wk)
  wkfac=2.d0/wktot


  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener,4,nspin) ,tmp_spec(nener,4))
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec(1:nener,1:4,1:nspin)=0.0d0
  tmp_spec(1:nener,1:4)=0.0d0


  ! allocate
  allocate( eigval(nbnd) )
  allocate( xas(nbnd,ncp,3) )

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'
  xmat_file=trim(filename)//'.xmat'

  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )
  call mp_file_open_dp( xmat_file, fhxmat, ionode_id, world_comm )


  iunstick=freeunit()
  stick_file=trim(filename)//'.stick'
  if( trim(nodenumber) /= '' ) then
    stick_file = trim(stick_file) // '.' // trim(nodenumber)
  endif
  open(iunstick,file=trim(stick_file),form='formatted')
  

  write(stdout,*) '    running...'
  if(.not. allocated(xas_xyz_sp)) allocate(xas_xyz_sp(3,nspin))
  do ispin=1,nspin
     do ik=1,nk
        if( mod(ik-1,nproc)/=mpime ) cycle
        
        write(stdout,*) ' reading ik= ', ik, ' of ', nk
        !write(*,*) kvec(:,ik), wk(ik)
        
        ! read eigenvalues
        offset = ((ispin-1)*nk + ik-1)*nbnd
        call mpi_file_read_at( fheigval, offset, &
             eigval, nbnd, &
             MPI_DOUBLE_PRECISION, status, ierr )
        ! read xas
        offset = ((ispin-1)*nk + ik-1)*nbnd*ncp*3*2
        call mpi_file_read_at( fhxmat, offset, &
             xas, nbnd*ncp*3*2, &
             MPI_DOUBLE_PRECISION, status, ierr )
        
        write(fmtstr,'(a,i6,a)') '(2i8,e20.12,e20.12,',3,'e20.12)'
        do i=1,nbnd
           if( have_efermi .and. eigval(i)<(efermi-sigma) ) cycle
           ! this line might need modification for L2/L3 branching ratio
           forall(j=1:3) xas_xyz_sp(j,ispin) = sum(conjg(xas(i,1:ncp,j))*xas(i,1:ncp,j))*wk(ik)
           ! I may want to split the spectrum over spin-channels
           if(have_stretch) eigval(i)=((eigval(i)-efermi)*estretch)+efermi
           tmp_spec(:,:)=0.0d0
           call add_gauss( nener, ener, eigval(i)+delta, &
                sigma, xas_xyz_sp(:,ispin), tmp_spec(1:nener,2:4) )
           spec(1:nener,2:4,ispin)=spec(1:nener,2:4,ispin)+tmp_spec(1:nener,2:4)
           write(iunstick,trim(fmtstr)) i, ik, (eigval(i)+delta)/evtory, sum(xas_xyz_sp(1:3,ispin))/3.d0, xas_xyz_sp(1:3,ispin)
        enddo
     enddo
     write(iunstick,*)
  enddo
  if(allocated(xas_xyz_sp)) deallocate(xas_xyz_sp)
  ! close 
  call mpi_file_close( fheigval, ierr )
  call mpi_file_close( fhxmat, ierr )

  close(iunstick)
  
  if( ionode ) write(stdout,*) '    waiting for other processors...'
  call mp_barrier

!  do i=1,nspin
  call mp_sum( spec ) !(:,:,ispin) )
!  enddo
  if( ionode ) then
    ! construct spatial average
     do ispin=1, nspin
        forall( i=1:nener) spec(i,1,ispin) = sum(spec(i,2:4,ispin))/3.d0
     enddo
    write(iunout,'(a,f12.5,a)') '# Applied a delta shift of ', delta/evtory, ' eV'
    do ispin=1, nspin
       do i=1,nener
          write(iunout,'(5e14.5e3)') ener(i)/evtory,spec(i,1:4,ispin)
       enddo
       write(iunout,*) 
    enddo
    close(iunout)
  endif

  call mp_barrier 
  call mp_end
  stop

  contains

    subroutine add_gauss( nener, ener, e, sigma, weight, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, weight(:), spec(:,:)

    real(dp) :: its2, pre(size(weight))
    real(dp) :: arg(nener)
    integer :: i, j

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2) * weight
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    forall( i=1:nener, j=1:size(weight) ) &
      spec(i,j) = spec(i,j) + pre(j) * exp( arg(i) )
     
    end subroutine add_gauss

  end program xas_para

