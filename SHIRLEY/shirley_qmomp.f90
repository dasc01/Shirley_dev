  program shirley_qmomp

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Feb 2007

  ! now parallelized
#include "f_defs.h"
  use hamq_shirley
  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_barrier, mp_end

  implicit none

!  integer,parameter :: stdout=6
  integer,parameter :: stdin =5
  integer,parameter :: dp = kind(1.d0)
  integer,parameter :: maxchar=255
  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  character(len=3) :: nodenumber
  character(maxchar) :: hamqfile
  character(maxchar) :: outfile
  character(maxchar) :: kpt_type
  character(maxchar) :: fmtstr
  integer :: nq, iq
  integer :: nqgrid(3), iqgrid(3)
  real(dp),allocatable :: qvec(:,:), wq(:)
  integer :: iqs
  integer :: nqtot
  integer,allocatable :: nq_local(:)
  integer,allocatable :: iqs_local(:)
  logical :: cartesian
  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: eigvec(:,:)
  complex(dp),allocatable :: mom(:,:,:)
  complex(dp),allocatable :: ztmp(:,:)
  integer :: iunhq, iunout, ierr

  integer :: nblock, lwork
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)

  integer,external :: freeunit
  integer,external :: ilaenv

  integer :: i,j,k
  integer :: ixyz

  namelist / input / hamqfile, outfile


  ! initialize mpi
  CALL start_shirley (nodenumber) 

! debugging
  if( .not. ionode ) stdout = 500+mpime
  call init_stdout( stdout )

  write(stdout,*) ' shirley_qmomp'
  write(stdout,*)

  if( ionode ) then

    read(stdin,nml=input,iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_qmomp','problem reading namelist &input',101)

    read(stdin,*) kpt_type ! type of coordinates
    kpt_type = adjustl(kpt_type)

    cartesian = .false.
    if( trim(kpt_type) == 'tpiba' .or. trim(kpt_type) == 'crystal' ) then
      read(stdin,*) nq
      allocate( qvec(3,nq), wq(nq) )
      do iq=1,nq
        read(stdin,*) qvec(1:3,iq), wq(iq)
      enddo
      if( trim(kpt_type) == 'tpiba' ) cartesian = .true.
    else if( trim(kpt_type) == 'automatic' ) then
      read(stdin,*) nqgrid, iqgrid
      if( any(iqgrid > 1) .or. any(iqgrid < 0) ) &
        call errore('shirley_qmomp','grid shift values should be 0 or 1',201)
      if( any(nqgrid < 1) ) &
        call errore('shirley_qmomp','grid values should be positive',202)
      nq = product(nqgrid)
      allocate( qvec(3,nq), wq(nq) )
      nq=0
      do i=1,nqgrid(1)
        do j=1,nqgrid(2)
          do k=1,nqgrid(3)
            nq=nq+1
            !  xk are the components of the complete grid in crystal axis
            qvec(1,nq) = DBLE(i-1)/nqgrid(1) - DBLE(iqgrid(1))/2/nqgrid(1)
            qvec(2,nq) = DBLE(j-1)/nqgrid(2) - DBLE(iqgrid(2))/2/nqgrid(2)
            qvec(3,nq) = DBLE(k-1)/nqgrid(3) - DBLE(iqgrid(3))/2/nqgrid(3)
          enddo
        enddo
      enddo
      wq = 2.d0/dble(nq)
    else
      write(stdout,*) ' kpoints flag unrecognized'
      write(stdout,*) ' should be: tpiba, crystal, or automatic'
      call errore('shirley_qmomp','stopping',1)
    endif
 
    ! now divide k-point set over processors
    nqtot = nq
    allocate( nq_local(nproc), iqs_local(nproc) )
    nq_local = nqtot / nproc
    do i=1,nqtot-nq_local(nproc)*nproc
      nq_local(i)=nq_local(i)+1
    enddo

    iqs_local(1)=0
    do i=2,nproc
      iqs_local(i)=iqs_local(i-1)+nq_local(i-1)
    enddo

  else  ! ionode

    allocate( nq_local(nproc), iqs_local(nproc) )

  endif ! ionode

  ! broadcast namelist input first
  call mp_bcast( hamqfile, ionode_id )
  call mp_bcast( outfile, ionode_id )

  ! broadcast k-point stuff next
  call mp_bcast( nqtot, ionode_id )
  call mp_bcast( cartesian, ionode_id )
  if( .not. ionode ) allocate( qvec(3,nqtot), wq(nqtot) )
  call mp_bcast( nq_local, ionode_id )
  call mp_bcast( iqs_local, ionode_id )
  call mp_bcast( qvec, ionode_id )
  nq = nq_local(mpime+1)
  iqs = iqs_local(mpime+1)
  
  write(stdout,*) ' # q-points         = ', nqtot
  write(stdout,*) ' # q-points locally = ', nq, ':', iqs+1, ' to ', iqs+nq


  if( ionode ) then
    iunhq = freeunit()
    open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) call errore('shirley_qmomp','problem opening file '//trim(hamqfile),102)
  endif

  iunout = freeunit()
  open(iunout,file=trim(outfile)//'.'//trim(nodenumber),form='formatted',iostat=ierr)
  if( ierr /= 0 ) call errore('shirley_qmomp','problem opening file '//trim(outfile)//'.'//trim(nodenumber),102)

  write(stdout,*) ' reading hamiltonian ...'
  if( ionode ) call read_hamq( iunhq )
  write(stdout,*) ' ... done'

  write(stdout,*) ' about to allocate eigenspace: ', nbasis
  allocate( eigval(nbasis), eigvec(nbasis,nbasis), &
            mom(nbasis,nbasis,3), ztmp(nbasis,nbasis), stat=ierr )
  if( ierr /= 0 ) call errore('shirley_qmomp','unable to allocate memory for eigen-problem',1)

  do iq=iqs+1,iqs+nq
    write(stdout,*)
    write(stdout,*) ' q-point ', iq-iqs, ' of ', nq

    call diag_hamq( qvec(1:3,iq), eigval, eigvec, cartesian )

    ! interpolate momentum components in shirley basis
    call mom_hamq( qvec(1:3,iq), cartesian, mom )
    
    ! expand in eigenbasis
    do ixyz=1,3
      CALL ZGEMM( 'N', 'N', nbasis, nbasis, nbasis, one, &
                  mom(1,1,ixyz), nbasis, eigvec, nbasis, zero, ztmp, nbasis )
      CALL ZGEMM( 'C', 'N', nbasis, nbasis, nbasis, one, &
                  eigvec, nbasis, ztmp, nbasis, zero, mom(1,1,ixyz), nbasis )
    enddo

    write(iunout,'(i,3f)') iq, qvec(1:3,iq)
    do j=1,nbasis
    do i=1,nbasis
      write(iunout,'(2i,2e,6e)') i, j, eigval(i), eigval(j), mom(i,j,1:3)
    enddo
    enddo
  enddo

  ! deallocate
  deallocate( qvec, wq )
  deallocate( nq_local, iqs_local )
  deallocate( eigvec, eigval, mom, ztmp )


  write(stdout,*) ' end shirley_qmomp'
  call mp_barrier
  call mp_end
  stop
  
  end program shirley_qmomp
