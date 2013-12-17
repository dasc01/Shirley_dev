  program shirley_proj

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and produce the eigensolutions
  ! dumping them in binary format

  ! David Prendergast, UCB, Feb 2008

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
!  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_end, mp_barrier
  use kpt_module
  use shirley_input_module
  use diag_module
  use splines_module

  implicit none

  type(spline_fit_type) :: proj_fit
  complex(dp),allocatable :: projk(:,:)

  character(255) :: filename

  complex(dp),allocatable :: zeigval(:)
  integer,external :: freeunit
  integer :: ik
  integer :: iunproj, n, m


  call shirley_input

  if( nspline_fit < 1 ) then
    call errore('shirley_proj','need to have at least one spline_fit for projections',1)
  endif

  ! load atomic projections
  if( ionode ) then
    iunproj = freeunit()
    open(iunproj,file=trim(spline_fit(1)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(1))
    call read_spline_fit( iunproj, proj_fit )
  endif
  call bcast_spline_fit( proj_fit, mpime, ionode_id )

  ! allocate projection space
  call getdims_spline_fit( n, m, proj_fit )
  write(stdout,*) ' proj size = ', n, m
  allocate( projk(n,m) )

  ! close iunout as formatted and reopen as unformatted
  inquire(unit=iunout,name=filename)
  close(iunout,status='delete')
  open(iunout,file=filename,form='unformatted')
  rewind(iunout)

  call diag_init

  write(stdout,*) ' shirley_proj'
  write(stdout,*)

  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

  allocate( zeigval(nbasis) )

! ======================================================================
  do ik=1,kpt%list%nk
! ======================================================================

    write(stdout,'(a,3f,i,a,i,a,a)') ' k-point ', &
      kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                  ' on node ', trim(nodenumber)

    write(iunout) kpt%list%kvec(1:3,ik), kpt%list%wk(ik)
    write(iunout) 2 ! eigenvalues and 1 set of eigenvectors plus projectors

    ! build the Hamiltonian for this q-point
    ! local contribution
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )

    ! diagonalize matrix to find eigenvalues and vectors
    call diag_ham

    ! interpolate projections
    call evaluate_spline_fit( kpt%list%kvec(1:3,ik), proj_fit, projk )

    ! expand in eigenbasis
    projk = matmul( projk, eigvec )

    ! write eigenvalues
    zeigval = cmplx( eigval )
    write(iunout) nbasis, 1
    write(iunout) zeigval

    write(iunout) n, m
    write(iunout) projk

    ! make sure iunout is written
    call flush(iunout)

! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! close binary dump
  close(iunout)

  deallocate( zeigval )

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_proj'
  stop
  
  end program shirley_proj
