  program shirley_opt

  ! Compute optical absorption using the Shirley basis

  ! David Prendergast, UCB, May 2007

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier
  use kpt_module
  use corerepair_module
  use shirley_input_module
  use diag_module

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  character(maxchar) :: fmtstr
  integer :: ik
  integer :: ierr

  integer :: iatom, it, np, npm
  complex(dp),allocatable :: mom(:,:,:), mtmp(:,:)

  integer,external :: freeunit

  integer :: i,j,k,ixyz


  call shirley_input

  call diag_init

! debugging
  if( .not. ionode ) stdout = 500+mpime
  call init_stdout( stdout )

  write(stdout,*) ' shirley_opt'
  write(stdout,*)

  call kpt_scatter( kpt%list, kpt%param, ionode_id )

  write(fmtstr,'(a,i,a)') '(i,4e,',nbasis,'e)'

  npm = maxval( nproj_type(1:ntype) )
  allocate( mom(nbasis,nbasis,3), mtmp(nbasis,nbasis) )

  do ik=1,kpt%list%nk

    write(stdout,*) ' k-point ', ik, ' of ', kpt%list%nk

    ! build the Hamiltonian for this k-point
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )

    ! diagonalize
    call diag_ham

    call mom_hamq( kpt%list%kvec(1:3,ik), kpt%param%cartesian, mom )

!    ! non-local contribution
!    if( nproj > 0 ) then
!      call build_hamq_projs( qvec(1:3,iq), cartesian, betaq )
!
!      do iatom=1,natom
!        it = type_atom(iatom)
!        np = nproj_type(it)
!        write(stdout,*) '     type_atom = ', it
!        write(stdout,*) '    nproj_type = ', np
!
!        forall( i=1:np, j=1:nbasis ) &
!          beta_block(i,j) = &
!            betaq(index_betaq(i,iatom),j)
!
!        forall( i=1:np, j=1:np ) &
!          cr(i,j) = cmplx( cr_atom(iatom)%matrix(i,j) )
!
!        CALL ZGEMM( 'N', 'N', np, nbasis, np, one, &
!                    cr, npm, &
!                    beta_block, npm, &
!                    zero, ztmp, npm )
!        CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, one, &
!                    beta_block, npm, &
!                    ztmp, npm, &
!                    one, mom, nbasis )
!      enddo
!
!      write(stdout,*) ' done with nonlocal'
!    endif

    ! now transform momentum matrix to this basis
    do ixyz=1,3
    CALL ZGEMM('N', 'N', nbasis, nbasis, nbasis, one, &
               mom(1,1,ixyz), nbasis, &
               eigvec, nbasis, &
               zero, mtmp, nbasis )
    CALL ZGEMM('C', 'N', nbasis, nbasis, nbasis, one, &
               eigvec, nbasis, &
               mtmp, nbasis, &
               zero, mom(1,1,ixyz), nbasis )
    enddo

    eigval = eigval * rytoev
    do i=1,nbasis
      do j=1,nbasis
        write(iunout,'(3i,f,6f)') ik, i, j, (eigval(j)-eigval(i)), mom(i,j,1:3)
      enddo
    enddo

    !  write(iunout,fmtstr) iq-1, matmul( trnlp2kin, qvec(1:3,iq)), eigval*rytoev
    
  enddo
!  write(stdout,*) ' done with k-points 1 to ', kpt%list%nk

  ! deallocate
  deallocate( mom, mtmp )

  ! close files
  close(iunout)

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_opt'
  stop

  end program shirley_opt
