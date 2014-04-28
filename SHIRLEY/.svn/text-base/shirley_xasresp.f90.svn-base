  program shirley_xasresp

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and produce x-ray absorption spectra
  ! by averaging over a k-point grid

  ! David Prendergast, UCB, Jul 2007

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
!  use diag_shirley
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


  character(255) :: fmtstr, filename

  complex(dp),allocatable :: zeigval(:), phi_omega(:), beta_block(:,:)

  complex(dp),allocatable :: posn(:,:,:)

  real(dp) :: pi
  complex(dp) :: prefac

  integer,external :: freeunit
  integer,external :: ilaenv

  integer :: ik, i, j
  integer :: icore, ibasis, ncm, ixyz
  integer :: iatom, it, np, npm

  pi = acos(-1.d0)

  call shirley_input

  call diag_init

  write(stdout,*) ' shirley_xasresp'
  write(stdout,*)
  write(stdout,*) '     omega = ', omega_resp
  write(stdout,*) '       eta = ', eta_resp
  write(stdout,*) 'plot style = ', plot_style
  ! convert to Ry
  omega_resp = omega_resp / rytoev
  eta_resp = eta_resp / rytoev

  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

  write(fmtstr,'(a,i,a)') '(i,4e,',nbasis,'e)'

!  npm = maxval( nproj_type_nl(1:ntype) )
  npm = maxval( nproj_type(1:ntype) )
  allocate( zeigval(nbasis), phi_omega(nbasis), beta_block(npm,nbasis) )
  ncm = maxval( corerep%core(1:corerep%ncore)%nproj2 )
  allocate( posn(nbasis,ncm,3) )

! ======================================================================
  do ik=1,kpt%list%nk
! ======================================================================

    write(stdout,'(a,3f,i,a,i,a,a)') ' k-point ', &
      kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                  ' on node ', trim(nodenumber)

    ! build the Hamiltonian for this q-point
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )

!    ! diagonalize
!    call diag_ham

    ! add the energy and broadening to the Hamiltonian
    eigvec = -1.d0 * eigvec
    forall(i=1:nbasis) eigvec(i,i) = cmplx(omega_resp, eta_resp) + eigvec(i,i)

    prefac = cmplx( 2.d0 * sqrt( pi * eta_resp ) )
    
    ! construct the RHS
    posn = zero
    do icore=1,corerep%ncore

      iatom = corerep%core(icore)%atom
      it = type_atom(iatom)
      np = nproj_type(corerep%core(icore)%species)

      forall( i=1:np, j=1:nbasis ) &
        beta_block(i,j) = &
          betaq(index_betaq(i,iatom),j)


      ! sum <nk|beta><psi|r|phi>
      do ixyz=1,3
        CALL ZGEMM( 'C', 'N', nbasis, corerep%core(icore)%nproj2, np, prefac, &
                    beta_block, npm, &
                    corerep%core(icore)%matrix(1,1,ixyz), &
                    corerep%core(icore)%nproj1, &
                    zero, posn(1,1,ixyz), nbasis )
 
        do i=1,corerep%core(icore)%nproj2
          ! now call the GMRES solver
          call gmres( nbasis, eigvec, posn(:,i,ixyz), phi_omega, 5 )

          write(iunout,'(2i)') ik, i+(ixyz-1)*corerep%core(icore)%nproj2
          write(iunout,'(2e)') omega_resp, eta_resp
          write(iunout,'(6e)') phi_omega
        enddo

      enddo

    enddo
    
! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! close dump
  close(iunout)


  ! deallocate
  deallocate( phi_omega, posn )

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_xasresp'
  stop
  
  end program shirley_xasresp
