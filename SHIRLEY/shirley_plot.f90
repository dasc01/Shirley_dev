  program shirley_plot

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list and plot some states

  ! David Prendergast, UCB, Dec 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end
  use kpt_module
  use corerepair_module
  use shirley_input_module
  use diag_module
  use plot_module

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  integer :: ik, ibnd, ierr

  integer :: i,j,k
  real(dp),allocatable :: eigvalk(:)
  complex(dp),allocatable :: eigveck(:,:)

  integer,external :: freeunit


  call shirley_input

  call diag_init

  write(stdout,*) kpt%list%nk

  write(stdout,*) plotspec%nplot
  write(stdout,'(2i)') plotspec%ink(1:2,1:plotspec%nplot)

  call plotspec_scatter( plotspec, ionode_id )

  allocate( eigvalk(plotspec%nplot), &
            eigveck(nbasis,plotspec%nplot), &
            stat=ierr )
  call errore( 'shirley_plot', &
    'unable to allocate space for plots - increase procs or decrease plots', &
    abs(ierr) )

  do i=1,plotspec%nplot
    write(stdout,*) ' plot ', i, ' of ', plotspec%nplot

    ibnd = plotspec%ink(1,i)
    ik   = plotspec%ink(2,i)

    if( i>1 .and. plotspec%ink(2,i-1) == ik ) then
      eigveck(:,i) = eigvec(:,ibnd)
      eigvalk(i) = eigval(i)
    else
      ! build the Hamiltonian for this k-point
      ! local contribution
      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )

      call diag_ham

      write(stdout,*) ik, eigval(ibnd)

      eigveck(:,i) = eigvec(:,ibnd)
      eigvalk(i) = eigval(i)
    endif

    ! dump
    write(iunout,'(2i)') ik, ibnd
    write(iunout,'(2e)') eigvalk(i), 0.d0
    write(iunout,'(6e)') eigveck(:,i)
  enddo

  write(stdout,*) ' done with plots ', 1, ' to ', plotspec%nplot
  close(iunout)

  write(stdout,*) ' end shirley_plot'
  call mp_end
  stop
  
  end program shirley_plot
