  program shirley_efermi

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Jan 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use shirley_constants, only : maxchar, rytoev, kelvin2rydberg
  use shirley_input_module
  use hamq_shirley
  use diag_module
  USE io_global,  ONLY : stdout
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_barrier, mp_end, mp_sum, mp_max, mp_min
  use fermi

  implicit none

  character(maxchar) :: fmtstr
  integer :: nq, iq, nqmax, nqmin
  real(dp),allocatable :: wq(:)
  integer :: iqs
  real(dp),allocatable :: eigvalk(:,:)
  integer,allocatable :: isk(:)
  integer :: ierr

  real(dp) :: ef, kT
  integer :: nspin

  integer,external :: freeunit

  integer :: i,j,k
  integer :: nbasis_subset

  real(dp),external :: efermit

  call shirley_input

  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis 
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

  call diagx_init( band_subset(1), band_subset(2) )

  write(stdout,*) ' shirley_efermi'
  write(stdout,*)

  ! divide up k-points among procs
  nq=kpt%list%nk/nproc
  iqs=nq*mpime
  if( mpime < kpt%list%nk-nq*nproc ) then
    iqs=iqs+mpime
  else
    iqs=iqs+kpt%list%nk-nq*nproc
  endif
  if( mpime < kpt%list%nk-nq*nproc ) nq=nq+1
  ! report
  nqmax = nq
  nqmin = nq
  call mp_max( nqmax )
  call mp_min( nqmin )
  write(stdout,*) ' total no. k-points = ', kpt%list%nk
  write(stdout,*) ' k-points/proc = ', dble(kpt%list%nk)/dble(nproc)
  write(stdout,*) ' max/proc = ', nqmax
  write(stdout,*) ' min/proc = ', nqmin


  ! attempt to allocate space for eigenvalues
  allocate( eigvalk(nbasis_subset,kpt%list%nk), wq(kpt%list%nk), stat=ierr )
  if( ierr/=0 ) call errore('shirley_efermi','unable to allocate space for eigenvalues and weights',abs(ierr))

  ! output formatting
  write(fmtstr,'(a,i,a)') '(i,3e,',nbasis_subset,'e)'

  ! loop over k-points
  eigvalk = 0.d0
  wq = 0.d0
  do iq=iqs+1,iqs+nq
    write(stdout,*) ' q-point ', iq-iqs, ' of ', nq, kpt%list%wk(iq)

    call diag_build_hamk( kpt%list%kvec(1:3,iq), kpt%param%cartesian )

    call diagx_ham

    eigvalk(1:nbasis_subset,iq) = eigval(band_subset(1):band_subset(2))

!    write(iunout,fmtstr) iq-1, kpt%list%kvec(1:3,iq), eigvalk(1:nbasis_subset,iq-iqs)*rytoev

    ! redefine local version of wq
    wq(iq) = kpt%list%wk(iq)
  enddo

  call mp_barrier

  if( kpt%param%grid_type == 'tetrahedra' ) then
    write(stdout,*) ' reduce to proc 0'
    call mp_sum( eigvalk )
    write(stdout,*) eigvalk(1,:)
  
    if( mpime==0 ) then
      nspin=1
      allocate( isk(kpt%list%nk) )
      isk=1
      ef = efermit (eigvalk, nbasis_subset, kpt%list%nk, nelec, nspin, kpt%tetra%ntetra, kpt%tetra%tetra, 0, isk)
    endif

  else

    ! convert input temperature to Ry
    kT = temperature * kelvin2rydberg

    ! determine fermi energy
    call fermi_energy( nelec, nq, nbasis_subset, wq(iqs+1), eigvalk(:,iqs+1), kT, smearing, ef )

  endif

  write(stdout,*) ' smearing = ', trim(smearing)
  write(stdout,*) ' Electron Temp = ', kT*rytoev, ' eV'
  write(stdout,*) ' Electron Temp = ', kT/kelvin2rydberg, ' K'
  write(stdout,*) ' Fermi Energy = ', ef*rytoev, ' eV'

  ! deallocate
  deallocate( wq, eigvalk )

  write(stdout,*) ' end shirley_efermi'
  call mp_barrier
  call mp_end
  stop
  
  end program shirley_efermi
