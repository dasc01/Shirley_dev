  program shirley_mcsamp

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and generate dH/dk, i.e. derivatives
  ! of the Hamiltonian in k-space - useful matrix elements for
  ! computing optical conductivity

  ! David Prendergast, UCB, Feb 2008

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
!  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use kpt_module
  use corerepair_module
  use shirley_input_module
  use diag_module
  use fermi

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg=6.333630d-6
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp) :: kT

  character(maxchar) :: filename

  integer :: ik
  integer :: ibasis, jbasis, ideriv
  integer :: nbasis_subset

  integer :: iter_mc 
  real(dp) :: mcrand, rand3(3)
  integer :: ik3(3), idk3(3), kstep3(3)
  integer :: ik3_trial(3)
  real(dp) :: kvec_trial(3), dkvec(3)
  real(dp),allocatable :: pkvec(:)
  real(dp) :: pkvec_trial, pratio, ptot, ptot_samp
  real(dp),allocatable :: ptot_proc(:), ptot_samp_proc(:)
  logical :: cartesian
  integer :: i, j

  integer :: itotk


  call shirley_input

  ! rescale temperature
  kT = temperature * kelvin2rydberg
  efermi = efermi / rytoev

  write(stdout,*)
  write(stdout,*) ' electron temp = ', temperature, ' K'
  write(stdout,*) ' electron temp = ', kT*rytoev, ' eV'
  write(stdout,*) ' fermi energy = ', efermi*rytoev, ' eV'
  write(stdout,*)

  ! rescale spectral info
  spec_min = spec_min / rytoev
  spec_max = spec_max / rytoev
  spec_sigma = spec_sigma / rytoev

  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

!  ! close iunout and delete
!  inquire(unit=iunout,name=filename)
!  close(iunout,status='delete')

  call diagx_init( band_subset(1), band_subset(2) )

  write(stdout,*) ' shirley_mcsamp'
  write(stdout,*)

  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 


  cartesian=.false.
  do ik=1,kpt%list%nk
  ! convert k-points to non-cartesian
    if( kpt%param%cartesian ) then
      kpt%list%kvec(1:3,ik) = matmul( trkin2nlp, kpt%list%kvec(1:3,ik) )
    endif
    kpt%list%kvec(1:3,ik) = kpt%list%kvec(1:3,ik) &
                          - floor( kpt%list%kvec(1:3,ik) )
    write(stdout,'(a,3f)') ' k-point ', kpt%list%kvec(1:3,ik)
  enddo
  write(stdout,*)
    
  ! allocate space for eigenvalues and weights
  allocate( pkvec(kpt%list%nk) )

  write(stdout,*)
  write(stdout,*) ' equilibrium MC steps = ', niter_equil_mc
  write(stdout,*) '    sampling MC steps = ', niter_samp_mc
  write(stdout,*) '   radius of MC steps = ', mc_step_radius
  write(stdout,*)

  ! initialize random number generator
  call random_seed

  do ik=1,kpt%list%nk
    call diag_build_hamk( kpt%list%kvec(1:3,ik), cartesian )
    call diagx_ham
    ! probabilities
    call probk( eigval(band_subset(1):band_subset(2)), pkvec(ik) )
    write(stdout,'(a,3f,f)') '# k-point ', kpt%list%kvec(1:3,ik), pkvec(ik)
  enddo

  ptot=0.d0
  ptot_samp=0.d0
  itotk=0
  do iter_mc=-niter_equil_mc+1,niter_samp_mc

    write(stdout,*) ' iteration = ', iter_mc

    ! MC step

    ! loop over independent walkers
    do ik=1,kpt%list%nk
      call random_number( rand3 )
      dkvec(:) = (rand3(:)*2.d0-1.d0) * mc_step_radius
 
      kvec_trial = kpt%list%kvec(:,ik) + dkvec
      kvec_trial = kvec_trial - floor( kvec_trial )
      call diag_build_hamk( kvec_trial, cartesian )
      call diagx_ham

      call probk( eigval(band_subset(1):band_subset(2)), pkvec_trial )

    
      ! compute relative probability
      if( pkvec(ik) > 1.d-12 ) then
        pratio = min( 1.d0, pkvec_trial / pkvec(ik) )
      else
        pratio = 1.d0
!        if( pkvec_trial > 1.d-12 ) then
!          pratio = 1.d0
!        else
!          pratio = 0.d0
!        endif
      endif
    
      ! test
      call random_number( mcrand )
      if( pratio >= mcrand ) then

        if( iter_mc > 0 ) then
          itotk=itotk+1
          write(iunout,'(3f,f,2i)') kvec_trial, &
            1.d0/dble(kpt%param%nktot)/dble(niter_samp_mc), iter_mc, ik
        endif

        kpt%list%kvec(:,ik) = kvec_trial
        pkvec(ik) = pkvec_trial

        !sample with probability pratio
        ptot = ptot + pratio
        if( iter_mc > 0 ) ptot_samp = ptot_samp + pratio

!        if( iter_mc > 0 ) then
!          write(stdout,'(i,a,2f)') iter_mc, ' acceptance = ', &
!            ptot/dble(iter_mc+niter_equil_mc)/dble(kpt%list%nk), &
!            ptot_samp/dble(iter_mc)/dble(kpt%list%nk)
!        else
!          write(stdout,'(i,a,2f)') iter_mc, ' acceptance = ', &
!            ptot/dble(iter_mc+niter_equil_mc)/dble(kpt%list%nk)
!        endif
  
      endif

    enddo

  enddo

  allocate( ptot_proc(nproc), ptot_samp_proc(nproc) )
  ptot_proc = 0.d0
  ptot_proc(mpime+1) = ptot/dble(kpt%list%nk)
  ptot_samp_proc = 0.d0
  ptot_samp_proc(mpime+1) = ptot_samp/dble(kpt%list%nk)

  ! report acceptance
  call mp_sum( ptot_proc )
  ptot_proc = ptot_proc / dble(niter_equil_mc+niter_samp_mc)
  call mp_sum( ptot_samp_proc )
  ptot_samp_proc = ptot_samp_proc / dble(niter_samp_mc)

  if( ionode ) then
    write(stdout,*)
    write(stdout,'(2x,a5,a20,a20)') &
      'proc', 'total acceptance', 'sampled acceptance'
    do i=1,nproc
      write(stdout,'(2x,i5,f,f)') i, ptot_proc(i), ptot_samp_proc(i)
    enddo
    write(stdout,*)
  endif

  call mp_sum( ptot )
  ptot = ptot / dble(kpt%param%nktot)
  call mp_sum( ptot_samp )
  ptot_samp = ptot_samp / dble(kpt%param%nktot)
  call mp_sum( itotk )
  write(stdout,'(a,2f)') ' final acceptance = ', &
     ptot/dble(niter_equil_mc+niter_samp_mc), &
     ptot_samp/dble(niter_samp_mc)
  write(stdout,'(a,i)') ' accepted k-points = ', itotk

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_mcsamp'
  stop
  
  contains

  subroutine probk( eigval, p )

  ! sample the negative of the fermideriv to stay near the fermi surface
  ! (at low temperatures and at equilibrium)
  use shirley_constants, only : rytoev
  real(dp) :: p, dw, de, pexp
  real(dp) :: eigval(nbasis_subset)
  integer :: i,j

  p=0.d0
  do i=1,nbasis_subset
    pexp=-1.d0 * fermideriv( eigval(i), efermi, kT )
    p = max( p, pexp )
  enddo

  end subroutine probk

  end program shirley_mcsamp
