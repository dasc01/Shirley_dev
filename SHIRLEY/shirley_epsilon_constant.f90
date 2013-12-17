  program shirley_epsilon_constant

  ! Generate the full dielectric matrix using Shirley interpolation
  ! but in the shirley basis (not Fourier space)

  ! David Prendergast, MF, Jan 2008

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum, mp_get, &
                 mp_scatter_size, mp_scatter
  use kpt_module
  use corerepair_module
  use shirley_input_module
  use diag_module
  use fermi
  use splines_module

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6

  COMPLEX(DP),PARAMETER :: ZERO=CMPLX(0.d0,0.d0)
  COMPLEX(DP),PARAMETER :: ONE =CMPLX(1.d0,0.d0)

  integer :: iuntri
  integer :: neig_tri
  complex(dp),allocatable :: tri(:,:,:)
  complex(dp),allocatable :: tri_mat(:,:), tri_tmp(:,:)

  integer :: nnz, inz
  integer,allocatable :: index_nnz(:,:)
  complex(dp),allocatable :: wx(:), xivec(:,:), xjvec(:)

  integer :: ierr
  integer :: neig2
  real(dp) :: qvec(3), qvec_cart(3)

  real(dp) :: kpqvec(3)
  complex(dp),allocatable :: eigvec_k(:,:)
  real(dp),allocatable :: eigval_k(:)
  real(dp),allocatable :: focc_k(:), focc_kpq(:)
  real(dp),allocatable :: dfocc(:)
  logical,allocatable :: wmask(:,:)
  real(dp) :: df
  integer :: ik, iq, i, j, k

  complex(dp),allocatable :: eps(:,:)
  complex(dp),allocatable :: zeigval(:)
  character(255) :: filename

  integer :: iuneps, n, m
  type(spline_dataset_type) :: eps_dataset
  type(spline_fit_type) :: eps_fit
  integer :: nq, ksplord(3)

  real(dp) :: kT

  real(dp) :: mem

  complex(dp),allocatable :: B0vector(:)
  real(dp) :: r2(2)
  complex(dp) :: eps00

  ! system variables
  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba

  complex(dp),external :: ZDOTC, ZDOTU

  integer,external :: freeunit


  call start_clock( 'shirley' )


  call shirley_input

  ! load the epsilon matrix
  if( ionode ) then
    iuneps = freeunit()
    open(iuneps,file=trim(spline_fit(1)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(1))
    call read_spline_fit( iuneps, eps_fit )
  endif
  call bcast_spline_fit( eps_fit, mpime, ionode_id )
  call mp_bcast( n, ionode_id )
  call mp_bcast( m, ionode_id )
  ! close iunout as formatted and reopen as unformatted

  inquire(unit=iunout,name=filename)
  close(iunout,status='delete')

  call diag_init

  ! B0vector
  allocate( B0vector(neig) )
  if( ionode ) then
    read(stdin,*)
    do i=1,neig
      read(stdin,*) j, r2
      B0vector(i) = cmplx( r2(1), r2(2) )
    enddo
  endif
  call mp_bcast( B0vector, ionode_id )

  write(stdout,*) ' B0vector'
  do i=1,neig
    write(stdout,*) i, B0vector(i)
  enddo

  ! get system details from hamq_shirley
  call dump_system( nelec_, alat, omega, at, bg, tpiba )

  write(stdout,*) '     cell volume = ', omega
  write(stdout,*) '          efermi = ', efermi
  write(stdout,*) '   electron temp = ', temperature

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

  ! convert efermi to Ha
  efermi = efermi * 0.5d0

  allocate( eps(neig,neig) )

  ! load the epsilon matrix
  if( ionode ) then
    iuneps = freeunit()
    open(iuneps,file=trim(spline_fit(1)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(1))
    call read_spline_fit( iuneps, eps_fit )
  endif
  call bcast_spline_fit( eps_fit, mpime, ionode_id )

! ======================================================================
! loop over qvec to compute epsilon q
  do iq=1,qpt%list%nk
! ======================================================================
    qvec = qpt%list%kvec(:,iq)

    write(stdout,'(70("-"))')
    write(stdout,*) ' epsilon for q-point ', &
                      iq, ' of ', qpt%list%nk
    write(stdout,'(3f)') qvec
    write(stdout,*)

    ! convert qvec to Cartesian (units of tpiba)
    qvec_cart = matmul( trnlp2kin, qvec )
  
    write(stdout,*) ' interpolate the dielectric matrix '
    ! now generate the interaction
    call evaluate_spline_fit( qvec, eps_fit, eps )

    ! make the dielectric constant
    write(stdout,*) ' eps11 = ', eps(1,1)
    do j=1,neig
      do i=1,neig
        write(500+iq,'(2i,2e)') i, j, real(eps(i,j)), dimag(eps(i,j))
      enddo
      write(500+iq,*)
    enddo

    eps00=zero
    do j=1,neig
    do i=1,neig
      eps00 = eps00 + conjg(B0vector(i))*eps(i,j)*(B0vector(j))
      write(stdout,*) i, j, eps00
    enddo
    enddo
    write(stdout,*) ' eps00 = ', eps00

! ======================================================================
  enddo ! loop over q-points iq
! ======================================================================


  close(iuntri)

  write(stdout,*) ' waiting for other nodes'
  call stop_clock( 'shirley' )
  call stop_shirley
  
  end program shirley_epsilon_constant
