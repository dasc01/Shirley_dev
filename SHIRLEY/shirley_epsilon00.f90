  program shirley_epsilon00

  ! Generate the head of the dielectric matrix using Shirley interpolation
  ! This is for q->0 only and involves velocity matrix elements

  ! David Prendergast, MF, Jun 2008

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum, mp_get, &
                 mp_scatter_size, mp_scatter
  use kpt_module
  use shirley_input_module
  use diag_module
  use fermi

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6

  COMPLEX(DP),PARAMETER :: ZERO=CMPLX(0.d0,0.d0)
  COMPLEX(DP),PARAMETER :: ONE =CMPLX(1.d0,0.d0)

  integer :: nnz, inz
  integer,allocatable :: index_nnz(:,:)
  integer :: neig2

  integer :: ik, i, j, k
  integer :: ipol, apol

  character(255) :: filename

  real(dp) :: fpi
  real(dp) :: kT
  complex(dp) :: fac

  complex(dp),allocatable :: pwk(:,:), ztmp(:,:)
  real(dp) :: dq, kpq(3)
  integer :: iq

  ! tolerance for determining zeros
  real(dp) :: tol

  real(dp) :: mem

  ! system variables
  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba

  integer :: iuneps
  integer,external :: freeunit


  call start_clock( 'shirley' )

  fpi = acos(-1.d0)*4.d0

  call shirley_input

  call diag_init

!  ! B0vector
!  ! we'll need this at the end to convert epsilon00 to Shirley basis
!  allocate( B0vector(neig) )
!  if( ionode ) then
!    read(stdin,*)
!    do i=1,neig
!      read(stdin,*) j, r2
!      B0vector(i) = cmplx( r2(1), r2(2) )
!    enddo
!  endif
!  call mp_bcast( B0vector, ionode_id )
!
!  write(stdout,*) ' B0vector'
!  do i=1,neig
!    write(stdout,*) i, B0vector(i)
!  enddo

  ! band_subset
  if( band_subset(1) > band_subset(2) ) band_subset = cshift( band_subset, 1 )
  if( band_subset(1) < 1 .or. band_subset(1) > neig ) band_subset(1)=1
  if( band_subset(2) < 1 .or. band_subset(2) > neig ) band_subset(2)=neig
  write(stdout,*) ' band_subset = ', band_subset

  ! scatter k-points
  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

  write(stdout,*) '          efermi [eV] = ', efermi
  write(stdout,*) '   electron temp  [K] = ', temperature

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

  write(stdout,*) ' Total space to be allocated: this is incorrect'
  write(stdout,*) ' nproc = ', nproc
  write(stdout,*) ' neig = ', neig
  mem=dble(neig*neig)+3.d0*dble(neig)+dble(neig*neig)+dble((neig/nproc)*neig)+dble(neig*neig*(neig/nproc))+dble(neig*neig)
  mem=mem*16.d0/(1024.d0**3.d0)
  write(stdout,*) 'Memory = ', mem, 'GB'

  ! allocate space
  neig2=neig*neig

  ! get system details from hamq_shirley
  call dump_system( nelec_, alat, omega, at, bg, tpiba )

  WRITE( stdout, '(5X, &
       &     "crystal axes: (cart. coord. in units of a_0)",/, &
       &       3(15x,"a(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,  &
       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  WRITE( stdout, '(5x, &
       &   "reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  write(stdout,*) '          volume = ', omega


! ======================================================================
! Begin calculation of epsilon00(q->0) tensor
! ======================================================================

  ! sampling dependent tolerance
  tol = 1.d-10/dble(kpt%list%nk)

  call head
!  call head_smallq

  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call stop_clock( 'shirley' )
  call stop_shirley
  
  contains


! ======================================================================
  subroutine head
! ======================================================================

  complex(dp),allocatable :: wx(:), de(:), dvec(:,:), fdvec(:,:)
  real(dp),allocatable :: focc(:)
  complex(dp) :: eps00(3,3)
  complex(dp),allocatable :: dhdk(:,:,:)
  integer :: ideriv

  allocate( index_nnz(2,neig2), focc(neig) )
  allocate( dhdk(neig,neig,3) )

  eps00 = zero

  ! loop over k
  do ik=1,kpt%list%nk
    write(stdout,*) ' construct contribution from k-point ', &
                    ik, ' of ', kpt%list%nk, kpt%list%wk(ik)

    ! ----------------------------------------------------------------------
    ! construct H(k) and diagonalize
    ! ----------------------------------------------------------------------
    write(stdout,*) ' build ham k'
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )
    write(stdout,*) ' diag ham k'
    call diag_ham
    write(stdout,*) ' done with k'

    ! ----------------------------------------------------------------------
    ! construct dH/dk(k) in the eigenbasis
    ! ----------------------------------------------------------------------
    do ideriv=1,3
      call diag_build_dhamk( ideriv, kpt%list%kvec(1:3,ik), &
                             kpt%param%cartesian, dhdk(:,:,ideriv) )
    enddo

    ! ----------------------------------------------------------------------
    ! fermi occupation factors
    ! efermi read from input in epsilon_module
    ! ----------------------------------------------------------------------
    do i=1,neig
      focc(i) = fermifunc( eigval(i), efermi, kT )
    enddo

    ! determin non-zero occupation factors
    nnz=0
    do j=band_subset(1),band_subset(2)
      do i=band_subset(1),j
        if( abs( focc(i) - focc(j) ) > tol ) then
          nnz=nnz+1
          index_nnz(:,nnz) = (/ i, j /)
          if( i/=j ) then
            nnz=nnz+1
            index_nnz(:,nnz) = (/ j, i /)
          endif
        endif
      enddo
    enddo

    write(stdout,*) ' non-zero occupation differences = ', nnz, &
                    ' out of possible ', neig2
    allocate( wx(nnz), de(nnz) )
    allocate( dvec(nnz,3), fdvec(nnz,3) )

    ! need a check here for de=zero
    do inz=1,nnz
      i = index_nnz(1,inz)
      j = index_nnz(2,inz)
      de(inz) = cmplx( eigval(j)-eigval(i) )
      wx(inz) = cmplx( focc(i) - focc(j) )/de(inz)* 2.0d0 ! Ry to Ha
    enddo

    ! make velocity matrix elements by dividing dhdk by de
    do ideriv=1,3
      do inz=1,nnz
        i = index_nnz(1,inz)
        j = index_nnz(2,inz)
        dvec(inz,ideriv) = conjg( dhdk(i,j,ideriv) )/de(inz)
        fdvec(inz,ideriv) = wx(inz)*dvec(inz,ideriv)
      enddo
    enddo

    ! factor of 2 for spin-degeneracy should be included in the kpt weights
    fac = cmplx( kpt%list%wk(ik)*fpi/omega )
    ! make the eps00 tensor incrementing the current estimate
    call ZGEMM( 'C', 'N', 3, 3, nnz, fac, &
                dvec, nnz, fdvec, nnz, ONE, eps00, 3 )

    write(stdout,*) ' esp00 = '
    write(stdout,'(6f)') eps00
    write(stdout,*)
        
    deallocate( wx, de, dvec, fdvec )
      
! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! sum eps00 over processors
  write(stdout,*) ' sum over processors'
  call mp_sum( eps00 )

  ! transform to dimensionless
  eps00 = matmul( transpose( trkin2nlp ), matmul( eps00, trkin2nlp ) )

  ! dump to file
  if( ionode ) then
    iuneps = freeunit()
    open(iuneps,file=trim(eps00file),form='unformatted')
    write(iuneps) eps00
    close(iuneps)
  endif

  ! add identity for reporting purposes
  do i=1,3
    eps00(i,i) = 1.d0 + eps00(i,i)
  enddo
  write(stdout,*) ' eps00 = '
  write(stdout,'(6f)') eps00

  deallocate( index_nnz, focc )
  deallocate( dhdk )

  end subroutine head


! ======================================================================
  subroutine head_smallq
! ======================================================================

  real(dp) :: eps0, qvec(3), qvec_cart(3), qvec2, kpq(3)
  complex(dp),allocatable :: eigvec_k(:,:), pwk(:,:)
  real(dp),allocatable ::  eigval_k(:)
  real(dp),allocatable :: focc_k(:), focc(:)
  complex(dp),allocatable :: wx(:), de(:), dvec(:), fdvec(:)
  complex(dp),external :: ZDOTC

  allocate( index_nnz(2,neig2), focc(neig), focc_k(neig) )
  allocate( eigval_k(neig), eigvec_k(neig,neig), pwk(neig,neig) )

  ! lattice coordinates for small q
  !qvec = (/ 0.001d0, 0.0d0, 0.0d0 /)
  qvec = (/ 0.0001d0, 0.0d0, 0.0d0 /)
  qvec_cart = matmul( bg, qvec ) * tpiba
  qvec2 = dot_product( qvec_cart, qvec_cart )

  eps0 = 0.d0

  ! loop over k
  do ik=1,kpt%list%nk
    write(stdout,*) ' construct contribution from k-point ', &
                    ik, ' of ', kpt%list%nk, kpt%list%wk(ik)

    ! ----------------------------------------------------------------------
    ! construct H(k) and diagonalize
    ! ----------------------------------------------------------------------
    write(stdout,*) ' build ham k'
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )
    write(stdout,*) ' diag ham k'
    call diag_ham
    write(stdout,*) ' done with k'
    eigvec_k = eigvec
    eigval_k = eigval

    ! ----------------------------------------------------------------------
    ! construct H(k+q) and diagonalize
    ! ----------------------------------------------------------------------
    kpq = kpt%list%kvec(1:3,ik)+qvec
    write(stdout,*) ' build ham k+q'
    call diag_build_hamk( kpq, kpt%param%cartesian )
    write(stdout,*) ' diag ham k+q'
    call diag_ham
    write(stdout,*) ' done with k+q'

    ! also construct the plane-wave matrix element
    call ZGEMM( 'C', 'N', neig, neig, neig, ONE, eigvec_k, neig, &
                eigvec, neig, ZERO, pwk, neig )

    ! ----------------------------------------------------------------------
    ! fermi occupation factors
    ! efermi read from input in epsilon_module
    ! ----------------------------------------------------------------------
    do i=1,neig
      focc_k(i) = fermifunc( eigval_k(i), efermi, kT )
      focc(i) = fermifunc( eigval(i), efermi, kT )
    enddo

    ! determin non-zero occupation factors
    nnz=0
    do j=band_subset(1),band_subset(2)
      do i=band_subset(1),j
        if( abs( focc_k(i) - focc(j) ) > tol ) then
          nnz=nnz+1
          index_nnz(:,nnz) = (/ i, j /)
          if( i/=j ) then
            nnz=nnz+1
            index_nnz(:,nnz) = (/ j, i /)
          endif
        endif
      enddo
    enddo

    write(stdout,*) ' non-zero occupation differences = ', nnz, &
                    ' out of possible ', neig2
    allocate( wx(nnz), de(nnz) )
    allocate( dvec(nnz), fdvec(nnz) )

    ! need a check here for de=zero
    do inz=1,nnz
      i = index_nnz(1,inz)
      j = index_nnz(2,inz)
      de(inz) = cmplx( eigval(j)-eigval_k(i) )
      wx(inz) = cmplx( focc_k(i) - focc(j) )/de(inz)*2.d0 ! Ry to Ha
    enddo

    ! make vectors
    do inz=1,nnz
      i = index_nnz(1,inz)
      j = index_nnz(2,inz)
      dvec(inz) = pwk(i,j)
      fdvec(inz) = pwk(i,j)*wx(inz)
    enddo

    ! factor of 2 for spin-degeneracy should be included in the kpt weights
    fac = cmplx( kpt%list%wk(ik)*fpi/omega )
    fac = fac / qvec2
    eps0 = eps0 + fac * ZDOTC( nnz, dvec, 1, fdvec, 1 )

    write(stdout,*) ' esp0_smallq = '
    write(stdout,'(2f)') eps0
    write(stdout,*)
        
    deallocate( wx, de, dvec, fdvec )
      
! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! sum eps00 over processors
  write(stdout,*) ' sum over processors'
  call mp_sum( eps0 )
  ! add identity
  eps0 = 1.d0 + eps0

  write(stdout,*) ' eps0 = '
  write(stdout,'(2f)') eps0

  deallocate( index_nnz, focc, focc_k )
  deallocate( eigval_k, eigvec_k, pwk )

  end subroutine head_smallq

  end program shirley_epsilon00
