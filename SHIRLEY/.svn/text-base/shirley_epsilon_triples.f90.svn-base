  program shirley_epsilon_triples

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

  complex(dp),allocatable :: chimat(:,:), epsmat(:,:)
  complex(dp),allocatable :: zeigval(:)
  character(255) :: filename

  type(spline_fit_type) :: root_v_fit
  integer :: iunrvf, n, m
  complex(dp),allocatable :: root_v(:,:)
  type(spline_fit_type) :: qmin_fit
  complex(dp) :: qmin_array(1,1)
  real(dp) :: qmin

  complex(dp),allocatable :: dhdk(:,:,:)
  complex(dp),allocatable :: tridhdk(:,:)
  integer :: ideriv

  type(spline_dataset_type) :: eps_dataset
  type(spline_fit_type) :: eps_fit
  integer :: nq, ksplord(3)

  real(dp) :: fpi
  real(dp) :: epsilon_EC, epsilon_LF
  real(dp) :: kT

  real(dp) :: mem

!  complex(dp),allocatable :: B0vector(:)
!  real(dp) :: r2(2)
!  complex(dp) :: eps00

  ! system variables
  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba

  complex(dp),external :: ZDOTC, ZDOTU

  integer :: iuneps
  integer,external :: freeunit

  call start_clock( 'shirley' )

  fpi = acos(-1.d0)*4.d0

  call shirley_input

  ! close iunout as formatted 
  inquire(unit=iunout,name=filename)
  close(iunout,status='delete')

  ! load the root interaction fit
  if( nspline_fit < 1 ) &
    call errore('shirley_epsilon_triples','need the root interaction spline',1)
  if( ionode ) then
    iunrvf = freeunit()
    open(iunrvf,file=trim(spline_fit(1)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(1))
    call read_spline_fit( iunrvf, root_v_fit )
    call read_spline_fit( iunrvf, qmin_fit )
  endif
  call bcast_spline_fit( root_v_fit, mpime, ionode_id )
  call bcast_spline_fit( qmin_fit, mpime, ionode_id )

  ! allocate interaction space
  call getdims_spline_fit( n, m, root_v_fit )
  write(stdout,*) ' root_v size = ', n, m
  allocate( root_v(n,m) )

  call diag_init

  ! check
  if( n/=neig .and. m/=neig ) &
    call errore('shirley_epsilon_triples','size of root_v is incorrect',1)

!  ! B0vector
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

  ! scatter k-points
  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

!  write(stdout,*) '     cell volume = ', omega
  write(stdout,*) '          efermi = ', efermi
  write(stdout,*) '   electron temp = ', temperature

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

  ! convert efermi to Ha
  efermi = efermi * 0.5d0

  write(stdout,*) ' Total space to be allocated: '
  write(stdout,*) ' nproc = ', nproc
  write(stdout,*) ' neig = ', neig
  mem=dble(neig*neig)+3.d0*dble(neig)+dble(neig*neig)+dble((neig/nproc)*neig)+dble(neig*neig*(neig/nproc))+dble(neig*neig)
  mem=mem*16.d0/(1024.d0**3.d0)
  write(stdout,*) 'Memory = ', mem, 'GB'
!  WRITE( stdout, '(5X, &
!       &     "crystal axes: (cart. coord. in units of a_0)",/, &
!       &       3(15x,"a(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,  &
!       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
!  !
!  WRITE( stdout, '(5x, &
!       &   "reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
!       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
!       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)

  neig2=neig*neig
  allocate( eigvec_k(neig,neig), eigval_k(neig), &
            focc_k(neig), focc_kpq(neig), zeigval(neig) )


  allocate( epsmat(neig,neig), chimat(neig,neig), tri(neig,neig,neig) )

! read three-center integrals from file
 ! if( ionode ) then
    iuntri=freeunit()
    open(iuntri,file=trim(trifile),form='unformatted',access='direct', &
         recl=neig2*2*DIRECT_IO_FACTOR)
    do i=1,neig
      read(iuntri,rec=i) tri(:,:,i)
    enddo
 ! endif

  ! get system details from hamq_shirley
  call dump_system( nelec_, alat, omega, at, bg, tpiba )

  if( qpt%param%grid_type == 'automatic' ) then
  ! allocate space for epsilon spline
  call create_spline_dataset( neig, neig, qpt%param%nkgrid, qpt%param%ikgrid, &
                              bg*tpiba, transpose(at)/tpiba, &
                              eps_dataset )
  ksplord = 0
  call create_spline_fit( neig, neig, ksplord, eps_dataset, eps_fit )
  endif


  allocate( index_nnz(2,neig2) )

! ======================================================================
! loop over qvec to compute epsilon q
  if( qpt%param%grid_type == 'automatic' ) then
    nq = numgridcoord_spline_dataset( eps_dataset )
  else
    nq = qpt%list%nk
  endif
  do iq=1,nq
! ======================================================================
    if( qpt%param%grid_type == 'automatic' ) then
      call getgridcoord_spline_dataset( iq, qvec, eps_dataset )
    else
      qvec = qpt%list%kvec(:,iq)
    endif

    write(stdout,'(70("-"))')
    write(stdout,*) ' epsilon for q-point ', &
                      iq, ' of ', nq
    write(stdout,'(3f)') qvec
    write(stdout,*)

    ! convert qvec to Cartesian (units of tpiba)
    qvec_cart = matmul( trnlp2kin, qvec )
  
    ! evaluate spline of qmin
    call evaluate_spline_fit( qvec, qmin_fit, qmin_array )
    qmin = abs(qmin_array(1,1))
    write(stdout,'(a,f)') ' splined qmin = ', qmin

    if( abs(qmin) < 1.d-12 ) then
      if( all(deltaq < 1.d-12) ) then
        deltaq = 0.001d0
      endif
      write(stdout,*) ' I am cheating here by using a small deltaq = '
      write(stdout,*) deltaq 

      qvec = deltaq

      ! evaluate spline of qmin
      call evaluate_spline_fit( qvec, qmin_fit, qmin_array )
      qmin = abs(qmin_array(1,1))
      write(stdout,'(a,f)') ' splined qmin = ', qmin

!      allocate( dhdk(neig,neig,3), tridhdk(neig,neig) )
!      do ideriv=1,3
!        call diag_build_dham( ideriv, kpt%list%kvec(1:3,ik), kpt%param%cartesian, dhdk(:,:,ideriv) )
!      enddo
    endif

    chimat = ZERO

    ! loop over k
    do ik=1,kpt%list%nk
      write(stdout,*) ' construct contribution from k-point ', &
                      ik, ' of ', kpt%list%nk

      if( kpt%param%cartesian ) then
        kpqvec = kpt%list%kvec(:,ik) + qvec_cart
      else
        kpqvec = kpt%list%kvec(:,ik) + qvec
      endif
    
      ! ----------------------------------------------------------------------
      ! construct H(k) and H(k+q) and diagonalize
      ! ----------------------------------------------------------------------

      ! build hamiltonian for this k ...
      write(stdout,*) ' build ham k'
      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )
      ! ... and diagonalize
      write(stdout,*) ' diag ham k'
      call diag_ham
      ! store
      eigvec_k = eigvec
      eigval_k = eigval * 0.5d0 ! Hartree
      write(stdout,*) ' done with k'

      ! build hamiltonian for k+q ...
      write(stdout,*) ' build ham k+q'
      call diag_build_hamk( kpqvec, kpt%param%cartesian )
      ! ... and diagonalize
      write(stdout,*) ' diag ham k+q'
      call diag_ham
      eigval = eigval * 0.5d0 ! Hartree

      ! fermi occupation factors
      ! efermi read from input in epsilon_module
      do i=1,neig
        focc_k(i) = fermifunc( eigval_k(i), efermi, kT )
        focc_kpq(i) = fermifunc( eigval(i), efermi, kT )
      enddo

      ! determin non-zero occupation factors
      nnz=0
      do j=band_subset(1),band_subset(2)
        do i=1,j
          if( abs( focc_k(i) - focc_kpq(j) ) > 1.d-12 ) then
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
      allocate( wx(nnz) )

      do inz=1,nnz
        i = index_nnz(1,inz)
        j = index_nnz(2,inz)
        wx(inz) = cmplx( (focc_k(i) - focc_kpq(j))/(eigval(j)-eigval_k(i)) )
      enddo

      ! construct X for each local component
      write(stdout,*) 'Construct X'
      allocate( xivec(nnz,neig), xjvec(max(neig,nnz)) )

!      ! check for q=0
!      if( abs(qmin) > 1.d-12 ) then
 
      do k=1,neig
        !read(iuntri,rec=k) tri
        do inz=1,nnz
          i = index_nnz(1,inz)
          j = index_nnz(2,inz)
          call ZGEMV( 'N', neig, neig, ONE, tri(1,1,k), neig, &
                      eigvec(1,j), 1, ZERO, xjvec(1), 1 )
          xivec(inz,k) = ZDOTC( neig, eigvec_k(1,i), 1, xjvec(1), 1 )
        enddo
      enddo
      xivec = xivec / abs(qmin)

!      else
!
!      write(stdout,*) ' Warning: I am deliberately only using the first &
!component of the derivative of H wrt k assuming that the system dispersion is &
!isotropic'
!      do k=1,neig
!        !read(iuntri,rec=k) tri
!        call ZGEMM( 'N', 'N', neig, neig, neig, ONE, tri(1,1,k), neig, &
!                    dhdk(1,1,1), neig, ZERO, tridhdk(1,1), neig )
!        do inz=1,nnz
!          i = index_nnz(1,inz)
!          j = index_nnz(2,inz)
!          call ZGEMV( 'N', neig, neig, ONE, tridhdk(1,1), neig, &
!                      eigvec(1,j), 1, ZERO, xjvec(1), 1 )
!          xivec(inz,k) = ZDOTC( neig, eigvec_k(1,i), 1, xjvec(1), 1 )
!          xivec(inz,k) = xivec(inz,k) / (eigval(j)-eigval_k(i))
!        enddo
!      enddo
!
!      endif  ! qmin==0

      do j=1,neig
        xjvec(:) = wx(:) * conjg( xivec(:,j) )
        do i=1,neig
          ! factor of 2 for spin-degeneracy should be included in the kpt weights
          chimat(i,j) = chimat(i,j) + &
            kpt%list%wk(ik) * ZDOTU( nnz, xjvec(1), 1, xivec(1,i), 1 )
        enddo
      enddo
          
      write(stdout,*) chimat(1,1)
        
      deallocate( wx, xivec, xjvec )
      
! ======================================================================
    enddo ! loop over k-points ik
! ======================================================================

!    if( abs(qmin) < 1.d-12 ) then
!      deallocate( dhdk, tridhdk )
!    endif

    ! sum chimat over processors
    write(stdout,*) ' sum over processors'
    call mp_sum( chimat )

    if( ionode ) then

    write(stdout,*) chimat(1,1)
        
    write(stdout,*) ' interpolate the interaction '
    ! now generate the interaction
    call evaluate_spline_fit( qvec, root_v_fit, root_v )

    write(stdout,*) ' make epsilon'
    epsmat = chimat
    call ZGEMM( 'N', 'N', neig, neig, neig, ONE, &
                epsmat, neig, root_v, neig, ZERO, chimat, neig )
    call ZGEMM( 'N', 'N', neig, neig, neig, ONE, &
                root_v, neig, chimat, neig, ZERO, epsmat, neig )

    ! add unity
    forall( i=1:neig ) epsmat(i,i) = epsmat(i,i) + ONE
    !forall( i=1:neig ) epsmat(i,i) = epsmat(i,i) + qmin**2.d0

    write(stdout,*) epsmat(1,1)
    ! make the dielectric constant
!    eps00=zero
!    do j=1,neig
!    do i=1,neig
!      eps00 = eps00 + B0vector(i)*epsmat(i,j)*conjg(B0vector(j))
!    enddo
!    enddo
!    write(stdout,*) ' eps00 = ', eps00
!
    ! now diagonalize
    eigvec = epsmat
    call diag_ham

    write(stdout,*) ' dielectric eigenvalues and inverses:'
    do i=1,neig
      write(stdout,'(i,2f)') i, eigval(i), 1.d0/eigval(i)
    enddo

    ! dump in complex datablock format
!    zeigval = cmplx( eigval )
!    write(iunout) neig, 1
!    write(iunout) zeigval
!
!    write(iunout) neig, neig
!    write(iunout) eigvec
!
    if( qpt%param%grid_type == 'automatic' ) then
    ! store the epsilon matrix
    do j=1,neig
    do i=1,neig
      call putelement_spline_dataset( i, j, iq, epsmat(i,j), eps_dataset )
    enddo
    enddo
    endif
    
    endif ! ionode

! ======================================================================
  enddo ! loop over q-points iq
! ======================================================================

  ! dump to disk
  if( ionode ) then
    if( qpt%param%grid_type == 'automatic' ) then
    ! open epsilon file output
      iuneps = freeunit()
      open(iuneps,file=trim(epsfile),form='unformatted')

      call fit_spline_to_dataset( eps_dataset, eps_fit )
      write(stdout,*) ' dumping the epsilon matrix spline '
      call dump_spline_dataset( eps_dataset )
      call write_spline_fit( iuneps, eps_fit )

      close(iuneps)
    endif
  endif

  deallocate( index_nnz )
  deallocate( root_v )
  close(iuntri)

  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call stop_clock( 'shirley' )
  call stop_shirley
  
  end program shirley_epsilon_triples
