  program shirley_selfenergy_triples

  ! Generate the full self-energy matrix using Shirley interpolation
  ! but in the shirley basis (not Fourier space)

  ! David Prendergast, MF, Mar 2008

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

  integer :: iunvxc
  complex(dp),allocatable :: vxc(:,:), vxck(:,:)

  integer :: nnz, inz
  integer,allocatable :: index_nnz(:)
  complex(dp),allocatable :: xivec(:,:), xjvec(:,:), xtmp(:,:)
  complex(dp),allocatable :: sivec(:,:), sjvec(:,:)

  integer :: ierr
  integer :: neig2
  real(dp) :: kvec(3), kvec_cart(3)

  real(dp) :: kmqvec(3)
  complex(dp),allocatable :: eigvec_k(:,:), eigvec_kmq(:,:)
  real(dp),allocatable :: eigval_k(:), eigval_kmq(:)
  real(dp),allocatable :: focc_kmq(:)
  integer :: ik, iq, i, j, k, l, n, m

  complex(dp),allocatable :: zeigval(:)
  character(255) :: filename, fmtstr

  ! splines
  type(spline_fit_type) :: rtv_fit
  integer :: iunrvf
  complex(dp),allocatable :: rtv(:,:)

  type(spline_fit_type) :: eps_fit
  integer :: iuneps
  complex(dp),allocatable :: eps(:,:)
  integer :: ieps_start

  type(spline_fit_type) :: pol_fit
  integer :: iunpol
  complex(dp),allocatable :: pol(:,:)

  complex(dp),allocatable :: sigma_x(:,:), sigma_sx(:,:), sigma_c(:,:)
  complex(dp) :: sigma_l, sigma_ct, sigma_xt, sigma_sxt
  complex(dp),allocatable :: zqvec(:), zq(:)
  real(dp),allocatable :: wq(:), denom(:)
  real(dp),allocatable :: deqp(:,:), eqp(:)

  real(dp) :: fpi
  real(dp) :: kT

  real(dp) :: mem

!  complex(dp),allocatable :: B0vector(:)
!  real(dp) :: r2(2)
!  complex(dp) :: eps00

  ! system variables
  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba

  complex(dp),external :: ZDOTC, ZDOTU

  integer,external :: freeunit


  call start_clock( 'shirley' )

  fpi = acos(-1.d0)*4.d0

  call shirley_input

  ! scatter k-points
  call kpt_scatter( qpt%list, qpt%param, ionode_id ) 

  ! output for band structure
  close(iunout,status='delete')
  if( ionode ) open(iunout,file=trim(outfile),form='formatted')

! ----------------------------------------------------------------------
! splines
! ----------------------------------------------------------------------
  if( nspline_fit < 3 ) &
    call errore('shirley_selfenergy_triples','need 3 splines',1)

  ! load the root interaction fit
  if( ionode ) then
    iunrvf = freeunit()
    open(iunrvf,file=trim(spline_fit(1)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(1))
    call read_spline_fit( iunrvf, rtv_fit )
    close(iunrvf)
  endif
  call bcast_spline_fit( rtv_fit, mpime, ionode_id )

  ! allocate interaction space
  call getdims_spline_fit( n, m, rtv_fit )
  write(stdout,*) ' rtv size = ', n, m
  allocate( rtv(n,m) )

  ! load the epsilon matrix
  if( ionode ) then
    iuneps = freeunit()
    open(iuneps,file=trim(spline_fit(2)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(2))
    call read_spline_fit( iuneps, eps_fit )
    close(iuneps)
  endif
  call bcast_spline_fit( eps_fit, mpime, ionode_id )

  ! allocate interaction space
  call getdims_spline_fit( n, m, eps_fit )
  write(stdout,*) ' eps size = ', n, m
  allocate( eps(n,m) )

  ! load the pole-strength matrix
  if( ionode ) then
    iunpol = freeunit()
    open(iunpol,file=trim(spline_fit(3)),form='unformatted')
    write(stdout,*) ' opening file ', trim(spline_fit(3))
    call read_spline_fit( iunpol, pol_fit )
    close(iunpol)
  endif
  call bcast_spline_fit( pol_fit, mpime, ionode_id )

  ! allocate interaction space
  call getdims_spline_fit( n, m, pol_fit )
  write(stdout,*) ' pol size = ', n, m
  allocate( pol(n,m) )
! ----------------------------------------------------------------------

  call mp_barrier

  call diag_init

  ! check
  if( n/=neig .and. m/=neig ) &
    call errore('shirley_selfenergy_triples','size of splines is incorrect',1)


!  write(stdout,*) '     cell volume = ', omega
  write(stdout,*) '          efermi = ', efermi
  write(stdout,*) '   electron temp = ', temperature

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

  ! convert efermi to Ha
  efermi = efermi * 0.5d0

  neig2=neig*neig
  allocate( eigvec_k(neig,neig), eigval_k(neig), &
            eigvec_kmq(neig,neig), eigval_kmq(neig), &
            focc_kmq(neig) )
  allocate( tri(neig,neig,neig), vxc(neig,neig) )

! read three-center integrals from file
 ! if( ionode ) then
    iuntri=freeunit()
    open(iuntri,file=trim(trifile),form='unformatted',access='direct', &
         recl=neig2*2*DIRECT_IO_FACTOR)
    do i=1,neig
      read(iuntri,rec=i) tri(:,:,i)
    enddo
 ! endif

! read VXC
    iunvxc=freeunit()
    open(iunvxc,file=trim(vxcfile),form='unformatted') 
    read(iunvxc) vxc
    close(iunvxc)

  ! get system details from hamq_shirley
  call dump_system( nelec_, alat, omega, at, bg, tpiba )

  allocate( index_nnz(neig) )
  allocate( sigma_x( neig, neig ), sigma_sx( neig, neig ), sigma_c( neig, neig ) )
  allocate( vxck(neig,neig), xtmp(neig,neig) )

  allocate( xjvec(neig,neig), xivec(neig,neig) )
  allocate( sjvec(neig,neig), sivec(neig,neig) )
  allocate( zq(neig), zqvec(neig), wq(neig), denom(neig) )

! ======================================================================
! loop over kvec to compute selfenergy k
  do ik=1,kpt%list%nk
! ======================================================================
    kvec = kpt%list%kvec(:,ik)

    write(stdout,'(70("-"))')
    write(stdout,*) ' selfenergy for k-point ', &
                      ik, ' of ', kpt%list%nk
    write(stdout,'(3f)') kvec
    write(stdout,*)

    ! convert qvec to Cartesian (units of tpiba)
    kvec_cart = matmul( trnlp2kin, kvec )
  
    ! ----------------------------------------------------------------------
    ! construct H(k) and diagonalize
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

    ! construct vxck
    call ZGEMM( 'N', 'N', neig, neig, neig, ONE, vxc, neig, &
                eigvec_k, neig, ZERO, xtmp, neig )
    call ZGEMM( 'C', 'N', neig, neig, neig, ONE, eigvec_k, neig, &
                xtmp, neig, ZERO, vxck, neig )
    ! convert to Hartree
    vxck = vxck * 0.5d0
    write(stdout,*) ' Vxc at this k-point:'
    do i=band_subset(1),band_subset(2)
      write(stdout,'(i,f)') i, real(vxck(i,i))*2.d0*rytoev
    enddo

    sigma_x = ZERO
    sigma_sx = ZERO
    sigma_c = ZERO

    ! loop over q
    do iq=1,qpt%list%nk
      write(stdout,*) ' construct contribution from q-point ', &
                      iq, ' of ', qpt%list%nk

      if( qpt%param%cartesian ) then
        kmqvec = kvec_cart - qpt%list%kvec(:,iq)
      else
        kmqvec = kvec      - qpt%list%kvec(:,iq)
      endif
    
      ! ----------------------------------------------------------------------
      ! construct H(k-q) and diagonalize
      ! ----------------------------------------------------------------------

      ! build hamiltonian for k-q ...
      write(stdout,*) ' build ham k-q'
      call diag_build_hamk( kmqvec, qpt%param%cartesian )
      ! ... and diagonalize
      write(stdout,*) ' diag ham k-q'
      call diag_ham
      ! store
      eigvec_kmq = eigvec
      eigval_kmq = eigval * 0.5d0 ! Hartree
      write(stdout,*) ' done with k-q'

      ! fermi occupation factors
      ! efermi read from input in selfenergy_module
      do i=1,neig
        focc_kmq(i) = fermifunc( eigval_kmq(i), efermi, kT )
      enddo

      ! determin non-zero occupation factors
      nnz=0
      do i=1,neig
        if( focc_kmq(i) > 1.d-12 ) then
          nnz=nnz+1
          index_nnz(nnz) = i
        endif
      enddo

      write(stdout,*) ' non-zero occupations = ', nnz, &
                      ' out of possible ', neig

      ! load epsilon matrix at q
      call evaluate_spline_fit( qpt%list%kvec(:,iq), eps_fit, eps )

      ! load pole strength matrix at q
      call evaluate_spline_fit( qpt%list%kvec(:,iq), pol_fit, pol )

      ! load root interaction matrix at q
      call evaluate_spline_fit( qpt%list%kvec(:,iq), rtv_fit, rtv )

      ! diagonalize epsilon
      eigvec = eps
      call diag_ham

      write(stdout,*) ' index   eps(q)  1/eps(q)'
      do i=1,neig
        write(stdout,*) i, eigval(i), 1.d0/eigval(i)
      enddo
      ! truncate epsilon eigenspace
      do i=1,neig
        if( (1.d0-1.d0/eigval(i)) > epstol ) exit
      enddo
      ieps_start=min(i,neig)
      write(stdout,*) ' truncating epsilon eigenvalues from ', neig, ' to ', neig-ieps_start+1
      write(stdout,*) ' using cut-off tolerance of ', epstol
      write(stdout,*)
      ! expand pole strength in this basis
      write(stdout,*) ' pole strength'
      write(stdout,*) ' index 1/eps(q)  z(q) w(q) [eV]'
      do i=ieps_start,neig
        call ZGEMV( 'N', neig, neig, ONE, pol, neig, &
                    eigvec(1,i), 1, ZERO, zqvec, 1 )
        zq(i) = ZDOTC( neig, eigvec(1,i),1, zqvec, 1 )
        wq(i) = sqrt(real(zq(i))/(1.d0-1.d0/eigval(i)))

        write(stdout,'(i,f,2f,f)') i, 1.d0/eigval(i), zq(i), wq(i)*2.d0*rytoev
      enddo

      ! loop over band indices for self-energy
      do j=band_subset(1),band_subset(2)

        ! construct the triple xjvec
        do k=1,neig
          call ZGEMM( 'C', 'N', neig, neig, neig, ONE, eigvec_kmq, neig, &
                      tri(1,1,k), neig, ZERO, xtmp, neig )
          call ZGEMV( 'N', neig, neig, ONE, xtmp, neig, &
                      eigvec_k(1,j), 1, ZERO, xjvec(1,k), 1 )
        enddo
        ! multiply by root interaction
        xtmp = xjvec
        call ZGEMM( 'N', 'N', neig, neig, neig, ONE, xtmp, neig, &
                    rtv, neig, ZERO, xjvec, neig )
        ! generalized overlap
        call ZGEMM( 'N', 'N', neig, neig-ieps_start+1, neig, ONE, xjvec, neig, &
                    eigvec(1,ieps_start), neig, ZERO, sjvec(1,ieps_start), neig )

!        do i=band_subset(1),band_subset(2)
        i=j

          if( i/=j ) then
          ! construct the triple xivec
          do k=1,neig
            call ZGEMM( 'C', 'N', neig, neig, neig, ONE, eigvec_kmq, neig, &
                        tri(1,1,k), neig, ZERO, xtmp, neig )
            call ZGEMV( 'N', neig, neig, ONE, xtmp, neig, &
                        eigvec_k(1,i), 1, ZERO, xivec(1,k), 1 )
          enddo
          ! multiply by root interaction
          xtmp = xivec
          call ZGEMM( 'N', 'N', neig, neig, neig, ONE, xtmp, neig, &
                      rtv, neig, ZERO, xivec, neig )
          ! generalized overlap
          call ZGEMM( 'N', 'N', neig, neig-ieps_start+1, neig, ONE, xivec, neig, &
                      eigvec(1,ieps_start), neig, ZERO, sivec(1,ieps_start), neig )
          else
          ! copy from j
          xivec=xjvec
          sivec=sjvec
          endif

          sigma_xt = ZERO
          sigma_sxt = ZERO
          sigma_ct = ZERO

          ! exchange contribution to self-energy
          do inz=1,nnz
            k = index_nnz(inz)
            sigma_xt = sigma_xt - focc_kmq(k) &
                         * ZDOTC( neig, xivec(k,1), neig, xjvec(k,1), neig )
          enddo

          ! screened exchange contribution to self-energy
          do l=ieps_start,neig
            do inz=1,nnz
              k = index_nnz(inz)
              denom(k) = ( eigval_k(i) - eigval_kmq(k) )**2.d0 + &
                                 wq(l)**2.d0 
            enddo
            sigma_l = sum( conjg(sivec(1:nnz,l))*sivec(1:nnz,l)/denom(1:nnz) )
            sigma_sxt = sigma_sxt + zq(l) * sigma_l
          enddo

          ! correlation contribution to self-energy
          do l=ieps_start,neig
            do k=1,neig
              denom(k) = wq(l)*( eigval_k(i) - eigval_kmq(k) + &
                                 sign(wq(l),efermi-eigval_kmq(k)) )
            enddo
            sigma_l = sum( conjg(sivec(:,l))*sivec(:,l)/denom(:) )
            !write(stdout,'(i,f,2f)') l, real(zq(l)), sigma_l
            sigma_ct = sigma_ct + 0.5d0 * zq(l) * sigma_l
          enddo

          write(stdout,'(a,2i,6f)') ' sigma ', i, j, &
            eigval_k(i)*2.d0*rytoev, &
            real(sigma_xt)*4.d0*rytoev, real(sigma_sxt)*4.d0*rytoev, &
            real(sigma_ct)*4.d0*rytoev, &
            real(vxck(i,j))*2.d0*rytoev, & 
            real(sigma_xt+sigma_sxt+sigma_ct-0.5d0*vxck(i,j))*4.d0*rytoev
            ! extra factors of 2 above come from a missing sum over spin which
            ! doesn't appear till the weigths wk are applied just next...
          
          ! accumulate BZ average
          sigma_x(i,j) = sigma_x(i,j) + qpt%list%wk(iq)*sigma_xt
          sigma_sx(i,j) = sigma_sx(i,j) + qpt%list%wk(iq)*sigma_sxt
          sigma_c(i,j) = sigma_c(i,j) + qpt%list%wk(iq)*sigma_ct

    !    enddo ! band i

      enddo ! band j
!      
! ======================================================================
    enddo ! loop over q-points iq
! ======================================================================

    ! sum sigma over processors
    write(stdout,*) ' sum over processors'
    call mp_barrier
    call mp_sum( sigma_x )
    call mp_sum( sigma_sx )
    call mp_sum( sigma_c )

    if( ionode ) then
      allocate( deqp(neig,neig), eqp(neig) )

      write(stdout,*) ' i   j   E_LDA   DeltaE_QP   E_QP'
      do j=band_subset(1),band_subset(2)
      !do i=band_subset(1),band_subset(2)
        i=j
        deqp(i,j) = real(sigma_x(i,j)+sigma_sx(i,j)+sigma_c(i,j)-vxck(i,j))
        write(stdout,'(2i,3f)') i, j, eigval_k(i)*2.d0*rytoev, &
          deqp(i,j)*2.d0*rytoev, (eigval_k(i)+deqp(i,j))*2.d0*rytoev
      !enddo
        eqp(j) = eigval_k(j) + deqp(j,j)
      enddo
        
      if( kpt%param%grid_type == 'bandstructure' ) then
        write(fmtstr,'(a,i,a)') '(i,4e,',band_subset(2)-band_subset(1)+1,'e)'
        write(iunout,fmtstr) ik-1, matmul( trnlp2kin, kpt%list%kvec(1:3,ik)), &
          kpt%list%kpathlen(ik), &
          ( eqp(band_subset(1):band_subset(2)) )*2.d0*rytoev

        call shirley_output_bandstructure( 400, trim(outfile), &
               band_subset(2)-band_subset(1)+1, kpt%bandstructure, 'no' )
      endif

      deallocate( deqp, eqp )
    endif ! ionode

! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  deallocate( eigvec_k, eigval_k, &
              eigvec_kmq, eigval_kmq, &
              focc_kmq )
  deallocate( tri )

  deallocate( index_nnz )
  deallocate( rtv, eps, pol )
  deallocate( xtmp, vxck, vxc )
  deallocate( sigma_x, sigma_sx, sigma_c )
  deallocate( xivec, xjvec )
  deallocate( sivec, sjvec )
  deallocate( zq, zqvec, wq, denom )

  close(iuntri)

  write(stdout,*) ' waiting for other nodes'
  call stop_clock( 'shirley' )
  call stop_shirley
  
  end program shirley_selfenergy_triples
