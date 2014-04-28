! ---------------------------------------------------------------------- 
  subroutine pole_strength( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode
  USE constants, ONLY : rytoev, fpi, e2, tpi
  USE ions_base, ONLY : ntyp=>nsp, nat, ityp
  USE cell_base
  USE gvect  
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE scf,      ONLY : rho
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use ener, only : ehart, etxc, vtxc
  use mp,        only : mp_max, mp_sum, mp_bcast
  use mp_global, ONLY : mpime, nproc
  USE mp_wave,   ONLY : mergewf
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  use fft_base, only : dffts
  use fft_interfaces, only : fwfft, invfft

  !use hamq_shirley
  use splines_module
  use shirley_ham_input, only : nkgrid, ikgrid, ksplord, &
                                local_channel, debug, ncpp, &
                                band_subset

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one=cmplx(1.d0,0.d0)

  real(dp),allocatable :: xk_x(:), xk_y(:), xk_z(:)
  real(dp) :: xk_cart(3)
  integer,allocatable :: kmap(:,:)
  integer :: kxord, kyord, kzord, nk3(3)
  integer :: nkr

  integer :: ik, ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  real(dp) :: v0, charge, fac

  complex(dp),allocatable :: ham(:), dnl(:,:)
  integer,allocatable :: ikb_ham(:)
  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: eigvec(:,:)

  integer :: ierr, ixyz, ij, ig, i, j, k, n
  real(dp) :: qvec(3), tqvec(3)

  real(dp) :: q2min
  real(dp),allocatable :: q2min_all(:)
  integer :: igmin
  integer,allocatable :: igmin_all(:)
  integer :: minproc(1)

  complex(dp),allocatable :: rhoG(:)
  complex(dp),allocatable :: rtmp(:,:,:), jtmp(:), gtmp(:,:,:)
  complex(dp),allocatable :: pole(:,:)

  real(dp),allocatable :: qpG(:,:)
  integer :: iunpol, igwx
  logical :: exst
  real(dp) :: qg, q2
  complex(dp),external :: ZDOTC
  integer,external :: freeunit

  type(spline_dataset_type) :: pole_dataset
  type(spline_fit_type) :: pole_fit

  WRITE( stdout, '(/5x,"Calling pole_strength .... ",/)')

  ! ======================================================================
  ! sort out which band subset we will work with
  ! ======================================================================
  if( band_subset(1) > band_subset(2) ) then
    i=band_subset(2)
    band_subset(2) = band_subset(1)
    band_subset(1) = i
  endif
  if( band_subset(1)>=nbnd .or. band_subset(1)<=0 ) band_subset(1) = 1
  if( band_subset(2)>=nbnd .or. band_subset(2)<=0 ) band_subset(2) = nbnd

  if( band_subset(2)-band_subset(1)+1 < nbnd ) then
    write(stdout,*) ' Requested band subset:', band_subset(1), &
                    ' ... ', band_subset(2)
    write(stdout,*) ' Reducing total number of bands from ', nbnd, &
                    ' to ', band_subset(2)-band_subset(1)+1
  endif

  nbnd = band_subset(2)-band_subset(1)+1
  allocate( ibnd_indx(nbnd) )
  do i=1,nbnd
    ibnd_indx(i) = band_subset(1)+i-1
  enddo

  ! allocate space for spline data 
  call create_spline_dataset( nbnd, nbnd, nkgrid, ikgrid, &
                              bg*tpiba, transpose(at)/tpiba, &
                              pole_dataset )
  ! allocate space for spline fit
  call create_spline_fit( nbnd, nbnd, ksplord, pole_dataset, pole_fit )

  call flush_unit( stdout )

  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
  ! ======================================================================
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  ! ======================================================================
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)
!  !call sum_band
  write(stdout,*) ' ecutwfc = ', ecutwfc
  write(stdout,*) ' tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nktot = ', nks, nkstot
  write(stdout,*) ' xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '     npw = ', ngk(1:nks)
  write(stdout,*) '    npwx = ', npwx


!  ! gamma-point only
!  qvec = 0.d0
!  current_k = 1
!  if( lsda ) current_spin = isk(1)

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! ======================================================================
  ! allocate space for matrix elements
  ! ======================================================================
  allocate( pole(nbnd,nbnd), rhoG(npw), rtmp(nrxxs,nbnd,3), jtmp(npw), &
            gtmp(npw,nbnd,3), stat=ierr )
  if( ierr/=0 ) then
    call errore('pole_strength','problem allocating realspace wavefunctions ',1)
  endif

  write(stdout,*)
  write(stdout,*) ' construct pole_strength:'
  write(stdout,*)

  ! ======================================================================
  ! pole_strength matrix elements
  ! ======================================================================
  write(stdout,*) ' pole_strength matrix elements'

  ! k-point grid - uniform only
  write(stdout,*) ' generate nonlocal projectors on the following k-point grid:'
  write(stdout,'(2(2x,a,3i6))') 'nkgrid =', nkgrid(1:3), 'ikgrid =', ikgrid(1:3)

  nkr = numgridcoord_spline_dataset( pole_dataset )

  allocate( qpG(npw,3) )
  if( allocated(psic) ) deallocate( psic )
  allocate( psic(nrxxs) )

  allocate( q2min_all(nproc), igmin_all(nproc) )
  
  do ik=1,nkr

    call getcartcoord_spline_dataset( ik, xk_cart, pole_dataset )
    xk_cart = xk_cart / tpiba

    write(stdout,*) ' load pole_strength for ik = ', ik, ' of ', nkr
    write(stdout,'(x,a,3f12.5)') ' xk_cart = ', xk_cart(1:3)
    call flush_unit(stdout)

    ! I think one should find the minimum length q
    ! to make this function periodic in the Brillouin zone
    q2min = xk_cart(1)**2 + xk_cart(2)**2 + xk_cart(3)**2
    write(stdout,'(x,a,3f12.5,2x,f12.5)') '        xk_cart =    ', xk_cart, q2min
    igmin=0
    do ig = 1, npw
      qpG (ig,1) = xk_cart(1) + g(1, igk(ig) )
      qpG (ig,2) = xk_cart(2) + g(2, igk(ig) )
      qpG (ig,3) = xk_cart(3) + g(3, igk(ig) )
      q2 = qpG(ig,1)**2 +  qpG(ig,2)**2 + qpG(ig,3)**2

      if( q2 < q2min ) then
        q2min = q2
        igmin = ig
      endif
    enddo

    ! now find minimum over all processes
    q2min_all = 0.d0
    q2min_all(mpime+1) = q2min
    igmin_all = 0
    igmin_all(mpime+1) = igmin
    call mp_sum( q2min_all )
    call mp_sum( igmin_all )
    q2min = minval( q2min_all )
    minproc=minloc( q2min_all )-1
    if( mpime == minproc(1) .and. igmin/=0 ) then
      xk_cart = xk_cart + g(:,igk(igmin))
    endif
    call mp_bcast( xk_cart, minproc(1) )

    write(stdout,'(x,a,3f12.5,2x,f12.5)') ' minimal length xk = ', xk_cart, q2min
    write(stdout,*) ' from process ', minproc(1), ' position ', igmin_all(minproc(1)+1)

    ! now construct the lengths again
    do ig = gstart, npw
      qpG (ig,1) = xk_cart(1) + g(1, igk(ig) )
      qpG (ig,2) = xk_cart(2) + g(2, igk(ig) )
      qpG (ig,3) = xk_cart(3) + g(3, igk(ig) )
      qg = qpG(ig,1)**2 +  qpG(ig,2)**2 + qpG(ig,3)**2
      qg = sqrt( qg )
      qpG (ig,:) = qpG(ig,:) / qg
    enddo
    ! This point is an issue. It should be treated as a tensor with each
    ! limit as q->0 treated separately in Cartesian space.
    ! for now I'm fudging it by assuming an isotropic dependence
    if( gstart==2 .and. q2min < 1.d-12 ) then
      qg = sqrt( 3.d0 )
      qpG(1,:) = (/ 1.d0, 1.d0, 1.d0 /)
      qpG(1,:) = qpG(1,:) / qg
    else if( gstart==2 ) then
      qpG(1,:) = xk_cart(:) / sqrt( q2min )
    endif
      
    write(stdout,*) ' check on zeroth element'
    write(stdout,*) gstart, qpG(1,:)

    ! multiply by sqrt(q2*v) the rescaled interaction
    do ig=1,npw
      qpG(ig,:) = sqrt( fpi ) * qpG(ig,:)
    enddo

    do ixyz=1,3
    do ibnd=1,nbnd

      gtmp(1:npw,ibnd,ixyz) = evc(1:npw,ibnd_indx(ibnd)) * qpG(1:npw,ixyz)
      rtmp(1:nrxxs,ibnd,ixyz) = zero
      rtmp(nls(igk(1:npw)),ibnd,ixyz) = gtmp(1:npw,ibnd,ixyz)

      !CALL cft3s( rtmp(:,ibnd,ixyz), nr1s, nr2s, nr3s, &
      !                               nr1sx, nr2sx, nr3sx, 2 )
      CALL invfft ('Wave', rtmp(:,ibnd,ixyz), dffts)

    enddo
    enddo

!    psic(1:nrxxs) = cmplx( rho%of_r(1:nrxxs,1) )
!    if( nspin == 2 ) psic(1:nrxxs) = psic(1:nrxxs) + cmplx( rho%of_r(1:nrxxs,2) )
!    CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
!    rhoG(1:npw) = psic(nls(igk(1:npw)))
!    
!    do jbnd=1,nbnd
!    do ibnd=1,nbnd
!      do i=1,nrxxs
!        ! dot-product
!        psic(i) = sum( rtmp(i,ibnd,1:3) * conjg(rtmp(i,jbnd,1:3)) )
!      enddo
!      ! Fourier transform
!      CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
!      jtmp(1:npw) = psic(nls(igk(1:npw)))
!
!      ! integral with density
!      pole(ibnd,jbnd) = ZDOTC( npw, jtmp, 1, rhoG, 1 )
!    enddo
!    enddo

    ! load density into complex variable
    psic(1:nrxxs) = cmplx( rho%of_r(1:nrxxs,1) )
    if( nspin == 2 ) psic(1:nrxxs) = psic(1:nrxxs) + cmplx( rho%of_r(1:nrxxs,2) )
    
    ! multiply by real-space wave functions
    do ixyz=1,3
    do ibnd=1,nbnd
      rtmp(:,ibnd,ixyz) = rtmp(:,ibnd,ixyz) * psic(:)
    enddo
    enddo

    pole = ZERO
    do jbnd=1,nbnd
      do ixyz=1,3
        psic(1:nrxxs) = rtmp(1:nrxxs,jbnd,ixyz)
        ! Fourier transform
        !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
        CALL fwfft ('Wave', psic, dffts)
        jtmp(1:npw) = psic(nls(igk(1:npw)))

        call ZGEMV( 'C', npw, nbnd, ONE, gtmp(1,1,ixyz), npwx, &
                    jtmp, 1, ONE, pole(1,jbnd), 1 )
      enddo
    enddo
    call mp_sum( pole, intra_pool_comm )

    write(stdout,*) ' putelement_spline_dataset '
    write(stdout,*) pole(1,1)
    do jbnd=1,nbnd
    do ibnd=1,nbnd
      call putelement_spline_dataset( ibnd, jbnd, ik, pole(ibnd,jbnd), pole_dataset )
    enddo
    enddo

  enddo

  call fit_spline_to_dataset( pole_dataset, pole_fit )

  ! dump to disk
  if( ionode ) then
    iunpol = freeunit()
    call seqopn(iunpol,'polstr','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping pole_strength to file ', trim(prefix) // '.polstr'

    call dump_spline_dataset( pole_dataset )

    call write_spline_fit( iunpol, pole_fit )
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( pole )
  deallocate( ibnd_indx )

  write(stdout,*) ' Done with pole_strength'
  return

  end subroutine pole_strength
! ---------------------------------------------------------------------- 

