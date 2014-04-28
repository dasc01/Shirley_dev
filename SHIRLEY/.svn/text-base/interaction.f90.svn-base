! ---------------------------------------------------------------------- 
  subroutine interaction( )
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
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use ener, only : ehart, etxc, vtxc
  use mp,        only : mp_max, mp_sum, mp_bcast
  use mp_global, ONLY : mpime, nproc
  USE mp_wave,   ONLY : mergewf
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool

  !use hamq_shirley
  use splines_module
  use shirley_ham_input, only : nkgrid, ikgrid, ksplord, &
                                local_channel, debug, ncpp, &
                                band_subset

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  real(dp),allocatable :: xk_x(:), xk_y(:), xk_z(:)
  real(dp) :: xk_cart(3)
  integer,allocatable :: kmap(:,:)
  integer :: kxord, kyord, kzord, nk3(3)
  integer :: nkr

  complex(dp),allocatable :: becp(:,:)

  integer :: ik, ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  real(dp) :: v0, charge, fac

  integer :: ierr, ixyz, ij, ig, i, j, k, n
  complex(dp),allocatable :: gtmp(:)
  real(dp) :: qvec(3), tqvec(3), qpG(3)

  real(dp) :: q2min
  real(dp),allocatable :: q2min_all(:)
  integer :: igmin
  integer,allocatable :: igmin_all(:)
  integer :: minproc(1)

  real(dp),allocatable :: qg(:)
  integer :: iunint, igwx
  logical :: exst
  real(dp) :: q2
  complex(dp),external :: zdotc
  REAL(DP), EXTERNAL :: erf
  integer,external :: freeunit

  type(spline_dataset_type) :: vint_dataset
  type(spline_fit_type) :: vint_fit

  WRITE( stdout, '(/5x,"Calling interaction .... ",/)')

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
                              vint_dataset )
  ! allocate space for spline fit
  call create_spline_fit( nbnd, nbnd, ksplord, vint_dataset, vint_fit )

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


  ! ======================================================================
  ! allocate space for matrix elements
  ! ======================================================================
  allocate( becp(nbnd,nbnd), gtmp(npwx), stat=ierr )
  if( ierr/=0 ) then
    call errore('hamq','problem allocating space for shirley hamiltonian',1)
  endif

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

  write(stdout,*)
  write(stdout,*) ' construct interaction:'
  write(stdout,*)

  ! ======================================================================
  ! interaction matrix elements
  ! ======================================================================
  write(stdout,*) ' interaction matrix elements'

  ! k-point grid - uniform only
  write(stdout,*) ' generate interaction matrices on the following k-point grid:'
  write(stdout,'(2(2x,a,3i6))') 'nkgrid =', nkgrid(1:3), 'ikgrid =', ikgrid(1:3)

  nkr = numgridcoord_spline_dataset( vint_dataset )


  allocate( q2min_all(nproc), igmin_all(nproc) )
  allocate( qg(npw) )

  do ik=1,nkr

    call getcartcoord_spline_dataset( ik, xk_cart, vint_dataset )
    xk_cart = xk_cart / tpiba

    write(stdout,*)
    write(stdout,*) ' load interaction for ik = ', ik, ' of ', nkr
    call flush_unit(stdout)

    ! I think one should find the minimum length q
    ! to make this function periodic in the Brillouin zone
    q2min = xk_cart(1)**2 + xk_cart(2)**2 + xk_cart(3)**2
    write(stdout,'(x,a,3f12.5,2x,f12.5)') '        xk_cart =    ', xk_cart, q2min
!    igmin=0
!    do ig = 1, npw
!      qpG (1) = xk_cart(1) + g(1, igk(ig) )
!      qpG (2) = xk_cart(2) + g(2, igk(ig) )
!      qpG (3) = xk_cart(3) + g(3, igk(ig) )
!      q2 = qpG(1)**2 +  qpG(2)**2 + qpG(3)**2
!      
!      if( q2 < q2min ) then
!        q2min = q2
!        igmin = ig
!      endif
!    enddo
!
!    ! now find minimum over all processes
!    q2min_all = 0.d0
!    q2min_all(mpime+1) = q2min
!    igmin_all = 0
!    igmin_all(mpime+1) = igmin
!    call mp_sum( q2min_all )
!    call mp_sum( igmin_all )
!    q2min = minval( q2min_all )
!    minproc=minloc( q2min_all )-1
!    if( mpime == minproc(1) .and. igmin/=0 ) then
!      xk_cart = xk_cart + g(:,igk(igmin))
!    endif
!    call mp_bcast( xk_cart, minproc(1) )
!
!    write(stdout,'(x,a,3f,2x,f)') ' minimal length xk = ', xk_cart, q2min
!    write(stdout,*) ' from process ', minproc(1), ' position ', igmin_all(minproc(1)+1)

    ! now construct the lengths again
    do ig = 1, npw
      qpG (1) = xk_cart(1) + g(1, igk(ig) )
      qpG (2) = xk_cart(2) + g(2, igk(ig) )
      qpG (3) = xk_cart(3) + g(3, igk(ig) )
      qg(ig) = qpG(1)**2 +  qpG(2)**2 + qpG(3)**2
    enddo

    qg = qg * (tpiba**2.d0)
    do ig=gstart,npw
      qg(ig) = sqrt( fpi / qg(ig) / omega )
    enddo
    ! do not include G=0 term
    if( gstart==2 ) then
      qg(1)=0.d0
    endif

    write(stdout,*) ' accumulate becp '
    do jbnd=1,nbnd
      gtmp(1:npw) = evc(1:npw,ibnd_indx(jbnd)) * qg(1:npw)
      do ibnd=1,nbnd
        becp(ibnd,jbnd) = zdotc( npw, evc(1,ibnd_indx(ibnd)), 1, gtmp, 1 )
      enddo
    enddo
    call mp_sum( becp, intra_pool_comm )

    write(stdout,*) ' putelement_spline_dataset '
    call getgridcoord_spline_dataset( ik, xk_cart, vint_dataset )
    write(stdout,'(3f12.5)') xk_cart
    write(stdout,*) ik, becp(1,1)
    do jbnd=1,nbnd
    do ibnd=1,nbnd
      call putelement_spline_dataset( ibnd, jbnd, ik, becp(ibnd,jbnd), vint_dataset )
    enddo
    enddo

  enddo
  deallocate( qg )

  write(stdout,*) ' splining the data...'
  call fit_spline_to_dataset( vint_dataset, vint_fit )
  write(stdout,*) ' ...done'

  ! dump to disk
  if( ionode ) then
    iunint = freeunit()
    call seqopn(iunint,'intrxn','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping interaction to file ', trim(prefix) // '.intrxn'

    call dump_spline_dataset( vint_dataset )

    call write_spline_fit( iunint, vint_fit )
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( becp, gtmp, stat=ierr )
  deallocate( ibnd_indx )

  write(stdout,*) ' Done with interaction'
  return

  end subroutine interaction
! ---------------------------------------------------------------------- 

! ---------------------------------------------------------------------- 
  subroutine interaction0( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode
  USE cell_base
  USE gvect  
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
! davegp
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm

  !use hamq_shirley
  use shirley_ham_input, only : nkgrid, ikgrid, ksplord, &
                                local_channel, debug, ncpp, &
                                band_subset

  implicit none

  integer :: ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)

  integer :: ierr, ig, i
  complex(dp),allocatable :: becp(:,:)
  complex(dp),allocatable :: gtmp(:)

  integer :: iunint
  logical :: exst
  complex(dp),external :: zdotc
  integer,external :: freeunit


  WRITE( stdout, '(/5x,"Calling interaction0 .... ",/)')

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


  ! ======================================================================
  ! allocate space for matrix elements
  ! ======================================================================
  allocate( becp(nbnd,nbnd), gtmp(npwx), stat=ierr )
  if( ierr/=0 ) then
    call errore('hamq','problem allocating space for shirley hamiltonian',1)
  endif

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

  write(stdout,*)
  write(stdout,*) ' construct interaction0:'
  write(stdout,*)

  ! ======================================================================
  ! interaction0 matrix elements
  ! ======================================================================
  write(stdout,*) ' interaction0 matrix elements'

  do ig=gstart,npw
    g2kin(ig) = sqrt( 1.d0 / g2kin(ig) )
  enddo
  ! do not include G=0 term
  if( gstart==2 ) then
    g2kin(1)=0.d0
  endif

  write(stdout,*) ' accumulate becp '
  do jbnd=1,nbnd
    forall (i=1:npw) gtmp(i) = evc(i,ibnd_indx(jbnd)) * g2kin(i)
    do ibnd=1,nbnd
      becp(ibnd,jbnd) = zdotc( npw, evc(1,ibnd_indx(ibnd)), 1, gtmp, 1 )
    enddo
  enddo
  call mp_sum( becp, intra_pool_comm )

  ! dump to disk
  if( ionode ) then
    iunint = freeunit()
    call seqopn(iunint,'intrxn0','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping interaction0 to file ', trim(prefix) // '.intrxn0'

    write(iunint) becp
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( becp, gtmp )
  deallocate( ibnd_indx )

  write(stdout,*) ' Done with interaction0'
  return

  end subroutine interaction0
! ---------------------------------------------------------------------- 
