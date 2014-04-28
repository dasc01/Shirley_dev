! ---------------------------------------------------------------------- 
  subroutine interaction_test( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode, ionode_id
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
  use kpt_module

  implicit none

  integer,parameter :: stdin =5
  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  real(dp),allocatable :: xk_x(:), xk_y(:), xk_z(:)
  real(dp) :: xk_cart(3)
  integer,allocatable :: kmap(:,:)
  integer :: kxord, kyord, kzord, nk3(3)
  integer :: nkr

  complex(dp),allocatable :: becp(:,:)

  integer :: ik, jk, ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  real(dp) :: charge, fac

  integer :: ierr, ixyz, ij, ig, i, j, k, n
  complex(dp),allocatable :: gtmp(:)
  real(dp) :: qvec(3), tqvec(3), qpG(3)
  integer :: ig0

  real(dp),allocatable :: qg(:)
  integer :: iunint, igwx
  logical :: exst
  real(dp) :: q2
  complex(dp),external :: zdotc
  integer,external :: freeunit
  character(255) :: fmtstr

  type(spline_dataset_type) :: vint_dataset
  type(spline_fit_type) :: vint_fit
  type(kpt_type) :: qpt

  WRITE( stdout, '(/5x,"Calling interaction_test .... ",/)')

  if( ionode ) then
    call kpt_read( stdin, qpt )
  endif
  call kpt_bcast( qpt, ionode_id )

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

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, nwordwfc, iunwfc, 1, - 1 )

  write(stdout,*)
  write(stdout,*) ' construct interaction:'
  write(stdout,*)

  ! ======================================================================
  ! interaction matrix elements
  ! ======================================================================
  write(stdout,*) ' interaction matrix elements'

  if( ionode ) then
    iunint = freeunit()
    call seqopn(iunint,'intrxn','unformatted',exst)

    write(iunint) qpt%list%nk, nbnd
  endif

  allocate( qg(npw) )

  do ik=1,qpt%list%nk
    write(stdout,'(i8,3f8.4,f8.4)') ik, qpt%list%kvec(1:3,ik), qpt%list%wk(ik)

    xk_cart = qpt%list%kvec(1:3,ik)
    if( .not. qpt%param%cartesian ) then
      xk_cart = matmul( bg, xk_cart )
    endif

    if( ionode ) then
      write(iunint) xk_cart*tpiba2
    endif

    ! now construct the lengths again
    do ig = 1, npw
      qpG (1) = xk_cart(1) + g(1, igk(ig) )
      qpG (2) = xk_cart(2) + g(2, igk(ig) )
      qpG (3) = xk_cart(3) + g(3, igk(ig) )
      qg(ig) = qpG(1)**2 +  qpG(2)**2 + qpG(3)**2
      if( qg(ig) < 1.d-12 ) ig0 = ig
    enddo
    qg = qg * (tpiba**2.d0)

    do ig=1,npw
      if( ig==ig0 ) cycle
      qg(ig) = sqrt( fpi / qg(ig) )
    enddo
    qg(ig0) = 0.d0

    do jbnd=1,nbnd
      gtmp(1:npw) = evc(1:npw,ibnd_indx(jbnd)) * qg(1:npw)
      do ibnd=1,nbnd
        becp(ibnd,jbnd) = zdotc( npw, evc(1,ibnd_indx(ibnd)), 1, gtmp, 1 )
      enddo
    enddo
    call mp_sum( becp, intra_pool_comm )

    if( ionode ) then
      write(stdout,'(a,2e12.5)') ' v(1,1) = ', becp(1,1)
      write(iunint) becp
    endif

  enddo

  deallocate( qg )

  write(stdout,*) ' ...done'

  if( ionode ) then
    close(iunint)
    write(stdout,*)
    write(stdout,*) ' interaction dumped to file ', trim(prefix) // '.intrxn'
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( becp, gtmp, stat=ierr )
  deallocate( ibnd_indx )

  write(stdout,*) ' Done with interaction_test'
  return

  end subroutine interaction_test
! ---------------------------------------------------------------------- 

