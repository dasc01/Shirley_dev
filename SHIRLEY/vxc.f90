! ---------------------------------------------------------------------- 
  subroutine vxc( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode
  USE constants, ONLY : rytoev, fpi, e2, tpi
  USE ions_base, ONLY : ntyp=>nsp, nat, ityp
  USE cell_base
  USE gvect  
  USE grid_dimensions, ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions, ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE scf,      ONLY : vrs, vltot, v, rho, rho_core
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use ener, only : ehart, etxc, vtxc
  use mp,        only : mp_max
  use mp_global, ONLY : mpime
  USE mp_wave,   ONLY : mergewf
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool

  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)

  integer :: ik, ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  real(dp) :: vxc0, charge, fac

  complex(dp),allocatable :: ham(:)
  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: eigvec(:,:)

  complex(dp),allocatable :: vxc_mat(:,:)

  integer :: ierr, ixyz, ij, ig, i, j, k, n
  complex(dp),allocatable :: jtmp(:), wtmp(:), gtmp(:)
  real(dp) :: qvec(3), tqvec(3)
  integer :: iunvxc, igwx
  logical :: exst
  REAL(DP), EXTERNAL :: erf
  integer,external :: freeunit


  WRITE( stdout, '(/5x,"Calling vxc .... ",/)')
  write(stdout,*) ' npw = ', npw
  write(stdout,*) ' nbnd = ', nbnd, nbndx

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
  ! allocate space for Hamiltonian
  ! ======================================================================
  allocate( ham((nbnd*(nbnd+1))/2), &
            jtmp(npwx), gtmp(npwx), stat=ierr )
  if( ierr/=0 ) then
    call errore('hamq','problem allocating space for shirley hamiltonian',1)
  endif

  ! set up output for Hamiltonian module hamq_shirley
  call init_stdout( stdout )

  jtmp = zero
  gtmp = zero


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

!  ! I'm not sure if this is even necessary - it's just to be consistent
!  IF ( qcutz > 0.D0 ) THEN
!     !
!     DO ig = 1, npw
!        !
!        g2kin(ig) = g2kin(ig) + qcutz * &
!                    ( 1.D0 + erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
!        !
!     END DO
!     !
!  END IF
  !

  write(stdout,*)
  write(stdout,*) ' construct exchange-correlation potential:'
  write(stdout,*)

  ! ======================================================================
  ! exchange-correlation potential
  ! ======================================================================

  CALL v_xc( rho, rho_core, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
             nrxx, nl, ngm, g, nspin, alat, omega, etxc, vtxc, v%of_r )

  ham = zero
  vxc0 = sum( v%of_r )/dble(nrxxs) /omega
  call reduce( 1, vxc0 )
  write(stdout,*) ' vxc0 [Ry] = ', vxc0
  ij=0
  do jbnd=1,nbnd
    !
    CALL start_clock( 'firstfft' )
    !
    psic(1:nrxxs) = ( 0.D0, 0.D0 )
    psic(nls(igk(1:npw))) = evc(1:npw,ibnd_indx(jbnd))
    CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
    !
    CALL stop_clock( 'firstfft' )
    !
    ! ... product with the potential vrs = (vltot+v%of_r) on the smooth grid
    !
    psic(1:nrxxs) = psic(1:nrxxs) * v%of_r(1:nrxxs,current_spin)
    !
    ! ... back to reciprocal space
    !
    CALL start_clock( 'secondfft' )
    !
    CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
    !
    ! ... store with correct ordering
    !
    jtmp(1:npw) = psic(nls(igk(1:npw)))
    jtmp = conjg(jtmp)
    !
    CALL stop_clock( 'secondfft' )
    !
    do ibnd=1,jbnd
      ij=ij+1
      ham(ij) = sum( evc(1:npw,ibnd_indx(ibnd)) * jtmp(1:npw) )
    enddo
  enddo
  call reduce( 2*ij, ham )
  ham = conjg(ham)
  !! call dump_hamq( 304, nbnd, ham )

  call flush_unit(stdout)

  ! dump the exchange-correlation potential to disk
  if( ionode ) then
    iunvxc = freeunit()
    call seqopn(iunvxc,'vxc','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping XC potential to file ', trim(prefix) // '.vxc'

    allocate( vxc_mat(nbnd,nbnd) )
    ij=0
    do jbnd=1,nbnd
      do ibnd=1,jbnd
        ij=ij+1
        vxc_mat(ibnd,jbnd) = ham(ij)
        if( ibnd/=jbnd ) vxc_mat(jbnd,ibnd) = conjg(vxc_mat(ibnd,jbnd))
      enddo
    enddo

    ! probably modify this
    write(iunvxc) vxc_mat
    close(iunvxc)
    deallocate( vxc_mat )
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( ham, gtmp, stat=ierr )


  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( ibnd_indx, jtmp )

  write(stdout,*) ' Done with vxc'
  return

  end subroutine vxc
! ---------------------------------------------------------------------- 

