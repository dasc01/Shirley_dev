! ---------------------------------------------------------------------- 
  subroutine pwmat( )
! ---------------------------------------------------------------------- 

  use constants, only : tpi
  USE io_global,  ONLY : stdout, ionode
  USE cell_base
  USE gvect  
  USE grid_dimensions,  ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk, diropn
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use mp,        only : mp_max, mp_sum
  USE mp_wave,   ONLY : mergewf
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  use fft_base, only : dffts
  use fft_interfaces, only : fwfft, invfft

  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset, &
                                pwmtxel_method, pwmtxel_cutoff

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  integer :: ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)

  integer :: ij, i, j, k, n
  complex(dp),allocatable :: jtmp(:), wtmp(:), psid(:)
  integer :: iunpwm, iunpwi, igwx, ixyz
  character(3) :: nd_nmbr_tmp
  integer :: nwordpwm
  logical :: exst

  integer :: ngk_pw, npwx_pw, ngk_pwg, igwx_pw
  integer,allocatable :: igk_l2g(:,:)
  real(dp),allocatable :: g_pw(:,:), g2kin_pw(:)
  complex(dp),allocatable :: evcG(:,:), pw_matrix(:,:)
  real(dp),allocatable :: r(:,:)
  complex(dp),allocatable :: expikr(:)

  integer,external :: freeunit


  WRITE( stdout, '(/5x,"Calling pwmat .... ",/)')
  write(stdout,*) ' npwx  = ', npwx
  write(stdout,*) ' npw   = ', npw
  write(stdout,*) ' nbnd  = ', nbnd
  write(stdout,*) ' nbndx = ', nbndx

  ! sort out which band subset we will work with
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

  if( ionode ) then
    write(stdout,*) ' generating plane-wave matrix elements in this basis'
    write(stdout,*) '   < B_i | exp(-iG.r) | B_j > '
  endif


  call flush_unit( stdout )

  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)

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

  ! I no longer have access to igk_l2g - need to regenerate
  if( allocated(igk_l2g) ) deallocate( igk_l2g )
  allocate( igk_l2g( npwx, nks ) )
  igk_l2g = 0

  call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! report norms
  allocate( norm(nbnd) )
  do ibnd=1,nbnd
    norm(ibnd) = dot_product( evc(:,ibnd_indx(ibnd)), evc(:,ibnd_indx(ibnd)) )
  enddo
  call mp_sum( norm, intra_pool_comm )
  do ibnd=1,nbnd
    if( abs(norm(ibnd)-1.d0) > eps ) then
      write(stdout,'(a,i6,a,f14.10)') ' band ', ibnd, ' norm = ', norm(ibnd)
      call errore('hamq','wave function norm is not 1',1)
    endif
  enddo
  deallocate( norm )
 

  ! ======================================================================
  ! plane-wave matrix elements
  ! ======================================================================
  write(stdout,*) ' plane-wave matrix elements'

  if( pwmtxel_method == 1 ) then

  ! ======================================================================
  write(stdout,*) ' pwmtxel_method I'
  ! ======================================================================

  if( ionode ) then
    iunpwm = freeunit()
    open(iunpwm,file=trim(prefix)//'.pwm',form='unformatted')
    write(stdout,*) ' will be saved to file: ', trim(prefix)//'.pwm'
    !
    iunpwi = freeunit()
    open(iunpwi,file=trim(prefix)//'.pwi',form='unformatted')
    write(stdout,*) ' useful info also saved to file: ', trim(prefix)//'.pwi'
  endif

  ! find total number of G vectors over processes
  igwx = maxval( igk_l2g(1:npw,1) )
  call mp_max( igwx )
  allocate( psid(nrxxs), wtmp(igwx), jtmp(npw) )


  ! ======================================================================
  ! dump info
  ! ======================================================================
  if( ionode ) then
    write(stdout,*) ' dumping info file'
    write(iunpwi) pwmtxel_cutoff ! cut-off
    write(iunpwi) igwx  ! number of G vectors
    write(iunpwi) nbnd  ! number of bands
    write(iunpwi) at    ! transformation matrix of bravais lattice
    write(iunpwi) bg    ! transformation matrix of reciprocal lattice
    write(iunpwi) omega ! cell volume
    write(iunpwi) tpiba ! 2 * pi / alat
  endif
  ! merge the distributed G vectors for writing
  do ixyz=1,3
    jtmp(1:npw) = cmplx(g(ixyz,igk(1:npw)),0.d0)
    wtmp = zero
    call mergewf(jtmp, wtmp, npw, igk_l2g(:,1), &
                 me_pool, nproc_pool, root_pool, intra_pool_comm)
    if( ionode ) write(iunpwi) real(wtmp) ! write ixyz terms from g
  enddo
  if( ionode ) close(iunpwi)


  ij=0
  do ibnd=1,nbnd
    !
    CALL start_clock( 'firstfft' )
    !
    psic(1:nrxxs) = ( 0.D0, 0.D0 )
    psic(nls(igk(1:npw))) = conjg( evc(1:npw,ibnd_indx(ibnd)) )
    !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
    CALL invfft ('Wave', psic, dffts)
    !
    CALL stop_clock( 'firstfft' )
    !
    ! ... product with another state 
    !
    do jbnd=1,nbnd
      !
      CALL start_clock( 'secondfft' )
      !
      psid(1:nrxxs) = ( 0.D0, 0.D0 )
      psid(nls(igk(1:npw))) = evc(1:npw,ibnd_indx(jbnd))
      !CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
      CALL invfft ('Wave', psid, dffts)
      !
      CALL stop_clock( 'secondfft' )
      !
      psid(1:nrxxs) = psic(1:nrxxs) * psid(1:nrxxs)
      !
      ! ... back to reciprocal space
      !
      CALL start_clock( 'secondfft' )
      !
      !CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
      CALL fwfft ('Wave', psid, dffts)
      !
      ! ... store with correct ordering
      !
      jtmp(1:npw) = psic(nls(igk(1:npw)))
      !
      CALL stop_clock( 'secondfft' )
      !
      ! now merge this set of plane-wave matrix elements for (ibnd,jbnd)
      !
      wtmp = zero
      call mergewf(jtmp, wtmp, npw, igk_l2g(:,1), &
                   me_pool, nproc_pool, root_pool, intra_pool_comm)
      !
      ! and dump to file
      !
      if( ionode ) then
        write(iunpwm) wtmp
        if( debug ) write(stdout,*) ibnd, jbnd, ' plane-wave matrix elements'
      endif
      !
    enddo
  enddo

  close(iunpwm)
  deallocate( wtmp, psid )

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( ibnd_indx, jtmp )


  else


  ! ======================================================================
  write(stdout,*) ' pwmtxel_method II'
  ! ======================================================================

  if( ionode ) then
    ! open output files
    nwordpwm = nbnd*nbnd*2
    iunpwm = freeunit()
    ! since only the ionode is going to write to this file
    ! there is no need for the nd_nmbr suffix
    nd_nmbr_tmp = nd_nmbr
    nd_nmbr=''
    ! open for direct access
    CALL diropn( iunpwm, 'pwm', nwordpwm, exst )
    ! restore in case needed again
    nd_nmbr = nd_nmbr_tmp
    write(stdout,*) ' will be saved to file: ', trim(tmp_dir)//trim(prefix)//'.pwm'

    iunpwi = freeunit()
    open(iunpwi,file=trim(tmp_dir)//trim(prefix)//'.pwi',form='unformatted')
    write(stdout,*) ' useful info also saved to file: ', trim(tmp_dir)//trim(prefix)//'.pwi'
  endif

  write(stdout,*) ' pwmtxel_cutoff = ', pwmtxel_cutoff
  
  ! determine the local number of plane-waves less then pwmtxel_cutoff
  call n_plane_waves (pwmtxel_cutoff, tpiba2, 1, xk, g, ngm, npwx_pw, ngk_pw)

  write(100+me_pool,*) ' ngk = ', ngk(1), ' ngk_pw = ', ngk_pw
  write(100+me_pool,*) ' npwx = ', npwx, ' npwx_pw = ', npwx_pw

  ! find the global number
  ngk_pwg = ngk_pw
  call mp_sum( ngk_pwg )

  write(stdout,*) ' nkg_pwg = ', ngk_pwg
  allocate( g_pw(3,ngk_pwg), g2kin_pw(ngk_pwg) )

  ! generate sorted local list
  CALL gk_sort( xk(1,1), ngm, g, pwmtxel_cutoff / tpiba2, ngk_pw, igk, g2kin )
  call gk_l2gmap (ngm, ig_l2g(1), ngk_pw, igk, igk_l2g(1,1))
  g2kin = g2kin * tpiba2

  do i=1,ngk_pw
    write(100+me_pool,'(i6,3f12.5,2x,f12.5)') i, g(1:3,igk(i)), g2kin(i)
  enddo

  ! merge vectors
  igwx_pw = maxval( igk_l2g(1:ngk_pw,1) )
  call mp_max( igwx_pw )
  allocate( jtmp(npwx_pw), wtmp(igwx_pw) )
  do ixyz=1,3
    jtmp(1:ngk_pw) = cmplx(g(ixyz,igk(1:ngk_pw)),0.d0)
    wtmp = zero
    call mergewf(jtmp, wtmp, ngk_pw, igk_l2g(:,1), &
                 me_pool, nproc_pool, root_pool, intra_pool_comm)
    g_pw(ixyz,1:ngk_pwg) = wtmp(1:ngk_pwg)
  enddo
  ! merge sq magn
  jtmp(1:ngk_pw) = cmplx(g2kin(1:ngk_pw),0.d0)
  wtmp = zero
  call mergewf(jtmp, wtmp, ngk_pw, igk_l2g(:,1), &
               me_pool, nproc_pool, root_pool, intra_pool_comm)
  g2kin_pw(1:ngk_pwg) = wtmp(1:ngk_pwg)

  if( ionode ) then
    write(iunpwi) pwmtxel_cutoff ! cut-off
    write(iunpwi) ngk_pwg  ! number of G vectors
    write(iunpwi) nbnd     ! number of bands
    write(iunpwi) at       ! transformation matrix of bravais lattice
    write(iunpwi) bg       ! transformation matrix of reciprocal lattice
    write(iunpwi) omega    ! cell volume
    write(iunpwi) tpiba    ! 2 * pi / alat

    write(iunpwi) g_pw  ! the G-vectors in Cartesian coords
    close(iunpwi)
  endif

  ! convert G-vectors back to lattice coords
  g_pw = matmul( transpose(at), g_pw )
  ! report
  write(stdout,*) ' selected G-vectors'
  do i=1,ngk_pwg
    write(stdout,'(i6,3f12.5,2x,f12.5)') i, g_pw(1:3,i), g2kin_pw(i)
  enddo
  
  ! allocate space
  allocate( pw_matrix(nbnd,nbnd), &
            evcG(npwx,nbnd) )

  ! make real-space grid for r
  write(stdout,*) ' rgrid '
  allocate( r(nrxx,3), expikr(nrxx) )
  call rgrid( r )
!  if( doublegrid ) then
!    do ixyz=1,3
!      call interpolate( r(1,ixyz), r(1,ixyz), -1 )
!    enddo
!  endif

  deallocate( jtmp, wtmp )
  allocate( wtmp(nrxx) )

  ! determine igk again
  CALL gk_sort (xk(1,1), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
  !call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
 
  ! reorder/resize evc in the correct way
  write(stdout,*) ' resizing evc '
  do ibnd=1,nbnd
    evcG(1:npw,ibnd) = evc(1:npw,ibnd_indx(ibnd))
    evcG(npw+1:npwx,ibnd) = zero
  enddo
  deallocate( evc )
  allocate( evc(npwx,nbnd) )
  evc = evcG

  ! loop over G-vectors
  do i=1,ngk_pwg

    ! plane-wave matrix element for G

    ! construct phase
    ! e^(i 2 pi xk_map.r )
    forall( j=1:nrxx ) expikr(j) = exp( iota &
      * cmplx( tpi * sum( r(j,:) * g_pw(:,i) ) ) )
!    expikr(:) = matmul( r, g_pw(:,i) )
!    expikr(:)= exp( ( iota * tpi ) * expikr(:) )

    ! loop over bands
    do ibnd=1,nbnd
      !
      CALL start_clock( 'firstfft' )
      !
      ! fft to real-space
      psic(1:nrxxs) = zero
      psic(nls(igk(1:npw))) = evc(1:npw,ibnd)
      !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
      CALL invfft ('Wave', psic, dffts)
      !
      CALL stop_clock( 'firstfft' )
      !
      ! ... product with the phase factor expikr =  exp(i xmap.r) 
      ! ... make the realspace product at higher resolution
      call cinterpolate( wtmp, psic, +1 )
      wtmp(1:nrxx) = wtmp(1:nrxx) * expikr(1:nrxx)
      call cinterpolate( wtmp, psic, -1 )
      !
      ! ... back to reciprocal space
      !
      CALL start_clock( 'secondfft' )
      !
      !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
      CALL fwfft ('Wave', psic, dffts)
      !
      ! ... store with correct ordering
      ! ... as evcG
      !
      evcG(1:npw,ibnd) = psic(nls(igk(1:npw)))
      evcG(npw+1:npwx,ibnd) = zero
      !
      CALL stop_clock( 'secondfft' )
      !
    enddo
    ! end loop over bands

    ! find matrix product of conjg(evcG) and evc = pw_matrix
    call ZGEMM( 'C', 'N', nbnd, nbnd, npw, one, &
                evcG, npwx, evc, npwx, zero, pw_matrix, nbnd ) 
    call mp_sum( pw_matrix, intra_pool_comm )

    ! write to file
    write(stdout,*) '     dumping pw_matrix for G-vector ', i
!    do jbnd=1,nbnd
!      do ibnd=1,nbnd
!        write(stdout,*) ibnd,jbnd, pw_matrix(ibnd,jbnd)
!      enddo
!    enddo
    if( ionode ) call davcio( pw_matrix, nwordpwm, iunpwm, i, +1 )

  enddo

  close(iunpwm)

  deallocate( wtmp )
  deallocate( r, expikr )
  deallocate( evcG )
  deallocate( igk_l2g )

  endif

  return

  end subroutine pwmat
! ---------------------------------------------------------------------- 

