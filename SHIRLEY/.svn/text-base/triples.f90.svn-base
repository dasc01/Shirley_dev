! ---------------------------------------------------------------------- 
  subroutine triples( )
! ---------------------------------------------------------------------- 

  use constants, only : tpi
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base
  USE gvect  
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk, diropn
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use mp,        only : mp_max, mp_sum, mp_scatter_size, mp_scatter
  USE mp_wave,   ONLY : mergewf
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  use fft_base, only : dffts
  use fft_interfaces, only : fwfft, invfft

  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset, triples_subset, reduced_io

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  integer :: ibnd, jbnd, kbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)

  integer :: ij, i, j, k, n
  complex(dp),allocatable :: jtmp(:), psid(:)
  complex(dp),allocatable :: psic_all(:,:)
  complex(dp),allocatable :: tri_2(:,:)
  complex(dp),allocatable :: tri_3(:,:,:)
  integer :: iuntri, iunwfr
  integer :: nwordtri, nwordwfr
  logical :: exst
  character(3) :: nd_nmbr_tmp

  integer :: nbnd_proc
  integer,allocatable :: band(:), band_proc(:)

  integer :: ierr
  logical :: realspace_incore

  integer,external :: freeunit
  complex(dp),external :: ZDOTC


  WRITE( stdout, '(/5x,"Calling triples .... ",/)')
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
    write(stdout,*) ' generating triple integrals in this basis'
    write(stdout,*) '   < B_i | B_j* B_k > '
  endif

  if( triples_subset(1) < 1 .or. triples_subset(1) > nbnd ) &
    triples_subset(1)=1
  if( triples_subset(2) < 1 .or. triples_subset(2) > nbnd ) &
    triples_subset(2)=nbnd
  if( triples_subset(1) > triples_subset(2) ) then
    i=triples_subset(2)
    triples_subset(2) = triples_subset(1)
    triples_subset(1) = i
  endif
  if( ionode ) write(stdout,*) ' subset ', triples_subset

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
  ! I don't think I need this
  !call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  if( .not. reduced_io ) then
    ! scatter bands
    allocate( band(nbnd) )
    forall ( i=1:nbnd ) band(i)=i
    call mp_scatter_size( nbnd, nbnd_proc, ionode_id )
    allocate( band_proc(nbnd_proc) )
    call mp_scatter( band, band_proc, ionode_id )
  endif

  ! ======================================================================
  ! triple integrals
  ! ======================================================================
  write(stdout,*) ' triple integrals'

  if( reduced_io .and. ionode ) then
    nwordtri = nbnd*nbnd*2
    iuntri = freeunit()
    ! open for direct access
    nd_nmbr_tmp = nd_nmbr
    nd_nmbr=''
    CALL diropn( iuntri, 'tri', nwordtri, exst )
    nd_nmbr = nd_nmbr_tmp
    write(stdout,*) ' will be saved to file: ', trim(prefix)//'.tri'
  endif

  if( .not. reduced_io ) then
    nwordtri = nbnd*nbnd*2
    iuntri = freeunit()
    ! open for direct access
    nd_nmbr_tmp = nd_nmbr
    nd_nmbr=''
    CALL diropn( iuntri, 'tri', nwordtri, exst )
    nd_nmbr = nd_nmbr_tmp
    write(stdout,*) ' will be saved to file: ', trim(prefix)//'.tri'
  endif

  ! space
  allocate( psid(nrxxs), tri_2(nbnd,nbnd), jtmp(npw) )


  if( reduced_io ) then

    write(stdout,*)
    write(stdout,*) ' local memory:'
    write(stdout,*) ' memory required for real-space wave functions = ', &
      dble(nrxxs)*dble(nbnd)*16.d0/1024.d0**2.d0, 'MB'
    write(stdout,*) ' memory required for                   triples = ', &
      dble(nbnd)**3.d0*16.d0/1024.d0**2.d0, 'MB'
    write(stdout,*)
  
    allocate( psic_all(nrxxs,nbnd), &
              tri_3(nbnd,nbnd,nbnd), stat=ierr )
    if( ierr /= 0 ) then
      write(stdout,*) ' unable to allocate space for all real-space wave functions and triples'
      write(stdout,*) ' you might consider setting reduced_io to .false.'
      call errore('triples','Allocation error',abs(ierr))
    endif
  
    if( ionode ) write(stdout,*) ' loading real-space basis functions '
      !
      do ibnd=1,nbnd
        !
        CALL start_clock( 'firstfft' )
        !
        psic_all(1:nrxxs,ibnd) = ( 0.D0, 0.D0 )
        psic_all(nls(igk(1:npw)),ibnd) = conjg( evc(1:npw,ibnd_indx(ibnd)) )
!        CALL cft3s( psic_all(:,ibnd), nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
        CALL invfft ('Wave', psic_all(:,ibnd), dffts)
        !
        CALL stop_clock( 'firstfft' )
        !
      enddo
  
  else ! not reduced_io

    write(stdout,*)
    write(stdout,*) ' local memory:'
    write(stdout,*) ' memory required for                   triples = ', &
      dble(nbnd)**2.d0*dble(nbnd_proc)*16.d0/1024.d0**2.d0, 'MB'
    write(stdout,*) ' memory required for real-space wave functions = ', &
      dble(nrxxs)*dble(nbnd)*16.d0/1024.d0**2.d0, 'MB'
    write(stdout,*)
  
    ! try to allocate a portion of tri_3 locally
    allocate( psic_all(nrxxs,nbnd), &
              tri_3(nbnd,nbnd,nbnd_proc), stat=ierr )
    if( ierr/=0 ) then
      write(stdout,*) 'unable to allocate space for local copy of triples'
      call errore('triples','Allocation error',abs(ierr))
    endif
  
    realspace_incore = .true.
    allocate( psic_all(nrxxs,nbnd), stat=ierr )
    if( ierr/=0 ) then
      realspace_incore = .false.
      write(stdout,*) 'unable to allocate space for realspace wave functions'
  
      ! generate real-space basis functions and dump to file
      nwordwfr = nrxxs*2
      iunwfr = freeunit()
      ! open for direct access
      CALL diropn( iunwfr, 'wfr', nwordwfr, exst )
  
      if( .not. exst ) then
        if( ionode ) write(stdout,*) ' dumping real-space basis functions '
        !
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
          CALL davcio( psic, nwordwfr, iunwfr, ibnd, + 1 )
          !
        enddo
      endif
    else
      if( ionode ) write(stdout,*) ' loading real-space basis functions '
      !
      do ibnd=1,nbnd
        !
        CALL start_clock( 'firstfft' )
        !
        psic_all(1:nrxxs,ibnd) = ( 0.D0, 0.D0 )
        psic_all(nls(igk(1:npw)),ibnd) = conjg( evc(1:npw,ibnd_indx(ibnd)) )
        !CALL cft3s( psic_all(:,ibnd), nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
        CALL invfft ('Wave', psic_all(:,ibnd), dffts)
        !
        CALL stop_clock( 'firstfft' )
        !
      enddo
    endif

  endif ! reduced_io


  if( ionode ) write(stdout,*) ' now the integrals'

  do kbnd=1,nbnd
    !
    if( reduced_io .or. realspace_incore ) then
      psic = psic_all(:,kbnd)
    else
      call davcio( psic, nwordwfr, iunwfr, kbnd, - 1 )
    endif
    psic = conjg( psic )

    do jbnd=1,nbnd
      !
      if( reduced_io .or. realspace_incore ) then
        psid = psic_all(:,jbnd)
      else
        call davcio( psid, nwordwfr, iunwfr, jbnd, - 1 )
      endif
      psid = psid * psic
      !
      CALL start_clock( 'secondfft' )
      !
      !CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
      CALL fwfft ('Wave', psid, dffts)
      !
      ! ... store with correct ordering
      !
      jtmp(1:npw) = psid(nls(igk(1:npw)))
      !
      CALL stop_clock( 'secondfft' )
      !
      ! matrix product of pair fourier transform with eigenvectors
      !
      call ZGEMV( 'C', npw, nbnd, one, &
                  evc, npwx, jtmp, 1, zero, tri_2(1,jbnd), 1 ) 
      !
    enddo ! jbnd
    !
    call mp_sum( tri_2, intra_pool_comm )
    !
    if( reduced_io .and. ionode ) then
      forall( j=1:nbnd, i=1:nbnd ) tri_3(j,kbnd,i) = tri_2(i,j)
    endif
    !
    if( .not. reduced_io ) then
      forall( j=1:nbnd, i=1:nbnd_proc ) &
        tri_3(j,kbnd,i) = tri_2(band_proc(i),j)
    endif
    !
    if( ionode ) write(stdout,*) ' triple subset ', kbnd, ' of ', nbnd
    !
  enddo ! kbnd 

  if( reduced_io .and. ionode ) then
    write(stdout,*) ' dumping all of triples'
    do ibnd=1,nbnd
      !
      CALL davcio( tri_3(:,:,ibnd), nwordtri, iuntri, ibnd, + 1 )
      !
    enddo
    close(iuntri)
  endif

  if( .not. reduced_io ) then
    write(stdout,*) ' dumping local triples'
    do ibnd=1,nbnd_proc
      !
      CALL davcio( tri_3(:,:,ibnd), nwordtri, iuntri, band_proc(ibnd), + 1 )
      !
    enddo
    close(iuntri)
  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( tri_2, psid )
  deallocate( ibnd_indx, jtmp )
  deallocate( psic_all )
  deallocate( tri_3 )
  if( .not.reduced_io ) then
    deallocate( band, band_proc )
  endif


  return

  end subroutine triples
! ---------------------------------------------------------------------- 

