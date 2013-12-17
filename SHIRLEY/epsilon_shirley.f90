  module epsilon_shirley

  use kinds, only : dp
  USE io_global,  ONLY : stdout, ionode, ionode_id

  implicit none

  real(dp),parameter :: real2MB = 8.d0 / 2.d0**20
  real(dp),parameter :: cmplx2MB = real2MB*2.d0
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  ! file info
  integer :: iunwfcr, nwordwfcr

  ! k-points
  integer :: nkpt
  real(dp),allocatable :: xkpt(:,:), wkpt(:)
  logical :: kcartesian

  ! q-points
  integer :: nqpt
  real(dp),allocatable :: xqpt(:,:), wqpt(:)
  logical :: qcartesian

  ! eigenspace
  complex(dp),allocatable :: eigvec_k(:,:)
  real(dp),allocatable :: eigval_k(:)
  complex(dp),allocatable :: eigvec_kpq(:,:)
  real(dp),allocatable :: eigval_kpq(:)
  real(dp),allocatable :: f_ek(:), f_ekpq(:)
  complex(dp),allocatable :: evc_k(:,:), evc_kpq(:,:)

  ! R space
  integer :: n1, n2, n3, nx
  logical,allocatable :: nonzero(:)

  ! G space
  integer :: nqpgx
  integer,allocatable :: nqpg(:)
  integer :: ngx
  real(dp), allocatable :: g2kinq(:)
  integer,allocatable :: igkq(:)
  integer,allocatable :: igk_l2gq(:,:)

  ! chiq
  integer :: nchiq
  integer,allocatable :: indx_chiq(:,:)
  complex(dp),allocatable :: chiq(:)

  ! tmp space
  complex(dp),allocatable :: chi_r(:), x_g(:), x_gloc(:)
  complex(dp),allocatable :: chi_gloc(:), chi_g(:)

  contains

! ---------------------------------------------------------------------- 
  subroutine basis_to_rspace( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode
  USE cell_base
  USE gvect  
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use scf, only : rho
  use mp, only : mp_bcast, mp_barrier

  use shirley_epsilon_input, only : debug, band_subset, nstride

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)

  integer :: ibnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  complex(dp),allocatable :: psir(:,:)
  complex(dp),allocatable :: psid(:,:,:)
  real(dp),allocatable :: rhor(:), rhoave(:)
  integer :: nrho
  real(dp) :: minrho, maxrho

  integer :: i, ierr
  logical :: exst
  integer,external :: freeunit


  WRITE( stdout, '(/5x,"Calling basis_to_rspace .... ",/)')
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

  write(stdout,*) ' transforming basis functions from Fourier to Real space'
  write(stdout,*) ' nrxxs = ', nrxxs
  write(stdout,*) ' nrxx  = ', nrxx
  write(stdout,*) ' npwx  = ', npwx
  write(stdout,*) ' nr1x  = ', nr1x 
  write(stdout,*) ' nr2x  = ', nr2x 
  write(stdout,*) ' nr3x  = ', nr3x 
  write(stdout,*) ' nr1sx = ', nr1sx
  write(stdout,*) ' nr2sx = ', nr2sx
  write(stdout,*) ' nr3sx = ', nr3sx
  write(stdout,*) ' prod  = ', nr1sx*nr2sx*nr3sx

!  ! estimate reduction in real-space points using density
!  allocate( rhor(nr1*nr2*nr3), rhoave(nrxxs), nonzero(nr1*nr2*nr3) )
!  rhoave= rho(:,1)
!  if( nspin == 2 ) rhoave = rhoave + rho(:,2)
!  call gather( rhoave, rhor )
!  ! search for points with values less than some epsilon
!  nrho=0
!  maxrho = maxval(rhor)
!  minrho = maxrho*1.d-3
!  nonzero = ( abs(rhor(:)) > minrho )
!  do i=1,nr1*nr2*nr3
!    if( nonzero(i) ) then
!      nrho=nrho+1
!    endif
!  enddo
!  write(stdout,*) '         density maximum,norm : ', maxrho, sum(rhor)
!  write(stdout,*) ' cut-off density minimum,norm : ', minrho, sum(rhor,nonzero)
!  write(stdout,*) ' number of    total real-space grid points: ', nr1*nr2*nr3
!  write(stdout,*) ' number of non-zero real-space grid points: ', nrho
!  deallocate( rhor, rhoave )
      
  
 
  ! try to allocate space for all real wave functions
  ! but with a smaller mesh size - twice as small
  allocate( psid(nr1,nr2,nr3) )
  n1 = size(psid(1:nr1:nstride(1),:,:),1)
  n2 = size(psid(:,1:nr2:nstride(2),:),2)
  n3 = size(psid(:,:,1:nr3:nstride(3)),3)
  deallocate( psid )

  nx=n1*n2*n3
  write(stdout,*) ' new reduced size: ', nx
  write(stdout,*) n1, n2, n3
  write(stdout,*) ' attempt to allocate ', dble(nx)*dble(n1*n2*n3)*16.d0/2.d0**20.d0, ' MB'
  allocate( psir(nbnd,nx), stat=ierr )
  if( ierr/=0 ) then
    write(stdout,*) ' unable to allocate space for all states'
    call errore('basis_to_rspace','error allocating',1)
  endif

  nwordwfcr = 2 * nbnd
!  if( ionode ) then
    iunwfcr = freeunit()
    call diropn( iunwfcr, 'wfcr', nwordwfcr, exst )

    if( exst ) then
      write(stdout,*) ' Real space file exists '
      write(stdout,*) ' Testing to see if correct size'
      CALL davcio( psir(:,1), nwordwfcr, iunwfcr, 1, - 1 )
      write(stdout,*) ' first basis function ok'
      CALL davcio( psir(:,nx), nwordwfcr, iunwfcr, nx, - 1 )
      write(stdout,*) ' last  basis function ok'
      goto 999
    endif
!  endif

  call mp_barrier

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

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )


  ! ======================================================================
  ! transform
  ! ======================================================================
  write(stdout,*) ' transforming to Real space'
  allocate( psid(nr1,nr2,nr3) )
  do ibnd=1,nbnd
    !
    write(stdout,*) ' band ', ibnd
    CALL start_clock( 'firstfft' )
    !
    psic(1:nrxxs) = ( 0.D0, 0.D0 )
    psic(nls(igk(1:npw))) = conjg( evc(1:npw,ibnd_indx(ibnd)) )
    CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
    !
    CALL stop_clock( 'firstfft' )
    !
    ! note that I need a complex version of this
    call cgather_sym( psic, psid )
    ! now reduce the resolution
    psir(ibnd,:) = reshape( &
                     psid( 1:nr1:nstride(1), &
                           1:nr2:nstride(2), &
                           1:nr3:nstride(3) ) , shape=(/ nx /) )
  enddo
  deallocate( psid )

!  if( ionode ) then
    write(stdout,*) ' writing to file'
    do i=1,nx
      CALL davcio( psir(:,i), nwordwfcr, iunwfcr, i, + 1 )
    enddo
!  endif

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  999 continue
  deallocate( ibnd_indx )
  if( allocated(evc) ) deallocate(evc)

  return
  end subroutine basis_to_rspace
! ---------------------------------------------------------------------- 


! ---------------------------------------------------------------------- 
  subroutine make_chir
! ---------------------------------------------------------------------- 

  use wvfct, only : nbnd
  use mp, only : mp_barrier
  use mp_global, only : nproc, mpime

  use shirley_epsilon_input, only : block_size, debug
  use diag_shirley

  real(dp),parameter :: eps=1.d-6

  complex(dp),allocatable :: bi(:,:), bj(:,:)
  complex(dp),allocatable :: unk(:,:), unkpq(:,:), btmp(:,:)
  complex(dp),allocatable :: vk(:), vkpq(:), fur(:), ful(:), chiq(:,:)
  complex(dp),allocatable :: fr(:,:), er(:,:)
  complex(dp),allocatable :: fl(:,:), el(:,:)

  integer :: ierr
  integer :: nr, nri, nrj, nrij
  integer :: i, j, ir, jr
  integer :: iq, ik
  real(dp) :: kpq(3)
  integer :: nocc, nemp, nemp0
  integer :: nocc_k, nemp_k
  integer :: nocc_kpq, nemp_kpq
  integer :: ibnd, jbnd

  complex(dp),external :: ZDOTC
  real,external :: etime
  real :: elapsed(2), total
  logical,allocatable :: block_todo(:,:)
  integer :: nblock_g, nblock, iblock
  complex(dp) :: cr, cl


  allocate( eigvec_k(nbasis,nbasis), eigval_k(nbasis), f_ek(nbasis), &
            eigvec_kpq(nbasis,nbasis), eigval_kpq(nbasis), f_ekpq(nbasis), &
            stat=ierr )
  if( ierr/= 0 ) call errore('make_chir','unable to allocate eigenspace',abs(ierr))

  ! block size that can be stored
  write(stdout,*) ' block_size = ', block_size
  nr=2*block_size
  allocate( bi(nbnd,nr), bj(nbnd,nr), & 
            unk(nbasis,nr), unkpq(nbasis,nr), &
            btmp(nbasis,nr), &
            vk(nbasis), vkpq(nbasis), & 
            fur(nbasis), ful(nbasis), &
            chiq(block_size,block_size), & 
            fr(nbasis,nbasis), er(nbasis,nbasis), &
            fl(nbasis,nbasis), el(nbasis,nbasis), &
            block_todo(nx,nx), &
            stat=ierr )
  if( ierr/=0 ) call errore('make_chir','unable to allocate workspace',abs(ierr))

  call mp_barrier
  block_todo = .false.
  nblock_g=0
  nblock=0
  do j=1,nx,block_size
   do i=1,j,block_size
     if( mod(nblock_g,nproc) == mpime ) then
       block_todo(i,j) = .true.
       nblock=nblock+1
     endif
     nblock_g=nblock_g+1
   enddo
  enddo
  write(stdout,*) ' number of r,rp blocks: ', nblock_g
  write(stdout,*) ' number on this proc:   ', nblock
  write(mpime+500,*) ' number of r,rp blocks: ', nblock_g
  write(mpime+500,*) ' number on this proc:   ', nblock
  

  iblock=0
  do j=1,nx,block_size
    nrj = min( block_size, nx-j )
    do ir=1,nrj
      CALL davcio( bj(:,ir), nwordwfcr, iunwfcr, j+ir-1, - 1 )
    enddo
    btmp(:,1:nrj) = bj(:,1:nrj)

    do i=1,j,block_size

      if( .not. block_todo(i,j) ) cycle

      iblock=iblock+1

      nri=0
      if( i/=j ) then
        nri = min( block_size, nx-i )
        do ir=1,nri
          CALL davcio( bi(:,ir), nwordwfcr, iunwfcr, i+ir-1, - 1 )
        enddo
        btmp(:,nrj+1:nrj+nri) = bi(:,1:nri)
      endif
      
      nrij = nri+nrj

      total = etime(elapsed)
!      write(stdout,*) ' (b,bp) = ', i, j, nri, nrj, nrij
      write(stdout,'(a,i,a,i,a,f)') ' block ', iblock, ' of ', nblock, ' t = ', total

      qpt_loop: do iq=1,nqpt

        if( debug ) write(stdout,*) ' q-point ', iq

        chiq = 0.d0

        kpt_loop: do ik=1,nkpt

          if( debug ) write(stdout,*) ' k-point ', ik

          ! find eigenstates and values at k and k+q
          call diag_hamq( xkpt(:,ik), eigval_k, eigvec_k, kcartesian )
          kpq = xkpt(:,ik) + xqpt(:,iq)
          call diag_hamq( kpq, eigval_kpq, eigvec_kpq, kcartesian )
      
          ! occupancies
          f_ek = fermifunc( eigval_k )
          f_ekpq = fermifunc( eigval_kpq )

          do nocc_k=1,nbasis
            if( abs(f_ek(nocc_k)) < eps ) exit
          enddo
          nemp_k = nocc_k
          nocc_k = nocc_k - 1
    
          do nocc_kpq=1,nbasis
            if( abs(f_ek(nocc_kpq)) < eps ) exit
          enddo
          nemp_kpq = nocc_kpq
          nocc_kpq = nocc_kpq - 1
    
          forall(ibnd=1:nocc_k, jbnd=nemp_kpq:nbasis) &
            fr(ibnd,jbnd-nemp_kpq+1) = cmplx(f_ekpq(jbnd) - f_ek(ibnd))
          forall(ibnd=nemp_k:nbasis, jbnd=1:nocc_kpq) &
            fl(ibnd-nemp_k+1,jbnd) = cmplx(f_ekpq(jbnd) - f_ek(ibnd))

          forall(ibnd=1:nocc_k, jbnd=nemp_kpq:nbasis) &
            er(ibnd,jbnd-nemp_kpq+1) = cmplx(eigval_kpq(jbnd) - eigval_k(ibnd))
          forall(ibnd=nemp_k:nbasis, jbnd=1:nocc_kpq) &
            el(ibnd-nemp_k+1,jbnd) = cmplx(eigval_kpq(jbnd) - eigval_k(ibnd))

          ! divide by the energy
          where( abs(fr) >= eps ) fr = fr / er
          where( abs(fr) <  eps ) fr = zero
          where( abs(fl) >= eps ) fl = fl / el
          where( abs(fl) <  eps ) fl = zero
          
!          write(stdout,*) ' occupied k   = ', nocc_k
!          write(stdout,*) '    empty k   = ', nemp_k
!          write(stdout,*) ' occupied k+q = ', nocc_kpq
!          write(stdout,*) '    empty k+q = ', nemp_kpq

!          write(400+mpime,'(4e)') (eigval_k(ibnd), f_ek(ibnd), &
!                                  eigval_kpq(ibnd), f_ekpq(ibnd), &
!                                  ibnd=1,nbasis)
!           write(300+mpime,'(2i,2e)') ((ibnd, jbnd, f(ibnd,jbnd), ibnd=1,nbasis), jbnd=1,nbasis)


          ! expand unk(r)
          call ZGEMM( 'C', 'N', nbasis, nrij, nbasis, one, &
                      eigvec_k, nbasis, btmp, nbasis, zero, unk, nbasis )

          ! expand unk+q(r)
          call ZGEMM( 'C', 'N', nbasis, nrij, nbasis, one, &
                      eigvec_kpq, nbasis, btmp, nbasis, zero, unkpq, nbasis )

          if( i==j ) then
            do jr=1,nrj
              do ir=1,jr
                vk(:)   = unk(:,ir)   * conjg( unk(:,jr) )
                vkpq(:) = unkpq(:,ir) * conjg( unkpq(:,jr) )

!                cr = dot_product( vk(1:nocc_k), &
!                       matmul( fr, vkpq(nemp_kpq:nbasis) ) )
                call zgemv( 'N', nocc_k, nbasis-nemp_kpq+1, one, &
                            fr, nbasis, vkpq(nemp_kpq), 1, zero, fur, 1 )
                cr = zdotc( nocc_k, vk, 1, fur, 1 )

!                cl = dot_product( vk(nemp_k:nbasis), &
!                       matmul( fl, vkpq(1:nocc_kpq) ) )
                call zgemv( 'N', nbasis-nemp_k+1, nocc_kpq, one, &
                            fl, nbasis, vkpq, 1, zero, ful, 1 )
                cr = zdotc( nbasis-nemp_k+1, vk(nemp_k), 1, ful, 1 )

                chiq(ir,jr) = chiq(ir,jr) + wkpt(ik) * (cr+cl)
              enddo
            enddo
          else
            do jr=1,nrj
              do ir=1,nri
                vk(:)   = unk(:,nrj+ir)   * conjg( unk(:,jr) )
                vkpq(:) = unkpq(:,nrj+ir) * conjg( unkpq(:,jr) )

!                cr = dot_product( vk(1:nocc_k), &
!                       matmul( fr, vkpq(nemp_kpq:nbasis) ) )
                call zgemv( 'N', nocc_k, nbasis-nemp_kpq+1, one, &
                            fr, nbasis, vkpq(nemp_kpq), 1, zero, fur, 1 )
                cr = zdotc( nocc_k, vk, 1, fur, 1 )

!                cl = dot_product( vk(nemp_k:nbasis), &
!                       matmul( fl, vkpq(1:nocc_kpq) ) )
                call zgemv( 'N', nbasis-nemp_k+1, nocc_kpq, one, &
                            fl, nbasis, vkpq, 1, zero, ful, 1 )
                cr = zdotc( nbasis-nemp_k+1, vk(nemp_k), 1, ful, 1 )

                chiq(ir,jr) = chiq(ir,jr) + wkpt(ik) * (cr+cl)
              enddo
            enddo
          endif

        enddo kpt_loop

        ! report
        if( i==j ) then
          do jr=1,nrj
            do ir=1,jr
              write(600+mpime,'(i4,a,2i,2e)') iq, ' chiq ', i+ir-1, j+jr-1, chiq(ir,jr)
            enddo
          enddo
        else
          do jr=1,nrj
            do ir=1,nri
              write(600+mpime,'(i4,a,2i,2e)') iq, ' chiq ', i+ir-1, j+jr-1, chiq(ir,jr)
            enddo
          enddo
        endif
        
      enddo qpt_loop
      
    enddo
  enddo
  deallocate( bi, bj )

  deallocate( bi, bj,  unk, unkpq, btmp, vk, vkpq, & 
              fur, ful, chiq,  fr, er, fl, el )
  deallocate( eigvec_k, eigval_k, f_ek, &
              eigvec_kpq, eigval_kpq, f_ekpq )

  end subroutine make_chir
! ---------------------------------------------------------------------- 


! ---------------------------------------------------------------------- 
  subroutine make_kpoints
! ---------------------------------------------------------------------- 
  use shirley_epsilon_input, only : nkgrid, ikgrid
  use cell_base, only : bg

  integer :: ik, i1, i2, i3

  nkpt = product( nkgrid )
  allocate( xkpt(3,nkpt), wkpt(nkpt) )
  ik=0
  do i1=1,nkgrid(1) 
  do i2=1,nkgrid(2) 
  do i3=1,nkgrid(3) 
    ik=ik+1
    xkpt(1,ik) = dble(i1-1)/dble(nkgrid(1))+ikgrid(1)/2.d0/dble(nkgrid(1))
    xkpt(2,ik) = dble(i2-1)/dble(nkgrid(2))+ikgrid(2)/2.d0/dble(nkgrid(2))
    xkpt(3,ik) = dble(i3-1)/dble(nkgrid(3))+ikgrid(3)/2.d0/dble(nkgrid(3))
  enddo
  enddo
  enddo
  kcartesian = .true.
  xkpt = matmul( bg, xkpt )
  wkpt = 2.d0 / dble(nkpt)

  write(stdout,'(/,a,/)') ' k-points: '
  do ik=1,nkpt
    write(stdout,'(i,3f)') ik, xkpt(1:3,ik)
  enddo
  write(stdout,'(/)')

  return
  end subroutine make_kpoints

! ---------------------------------------------------------------------- 
  subroutine make_qpoints
! ---------------------------------------------------------------------- 
  use shirley_epsilon_input, only : nqgrid, iqgrid
  use cell_base, only : bg

  integer :: iq, i1, i2, i3

  nqpt = product( nqgrid )
  allocate( xqpt(3,nqpt), wqpt(nqpt) )
  iq=0
  do i1=1,nqgrid(1) 
  do i2=1,nqgrid(2) 
  do i3=1,nqgrid(3) 
    iq=iq+1
    xqpt(1,iq) = dble(i1-1)/dble(nqgrid(1))+iqgrid(1)/2.d0/dble(nqgrid(1))
    xqpt(2,iq) = dble(i2-1)/dble(nqgrid(2))+iqgrid(2)/2.d0/dble(nqgrid(2))
    xqpt(3,iq) = dble(i3-1)/dble(nqgrid(3))+iqgrid(3)/2.d0/dble(nqgrid(3))
  enddo
  enddo
  enddo
  qcartesian = .true.
  xqpt = matmul( bg, xqpt )
  wqpt = 2.d0 / dble(nqpt)

  write(stdout,'(/,a,/)') ' q-points: '
  do iq=1,nqpt
    write(stdout,'(i,3f)') iq, xqpt(1:3,iq)
  enddo
  write(stdout,'(/)')

  return
  end subroutine make_qpoints


! ---------------------------------------------------------------------- 
  subroutine make_epsilon_lowmem
! ---------------------------------------------------------------------- 

  use cell_base
  use gvect
  use gvecs
  use smooth_grid_dimensions
  use wvfct
  use wavefunctions_module
  use io_files
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool, &
                        mpime
  USE mp_wave,   ONLY : mergewf
  use diag_shirley
  use shirley_epsilon_input, only : ecuteps, debug

  real(dp),parameter :: eps=1.d-10

  real(dp) :: kpq(3)
  integer :: ik, iq, ierr
  integer :: nocc, nemp, nemp0
  integer :: ibasis, jbasis
  integer :: i, j, ig, jg, n, ij
  complex(dp) :: fac
  complex(dp),allocatable :: psid(:)


  write(500+me_pool,*) ' make_epsilon '
  call gspace_dim

  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  write(stdout,*) ' attempt to allocate space: '
  write(stdout,*) dble(nbasis*nbasis)*cmplx2MB + dble(nbasis*4)*real2MB
  write(stdout,*) dble(nrxxs)*cmplx2MB

  write(500+me_pool,*) nqpgx, ngx

  allocate( eigvec_k(nbasis,nbasis), eigval_k(nbasis), &
            eigvec_kpq(nbasis,nbasis), eigval_kpq(nbasis), &
            f_ek(nbasis), f_ekpq(nbasis), &
            psid(nrxxs), &
            stat=ierr )
  if( ierr/= 0 ) then
    call errore('make_epsilon','error allocating eigenspace',abs(ierr))
  endif
  !          psid(nrxxs), &
  !          chi_gloc(nqpgx), chi_g(ngx), &
  !          stat=ierr )

  qpt_loop: do iq=1,nqpt

    CALL gk_sort( xqpt(1,iq), ngm, g, ecuteps / tpiba2, nqpg(iq), igkq, g2kinq )
    where( g2kinq > 0.d0 ) g2kinq = 1.d0 / ( sqrt(g2kinq) * tpiba )

    write(500+me_pool,*) ' iq = ', iq

    chiq = zero

    kpt_loop: do ik=1,nkpt

      write(500+me_pool,*) ' ik = ', ik

      ! find eigenstates and values at k and k+q
      call diag_hamq( xkpt(:,ik), eigval_k, eigvec_k, kcartesian )
      kpq = xkpt(:,ik) + xqpt(:,iq)
      call diag_hamq( kpq, eigval_kpq, eigvec_kpq, kcartesian )
      
      ! occupancies
      f_ek=0.d0
      do nocc=1,nbasis
        f_ek(nocc) = fermifunc( eigval_k(nocc) )
        if( f_ek(nocc) < eps ) exit
      enddo
      nocc = nocc - 1

      f_ekpq=0.d0
      do nemp0=nbasis,1,-1
        f_ekpq(nemp0) = fermifunc( eigval_kpq(nemp0) )
        if( f_ekpq(nemp0) >= eps ) exit
      enddo
      nemp = nbasis - nemp0
        
      write(stdout,*) ' occupied = ', nocc
      write(stdout,*) '    empty = ', nemp

      write(500+me_pool,*) 1, nocc, nemp0, nbasis

  write(stdout,*) ' attempt to allocate space: '
  write(stdout,*) (dble(npwx)*dble(nocc+nemp)+dble(nqpgx)+dble(ngx))*cmplx2MB

      allocate( evc_k(npwx,nocc), evc_kpq(npwx,nemp), &
                x_gloc(nqpgx), x_g(ngx), &
                stat=ierr )
      if( ierr/= 0 ) &
        call errore('make_epsilon','error allocating evc',abs(ierr))

      call ZGEMM( 'N', 'N', npw, nocc, nbasis, one, &
                  evc, npwx, eigvec_k, nbasis, zero, evc_k, npwx )
      
      call ZGEMM( 'N', 'N', npw, nemp, nbasis, one, &
                  evc, npwx, eigvec_kpq(:,nemp0+1:nbasis), nbasis, zero, evc_kpq, npwx )
      
      ! ... to real space
      !
      CALL start_clock( 'firstfft' )
      do i=1,nocc
        if( debug ) write(stdout,*) ' occupied ', i

        psic(nls(igk(1:npw))) = evc_k(1:npw,i)
        CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )
        psic = conjg( psic )
        do j=1,nemp
          fac = sqrt( (f_ekpq(nemp0+j)-f_ek(i))/(eigval_kpq(nemp0+j)-eigval_k(i)) )
          psid(nls(igk(1:npw))) = evc_kpq(1:npw,j)
          CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )
          psid = psid * psic
          CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )

          x_gloc(1:nqpg(iq)) = &
            fac * psid(nls(igkq(1:nqpg(iq)))) / g2kinq(1:nqpg(iq))

!          call mergewf(x_gloc, x_g, nqpg(iq), igk_l2gq(:,iq), &
!                       me_pool, nproc_pool, root_pool, intra_pool_comm)

!          do ij=1,nchiq
!            chiq(ij) = chiq(ij) + &
!              x_g(indx_chiq(1,ij)) * conjg( x_g(indx_chiq(2,ij)) )
!          enddo

        enddo
      enddo
      deallocate( x_g, x_gloc, evc_k, evc_kpq )
      !
    enddo kpt_loop

    write(600+mpime,*) iq
    write(600+mpime,*) chiq

  enddo qpt_loop

  return
  end subroutine make_epsilon_lowmem

! ---------------------------------------------------------------------- 
  subroutine make_epsilon_highmem
! ---------------------------------------------------------------------- 

  use cell_base
  use gvect
  use gvecs
  use smooth_grid_dimensions
  use klist, only : xk
  use wvfct
  use wavefunctions_module
  use io_files
  use mp, only : mp_sum, mp_bcast, mp_barrier, mp_max
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool, &
                        mpime
  USE mp_wave,   ONLY : mergewf, splitwf
  use diag_shirley
  use shirley_epsilon_input, only : ecuteps, debug

  real(dp),parameter :: eps=1.d-10

  integer :: npweps, igwxeps
  integer,allocatable :: igkeps(:), igk_l2geps(:)
  real(dp),allocatable :: g2kineps(:)
  real(dp) :: kpq(3)
  integer :: ik, iq, ierr
  integer :: nocc, nemp, nemp0
  integer :: ibasis, jbasis, kbasis, lbasis
  integer :: npair, ipair, jpair
  integer,allocatable :: indx(:,:)
  integer :: i, j, ig, jg, n, ij, ikq
  real(dp) :: fac
  logical,allocatable :: pw_written(:)
  complex(dp),allocatable :: pwelti(:), pweltj(:), pwelt_g(:)
  complex(dp),allocatable :: phi_ij(:)
  logical,allocatable :: kq_todo(:,:)
  integer :: nwordpw, iunpw
  logical :: exst
  integer,external :: freeunit
  complex(dp),allocatable :: psid(:)


  write(500+me_pool,*) ' make_epsilon '
!  call gspace_dim
  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
  allocate( igkeps(npw), g2kineps(npw) )
  CALL gk_sort( xk(1,1), ngm, g, ecuteps / tpiba2, npweps, igkeps, g2kineps )
  deallocate( igkeps, g2kineps )
  allocate( igkeps(npweps), g2kineps(npweps), igk_l2geps(npweps) )
  CALL gk_sort( xk(1,1), ngm, g, ecuteps / tpiba2, npweps, igkeps, g2kineps )
  CALL gk_l2gmap (ngm, ig_l2g(1), npweps, igkeps, igk_l2geps(1))
  igwxeps = maxval( igk_l2geps )
  call mp_max( igwxeps )
  

  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! divide (k,q) among processes
  allocate( kq_todo(nkpt,nqpt) )
  kq_todo = .false.
  ikq=0
  do iq=1,nqpt
    do ik=1,nkpt
      if( mod(ikq,nproc_pool) == me_pool ) kq_todo(ik,iq)=.true.
      ikq=ikq+1
    enddo
  enddo

  npair=nbnd*nbnd
  allocate( indx(2,npair) )
  npair=0
  do jbasis=1,nbnd
    do ibasis=1,nbnd
      npair=npair+1
      indx(:,npair) = (/ ibasis, jbasis /)
    enddo
  enddo

  allocate( pw_written(npair) )
  pw_written = .false.

  if( ionode ) then
    nwordpw = 2 * igwxeps
    iunpw = freeunit()
    call diropn( iunpw, 'pw', nwordpw, exst )
  endif
  call mp_bcast( exst, ionode_id )
  if( exst ) then
    pw_written = .true.
    write(stdout,*) ' Found plane-wave matrix file'
  endif


  ! big allocations
  if( allocated(psid) ) deallocate(psid)
  allocate( psid(nrxxs), pwelti(npweps), pweltj(npweps), pwelt_g(igwxeps), &
            eigvec_k(nbasis,nbasis), eigval_k(nbasis), &
            eigvec_kpq(nbasis,nbasis), eigval_kpq(nbasis), &
            f_ek(nbasis), f_ekpq(nbasis), &
            stat=ierr )
  if( ierr/=0 ) &
    call errore('make_epsilson','unable to allocate space',abs(ierr))


  ! new code
  ! outer loop over plane-wave matrix elements
  outer_pair: do ipair=1,npair

    ibasis=indx(1,ipair)
    jbasis=indx(2,ipair)
    write(stdout,*) ' outer_pair ', ibasis,jbasis

    if( .not. pw_written(ipair) ) then
      psic(nls(igk(1:npw))) = evc(1:npw,ibasis)
      CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )
      psic = conjg( psic )
    
      psid(nls(igk(1:npw))) = evc(1:npw,jbasis)
      CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )

      ! construct plane-wave matrix element
      psid = psic * psid
      CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )

      pwelti(1:npweps) = psid(nls(igkeps(1:npweps)))
      call mergewf(pwelti, pwelt_g, npweps, igk_l2geps, &
                   me_pool, nproc_pool, root_pool, intra_pool_comm)

      if( ionode ) call davcio( pwelt_g, nwordpw, iunpw, ipair, + 1 )
      call mp_barrier
      pw_written(ipair) = .true.
    else
      if( ionode ) call davcio( pwelt_g, nwordpw, iunpw, ipair, - 1 )
      call mp_barrier
      
      call splitwf(pwelti, pwelt_g, npweps, igk_l2geps, &
                   me_pool, nproc_pool, root_pool, intra_pool_comm)
    endif


    inner_pair: do jpair=1,ipair

      kbasis=indx(1,jpair)
      lbasis=indx(2,jpair)
      write(stdout,*) ' inner_pair ', kbasis,lbasis

      ! I think this should never happen
      if( .not. pw_written(jpair) ) then
        write(stdout,*) ' Warning: generating a pw-matrix element in the inner loop'

        psic(nls(igk(1:npw))) = evc(1:npw,kbasis)
        CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )
        psic = conjg( psic )

        psid(nls(igk(1:npw))) = evc(1:npw,lbasis)
        CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, +2 )

        ! construct plane-wave matrix element
        psid = psic * psid
        CALL cft3s( psid, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )

        pweltj(1:npweps) = psid(nls(igkeps(1:npweps)))
        call mergewf(pweltj, pwelt_g, npweps, igk_l2geps, &
                     me_pool, nproc_pool, root_pool, intra_pool_comm)

        if( ionode ) call davcio( pwelt_g, nwordpw, iunpw, jpair, + 1 )
        call mp_barrier
        pw_written(jpair) = .true.
      else
        if( ionode ) call davcio( pwelt_g, nwordpw, iunpw, jpair, - 1 )
        call mp_barrier

        call splitwf(pweltj, pwelt_g, npweps, igk_l2geps, &
                     me_pool, nproc_pool, root_pool, intra_pool_comm)
      endif

      ! now begin summing elements of chi

      allocate( phi_ij(nqpt) )
      phi_ij = zero

      qpt_loop: do iq=1,nqpt

!        write(stdout,*) ' q-point ', iq

!        CALL gk_sort( xqpt(1,iq), ngm, g, ecuteps / tpiba2, &
!                      nqpg(iq), igkq, g2kinq )
!        where( g2kinq > 0.d0 ) g2kinq = 1.d0 / ( sqrt(g2kinq) * tpiba )


        ! this part should be split over processors
        kpt_loop: do ik=1,nkpt

!          write(stdout,*) '   k-point ', ik

          if( .not. kq_todo(ik,iq) ) cycle kpt_loop

!          write(500+mpime,*) ' processor ', mpime, ' working on (ik,iq)=', ik,iq

          ! find eigenstates and values at k and k+q
          call diag_hamq( xkpt(:,ik), eigval_k, eigvec_k, kcartesian )
          kpq = xkpt(:,ik) + xqpt(:,iq)
          call diag_hamq( kpq, eigval_kpq, eigvec_kpq, kcartesian )
      
          ! occupancies
          f_ek=0.d0
          do nocc=1,nbasis
            f_ek(nocc) = fermifunc( eigval_k(nocc) )
            if( f_ek(nocc) < eps ) exit
          enddo
          nocc = nocc - 1
    
          f_ekpq=0.d0
          do nemp0=nbasis,1,-1
            f_ekpq(nemp0) = fermifunc( eigval_kpq(nemp0) )
            if( f_ekpq(nemp0) >= eps ) exit
          enddo
          nemp = nbasis - nemp0
        
!          write(stdout,*) ' occupied = ', nocc
!          write(stdout,*) '    empty = ', nemp

          do i=1,nocc
            do j=nemp0,nbasis
              
              fac = f_ekpq(j)-f_ek(i)
              if( fac < eps ) then
                fac = 0.d0
              else
                fac = (f_ekpq(j)-f_ek(i))/(eigval_kpq(j)-eigval_k(i))
              endif

              phi_ij(iq) = phi_ij(iq) &
                + fac * conjg(eigvec_k(ibasis,i))*eigvec_kpq(jbasis,j) &
                      * conjg(eigvec_kpq(kbasis,j))*eigvec_k(lbasis,i)

            enddo
          enddo

        enddo kpt_loop

      enddo qpt_loop

      call mp_sum( phi_ij, intra_pool_comm )

!      ! combine with plane-wave matrix elements
!      do iq=1,nqpt
!
!        do ig=1,igwxeps
!
!        ! load chi_q into memory
!        call dav
!
!        chi_q
!      enddo



      deallocate( phi_ij )

     enddo inner_pair

  enddo outer_pair

  deallocate( psid, pwelti, pweltj, pwelt_g )
  deallocate( kq_todo, indx, pw_written )
  deallocate( igkeps, g2kineps )

  return
  end subroutine make_epsilon_highmem

  subroutine gspace_dim

  use cell_base
  use gvect
  use wvfct
  use mp,        only : mp_max
  use mp_global, only : mpime, nproc
  use shirley_epsilon_input, only : ecuteps

  integer :: iq, ij, ig, jg, ierr

  if( ecuteps <= 0.d0 ) ecuteps = ecutwfc

  if( allocated(igkq) ) deallocate(igkq)
  allocate( igkq(ngm) )
  if( allocated(g2kinq) ) deallocate(g2kinq)
  allocate( g2kinq(ngm) )

  if( allocated(igk_l2gq) ) deallocate(igk_l2gq)
  allocate( igk_l2gq(ngm,nqpt) )

  allocate( nqpg(nqpt) )
  igk_l2gq=0
  do iq=1,nqpt
    CALL gk_sort( xqpt(1,iq), ngm, g, ecuteps / tpiba2, nqpg(iq), igkq, g2kinq )
    CALL gk_l2gmap (ngm, ig_l2g(1), nqpg(iq), igkq, igk_l2gq(1,iq))
  enddo
  nqpgx = maxval( nqpg )
  ngx = maxval( igk_l2gq )
  call mp_max(ngx)

  write(stdout,*) ' nqpgx = ', nqpgx
  write(stdout,*) ' nqpg = ', nqpg
  write(stdout,*) ' ngx   = ', ngx

  deallocate( igkq, g2kinq )
  allocate( igkq(nqpgx), g2kinq(nqpgx) )

  ij=0
  nchiq=0
  do jg=1,ngx
    do ig=1,jg
      if( mod(ij,nproc)==mpime ) nchiq=nchiq+1
      ij=ij+1
    enddo
  enddo
  write(stdout,*) ' attempt to allocate ', nchiq*cmplx2MB
  allocate( chiq(nchiq), stat=ierr )
  if( ierr/= 0 ) &
    call errore('gspace_dim','unable to allocate chiq',abs(ierr))

  allocate( indx_chiq(2,nchiq) )
  ij=0
  nchiq=0
  do jg=1,ngx
    do ig=1,jg
      if( mod(ij,nproc)==mpime ) then
        nchiq=nchiq+1
        indx_chiq(:,nchiq) = (/ ig, jg /)
      endif
      ij=ij+1
    enddo
  enddo
      
  return
  end subroutine gspace_dim


  elemental function fermifunc( e )
  use shirley_epsilon_input, only : efermi, etemp
  real(dp),intent(in) :: e
  real(dp) :: fermifunc

  if( etemp < 0.d0 ) then
    if( e < efermi ) then
      fermifunc=1.d0
      return
    else
      fermifunc=0.d0
      return
    endif
  else
    if( e < efermi - 6.d0*etemp ) then
      fermifunc = 1.d0
      return
    elseif( e > efermi + 6.d0*etemp ) then
      fermifunc = 0.d0
      return
    else
      fermifunc = 1.d0 / ( exp( (e-efermi)/etemp ) + 1.d0 )
      return
    endif
  endif
  end function

  end module epsilon_shirley
