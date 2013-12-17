  program shirley_epstci

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
  ! use pwmat_module
  ! use scalapack_module

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6

  COMPLEX(DP),PARAMETER :: ZERO=CMPLX(0.d0,0.d0)
  COMPLEX(DP),PARAMETER :: ONE =CMPLX(1.d0,0.d0)

  integer :: iuntci
  integer :: neig_tci
  complex(dp),allocatable :: tci(:,:)
  complex(dp),allocatable :: tci_mat(:,:), tci_tmp(:,:)

  complex(dp),allocatable :: wx(:), xivec(:,:), xjvec(:)
  integer :: ij

  ! scalapack 
  integer,allocatable :: gridmap(:,:)
  integer :: context, iam
  integer :: nprocs, nprow, npcol, myrow, mycol
  integer :: li, lj, prow, pcol
  integer :: nblock, nr_l, info
  integer :: desc_eps(9)
  integer,external :: numroc


  integer :: ierr
  integer :: neig2
  real(dp) :: qvec(3), qvec_cart(3)

  real(dp) :: kpqvec(3)
  complex(dp),allocatable :: eigvec_k(:,:)
  real(dp),allocatable :: eigval_k(:)
  real(dp),allocatable :: focc_k(:), focc_kpq(:)
  real(dp),allocatable :: dfocc(:)
  logical,allocatable :: wmask(:,:)
  integer :: nnz
  integer,allocatable :: nzlist(:,:)
  integer :: ikoccmax, ikpqoccmax
  real(dp) :: df
  integer :: ik, iq, i, j, inz, ip, jp

  integer :: nqpg_lr, nqpg_lc, ldaqpg_l
  complex(dp),allocatable :: epsmat(:,:)
  integer :: ig, jg, ig_l, jg_l, jg_last
  logical :: islocal
  integer :: ipol, apol

  integer :: ipsour
  integer :: iunpwmtmp

  real(dp) :: fpi
  real(dp) :: epsilon_EC, epsilon_LF
  real(dp) :: kT

  real(dp) :: spindeg

  real(dp) :: mem

  complex(dp),external :: zdotc
  interface fermifunc
    pure function fermifunc( e, ef, kT )
    use kinds, only : dp
    real(dp),intent(in) :: e, ef, kT
    real(dp) :: fermifunc
    end function fermifunc
  end interface fermifunc

  integer,external :: freeunit


  call start_clock( 'shirley' )

  fpi = acos(-1.d0)*4.d0

  call shirley_input

  call diag_init

  ! set spin degeneracy
  spindeg = 2.d0

!  write(stdout,*) '     cell volume = ', omega
  write(stdout,*) '          efermi = ', efermi
  write(stdout,*) '   electron temp = ', temperature
  write(stdout,*) ' spin degeneracy = ', spindeg

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

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
!  allocate( eigvec_k(neig,neig), eigval_k(neig), &
!            focc_k(neig), focc_kpq(neig), &
!            xmat(neig,neig), wmask(neig,neig) )
  allocate( eigvec_k(neig,neig), eigval_k(neig), &
            focc_k(neig), focc_kpq(neig) )


  ! initialize scalapack
  !call scalapack_init

  ! make a horizontal striped grid (not block cyclic)
  call blacs_pinfo( iam, nprocs )
  call blacs_get( -1, 0, context )
  allocate( gridmap(nprocs,1) )
  do i=1,nprocs
    gridmap(i,1)=i-1
  enddo
  call blacs_gridmap( context, gridmap, nprocs, nprocs, 1 )
  call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

  ! make the descriptor
  nblock = max( neig / nprow, 1 )
  nr_l = numroc( neig, nblock, myrow, 0, nprow )
  call descinit( desc_eps, neig, neig, nblock, 1, 0, 0, context, &
                 max(1,nr_l), info )

  allocate( epsmat(nr_l,neig), tci(neig2,nr_l), &
            tci_mat(neig,neig) )

  write(200+iam,*) ' nprocs = ', nprocs
  write(200+iam,*) ' nprow = ', nprow
  write(200+iam,*) ' npcol = ', npcol
  write(200+iam,*) ' myrow = ', myrow
  write(200+iam,*) ' mycol = ', mycol

  do j=1,neig
  write(200+iam,*) 
  write(200+iam,'(7a6)') 'i', 'j', 'li', 'lj', 'iam', 'myrow', 'mycol'
    do i=1,neig
      call infog2l( i,j,desc_eps,nprow,npcol,myrow,mycol,li,lj,prow,pcol )
      if( myrow==prow .and. mycol==pcol ) &
        write(200+iam,'(7i6)') i,j,li,lj,iam,myrow,mycol
    enddo
  enddo

! read three-center integrals from file
  if( ionode ) then
    iuntci=freeunit()
    open(iuntci,file=trim(tcifile),form='unformatted')
    read(iuntci) neig_tci
    if( neig_tci /= neig ) &
      call errore('shirley_epstci', &
        'number of basis functions in tci file disagrees',1)
  endif
  allocate( xjvec(neig) )
  do j=1,neig2
    if( ionode ) read(iuntci) xjvec
    call mp_bcast( xjvec, ionode_id )
    do i=1,neig
      call infog2l( i,1,desc_eps,nprow,npcol,myrow,mycol,li,lj,prow,pcol )
      if( myrow==prow .and. mycol==pcol ) &
        tci(j,li) = xjvec(i)
    enddo  
  enddo
  if( ionode ) close(iuntci)
  deallocate( xjvec )
  if( ionode ) then
    write(stdout,*) ' basis functions = ', neig
    write(stdout,*) ' number of tci = ', nr_l
  endif
  

  epsmat = ZERO
  ! should make epsilon unity here
  do i=1,neig
    call infog2l( i,i,desc_eps,nprow,npcol,myrow,mycol,li,lj,prow,pcol )
    if( myrow==prow .and. mycol==pcol ) &
      epsmat(li,lj) = ONE
  enddo


! ======================================================================
! loop over qvec to compute epsilon q
  do iq=1,qpt%list%nk
! ======================================================================
    qvec = qpt%list%kvec(:,iq)

    ! convert qvec to Cartesian (units of tpiba)
    qvec_cart = matmul( trnlp2kin, qvec )
  
    ! loop over k
    do ik=1,kpt%list%nk
      write(stdout,*) ' construct contribution from k-point ', ik

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
      eigval_k = eigval
      write(stdout,*) ' done with k'

      ! build hamiltonian for k+q ...
      write(stdout,*) ' build ham k+q'
      call diag_build_hamk( kpqvec, kpt%param%cartesian )
      ! ... and diagonalize
      write(stdout,*) ' diag ham k+q'
      call diag_ham

      ! fermi occupation factors
      ! efermi read from input in epsilon_module
      do i=1,neig
        focc_k(i) = spindeg * fermifunc( eigval_k(i), efermi, kT )
        focc_kpq(i) = spindeg * fermifunc( eigval(i), efermi, kT )
      enddo

      !write(stdout,'(i,2f)') (i, focc_k(i), focc_kpq(i), i=1,neig)

      do i=neig,1,-1
        if( abs(focc_k(i)) > 1.d-12 ) exit
      enddo
      ikoccmax=i
      do j=neig,1,-1
        if( abs(focc_kpq(j)) > 1.d-12 ) exit
      enddo
      ikpqoccmax=j

      write(stdout,*) ' index of highest occupied band for   k = ', ikoccmax
      write(stdout,*) ' index of highest occupied band for k+q = ', ikpqoccmax

      ! redefine nnz based on two matrix blocks
      nnz = ikpqoccmax*(neig-ikoccmax)+ikoccmax*(neig-ikpqoccmax)
      write(stdout,*) ' non-zero occupation differences = ', nnz, &
                      ' out of possible ', neig*neig
      allocate( wx(nnz) )

      nnz=0
      do j=1,ikpqoccmax
        do i=ikoccmax+1,neig
          nnz=nnz+1
          wx(nnz) = cmplx( (focc_kpq(j) - focc_k(i))/(eigval(j)-eigval_k(i)) )
        enddo
      enddo
      do j=ikpqoccmax+1,neig
        do i=1,ikoccmax
          nnz=nnz+1
          wx(nnz) = cmplx( (focc_kpq(j) - focc_k(i))/(eigval(j)-eigval_k(i)) )
        enddo
      enddo

      ! construct X for each local component
      write(stdout,*) 'Construct X'
      allocate( xivec(nnz,nr_l), xjvec(nnz) )
      do i=1,nr_l
        ! write(stdout,*) 'X', i
        tci_mat = reshape( tci(:,i), (/ neig, neig /) )
        ! block 1 lower left
        allocate( tci_tmp(neig,ikpqoccmax) )
        call zgemm('N','N',neig,ikpqoccmax,neig,ONE,tci_mat,neig, &
                   eigvec,neig,ZERO,tci_tmp,neig)
        call zgemm('C','N',neig-ikoccmax,ikpqoccmax,neig,ONE,eigvec_k(1,ikoccmax+1),neig, &
                   tci_tmp,neig,ZERO,xivec(1:(neig-ikoccmax)*ikpqoccmax,i),neig-ikoccmax)

        ! block 2 upper right
        deallocate( tci_tmp )
        allocate( tci_tmp(neig,neig-ikpqoccmax) )
        call zgemm('N','N',neig,neig-ikpqoccmax,neig,ONE,tci_mat,neig, &
                   eigvec(1,ikpqoccmax+1),neig,ZERO,tci_tmp,neig)
        call zgemm('C','N',ikoccmax,neig-ikpqoccmax,neig,ONE,eigvec_k,neig, &
                   tci_tmp,neig,ZERO,xivec((neig-ikoccmax)*ikpqoccmax+1:nnz,i),ikoccmax)
        deallocate( tci_tmp )
      enddo

      ! now loop over indices of polarizability
      do j=1,neig
        call infog2l( j,1,desc_eps,nprow,npcol,myrow,mycol,lj,li,prow,pcol )
        if( myrow==prow .and. mycol==pcol ) xjvec = xivec(:,lj)
        !write(stdout,*) j, prow, pcol, gridmap(prow+1,pcol+1)
        call mp_bcast( xjvec, gridmap(prow+1,pcol+1) )
        xjvec(:) = kpt%list%wk(ik) * wx(:) * xjvec(:)
        do i=1,neig
          call infog2l( i,j,desc_eps,nprow,npcol,myrow,mycol,li,lj,prow,pcol )
          if( myrow==prow .and. mycol==pcol ) then
            epsmat(li,lj) = epsmat(li,lj) + &
              zdotc( nnz, xivec(1,li), 1, xjvec, 1 )
          endif
        enddo
      enddo

      write(stdout,*) epsmat(1,1)
        
      deallocate( wx, xivec, xjvec )
      
! ======================================================================
    enddo ! loop over k-points ik
! ======================================================================

! ======================================================================
  enddo ! loop over q-points iq
! ======================================================================

  goto 999

!! ======================================================================
!! loop over qvec to compute epsilon q
!  do iq=1,qpt%list%nk
!! ======================================================================
!    qvec = qpt%list%kvec(:,iq)
!
!    ! convert qvec to Cartesian (units of tpiba)
!    qvec_cart = matmul( bg, qvec )
!
!    ! generate q+G vectors inside the cutoff
!    call pw_qpgcut( qvec_cart, epsilon_cutoff )
!
!    write(stdout,*) ' q+G vectors : ', nqpg
!    do ig=1,nqpg
!      write(stdout,'(i,4f)') ig, qpgvec(1:3,ig), qpgvec_len(ig)
!    enddo
!    
!    ! establish size of epsilon matrix
!    write(stdout,*) ' size of epsilon matrix for this q = ', nqpg
!    ! report total size in MB
!    write(stdout,*) dble(nqpg*nqpg)*16.d0/(1024.d0**2.d0), ' MB'
!
!
!    ! distribute using scalapack
!    call scalapack_distrib( nqpg, nqpg, nqpg_lr, nqpg_lc )
!    write(200+mpime,*) ' allocating local part of epsilon matrix '
!    write(200+mpime,*) nqpg_lr, ' x ', nqpg_lc
!    allocate( epsmat(nqpg_lr, nqpg_lc) )
!
!    ! check inversion of identity
!    do jg=1,nqpg
!    do ig=1,nqpg
!      call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!      if( .not. islocal ) cycle
!
!      if( ig==jg ) then
!        epsmat(ig_l,jg_l) = ONE
!      else
!        epsmat(ig_l,jg_l) = ZERO
!      endif
!    enddo
!    enddo
!
!    call scalapack_invert( nqpg, epsmat )
!
!    do jg=1,nqpg
!    do ig=1,nqpg
!      call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!      if( .not. islocal ) cycle
!
!      if( ig==jg ) then
!        if( abs( epsmat(ig_l,jg_l) - ONE ) > 1.d-12 ) &
!          call errore('shirley_epsilon','failure to invert identity',1)
!      else
!        if( abs( epsmat(ig_l,jg_l) ) > 1.d-12 ) &
!          call errore('shirley_epsilon','failure to invert identity',1)
!      endif
!
!    enddo
!    enddo
!
!
!    ! loop over k
!    do ik=1,kpt%list%nk
!      write(stdout,*) ' construct contribution from k-point ', ik
!
!      if( kpt%param%cartesian ) then
!        kpqvec = kpt%list%kvec(:,ik) + qvec_cart
!      else
!        kpqvec = kpt%list%kvec(:,ik) + qvec
!      endif
!    
!      ! ----------------------------------------------------------------------
!      ! construct H(k) and H(k+q) and diagonalize
!      ! ----------------------------------------------------------------------
!
!      ! build hamiltonian for this k ...
!      write(stdout,*) ' build ham k'
!      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )
!      ! ... and diagonalize
!      write(stdout,*) ' diag ham k'
!      call diag_ham
!      ! store
!      eigvec_k = eigvec
!      eigval_k = eigval
!      write(stdout,*) ' done with k'
!
!      ! build hamiltonian for k+q ...
!      write(stdout,*) ' build ham k+q'
!      call diag_build_hamk( kpqvec, kpt%param%cartesian )
!      ! ... and diagonalize
!      write(stdout,*) ' diag ham k+q'
!      call diag_ham
!
!      ! fermi occupation factors
!      ! efermi read from input in epsilon_module
!      do i=1,neig
!        focc_k(i) = spindeg * fermifunc( eigval_k(i), efermi, kT )
!        focc_kpq(i) = spindeg * fermifunc( eigval(i), efermi, kT )
!      enddo
!
!      do i=neig,1,-1
!        if( abs(focc_k(i)) > 1.d-12 ) exit
!      enddo
!      ikoccmax=i
!      do j=neig,1,-1
!        if( abs(focc_kpq(j)) > 1.d-12 ) exit
!      enddo
!      ikpqoccmax=j
!
!      write(stdout,*) ' index of highest occupied band for   k = ', ikoccmax
!      write(stdout,*) ' index of highest occupied band for k+q = ', ikpqoccmax
!
!      ! redefine nnz based on two matrix blocks
!      nnz = ikpqoccmax*(neig-ikoccmax)+ikoccmax*(neig-ikpqoccmax)
!      write(stdout,*) ' non-zero occupation differences = ', nnz, &
!                      ' out of possible ', neig*neig
!      allocate( wx(nnz), xivec(nnz), xjvec(nnz) )
!
!      nnz=0
!      do j=1,ikpqoccmax
!        do i=ikoccmax+1,neig
!          nnz=nnz+1
!          wx(nnz) = cmplx( (focc_kpq(j) - focc_k(i))/(eigval(j)-eigval_k(i)) )
!        enddo
!      enddo
!      do j=ikpqoccmax+1,neig
!        do i=1,ikoccmax
!          nnz=nnz+1
!          wx(nnz) = cmplx( (focc_kpq(j) - focc_k(i))/(eigval(j)-eigval_k(i)) )
!        enddo
!      enddo
!
!      ! precalculate the plane-wave matrix elements 
!      call mp_barrier
!      write(stdout,*) ' dumping pw matrix to tmpfile'
!      if( ionode ) then
!        open(iunpwmtmp,file=trim(tmpfile),form='unformatted', &
!             access='direct',action='write',recl=(16*size(xivec)))
!      endif
!      do ig=1,nqpg
!        ! distribute
!        ipsour = mod(ig,nproc)
!        if( ipsour == mpime ) then
!          ! load plane-wave matrix elements from file for records G and G'
!          ! interpolating to k and k+q
!          call pwmtxel_int( ig, eigvec_k, eigvec, ikoccmax, ikpqoccmax, &
!                            xivec(1:(neig-ikoccmax)*ikpqoccmax), &
!                            xivec((neig-ikoccmax)*ikpqoccmax+1:nnz) )
!          ! send to ionode
!        endif
!        if( ipsour /= ionode_id ) then
!          CALL mp_get( xivec, xivec, mpime, ionode_id, ipsour, ig )
!        endif
!        if( ionode ) then
!          ! write
!          write(stdout,*) ' dumping G-vector ', ig
!          write(iunpwmtmp,rec=ig) xivec
!        endif
!      enddo
!      if( ionode ) close(iunpwmtmp)
!
!      call mp_barrier
!      write(stdout,*) ' opening tmpfile for read access only'
!      open(iunpwmtmp,file=trim(tmpfile),form='unformatted',status='old', &
!           access='direct',action='read',recl=(16*size(xivec)))
!          
!      ! ----------------------------------------------------------------------
!      ! loop over G' and G
!      ! ----------------------------------------------------------------------
!      jg_last = 0
!      do jg=1,nqpg
!      do ig=1,nqpg
!        call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!        if( .not. islocal ) cycle
!
!!        ! load plane-wave matrix elements from file for records G and G'
!!        ! interpolating to k and k+q
!!        call pwmtxel_int( ig, eigvec_k, eigvec, ikoccmax, ikpqoccmax, &
!!                          xivec(1:(neig-ikoccmax)*ikpqoccmax), &
!!                          xivec((neig-ikoccmax)*ikpqoccmax+1:nnz) )
!!!        xivec = pack( xmat, wmask )
!!        !xivec = wx * xivec 
!!        if( jg /= jg_last ) then
!!          call pwmtxel_int( jg, eigvec_k, eigvec, ikoccmax, ikpqoccmax, &
!!                            xjvec(1:(neig-ikoccmax)*ikpqoccmax), &
!!                            xjvec((neig-ikoccmax)*ikpqoccmax+1:nnz) )
!!!          xjvec = pack( xmat, wmask )
!!          xjvec = wx * xjvec 
!!          jg_last = jg
!!        endif
!
!        ! read from file
!        read(iunpwmtmp,rec=ig) xivec
!        if( jg /= jg_last ) then
!          read(iunpwmtmp,rec=jg) xjvec
!          xjvec = wx * xjvec
!          jg_last = jg
!        endif
!
!        ! construct contribution to epsmat as dot_product
!        ! the factor of two converts Rydberg to Hartree
!        epsmat(ig_l,jg_l) = epsmat(ig_l, jg_l) &
!          - kpt%list%wk(ik) * zdotc( nnz, xivec, 1, xjvec, 1 ) &
!            * cmplx( (2.d0*fpi)/(omega*qpgvec_len(ig)*qpgvec_len(jg)) )
!        write(stdout,'(2i,4e)') ig_l, jg_l, epsmat(ig_l,jg_l), &
!                                qpgvec_len(ig), qpgvec_len(jg)
!      enddo
!      enddo
!
!      deallocate( wx, xivec, xjvec )
!
!      close(iunpwmtmp)
!      call mp_barrier
!          
!! ======================================================================
!    enddo ! loop over k-points ik
!! ======================================================================
!
!    ! report epsilon without local fields
!    ig=1 ; jg=1
!    call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!    epsilon_EC = 0.d0
!    if( islocal ) then
!      epsilon_EC = real( epsmat(ig_l,jg_l) )
!    endif
!    call mp_sum( epsilon_EC )
!    write(stdout,*) ' epsilon (Ehrenreich-Cohen) = ', epsilon_EC
!
!  
!    write(stdout,*) ' inverting epsilon for this q using scalapack'
!    call scalapack_invert( nqpg, epsmat )
!
!    ! report epsilon with local fields
!    ig=1 ; jg=1
!    call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!    epsilon_LF = 0.d0
!    if( islocal ) then
!      epsilon_LF = 1.d0 / real( epsmat(ig_l,jg_l) )
!    endif
!    call mp_sum( epsilon_LF )
!    write(stdout,*) ' epsilon     (Local Fields) = ', epsilon_LF
!
!
!    do jg=1,nqpg
!    do ig=1,nqpg
!      call scalapack_localindex( ig, jg, ig_l, jg_l, islocal )
!      if( .not. islocal ) cycle
!
!      write(iunout,'(2i,2e)') ig, jg, epsmat(ig_l,jg_l)
!    enddo
!    enddo
!
!! ======================================================================
!  enddo ! loop over q-points iq
!! ======================================================================
!
!  999 continue
!
!  ! close binary dump
!  close(iunout)

  999 write(stdout,*) ' waiting for other nodes'
  call stop_clock( 'shirley' )
  call stop_shirley
  
  end program shirley_epstci
