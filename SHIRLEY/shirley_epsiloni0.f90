  program shirley_epsiloni0

  ! Generate the wing of the dielectric matrix using Shirley interpolation
  ! This is for q->0 only and involves velocity matrix elements

  ! David Prendergast, MF, Jun 2i08

#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum, mp_get, &
                 mp_scatter_size, mp_scatter
  use kpt_module
  use shirley_input_module
  use diag_module
  use fermi

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6

  COMPLEX(DP),PARAMETER :: ZERO=CMPLX(0.d0,0.d0)
  COMPLEX(DP),PARAMETER :: ONE =CMPLX(1.d0,0.d0)

  integer :: neig2, nnz, inz

  integer :: ik, i, j, k
  integer :: ipol, apol

  integer :: iunint, iuntri, iuneps, ierr
  logical :: exst

  real(dp) :: fpi
  real(dp) :: kT
  complex(dp) :: fac

  real(dp) :: dq, kpq(3)
  integer :: iq

  ! tolerance for determining zeros
  real(dp) :: tol

  real(dp) :: mem

  complex(dp),allocatable :: vint0(:,:)
  integer :: ntriple
  complex(dp),allocatable :: triple_i(:,:,:)
  real(dp) :: onegig, triple_size

  ! system variables
  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba

  ! external
  integer,external :: freeunit
  complex(dp),external :: ZDOTU


  call start_clock( 'shirley' )
 
  fpi = acos(-1.d0)*4.d0

  call shirley_input

  call diag_init

  ! allocations
  allocate( vint0(neig,neig) )

  ! check for size of triples
  if( ionode ) then
    onegig=1024.d0**3.d0 ! 1 GB
    triple_size=dble(neig)**2.d0*16.d0  ! complex numbers
    write(stdout,*) ' memory factor (of 1 GB/proc) = ', memfac
    ! try to have neig evenly divided between processes
    ! but if that is too much then it needs to be reduced
    ntriple=min( 2*neig/nproc, int(onegig*memfac/triple_size)*2 )
    ierr=1
    do while( ierr/=0 )
      if( ntriple==1 ) call errore('shirley_epsiloni0','unable to allocate space for triples',1)
      ntriple = max( ntriple / 2, 1 )
      allocate( triple_i(neig,neig,ntriple), stat=ierr )
    enddo
    write(stdout,*) ' number of triples allocated = ', ntriple
    write(stdout,*) ' memory required (GB) = ', dble(ntriple)*triple_size/onegig
  endif

  call mp_bcast( ntriple, ionode_id )
  if( .not. ionode ) then
    allocate( triple_i(neig,neig,ntriple), stat=ierr )
    if( ierr/=0 ) call errore('shirley_epsiloni0','unable to allocate space for triples off ionode',1)
  endif
  
  ! load int0
  if( ionode ) then
    iunint = freeunit()
    open(iunint,file=trim(int0file),form='unformatted')
    read(iunint) vint0
    close(iunint)
  endif
  call mp_bcast( vint0, ionode_id )

  ! open triples file
  neig2=neig*neig
  iuntri=freeunit()
  open(iuntri,file=trim(trifile),form='unformatted',access='direct', &
       recl=neig2*2*DIRECT_IO_FACTOR)

  ! band_subset
  if( band_subset(1) > band_subset(2) ) band_subset = cshift( band_subset, 1 )
  if( band_subset(1) < 1 .or. band_subset(1) > neig ) band_subset(1)=1
  if( band_subset(2) < 1 .or. band_subset(2) > neig ) band_subset(2)=neig
  write(stdout,*) ' band_subset = ', band_subset

  write(stdout,*) '          efermi [eV] = ', efermi
  write(stdout,*) '   electron temp  [K] = ', temperature

  ! convert efermi to Ry
  efermi = efermi / rytoev
  kT = temperature * kelvin2rydberg

  ! get system details from hamq_shirley
  call dump_system( nelec_, alat, omega, at, bg, tpiba )

  WRITE( stdout, '(5X, &
       &     "crystal axes: (cart. coord. in units of a_0)",/, &
       &       3(15x,"a(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,  &
       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  WRITE( stdout, '(5x, &
       &   "reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  write(stdout,*) '          volume = ', omega


! ======================================================================
! Begin calculation of wing epsiloni0(q->0) 
! ======================================================================

  ! sampling dependent tolerance
  tol = 1.d-10/dble(kpt%list%nk)

  call wing

  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call stop_clock( 'shirley' )
  call stop_shirley
  
  contains


! ======================================================================
  subroutine wing
! ======================================================================

  complex(dp),allocatable :: wx(:), de(:)
  real(dp),allocatable :: focc(:)
  complex(dp),allocatable :: epsi0(:,:)
  complex(dp),allocatable :: ztmp(:,:)
  complex(dp),allocatable :: dhdk(:,:,:)
  integer,allocatable :: index_nnz(:,:)
  integer :: nocc, npar, nemp
  integer,allocatable :: index_occ(:), index_par(:), index_emp(:), index_ope(:,:)
  integer :: n, m
  integer :: nf, nt, mf, mt
  complex(dp),allocatable :: eigvec_occ(:,:), eigvec_par(:,:), eigvec_emp(:,:)
  complex(dp),allocatable :: triple_op(:,:), triple_oe(:,:)
  complex(dp),allocatable :: triple_po(:,:), triple_pp(:,:), triple_pe(:,:)
  complex(dp),allocatable :: triple_eo(:,:), triple_ep(:,:)

  integer :: nprocs, iam
  integer :: ictxt(2)
  integer :: myrow(2), mycol(2)
  integer :: nprow(2), npcol(2)
  complex(dp),allocatable :: xvec1(:,:), xvec2(:,:)
  complex(dp),allocatable :: fdvec2(:,:)
  integer :: nx1(2), nx2(2)
  integer,dimension(9) :: desc1, desc2
  integer :: info
  integer :: ig, icol, irow, inzg, ideriv

  complex(dp) :: ztrkin2nlp(3,3)

  integer,external :: numroc


  ! set up 2 contexts 1xN and Nx1
  nprow(1) = 1
  npcol(1) = nproc
  call sl_init( ictxt(1), nprow(1), npcol(1) )

  nprow(2) = nproc
  npcol(2) = 1
  call sl_init( ictxt(2), nprow(2), npcol(2) )

  allocate( epsi0(3,neig) )
  allocate( ztmp(neig,neig) )
  allocate( dhdk(neig,neig,3) )
  allocate( focc(neig) )

  epsi0 = zero

  ! loop over k
  do ik=1,kpt%list%nk
    write(stdout,*) ' construct contribution from k-point ', &
                    ik, ' of ', kpt%list%nk, kpt%list%wk(ik)

    ! factor of 2 for spin-degeneracy should be included in the kpt weights
    fac = cmplx( kpt%list%wk(ik)*fpi/omega )

    ! ----------------------------------------------------------------------
    ! construct H(k) and diagonalize
    ! ----------------------------------------------------------------------
    write(stdout,*) ' build ham k'
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian )
    write(stdout,*) ' diag ham k'
    call diag_ham
    write(stdout,*) ' done with k'

    ! ----------------------------------------------------------------------
    ! construct dH/dk(k) in the Shirley basis (not eigenbasis)
    ! ----------------------------------------------------------------------
    do ideriv=1,3
      call diag_build_dham( ideriv, kpt%list%kvec(1:3,ik), &
                            kpt%param%cartesian, dhdk(:,:,ideriv) )
    enddo

    ! ----------------------------------------------------------------------
    ! fermi occupation factors
    ! efermi read from input in epsilon_module
    ! ----------------------------------------------------------------------
    do i=1,neig
      focc(i) = fermifunc( eigval(i), efermi, kT )
    enddo

    ! store occupied, partial, and empty indices
    allocate( index_occ(neig), index_par(neig), index_emp(neig), index_ope(2,neig) )
    nocc=0
    npar=0
    nemp=0
    do i=1,neig
      if( focc(i) > 1.d0 - tol ) then
        nocc=nocc+1
        index_occ(nocc) = i
        index_ope(:,i) = (/ 1, nocc /)
      elseif( focc(i) > tol ) then
        npar=npar+1
        index_par(npar) = i
        index_ope(:,i) = (/ 2, npar /)
      else
        nemp=nemp+1
        index_emp(nemp) = i
        index_ope(:,i) = (/ 3, nemp /)
      endif
    enddo

    write(stdout,*) ' nocc = ', nocc
    write(stdout,*) ' npar = ', npar
    write(stdout,*) ' nemp = ', nemp
    if( nocc+npar+nemp /= neig ) call errore('shirley_epsiloni0','nocc+npar+nemp /= neig',1)

    ! now set the eigenvector subsets
    allocate( eigvec_occ(neig,nocc), eigvec_par(neig,npar), eigvec_emp(neig,nemp) )
    forall(i=1:neig,j=1:nocc) eigvec_occ(i,j) = eigvec(i,index_occ(j))
    forall(i=1:neig,j=1:npar) eigvec_par(i,j) = eigvec(i,index_par(j))
    forall(i=1:neig,j=1:nemp) eigvec_emp(i,j) = eigvec(i,index_emp(j))

    allocate( triple_op(nocc,npar), triple_oe(nocc,nemp), &
              triple_po(npar,nocc), triple_pp(npar,npar), triple_pe(npar,nemp), &
              triple_eo(nemp,nocc), triple_ep(nemp,npar) )

    ! determine non-zero occupation factors
    allocate( index_nnz(2,neig2) )
    nnz=0
    do j=1,neig
      do i=1,neig
        if( abs(focc(i)-focc(j)) < tol ) cycle
        nnz=nnz+1
        index_nnz(:,nnz) = (/ i, j /)
      enddo
    enddo
    write(stdout,*) ' non-zero occupation differences = ', nnz, &
                    ' out of possible ', neig2

    allocate( wx(nnz), de(nnz) )
    ! need a check here for de=zero
    do inz=1,nnz
      i = index_nnz(1,inz)
      j = index_nnz(2,inz)
      de(inz) = cmplx( eigval(j)-eigval(i) )
      wx(inz) = cmplx( focc(i) - focc(j) )/de(inz)* 2.0d0 ! Ry to Ha
    enddo


    ! establish position in 1xN process grid(1), local size of xvec, and allocate
    call blacs_gridinfo( ictxt(1), nprow(1), npcol(1), myrow(1), mycol(1) )
    nx1(1)=numroc( nnz, ntriple, myrow(1), 0, nprow(1))
    nx2(1)=numroc(neig, ntriple, mycol(1), 0, npcol(1))
    allocate( xvec1(nx1(1),nx2(1)) ) 
    CALL DESCINIT( desc1, nnz, neig, ntriple, ntriple, 0, 0, &
                   ictxt(1), nx1(1), info )
    write(stdout,*) ' process grid 1: xvec dimensions ', size(xvec1,1), size(xvec1,2)

    ! establish position in Nx1 process grid(2), local size of xvec, and allocate
    call blacs_gridinfo( ictxt(2), nprow(2), npcol(2), myrow(2), mycol(2) )
    nx1(2)=numroc( nnz, ntriple, myrow(2), 0, nprow(2))
    nx2(2)=numroc(neig, ntriple, mycol(2), 0, npcol(2))
    allocate( xvec2(nx1(2),nx2(2)), fdvec2(nx1(2),3) ) 
    CALL DESCINIT( desc2, nnz, neig, ntriple, ntriple, 0, 0, &
                   ictxt(2), nx1(2), info )
    write(stdout,*) ' process grid 2: xvec dimensions ', size(xvec2,1), size(xvec2,2)
    write(stdout,*) ' process grid 2: fdvec dimensions ', size(fdvec2,1), size(fdvec2,2)


    ! compute xvec1 in distributed fashion reading every ntriple records
    write(stdout,*) ' compute xvec'
    do i=1,neig
      ! is this index to be handled locally?
      call INFOG1L( i, ntriple, nproc, mycol(1), 0, ig, icol )
      if( mycol(1) /= icol ) cycle
      
      write(stdout,*) ' triple_i=', i, ig, ' of ', nx2(1)
      if( mod(ig-1,ntriple)==0 ) then
        write(stdout,*) ' read ', min(ntriple,nx2(1)-ig+1), ' triples'
        do j=1,min(ntriple,nx2(1)-ig+1)
          ! load triple i
          read(iuntri,rec=i+j-1) triple_i(:,:,j)
        enddo
      endif

      j=mod(ig-1,ntriple)+1

      ! transform to eigenbasis 
      call matmul_CNN( nocc, npar, neig, eigvec_occ, triple_i(1,1,j), eigvec_par, ztmp, triple_op )
      call matmul_CNN( nocc, nemp, neig, eigvec_occ, triple_i(1,1,j), eigvec_emp, ztmp, triple_oe )
      call matmul_CNN( npar, nocc, neig, eigvec_par, triple_i(1,1,j), eigvec_occ, ztmp, triple_po )
      call matmul_CNN( npar, npar, neig, eigvec_par, triple_i(1,1,j), eigvec_par, ztmp, triple_pp )
      call matmul_CNN( npar, nemp, neig, eigvec_par, triple_i(1,1,j), eigvec_emp, ztmp, triple_pe )
      call matmul_CNN( nemp, nocc, neig, eigvec_emp, triple_i(1,1,j), eigvec_occ, ztmp, triple_eo )
      call matmul_CNN( nemp, npar, neig, eigvec_emp, triple_i(1,1,j), eigvec_par, ztmp, triple_ep )

      do inz=1,nnz
        n = index_nnz(1,inz)
        nf = index_ope(1,n)
        nt = index_ope(2,n)

        m = index_nnz(2,inz)
        mf = index_ope(1,m)
        mt = index_ope(2,m)

        if( nf==1 .and. mf==2 ) then
          xvec1(inz,ig) = triple_op(nt,mt)
        elseif( nf==1 .and. mf==3 ) then
          xvec1(inz,ig) = triple_oe(nt,mt)
        elseif( nf==2 .and. mf==1 ) then
          xvec1(inz,ig) = triple_po(nt,mt)
        elseif( nf==2 .and. mf==2 ) then
          xvec1(inz,ig) = triple_pp(nt,mt)
        elseif( nf==2 .and. mf==3 ) then
          xvec1(inz,ig) = triple_pe(nt,mt)
        elseif( nf==3 .and. mf==1 ) then
          xvec1(inz,ig) = triple_eo(nt,mt)
        elseif( nf==3 .and. mf==2 ) then
          xvec1(inz,ig) = triple_ep(nt,mt)
        endif
      enddo
    enddo  ! i=1:neig

    ! redistribute xvec col -> row
    write(stdout,*) ' redistribute'
    call PZGEMR2D( nnz, neig, xvec1, 1, 1, desc1, &
                   xvec2, 1, 1, desc2, ictxt(2) )


    ! compute dvec1 from dH/dk
    write(stdout,*) ' compute fdvec'
    do ideriv=1,3
      write(stdout,*) ' deriv=', ideriv, ' of ', 3

      ! transform to eigenbasis 
      call matmul_CNN( nocc, npar, neig, eigvec_occ, dhdk(1,1,ideriv), eigvec_par, ztmp, triple_op )
      call matmul_CNN( nocc, nemp, neig, eigvec_occ, dhdk(1,1,ideriv), eigvec_emp, ztmp, triple_oe )
      call matmul_CNN( npar, nocc, neig, eigvec_par, dhdk(1,1,ideriv), eigvec_occ, ztmp, triple_po )
      call matmul_CNN( npar, npar, neig, eigvec_par, dhdk(1,1,ideriv), eigvec_par, ztmp, triple_pp )
      call matmul_CNN( npar, nemp, neig, eigvec_par, dhdk(1,1,ideriv), eigvec_emp, ztmp, triple_pe )
      call matmul_CNN( nemp, nocc, neig, eigvec_emp, dhdk(1,1,ideriv), eigvec_occ, ztmp, triple_eo )
      call matmul_CNN( nemp, npar, neig, eigvec_emp, dhdk(1,1,ideriv), eigvec_par, ztmp, triple_ep )

      do inz=1,nnz
        call INFOG1L( inz, ntriple, nproc, myrow(2), 0, inzg, irow )
        if( myrow(2) /= irow ) cycle

        n = index_nnz(1,inz)
        nf = index_ope(1,n)
        nt = index_ope(2,n)

        m = index_nnz(2,inz)
        mf = index_ope(1,m)
        mt = index_ope(2,m)

        if( nf==1 .and. mf==2 ) then
          fdvec2(inzg,ideriv) = triple_op(nt,mt)
        elseif( nf==1 .and. mf==3 ) then
          fdvec2(inzg,ideriv) = triple_oe(nt,mt)
        elseif( nf==2 .and. mf==1 ) then
          fdvec2(inzg,ideriv) = triple_po(nt,mt)
        elseif( nf==2 .and. mf==2 ) then
          fdvec2(inzg,ideriv) = triple_pp(nt,mt)
        elseif( nf==2 .and. mf==3 ) then
          fdvec2(inzg,ideriv) = triple_pe(nt,mt)
        elseif( nf==3 .and. mf==1 ) then
          fdvec2(inzg,ideriv) = triple_eo(nt,mt)
        elseif( nf==3 .and. mf==2 ) then
          fdvec2(inzg,ideriv) = triple_ep(nt,mt)
        endif

        fdvec2(inzg,ideriv) = wx(inz) * fdvec2(inzg,ideriv)
      enddo  ! inz=1:nnz
    enddo  ! ideriv=1:3

        
    ! accumulate epsi0
    write(stdout,*) ' accumulate epsi0'
    call ZGEMM( 'C', 'N', neig, 3, nx1(2), fac, xvec2, nx1(2), fdvec2, nx1(2), ONE, epsi0, neig )
    write(stdout,*)
        

    deallocate( wx, de )
    deallocate( index_nnz )
    deallocate( index_occ, index_par, index_emp, index_ope )
    deallocate( eigvec_occ, eigvec_par, eigvec_emp )

    deallocate( triple_op, triple_oe, &
                triple_po, triple_pp, triple_pe, &
                triple_eo, triple_ep )

    deallocate( xvec1 )
    deallocate( xvec2, fdvec2 )
      
! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! sum epsi0 over processors
  write(stdout,*) ' sum over processors'
  call mp_sum( epsi0 )

  ! multiply by interaction and apply correct dimensions to dH/dk
  ztrkin2nlp = trkin2nlp
  call ZGEMM( 'N', 'N', neig, 3,    3, ONE, epsi0, neig, ztrkin2nlp,    3, ZERO,  ztmp, neig )
  call ZGEMM( 'C', 'N', neig, 3, neig, ONE, vint0, neig,       ztmp, neig, ZERO, epsi0, neig )

  ! dump to file
  if( ionode ) then
    write(stdout,*) 'epsi0'
    do j=1,neig
      write(stdout,'(i,6f)') j, epsi0(j,1:3)
    enddo
    iuneps = freeunit()
    open(iuneps,file=trim(epsi0file),form='unformatted')
    write(iuneps) epsi0
    close(iuneps)
  endif

  deallocate( focc )
  deallocate( epsi0, ztmp )

  end subroutine wing


  end program shirley_epsiloni0

! ======================================================================
  subroutine matmul_CNN( n, m, l, a, b, c, w, abc )
! ======================================================================

  use kinds, only : dp
  implicit none

  COMPLEX(DP),PARAMETER :: ZERO=CMPLX(0.d0,0.d0)
  COMPLEX(DP),PARAMETER :: ONE =CMPLX(1.d0,0.d0)

  integer,intent(in) :: n, m, l
  complex(dp),intent(in) :: a(l,n), b(l,l), c(l,m)
  complex(dp) :: w(l,l)
  complex(dp),intent(out) :: abc(n,m)

  if( n==0 .OR. m==0 .OR. l==0 ) return

  if( n < m ) then
    ! (A*B)*C
    call ZGEMM( 'C', 'N', n, l, l, ONE, a, l, b, l, ZERO, w, l )
    call ZGEMM( 'N', 'N', n, m, l, ONE, w, l, c, l, ZERO, abc, n )
  else
    ! A*(B*C)
    call ZGEMM( 'N', 'N', l, m, l, ONE, b, l, c, l, ZERO, w, l )
    call ZGEMM( 'C', 'N', n, m, l, ONE, a, l, w, l, ZERO, abc, n )
  endif

  end subroutine matmul_CNN

