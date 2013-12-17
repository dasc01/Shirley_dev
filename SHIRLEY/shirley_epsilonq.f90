  program shirley_epsilonq

  ! Generate the dielectric matrix using Shirley interpolation
  ! This is for q/=0 only 

  ! David Prendergast, MF, Jun 2008

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

  if( ionode ) then
    onegig=1024.d0**3.d0 ! 1 GB
    triple_size=dble(neig)**2.d0*16.d0  ! complex numbers
    write(stdout,*) ' memory factor (of 1 GB/proc) = ', memfac
    ! try to have neig evenly divided between processes
    ! but if that is too much then it needs to be reduced
    ntriple=min( 2*neig/nproc, int(onegig*memfac/triple_size)*2 )
    ierr=1
    do while( ierr/=0 )
      if( ntriple==1 ) call errore('shirley_epsilonq','unable to allocate space for triples',1)
      ntriple = max( ntriple / 2, 1 )
      allocate( triple_i(neig,neig,ntriple), stat=ierr )
    enddo
    write(stdout,*) ' number of triples allocated = ', ntriple
    write(stdout,*) ' memory required (GB) = ', dble(ntriple)*triple_size/onegig
  endif

  call mp_bcast( ntriple, ionode_id )
  if( .not. ionode ) then
    allocate( triple_i(neig,neig,ntriple), stat=ierr )
    if( ierr/=0 ) call errore('shirley_epsilonq','unable to allocate space for triples off ionode',1)
  endif
  

  ! load int0
  if( ionode ) then
    iunint = freeunit()
    open(iunint,file=trim(int0file),form='unformatted')
    read(iunint) vint0
    close(iunint)
!    do j=1,neig
!      do i=1,neig
!        write(200,*) vint0(i,j)
!      enddo
!      write(200,*)
!    enddo
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

!  ! scatter k-points
!  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

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
! Begin calculation of body epsilonq(q->0) 
! ======================================================================

  ! sampling dependent tolerance
  tol = 1.d-10/dble(kpt%list%nk)

  call body

  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call stop_clock( 'shirley' )
  call stop_shirley
  
  contains


! ======================================================================
  subroutine body
! ======================================================================

  use paulipairs

  complex(dp),allocatable :: wx(:)
  real(dp),allocatable :: focc(:)

  real(dp) :: qvec(3), kpq(3), qveclen
  character(255) :: fileq

  complex(dp),allocatable :: epsq(:,:)
  complex(dp),allocatable :: ztmp(:,:)
  integer :: n, m

  integer :: nprocs, iam
  integer :: ictxt(2)
  integer :: myrow(2), mycol(2)
  integer :: nprow(2), npcol(2)
  complex(dp),allocatable :: xvec1(:,:), xvec2(:,:)
  complex(dp),allocatable :: fxvec2(:,:)
  integer :: nx1(2), nx2(2)
  integer,dimension(9) :: desc1, desc2
  integer :: info
  integer :: ig, icol, irow, inzg

  integer,external :: numroc


  ! set up 2 contexts 1xN and Nx1
  nprow(1) = 1
  npcol(1) = nproc
  call sl_init( ictxt(1), nprow(1), npcol(1) )

  nprow(2) = nproc
  npcol(2) = 1
  call sl_init( ictxt(2), nprow(2), npcol(2) )

  allocate( epsq(neig,neig) )
  allocate( ztmp(neig,neig) )
  allocate( focc(neig) )

  ! this must be called before repeated calls to any of the paulipairs routines
  call init_paulipairs

! ======================================================================
! loop over q
! ======================================================================
  do iq=1,qpt%list%nk
    write(stdout,*) ' construct epsilon for q-point ', &
                    iq, ' of ', qpt%list%nk, qpt%list%wk(iq)
    write(stdout,'(3f)') qpt%list%kvec(1:3,iq)

    qvec = qpt%list%kvec(1:3,iq)
    qveclen = sum( qvec(1:3)**2 )
    if( qveclen < 1.d-12 ) then
      write(stdout,*) ' skipping this q-point as zero'
      cycle
    endif

    epsq = zero

  ! ======================================================================
  ! loop over k
  ! ======================================================================
    do ik=1,kpt%list%nk
      write(stdout,*) ' construct contribution from k-point ', &
                      ik, ' of ', kpt%list%nk, kpt%list%wk(ik)
      write(stdout,'(a,3f)') ' kvec = ', kpt%list%kvec(1:3,ik)

      kpq = qvec
      if( kpt%param%cartesian ) then
        if( .not. qpt%param%cartesian ) then
          kpq = matmul( trnlp2kin, kpq )
        endif
      else
        if( qpt%param%cartesian ) then
          kpq = matmul( trkin2nlp, kpq )
        endif
      endif
      write(stdout,'(a,3f)') ' qvec = ', kpq

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
      ! fermi occupation factors
      ! efermi read from input in epsilon_module
      ! ----------------------------------------------------------------------
      do i=1,neig
        focc(i) = fermifunc( eigval(i), efermi, kT )
      enddo
      call load_occ( focc, tol, 1 )
      call load_eigvec( eigvec, eigval, 1 )

      ! ----------------------------------------------------------------------
      ! construct H(k+q) and diagonalize
      ! ----------------------------------------------------------------------
      kpq = kpt%list%kvec(1:3,ik)+kpq
      write(stdout,'(a,3f)') ' k+q = ', kpq
      write(stdout,*) ' build ham k+q'
      call diag_build_hamk( kpq, kpt%param%cartesian )
      write(stdout,*) ' diag ham k+q'
      call diag_ham
      write(stdout,*) ' done with k+q'
      ! ----------------------------------------------------------------------
      ! fermi occupation factors
      ! ----------------------------------------------------------------------
      do i=1,neig
        focc(i) = fermifunc( eigval(i), efermi, kT )
      enddo
      call load_occ( focc, tol, 2 )
      call load_eigvec( eigvec, eigval, 2 )

      ! initialize space to store partitioned triple
      call init_triple( 2 )

      ! initialize non-zero pairs of occupation factors
      call init_nonzero_pairs( 2, tol )

      write(stdout,*) ' non-zero occupation differences = ', nnz, &
                      ' out of possible ', neig2

      allocate( wx(nnz) )
      do inz=1,nnz
        ! minus sign on energy denominator
        ! factor of 2 converts Ry to Ha
        wx(inz) = cmplx( df_nnz(inz) / (-1.d0 * de_nnz(inz)) * 2.0d0 )
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
      allocate( xvec2(nx1(2),nx2(2)), fxvec2(nx1(2),nx2(2)) ) 
      CALL DESCINIT( desc2, nnz, neig, ntriple, ntriple, 0, 0, &
                     ictxt(2), nx1(2), info )
      write(stdout,*) ' process grid 2: xvec dimensions ', size(xvec2,1), size(xvec2,2)


      ! compute xvec1 in distributed fashion reading every ntriple records
      write(stdout,*) ' compute xvec'
      xvec1=ZERO
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

        ! calculate the product: eigvec* triple eigvec = xvec
        call prod_triple( triple_i(:,:,j), 2, xvec1(:,ig) )
      enddo  ! i=1:neig

      ! redistribute xvec col -> row
      write(stdout,*) ' redistribute'
      xvec2=ZERO
      call PZGEMR2D( nnz, neig, xvec1, 1, 1, desc1, &
                     xvec2, 1, 1, desc2, ictxt(2) )

      ! make weighted fxvec
      write(stdout,*) ' fxvec'
      do i=1,neig
        do inz=1,nnz
          call INFOG1L( inz, ntriple, nproc, myrow(2), 0, inzg, irow )
          if( myrow(2) /= irow ) cycle

          fxvec2(inzg,i) = wx(inz) * xvec2(inzg,i)
        enddo
      enddo
      
        
      ! accumulate epsq
      write(stdout,*) ' accumulate epsq'
      call ZGEMM( 'C', 'N', neig, neig, nx1(2), fac, xvec2, nx1(2), fxvec2, nx1(2), ONE, epsq, neig )
      write(stdout,*)
        
      deallocate( wx )

      deallocate( xvec1 )
      deallocate( xvec2, fxvec2 )
      
  ! ======================================================================
    enddo ! loop over k-points ik
  ! ======================================================================

    ! sum epsq over processors
    write(stdout,*) ' sum over processors'
    call mp_sum( epsq )

!    ! multiply by interaction
!    call ZGEMM( 'N', 'N', neig, neig, neig, ONE, epsq, neig, vint0, neig, ZERO, ztmp, neig )
!    call ZGEMM( 'C', 'N', neig, neig, neig, ONE, vint0, neig, ztmp, neig, ZERO, epsq, neig )
  
    ! add unity to diagonal
    do j=1,neig
      epsq(j,j)=ONE+epsq(j,j)
    enddo

    ! dump to file
    if( ionode ) then
      write(stdout,*) 'epsq diagonals '
      do j=1,neig
        write(stdout,'(i,2f)') j, epsq(j,j)
      enddo
      iuneps = freeunit()
      write(fileq,*) iq
      write(fileq,*) trim(epsqfile) // '.' // trim(adjustl(fileq))
      open(iuneps,file=trim(fileq),form='unformatted')
      write(iuneps) epsq
      close(iuneps)
    endif


! ======================================================================
  enddo ! loop over q-points iq
! ======================================================================


  deallocate( focc )
  deallocate( epsq, ztmp )

  end subroutine body


  end program shirley_epsilonq
