  module paulipairs

  ! This module deals with expansions in the eigenbasis for pairs of
  ! eigenstates that are permitted by the Pauli exclusion principle.
  ! This requires that we know the occupations and pick only pairs
  ! that correspond to non-zero values of
  !   f(e_nk) - f(e_mk+q)
  ! This usually can reduce the workload for evaluating products like
  !   [ Z(k+q)* Xi Z(k) ]_nm

  use kinds, only : dp
  use io_global, only : stdout

  implicit none

  private
  public :: load_occ, load_eigvec, &
            init_triple, prod_triple, &
            init_nonzero_pairs, init_paulipairs
  public :: nnz, de_nnz, df_nnz

  integer :: nstate, nnz
  integer,allocatable :: index_nnz(:,:)
  real(dp),allocatable :: de_nnz(:), df_nnz(:)

  logical,save :: virgin=.true.

  type paulipair_type 
    real(dp),pointer :: focc(:)
    integer :: nocc, npar, nemp
    integer,pointer :: iocc(:), ipar(:), iemp(:)
    integer,pointer :: iope(:,:)
    complex(dp),pointer :: eocc(:,:), epar(:,:), eemp(:,:)
    real(dp),pointer :: eigval(:)
  end type paulipair_type
  type(paulipair_type) :: pp(2)

  complex(dp),allocatable :: triple_op(:,:), triple_oe(:,:), &
                             triple_po(:,:), triple_pp(:,:), triple_pe(:,:), &
                             triple_eo(:,:), triple_ep(:,:)

  contains

! ----------------------------------------------------------------------
  subroutine init_paulipairs
! ----------------------------------------------------------------------
  integer :: ik
  if( virgin ) then
    virgin = .false.

    do ik=1,size(pp)
      nullify( pp(ik)%focc, pp(ik)%iocc, pp(ik)%ipar, pp(ik)%iemp, &
               pp(ik)%iope, pp(ik)%eocc, pp(ik)%epar, pp(ik)%eemp, &
               pp(ik)%eigval )
    enddo
  endif
  end subroutine init_paulipairs

! ----------------------------------------------------------------------
  subroutine init_indices( focc, tol, nocc, npar, nemp, &
                           iocc, ipar, iemp, iope )
! ----------------------------------------------------------------------
  real(dp) :: focc(nstate), tol
  integer :: nocc, npar, nemp
  integer :: iocc(nstate), ipar(nstate), iemp(nstate), iope(2,nstate)
  integer :: i

  nocc=0
  npar=0
  nemp=0
  do i=1,nstate
    if( focc(i) > 1.d0 - tol ) then
      nocc=nocc+1
      iocc(nocc)=i
      iope(:,i)=(/ 1, nocc /)
    elseif( focc(i) > tol ) then
      npar=npar+1
      ipar(npar)=i
      iope(:,i)=(/ 2, npar /)
    else
      nemp=nemp+1
      iemp(nemp)=i
      iope(:,i)=(/ 3, nemp /)
    endif
  enddo
  if( nocc+npar+nemp /= nstate ) &
    call errore('init_indices','nocc+npar+nemp/=neig',1)

  end subroutine init_indices


! ----------------------------------------------------------------------
  subroutine load_occ( focc_, tol, ik )
! ----------------------------------------------------------------------
  real(dp),intent(in) :: focc_(:), tol
  integer,intent(in) :: ik

  if( ik < 1 .or. ik > 2 ) &
    call errore('load_occ','number of occupancies not 1 or 2',1)

  if( ik==1 ) then
    nstate=size(focc_)
  else
    if( nstate /= size(focc_) ) call errore('load_occ','number of states does not agree',1)
  endif

  if( associated( pp(ik)%focc ) ) deallocate( pp(ik)%focc )
  allocate( pp(ik)%focc(nstate) )
  pp(ik)%focc = focc_

  if( associated( pp(ik)%iocc ) ) deallocate( pp(ik)%iocc )
  if( associated( pp(ik)%ipar ) ) deallocate( pp(ik)%ipar )
  if( associated( pp(ik)%iemp ) ) deallocate( pp(ik)%iemp )
  if( associated( pp(ik)%iope ) ) deallocate( pp(ik)%iope )
  allocate( pp(ik)%iocc(nstate), pp(ik)%ipar(nstate), pp(ik)%iemp(nstate), &
            pp(ik)%iope(2,nstate) )

  call init_indices( pp(ik)%focc, tol, &
                     pp(ik)%nocc, pp(ik)%npar, pp(ik)%nemp, &
                     pp(ik)%iocc, pp(ik)%ipar, pp(ik)%iemp, &
                     pp(ik)%iope )
    
  write(stdout,*) ' for k ', ik
  write(stdout,*) ' nocc = ', pp(ik)%nocc
  write(stdout,*) ' npar = ', pp(ik)%npar
  write(stdout,*) ' nemp = ', pp(ik)%nemp

  end subroutine load_occ

! ----------------------------------------------------------------------
  subroutine init_nonzero_pairs( nk, tol )
! ----------------------------------------------------------------------
  ! determine which pairs of indices have nonzero pauli differences
  integer,intent(in) :: nk
  real(dp),intent(in) :: tol
  integer :: ik, jk
  integer :: i, j

  if( nk==1 ) then
    ! only one set of dimensions
    ik=1 ; jk=1
  elseif( nk==2 ) then
    ik=1 ; jk=2
  endif
  nnz=0
  do j=1,nstate
    do i=1,nstate
      if( abs(pp(ik)%focc(i)-pp(jk)%focc(j)) < tol ) cycle
      nnz=nnz+1
    enddo
  enddo

  if( allocated( index_nnz ) ) deallocate( index_nnz )
  if( allocated( de_nnz ) ) deallocate( de_nnz )
  if( allocated( df_nnz ) ) deallocate( df_nnz )
  allocate( index_nnz(2,nnz), de_nnz(nnz), df_nnz(nnz) )

  nnz=0
  do j=1,nstate
    do i=1,nstate
      if( abs(pp(ik)%focc(i)-pp(jk)%focc(j)) < tol ) cycle
      nnz=nnz+1
      index_nnz(:,nnz) = (/ i, j /)
      ! all these difference are (ik,i) - (jk,j)
      de_nnz(nnz) = pp(ik)%eigval(i) - pp(jk)%eigval(j)
      df_nnz(nnz) = pp(ik)%focc(i)-pp(jk)%focc(j)
    enddo
  enddo

  end subroutine init_nonzero_pairs

! ----------------------------------------------------------------------
  subroutine init_eigvec( eigvec_, ik )
! ----------------------------------------------------------------------
  complex(dp),intent(in) :: eigvec_(:,:)
  integer,intent(in) :: ik
  integer :: i, j

  if( associated( pp(ik)%eocc ) ) deallocate( pp(ik)%eocc )
  if( associated( pp(ik)%epar ) ) deallocate( pp(ik)%epar )
  if( associated( pp(ik)%eemp ) ) deallocate( pp(ik)%eemp )
  allocate( pp(ik)%eocc(nstate,pp(ik)%nocc), &
            pp(ik)%epar(nstate,pp(ik)%npar), & 
            pp(ik)%eemp(nstate,pp(ik)%nemp) )

  forall(i=1:nstate,j=1:pp(ik)%nocc) pp(ik)%eocc(i,j) = eigvec_(i,pp(ik)%iocc(j))
  forall(i=1:nstate,j=1:pp(ik)%npar) pp(ik)%epar(i,j) = eigvec_(i,pp(ik)%ipar(j))
  forall(i=1:nstate,j=1:pp(ik)%nemp) pp(ik)%eemp(i,j) = eigvec_(i,pp(ik)%iemp(j))

  end subroutine init_eigvec

! ----------------------------------------------------------------------
  subroutine load_eigvec( eigvec_, eigval_, ik )
! ----------------------------------------------------------------------
  ! load eigenvectors and distribute
  complex(dp),intent(in) :: eigvec_(:,:)
  real(dp),intent(in) :: eigval_(:)
  integer,intent(in) :: ik

  if( ik < 1 .or. ik > 2 ) &
    call errore('load_eigvec','number of occupancies not 1 or 2',1)

  call init_eigvec( eigvec_, ik )

  if( associated( pp(ik)%eigval ) ) deallocate( pp(ik)%eigval )
  allocate( pp(ik)%eigval(nstate) )
  pp(ik)%eigval = eigval_

  end subroutine load_eigvec

! ----------------------------------------------------------------------
  subroutine init_triple( nk )
! ----------------------------------------------------------------------
  ! initialize storage space for the partition of the triple
  integer,intent(in) :: nk
  integer :: ik, jk

  if( nk==1 ) then
    ! only one set of dimensions
    ik=1 ; jk=1
  elseif( nk==2 ) then
    ik=1 ; jk=2
  endif

  if( allocated( triple_op ) ) deallocate( triple_op )
  if( allocated( triple_oe ) ) deallocate( triple_oe )
  if( allocated( triple_po ) ) deallocate( triple_po )
  if( allocated( triple_pp ) ) deallocate( triple_pp )
  if( allocated( triple_pe ) ) deallocate( triple_pe )
  if( allocated( triple_eo ) ) deallocate( triple_eo )
  if( allocated( triple_ep ) ) deallocate( triple_ep )
  allocate( triple_op( pp(ik)%nocc, pp(jk)%npar ), &
            triple_oe( pp(ik)%nocc, pp(jk)%nemp ), &
            triple_po( pp(ik)%npar, pp(jk)%nocc ), &
            triple_pp( pp(ik)%npar, pp(jk)%npar ), &
            triple_pe( pp(ik)%npar, pp(jk)%nemp ), &
            triple_eo( pp(ik)%nemp, pp(jk)%nocc ), &
            triple_ep( pp(ik)%nemp, pp(jk)%npar ) )

  end subroutine init_triple

! ----------------------------------------------------------------------
  subroutine prod_triple( triple, nk, xvec )
! ----------------------------------------------------------------------
  ! generate the product: eigvec* triple eigvec
  complex(dp),intent(in) :: triple(nstate,nstate)
  integer,intent(in) :: nk
  complex(dp),intent(out) :: xvec(nnz)
  integer :: ik, jk
  integer :: inz, n, nf, nt, m, mf, mt
  complex(dp) :: ztmp(nstate,nstate)

  if( nk==1 ) then
    ! only one set of dimensions
    ik=1 ; jk=1
  elseif( nk==2 ) then
    ik=1 ; jk=2
  endif

  call matmul_CNN( pp(ik)%nocc, pp(jk)%npar, nstate, &
                   pp(ik)%eocc, triple, pp(jk)%epar, ztmp, triple_op )
  call matmul_CNN( pp(ik)%nocc, pp(jk)%nemp, nstate, &
                   pp(ik)%eocc, triple, pp(jk)%eemp, ztmp, triple_oe )
  call matmul_CNN( pp(ik)%npar, pp(jk)%nocc, nstate, &
                   pp(ik)%epar, triple, pp(jk)%eocc, ztmp, triple_po )
  call matmul_CNN( pp(ik)%npar, pp(jk)%npar, nstate, &
                   pp(ik)%epar, triple, pp(jk)%epar, ztmp, triple_pp )
  call matmul_CNN( pp(ik)%npar, pp(jk)%nemp, nstate, &
                   pp(ik)%epar, triple, pp(jk)%eemp, ztmp, triple_pe )
  call matmul_CNN( pp(ik)%nemp, pp(jk)%nocc, nstate, &
                   pp(ik)%eemp, triple, pp(jk)%eocc, ztmp, triple_eo )
  call matmul_CNN( pp(ik)%nemp, pp(jk)%npar, nstate, &
                   pp(ik)%eemp, triple, pp(jk)%epar, ztmp, triple_ep )

  do inz=1,nnz
    n = index_nnz(1,inz)
    nf = pp(ik)%iope(1,n)
    nt = pp(ik)%iope(2,n)

    m = index_nnz(2,inz)
    mf = pp(jk)%iope(1,m)
    mt = pp(jk)%iope(2,m)

    if( nf==1 .and. mf==2 ) then
      xvec(inz) = triple_op(nt,mt)
    elseif( nf==1 .and. mf==3 ) then
      xvec(inz) = triple_oe(nt,mt)
    elseif( nf==2 .and. mf==1 ) then
      xvec(inz) = triple_po(nt,mt)
    elseif( nf==2 .and. mf==2 ) then
      xvec(inz) = triple_pp(nt,mt)
    elseif( nf==2 .and. mf==3 ) then
      xvec(inz) = triple_pe(nt,mt)
    elseif( nf==3 .and. mf==1 ) then
      xvec(inz) = triple_eo(nt,mt)
    elseif( nf==3 .and. mf==2 ) then
      xvec(inz) = triple_ep(nt,mt)
    endif
  enddo

  end subroutine prod_triple

! ----------------------------------------------------------------------
  subroutine matmul_CNN( n, m, l, a, b, c, w, abc )
! ----------------------------------------------------------------------

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

  end module paulipairs
