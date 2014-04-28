  module diag_shirley

  use kinds, only : dp
  use hamq_shirley, only : nbasis, build_hamq
  use io_global, only : stdout

  implicit none

  private
  public :: diag_hamq

  interface diag_hamq
    module procedure diag_hamq_eigvals_only
    module procedure diag_hamq_eigvec
  end interface diag_hamq

  integer :: nbasis_current=0
  integer :: lwork
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  complex(dp),allocatable :: a(:,:)
  integer :: ierr

  contains

! ---------------------------------------------------------------------- 
  subroutine diag_hamq_eigvals_only( qvec, eigval, lkin )
! ---------------------------------------------------------------------- 

  ! diagonalize the hamiltonian for a given input q-vector
  ! this vector should be provided in both crystal and
  ! Cartesian coordinates including a factor of 2*pi/a
  ! It must be contained within the region of k-space for
  ! which the non-local potential is splined or an error is produced.
  ! It should be contained within the k-space region for over
  ! which the basis is chosen - the first Brillouin zone
  ! This Hamiltonian is not guaranteed to be periodic wrt q.

  real(dp),intent(in) :: qvec(3)
  real(dp),intent(out) :: eigval(nbasis)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis

  character,parameter :: EIGVALS_ONLY='N'

  allocate( a(nbasis,nbasis) )

  call  build_hamq( qvec, lkin, a )

  ! tmp space for diagonalization
  if( nbasis /= nbasis_current ) then

    call init_zheev( nbasis )

    nbasis_current = nbasis

  endif

  CALL ZHEEV( EIGVALS_ONLY, 'U', nbasis, a, nbasis, eigval, work, lwork, rwork, ierr )
  if( ierr /= 0 ) call errore('ZHEEV','problem diagonalizing',abs(ierr))

  deallocate( a )

  return
  end subroutine diag_hamq_eigvals_only

! ---------------------------------------------------------------------- 
  subroutine diag_hamq_eigvec( qvec, eigval, eigvec, lkin )
! ---------------------------------------------------------------------- 

  ! diagonalize the hamiltonian for a given input q-vector
  ! this vector should be provided in both crystal and
  ! Cartesian coordinates including a factor of 2*pi/a
  ! It must be contained within the region of k-space for
  ! which the non-local potential is splined or an error is produced.
  ! It should be contained within the k-space region for over
  ! which the basis is chosen - the first Brillouin zone
  ! This Hamiltonian is not guaranteed to be periodic wrt q.

  real(dp),intent(in) :: qvec(3)
  real(dp),intent(out) :: eigval(nbasis)
  complex(dp),intent(out) :: eigvec(nbasis,nbasis)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis

  character,parameter :: EIGVECS='V'


  call  build_hamq( qvec, lkin, eigvec )

  ! tmp space for diagonalization
  if( nbasis /= nbasis_current ) then

    call init_zheev( nbasis )

    nbasis_current = nbasis

  endif

  CALL ZHEEV( EIGVECS, 'U', nbasis, eigvec, nbasis, eigval, work, lwork, rwork, ierr )
  if( ierr /= 0 ) call errore('ZHEEV','problem diagonalizing',abs(ierr))

  return
  end subroutine diag_hamq_eigvec


  subroutine init_zheev( n )

  integer,intent(in) :: n
  integer :: nblock
  integer,external :: ilaenv

  if( allocated( rwork ) ) deallocate( rwork )
  if( allocated( work ) ) deallocate( work )

  nblock = ILAENV( 1, 'ZHETRD', 'U', n, - 1, - 1, - 1 )
  IF ( nblock < 1 ) nblock = MAX( 1, n )
  IF ( nblock == 1 .OR. nblock >= n ) THEN
     lwork = 2 * n - 1
  ELSE
     lwork = ( nblock + 1 ) * n
  END IF

  write(stdout,*) ' init_zheev allocating ', (dble(lwork)*2.d0 + dble(max(1,3*n-2)))*8.d0/2.d0**20.d0

  allocate( work(lwork), rwork(max(1,3*n-2)), stat=ierr )
  if( ierr/= 0 ) call errore('init_zheev','unable to allocate workspace',abs(ierr))
  end subroutine init_zheev

  end module diag_shirley
