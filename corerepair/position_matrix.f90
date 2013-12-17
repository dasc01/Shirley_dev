  subroutine position_matrix( n_f_, l_f, r_f_, f_, &
                              n_g_, l_g, r_g_, g_, &
                              rmax, posmat )

! Compute the position matrix elements between two atomic waves of
! particular angular momentum symmetry
! Makes use of some pre-evaluated integrals from Mathematica

  use kinds
  use splines

  implicit none

  integer,parameter :: ideg = 5
  real(dp),parameter :: zero = 0.d0
  real(dp),parameter :: eps12 = 1.d-12
  integer,parameter :: lmax = 3

  integer :: n_f, l_f
  integer :: n_f_
  real(dp) :: r_f_(n_f_), f_(n_f_)
  real(dp) :: r_f(n_f_+1), f(n_f_+1)
  integer :: n_g, l_g
  integer :: n_g_
  real(dp) :: r_g_(n_g_), g_(n_g_)
  real(dp) :: r_g(n_g_+1), g(n_g_+1)
  real(dp) :: rmax
  real(dp) :: posmat(3,2*l_f+1,2*l_g+1)

  integer :: i
  integer :: n_fmax
  integer :: ierr
  type(spline_struct) :: spl
  real(dp),allocatable :: gbr(:)
  real(dp) :: fgbr_int
  integer :: m_g, m_f, ixyz
  integer :: inonzero

  ! Angular integrals
  real(dp),allocatable :: Imat(:,:,:,:,:)


  ! set up angular integrals
  allocate( Imat(3,2*lmax+1,0:lmax,2*lmax+1,0:lmax) )
  Imat = zero
                                                                                
  Imat(3,1,0,1,1) =      0.5773502692
  Imat(1,1,0,2,1) =     -0.5773502692
  Imat(2,1,0,3,1) =     -0.5773502692
  Imat(3,1,1,1,0) =      0.5773502692
  Imat(3,1,1,1,2) =      0.5163977795
  Imat(1,1,1,2,2) =     -0.4472135955
  Imat(2,1,1,3,2) =     -0.4472135955
  Imat(1,2,1,1,0) =     -0.5773502692
  Imat(1,2,1,1,2) =      0.2581988897
  Imat(3,2,1,2,2) =      0.4472135955
  Imat(1,2,1,4,2) =     -0.4472135955
  Imat(2,2,1,5,2) =     -0.4472135955
  Imat(2,3,1,1,0) =     -0.5773502692
  Imat(2,3,1,1,2) =      0.2581988897
  Imat(3,3,1,3,2) =      0.4472135955
  Imat(2,3,1,4,2) =      0.4472135955
  Imat(1,3,1,5,2) =     -0.4472135955
  Imat(3,1,2,1,1) =      0.5163977795
  Imat(1,1,2,2,1) =      0.2581988897
  Imat(2,1,2,3,1) =      0.2581988897
  Imat(3,1,2,1,3) =      0.5070925528
  Imat(1,1,2,2,3) =     -0.4140393356
  Imat(2,1,2,3,3) =     -0.4140393356
  Imat(1,2,2,1,1) =     -0.4472135955
  Imat(3,2,2,2,1) =      0.4472135955
  Imat(1,2,2,1,3) =      0.2927700219
  Imat(3,2,2,2,3) =      0.4780914437
  Imat(1,2,2,4,3) =     -0.3779644730
  Imat(2,2,2,5,3) =     -0.3779644730
  Imat(2,3,2,1,1) =     -0.4472135955
  Imat(3,3,2,3,1) =      0.4472135955
  Imat(2,3,2,1,3) =      0.2927700219
  Imat(3,3,2,3,3) =      0.4780914437
  Imat(2,3,2,4,3) =      0.3779644730
  Imat(1,3,2,5,3) =     -0.3779644730
  Imat(1,4,2,2,1) =     -0.4472135955
  Imat(2,4,2,3,1) =      0.4472135955
  Imat(1,4,2,2,3) =      0.1195228609
  Imat(2,4,2,3,3) =     -0.1195228609
  Imat(3,4,2,4,3) =      0.3779644730
  Imat(1,4,2,6,3) =     -0.4629100499
  Imat(2,4,2,7,3) =     -0.4629100499
  Imat(2,5,2,2,1) =     -0.4472135955
  Imat(1,5,2,3,1) =     -0.4472135955
  Imat(2,5,2,2,3) =      0.1195228609
  Imat(1,5,2,3,3) =      0.1195228609
  Imat(3,5,2,5,3) =      0.3779644730
  Imat(2,5,2,6,3) =      0.4629100499
  Imat(1,5,2,7,3) =     -0.4629100499
  Imat(3,1,3,1,2) =      0.5070925528
  Imat(1,1,3,2,2) =      0.2927700219
  Imat(2,1,3,3,2) =      0.2927700219
  Imat(1,2,3,1,2) =     -0.4140393356
  Imat(3,2,3,2,2) =      0.4780914437
  Imat(1,2,3,4,2) =      0.1195228609
  Imat(2,2,3,5,2) =      0.1195228609
  Imat(2,3,3,1,2) =     -0.4140393356
  Imat(3,3,3,3,2) =      0.4780914437
  Imat(2,3,3,4,2) =     -0.1195228609
  Imat(1,3,3,5,2) =      0.1195228609
  Imat(1,4,3,2,2) =     -0.3779644730
  Imat(2,4,3,3,2) =      0.3779644730
  Imat(3,4,3,4,2) =      0.3779644730
  Imat(2,5,3,2,2) =     -0.3779644730
  Imat(1,5,3,3,2) =     -0.3779644730
  Imat(3,5,3,5,2) =      0.3779644730
  Imat(1,6,3,4,2) =     -0.4629100499
  Imat(2,6,3,5,2) =      0.4629100499
  Imat(2,7,3,4,2) =     -0.4629100499
  Imat(1,7,3,5,2) =     -0.4629100499

  ! quick check to see if integrals are zero by angular symmetry
  if( all( abs(Imat(:,:,l_g,:,l_f)) < eps12 ) ) then
    posmat = zero
    return
  endif

  ! add zeros
  if( r_f_(1) > zero ) then
    r_f = (/ zero, r_f_ /)
    f = (/ zero, f_ /)
    n_f = n_f_ + 1
  else
    r_f = r_f_
    f = f_
    n_f = n_f_ 
  endif

  if( r_g_(1) > zero ) then
    r_g = (/ zero, r_g_ /)
    g = (/ zero, g_ /)
    n_g = n_g_ + 1
  else
    r_g = r_g_
    g = g_
    n_g = n_g_ 
  endif


  ! check that rmax inside limits of radial grids
  if( rmax > r_f(n_f) .or. rmax > r_g(n_g) ) then
    write(0,*) ' position_matrix : rmax too large for these radial grids'
    write(0,*) ' rmax = ', rmax
    write(0,*) ' r_f(n_f) = ', r_f(n_f)
    write(0,*) ' r_g(n_g) = ', r_g(n_g)
    stop 1
  endif

  do i=1,n_f
    if( r_f(i) > rmax ) exit
  enddo
  if( i>=n_f ) then
    n_fmax = n_f
  else
    n_fmax = i+1
  endif

  ! should I modify rmax
  !if( r_g(n_g) < r_f(n_fmax) ) then
  !  n_fmax = n_fmax - 1
  !  rmax = r_f(n_fmax)
  !endif
  ! write(*,*) 'rmax = ', rmax

  allocate( gbr(n_fmax) )

  ! spline g
  call spline_fit( ideg, n_g, r_g, g, spl, ierr, zero, rmax )
  if( ierr > 0 ) then
    write(0,*) ' spline error 1', ierr
    stop
  endif

  ! evaluate g on f grid
  ! write(*,*) 'check :', n_fmax, r_f(1), r_f(n_fmax), zero, rmax
  call spline_eval( spl, n_fmax, r_f, gbr, ierr )
  if( ierr > 0 ) then
    write(0,*) ' spline error 2', ierr
    stop
  endif

  ! compute f g/r
  where (r_f > zero )
    gbr = f * gbr * r_f
  elsewhere
    gbr = zero
  endwhere
  
  call spline_fit( ideg, n_fmax, r_f, gbr, spl, ierr , zero, r_f(n_fmax) )
  if( ierr > 0 ) then
    write(0,*) ' spline error 5', ierr
    stop
  endif

  fgbr_int = spline_integral( spl, zero, rmax )
  ! write(*,*) 'fgbr_int = ', fgbr_int


  ! write(*,*) 'angular contributions'
  do m_g=1,2*l_g+1
  do m_f=1,2*l_f+1
    do ixyz=1,3
      posmat(ixyz,m_f,m_g) = fgbr_int  * Imat(ixyz,m_f,l_f,m_g,l_g)
    enddo
    !write(*,'(2x,4i4,3e20.10)') l_f, m_f, l_g, m_g, (Imat(ixyz,m_f,l_f,m_g,l_g), ixyz=1,3)
  enddo
  enddo

  return

  end subroutine position_matrix
