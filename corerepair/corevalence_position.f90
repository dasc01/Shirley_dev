  program corevalence_position

! Read the pseudo-waves
! Ask for corresponding ae wave
! Generate a position matrix for each pair of channels
!   for ae and ps
! Subtract to give corerepair term
!   <ae_l|p|ae_m> - <ps_l|p|ps_m>

  use kinds
  use fileio
  use atomic_waves
  use splines

  implicit none

  integer,external :: iargc

  real(dp),parameter :: zero  = 0.d0
  real(dp),parameter :: eps12 = 1.d-12

  integer :: narg, iunc, iunv
  character(255) :: filename, fmtstr
  
  type(atomic) :: core_wave, tmp_wave
  type(atomic),allocatable :: vlnc_wave(:)
  real(dp) :: rmax
  real(dp),allocatable :: posmat(:,:,:)

  integer,parameter :: ideg=5
  type(spline_struct) :: core_spl

  integer :: ngrid, l, igrid
  integer :: ixyz
  integer :: immax, jmmax, im, jm
  integer :: nmax, ierr, nnonzero
  real(dp) :: core_norm, crm

  complex(dp),allocatable :: corerep(:,:,:)
  integer :: iproj, jproj, nproj1, nproj2, i0, j0
  integer :: ivalence, nvalence


  ! Get user input
  narg = iargc()
  if( narg /= 2 ) then
    write(0,*) 'usage: corevalence_position corefile valencefile'
    stop
  endif

  call getarg( 1, filename )
  ! open file
  iunc = freeunit()
  open(iunc,file=trim(filename),form='formatted',err=911)

  call getarg( 2, filename )
  ! open file
  iunv = freeunit()
  open(iunv,file=trim(filename),form='formatted',err=911)

  ! Read the core file
#ifdef __DEBUG
  write(*,*) 'reading core file'
#endif
  read(iunc,*) 
  read(iunc,*) ngrid, l
  call init_atomic_wave( ngrid, l, core_wave )
  do igrid=1,ngrid
    read(iunc,*) core_wave%r(igrid), core_wave%f(igrid)
  enddo
  close(iunc)

  ! find extent of core
  do igrid=ngrid,1,-1
    if( abs(core_wave%f(igrid)) > eps12 ) exit
  enddo
  nmax = igrid
  rmax = core_wave%r(nmax)
#ifdef __DEBUG
  write(*,*) 'extent of core wave : rmax = ', rmax
#endif
  if( nmax==ngrid .and. abs(core_wave%f(nmax)) > 1e-6 ) then
    write(*,*) 'warning : core extends to limits of radial grid'
    write(*,*) '             rmax = ', rmax
    write(*,*) '          f(rmax) = ', core_wave%f(igrid)
    write(*,*) '          this may lead to inaccurate integrals'
  endif

  ! normalize core within this radius
#ifdef __DEBUG
  write(*,*) 'normalizing core'
#endif
  ngrid=core_wave%ngrid
  if( core_wave%r(1) > eps12 ) then
    ngrid=ngrid+1
    call init_atomic_wave( ngrid, core_wave%l, tmp_wave )
    tmp_wave%r = (/ zero, core_wave%r /)
    tmp_wave%f = (/ zero, core_wave%f**2 /)
  else
    call init_atomic_wave( ngrid, core_wave%l, tmp_wave )
    tmp_wave%r = core_wave%r
    tmp_wave%f = core_wave%f**2
  endif
  call spline_fit( ideg, tmp_wave%ngrid, tmp_wave%r, tmp_wave%f, core_spl, ierr )
  if ( ierr>0 ) then
    write(*,*) ' spline error ', ierr
    stop
  endif
  core_norm =  spline_integral( core_spl, zero, tmp_wave%r(tmp_wave%ngrid) )
  core_wave%f = core_wave%f / core_norm
#ifdef __DEBUG
  write(*,*) 'core_norm = ', core_norm
#endif

  ! Read the valence wave
#ifdef __DEBUG
  write(*,*) 'reading valence wave'
#endif
  read(iunv,*) nvalence

  write(*,'(2x,a)') 'position'
  write(*,'(2x,2i4,t24,a)') nvalence, 1, '! nwfc1, nwfc2'

  allocate( vlnc_wave(nvalence) )

  do ivalence=1,nvalence

  ! read the pseudo wave first
  read(iunv,*) ngrid, l
  do igrid=1,ngrid
    read(iunv,*)
  enddo

  ! read associated all-electron wave and keep
  read(iunv,*) ngrid, l
  call init_atomic_wave( ngrid, l, vlnc_wave(ivalence) )
  do igrid=1,ngrid
    read(iunv,*) vlnc_wave(ivalence)%r(igrid), vlnc_wave(ivalence)%f(igrid)
  enddo

  enddo ! ivalence

  write(fmtstr,*) '(2x,', nvalence, 'i2,t24,a)'
  write(*,fmtstr)       vlnc_wave(:)%l, '! lwfc1(1:nwfc1)'
  write(*,'(2x,i2,t24,a)') core_wave%l, '! lwfc2(1:nwfc2)'

  nproj1=0
  do ivalence=1,nvalence
    nproj1=nproj1+2*vlnc_wave(ivalence)%l+1
  enddo
  nproj2=2*core_wave%l+1
  allocate( corerep(3,nproj1,nproj2) )
  corerep = 0.d0


  j0=0
  do ivalence=1,nvalence

  ! check that valence extent large enough
  if( vlnc_wave(ivalence)%r(vlnc_wave(ivalence)%ngrid) < rmax ) then
    write(*,*) 'warning : not enough info'
    write(*,*) '          extent of valence wave : rmax = ', &
      vlnc_wave(ivalence)%r(vlnc_wave(ivalence)%ngrid)
    write(*,*) '          extent of    core wave : rmax = ', rmax
    write(*,*) 'warning : integrals will be inaccurate'
    write(*,*) '          using smaller rmax'
    rmax = vlnc_wave(ivalence)%r(vlnc_wave(ivalence)%ngrid)
  endif

  ! Compute position matrices
#ifdef __DEBUG
  write(*,*) 'generating position matrices'
#endif
  immax = 2 * core_wave%l + 1
  jmmax = 2 * vlnc_wave(ivalence)%l + 1
  ! write(*,*) vlnc_wave(ivalence)%l, core_wave%l, jmmax, immax
  allocate( posmat(3,jmmax,immax) )

  call position_matrix( &
    vlnc_wave(ivalence)%ngrid, vlnc_wave(ivalence)%l, vlnc_wave(ivalence)%r, vlnc_wave(ivalence)%f, &
    core_wave%ngrid, core_wave%l, core_wave%r, core_wave%f, &
    rmax, posmat )

  do im=1,immax
    do jm=1,jmmax
      do ixyz=1,3
        corerep(ixyz,j0+jm,im) = cmplx(posmat(ixyz,jm,im),0.d0)
      enddo
    enddo
  enddo

  deallocate( posmat )

  j0=j0+jmmax
  enddo ! ivalence

  close(iunv)

  ! determine number of nonzero elements
  nnonzero=0
  do iproj=1,nproj2
    do jproj=1,nproj1
      do ixyz=1,3
        crm = conjg(corerep(ixyz,jproj,iproj))*corerep(ixyz,jproj,iproj)
        if( crm > eps12 ) then
          nnonzero=nnonzero+1
        endif
      enddo
    enddo
  enddo
    
  write(*,'(2x,i,t24,a)') nnonzero, '! nonzero elements (i,j,ixyz,cR,cI)'

  do iproj=1,nproj2
    do jproj=1,nproj1
      do ixyz=1,3
        crm = conjg(corerep(ixyz,jproj,iproj))*corerep(ixyz,jproj,iproj)
        if( crm > eps12 ) then
          write(*,'(2x,3i3,x,2e)') jproj, iproj, ixyz, corerep(ixyz,jproj,iproj)
        endif
      enddo
    enddo
  enddo


  !do im=1,immax
  !do jm=1,jmmax
  !  write(*,'(x,4i4,x,3e20.10)') vlnc_wave(ivalence)%l, core_wave%l, jm, im, &
  !    (posmat(ixyz,jm,im), ixyz=1,3)
  !enddo
  !enddo

  stop

  911 write(0,*) 'corevalence_position : unable to open file ', trim(filename)
  stop

  end program corevalence_position
