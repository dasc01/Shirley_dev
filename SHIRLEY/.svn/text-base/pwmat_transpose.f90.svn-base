  program pwmat_transpose

  use kinds, only : dp
  use mp,    only : mp_start, mp_end, mp_env, mp_barrier

  implicit none

  integer,parameter :: maxchar=255

  integer :: nproc, mpime, gid, root
  logical :: amroot 

  integer :: stdin, stdout
  integer :: narg
  integer,external :: iargc

  character(maxchar) :: cecut, prefix
  real(dp) :: ecut

  integer :: iunpwi, iunpwm, iunpwt
  integer,external :: freeunit
  integer :: ierr

  ! pwi stuff
  integer :: igwx, nbnd
  real(dp) :: at(3,3), bg(3,3), tpiba
  real(dp),allocatable :: gvec(:,:)

  namelist / pwilist / igwx, nbnd, at, bg, tpiba

  real(dp) :: tpiba2
  real(dp),allocatable :: gvec2(:)
  integer :: igwx_trunc

  integer :: i, ixyz


  call mp_start()
  call mp_env( nproc, mpime, gid )

  root = 0
  amroot=.false.
  if( mpime==root ) amroot=.true.

  
  stdout=0
  stdin=0
  if( amroot ) then
    stdout=6
    stdin=5
  endif

  if( amroot ) then
    ! get arguments
    narg=iargc()
    if( narg /= 2 ) then
      write(stdout,*) ' usage: pwmat_transpose ecut prefix'
      goto 999
    endif

    call getarg( 1, cecut )
    call getarg( 2, prefix )

    read(cecut,*) ecut
    write(stdout,*) ' specified cut-off = ', ecut
    write(stdout,*) ' looking for files: '
    write(stdout,*) trim(prefix)//'.pwi'
    write(stdout,*) trim(prefix)//'.pwm'

    ! pwi - information on matrix elements
    iunpwi=freeunit()
    open(iunpwi,file=trim(prefix)//'.pwi',form='unformatted',status='old',iostat=ierr)
    if( ierr/=0 ) then
      write(stdout,*) ' problem opening ', trim(prefix)//'.pwi'
      goto 999
    endif

    ! pwm - matrix elements
    iunpwm=freeunit()
    open(iunpwm,file=trim(prefix)//'.pwm',form='unformatted',status='old',iostat=ierr)
    if( ierr/=0 ) then
      write(stdout,*) ' problem opening ', trim(prefix)//'.pwm'
      goto 999
    endif

    ! pwt - transposed matrix elements
    iunpwt=freeunit()
    open(iunpwt,file=trim(prefix)//'.pwt',form='unformatted',iostat=ierr)
    if( ierr/=0 ) then
      write(stdout,*) ' problem opening ', trim(prefix)//'.pwt'
      goto 999
    endif

    write(stdout,*) ' dumping transpose to file: '
    write(stdout,*) trim(prefix)//'.pwt'

  endif
    
  ! get info
  if( amroot ) then

    read(iunpwi,err=901) igwx  ! number of G vectors
    read(iunpwi,err=901) nbnd  ! number of bands
    read(iunpwi,err=901) at    ! transformation matrix of bravais lattice
    read(iunpwi,err=901) bg    ! transformation matrix of reciprocal lattice
    read(iunpwi,err=901) tpiba ! 2 * pi / alat

    allocate( gvec(igwx,3) )

    do ixyz=1,3
      read(iunpwi,err=901) gvec(1:igwx,ixyz)
    enddo

    close(iunpwi)

    write(stdout,nml=pwilist)

  endif

  ! reduce number of gvectors to inside cut-off
  if( amroot ) then
    allocate( gvec2(igwx) )
    tpiba2 = tpiba**2
    do i=1,igwx
      gvec2(i) = sum( gvec(i,1:3)**2 ) * tpiba2
    enddo
    ! check ordering
    do i=1,igwx-1
      if( gvec2(i+1) - gvec2(i) > 1.d-4 ) then
        write(stdout,*) ' error in G vector ordering'
        write(stdout,'(a,i,a,e,a,i,a,e)') &
          ' gvec2(', i+1, ')=', gvec2(i+1), &
          ' < gvec2(', i, ')=', gvec2(i)
      endif
    enddo
    ! cut-off
    do i=1,igwx
      if( gvec2(i) > ecut ) exit
    enddo
    igwx_trunc = i-1
    write(stdout,*) ' truncating number of G vectors to ', igwx_trunc

  endif

  call mp_barrier()



  999 call mp_end()
  stop

  901 write(stdout,*) ' error reading ', trim(prefix)//'.pwi'
  goto 999

  end program pwmat_transpose
