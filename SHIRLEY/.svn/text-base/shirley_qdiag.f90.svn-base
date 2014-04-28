  program shirley_qdiag

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Jan 2007

  use hamq_shirley

  implicit none

  integer,parameter :: stdout=6
  integer,parameter :: stdin =5
  integer,parameter :: dp = kind(1.d0)
  integer,parameter :: maxchar=255
  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  character(len=3)   :: nodenumber

  character(maxchar) :: hamqfile
  character(maxchar) :: outfile
  character(maxchar) :: vnlfile
  character(maxchar) :: kpt_type
  character(maxchar) :: fmtstr
  integer :: nq, iq
  integer :: nqgrid(3), iqgrid(3)
  real(dp),allocatable :: qvec(:,:), wq(:)
  logical :: cartesian
  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: eigvec(:,:)
  complex(dp),allocatable :: betaq(:,:)
  integer :: iunhq, iunout, iunvnl, ierr

  integer :: nksp, kres
  real(dp),allocatable :: xksp(:,:), qpathlensp(:)
  character(255),allocatable :: labelsp(:)
  real(dp) :: lambda

  integer :: nblock, lwork
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  real(dp) :: qpathlen, cqvec(3), dqvec(3)
  character(255) :: ic

  type(matrix_list),allocatable :: vnl_atom(:)
  integer :: iatom, it, np, npm
  complex(dp),allocatable :: beta_block(:,:), vnl(:,:), ztmp(:,:)

  integer,external :: freeunit
  integer,external :: ilaenv

  integer :: i,j,k

  namelist / input / hamqfile, outfile, vnlfile


  write(stdout,*) ' shirley_qdiag'
  write(stdout,*)

  ! initialize mpi
  CALL start_shirley (nodenumber) 

  read(stdin,nml=input,iostat=ierr)
  if( ierr /= 0 ) stop 101

  read(stdin,*) kpt_type ! type of coordinates
  kpt_type = adjustl(kpt_type)

  cartesian = .false.
  if( trim(kpt_type) == 'tpiba' .or. trim(kpt_type) == 'crystal' ) then
    read(stdin,*) nq
    allocate( qvec(3,nq), wq(nq) )
    do iq=1,nq
      read(stdin,*) qvec(1:3,iq), wq(iq)
    enddo
    if( trim(kpt_type) == 'tpiba' ) cartesian = .true.
  else if( trim(kpt_type) == 'automatic' ) then
    read(stdin,*) nqgrid, iqgrid
    if( any(iqgrid > 1) .or. any(iqgrid < 0) ) stop 201
    if( any(nqgrid < 1) ) stop 202
    nq = product(nqgrid)
    allocate( qvec(3,nq), wq(nq) )
    nq=0
    do i=1,nqgrid(1)
      do j=1,nqgrid(2)
        do k=1,nqgrid(3)
          nq=nq+1
          !  xk are the components of the complete grid in crystal axis
          qvec(1,nq) = DBLE(i-1)/nqgrid(1) - DBLE(iqgrid(1))/2/nqgrid(1)
          qvec(2,nq) = DBLE(j-1)/nqgrid(2) - DBLE(iqgrid(2))/2/nqgrid(2)
          qvec(3,nq) = DBLE(k-1)/nqgrid(3) - DBLE(iqgrid(3))/2/nqgrid(3)
        enddo
      enddo
    enddo
    wq = 2.d0/dble(nq)
  elseif( trim(kpt_type) == 'bandstructure' ) then
    read(stdin,*) nksp, kres
    allocate( xksp(3,nksp), labelsp(nksp), qpathlensp(nksp) )
    do i=1,nksp
      read(stdin,*) xksp(:,i), labelsp(i)
    enddo
    nq=1
    do i=2,nksp
      nq=nq+kres
    enddo
    allocate( qvec(3,nq), wq(nq) )
    nq=1
    qvec(:,nq) = xksp(:,1)
    do i=2,nksp
      do j=1,kres
        lambda=dble(j)/dble(kres)
        nq=nq+1
        qvec(:,nq) = (1.d0-lambda)*xksp(:,i-1) + lambda*xksp(:,i)
      enddo
    enddo
    wq = 1.d0
  else
    write(stdout,*) ' kpoints flag unrecognized'
    write(stdout,*) ' should be: tpiba, crystal, or automatic'
    stop 203
  endif

  iunhq = freeunit()
  open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
  if( ierr /= 0 ) stop 102

  iunout = freeunit()
  open(iunout,file=trim(outfile),form='formatted',iostat=ierr)
  if( ierr /= 0 ) stop 103

  iunvnl = freeunit()
  open(iunvnl,file=trim(vnlfile),form='unformatted',iostat=ierr)
  if( ierr /= 0 ) stop 104

  call read_hamq( iunhq )

  ! read atomic matrix elements for non-local potential
  allocate( vnl_atom(natom) )
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type(it)
    write(*,*) iatom, it, np
    allocate( vnl_atom(iatom)%matrix(np,np) )
    !read(iunvnl) vnl_atom(iatom)%matrix
  enddo
  call read_nloper( iunvnl, vnl_atom )
  write(*,*) ' nproj_nl = ', nproj_nl
  do iatom=1,natom
    write(*,*) iatom, vnl_atom(iatom)%matrix
  enddo
  do iatom=1,natom
    it=type_atom(iatom)
    write(*,*) ' index_betaq = ', iatom, it, index_betaq(1:nproj_type(it),iatom)
  enddo

  allocate( eigval(nbasis), eigvec(nbasis,nbasis) )
  allocate( betaq(nproj_nl,nbasis) )
  write(fmtstr,'(a,i,a)') '(i,3e,',nbasis,'e)'

  ! tmp space for diagonalization
  nblock = ILAENV( 1, 'ZHETRD', 'U', nbasis, - 1, - 1, - 1 )
  IF ( nblock < 1 ) nblock = MAX( 1, nbasis )
  IF ( nblock == 1 .OR. nblock >= nbasis ) THEN
     lwork = 2 * nbasis - 1
  ELSE
     lwork = ( nblock + 1 ) * nbasis
  END IF
  allocate( work(lwork), rwork(3*nbasis-2) )

  qpathlen=0.d0
  do iq=1,nq
    dqvec = cqvec
    cqvec = matmul( trnlp2kin, qvec(1:3,iq) )
    dqvec = cqvec - dqvec 
    if( iq==1 ) dqvec = 0.d0
    qpathlen = qpathlen + sqrt(dot_product( dqvec, dqvec ))
  enddo
  do iq=1,nq
    ! build the Hamiltonian for this q-point
    ! local contribution
    call build_hamq_local( qvec(1:3,iq), cartesian, eigvec )

    ! non-local contribution
    if( nproj_nl > 0 ) then
      call build_hamq_nlprojs( qvec(1:3,iq), cartesian, betaq )

      npm = maxval( nproj_type_nl(1:ntype) )
      allocate( beta_block(npm,nbasis), vnl(npm,npm), ztmp(npm,nbasis) )
      do iatom=1,natom
        it = type_atom(iatom)
        np = nproj_type_nl(it)
        write(*,*) '     type_atom = ', it
        write(*,*) ' nproj_type_nl = ', np

        forall( i=1:np, j=1:nbasis ) &
          beta_block(i,j) = &
            betaq(index_betaq(index_nlproj_type(i,it),iatom),j)

        forall( i=1:np, j=1:np ) &
          vnl(i,j) = cmplx( vnl_atom(iatom)%matrix( index_nlproj_type(i,it), &
                                                    index_nlproj_type(j,it) ) )

        CALL ZGEMM( 'N', 'N', np, nbasis, np, one, &
                    vnl, npm, &
                    beta_block, npm, &
                    zero, ztmp, npm )
        CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, one, &
                    beta_block, npm, &
                    ztmp, npm, &
                    one, eigvec, nbasis )
      enddo

      deallocate( beta_block, vnl, ztmp )
    endif

    CALL ZHEEV( 'V', 'U', nbasis, eigvec, nbasis, eigval, work, lwork, rwork, ierr )
    if( ierr /= 0 ) stop 105

    write(stdout,*) iq, eigval

    if( trim(kpt_type) == 'bandstructure' ) then
      dqvec = cqvec
      cqvec = matmul( trnlp2kin, qvec(1:3,iq) )
      dqvec = cqvec - dqvec 
      if( iq==1 ) dqvec = 0.d0
      qpathlen = qpathlen + sqrt(dot_product( dqvec, dqvec ))
      write(iunout,fmtstr) iq-1, matmul( trnlp2kin, qvec(1:3,iq)), qpathlen, eigval*rytoev
      if( mod(iq-1,kres)==0 ) qpathlensp((iq-1)/kres+1) = qpathlen
    else
      write(iunout,fmtstr) iq-1, matmul( trnlp2kin, qvec(1:3,iq)), eigval*rytoev
    endif

  enddo

  write(401,*) "set nokey"
  write(401,*) "set style data l"
  write(401,*) "bandcolor=-1"
  write(401,*) "efermi=0.0"
!  write(401,*) "set term postscript eps enhanced"
!  write(401,*) "set output 'qdiag.eps'"
  
  if( trim(kpt_type) == 'bandstructure' ) then
  write(401,*) "set xtics ( \"
  do i=1,nksp-1
    write(401,'(a,f,a)') '"'//trim(labelsp(i))//'"', qpathlensp(i), ", \"
  enddo
  write(401,'(a,f,a)') '"'//trim(labelsp(nksp))//'"', qpathlensp(nksp), ")"
  do i=2,nksp-1
    write(401,'(a,i,a,f,a,f,a)') "set arrow ", i, " from first ", qpathlensp(i), ", graph 0 to ", qpathlensp(i), ", graph 1 lt bandcolor lw 1 nohead" 
  enddo
  endif

  write(401,'(a, a, a)') "plot '", trim(outfile), "' u 5:($6-efermi) lt bandcolor \"
  do i=2,nbasis-1
    write(ic,'(i)') i+5 ; ic=adjustl(ic)
    write(401,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor \"
  enddo
  write(ic,'(i)') nbasis+5 ; ic=adjustl(ic)
  write(401,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor"

!  write(401,*) "set nokey"
!  write(401,'(a, a, a)',advance='no') "plot '", trim(outfile), "' u 1:5 lt 1"
!  do i=6,nbasis+4
!    if( i > 0   .and. i < 10 )   write(401,'(a,i1,a)',advance='no') ", '' u 1:", i, " lt 1"
!    if( i > 10  .and. i < 100 )  write(401,'(a,i2,a)',advance='no') ", '' u 1:", i, " lt 1" 
!    if( i > 100 .and. i < 1000 ) write(401,'(a,i3,a)',advance='no') ", '' u 1:", i, " lt 1"
!  enddo
!  write(401,*)
  
  write(stdout,*) ' end shirley_qdiag'
  
  end program shirley_qdiag
