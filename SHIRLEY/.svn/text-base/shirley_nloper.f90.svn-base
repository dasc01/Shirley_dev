  program shirley_nloper

  ! use this program to generate binary files containing
  ! information on the atomic matrix kernel of non-local operators
  ! these might include:
  !   1. the non-local potential
  !   2. core-repair matrix elements for momentum
  !   3. core-repair matrix elements for position

  ! David Prendergast, UCB, May 2007

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
  character(maxchar) :: vnlfile
  character(maxchar) :: nloperfile
  character(maxchar) :: heading

  integer :: iunhq, iunin, iunout, ierr

  type(matrix_list),allocatable :: nloper_atom(:)
  type(matrix_list),allocatable :: nloper_type(:)

  integer :: iatom, it, np, npm

  integer,external :: freeunit

  integer :: i,j,k

  namelist / input / hamqfile, nloper_inp, nloper_out


  write(stdout,*) ' shirley_nloper'
  write(stdout,*)

  ! initialize mpi - if necessary
  CALL start_shirley (nodenumber) 

  read(stdin,nml=input,iostat=ierr)
  if( ierr /= 0 ) stop 101

  iunhq = freeunit()
  open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
  if( ierr /= 0 ) stop 102

  iuninp = freeunit()
  open(iuninp,file=trim(nloper_inp),form='formatted',iostat=ierr)
  if( ierr /= 0 ) stop 103

  iunout = freeunit()
  open(iunout,file=trim(nloper_out),form='unformatted',iostat=ierr)
  if( ierr /= 0 ) stop 104

  ! read the Hamiltonian (for dimensions mostly)
  call read_hamq( iunhq )

  ! allocate space for nloper
  write(stdout,*) ' allocate space for non-local operator'
  allocate( nloper_atom(natom) )
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type(it)
    write(stdout,*) iatom, it, np
    allocate( nloper_atom(iatom)%matrix(np,np) )
  enddo

  ! read atomic matrix elements for non-local potential
  ! assuming that information does not vary within atoms of a given type
  ! (note that this may not be true for non-norm-conserving potentials)
  do it=1,ntype
    ! read matrix for this type
    read(iuninp,*) heading
    write(stdout,*) trim(heading)
    read(iuninp,*) np1, np2
    
  enddo

  ! transfer to atom-based form
  do it=1,ntype
    do iatom=1,natom
      if( type_atom(iatom) == it ) then
        nloper_atom(iatom)%matrix = nloper_type(it)%matrix
      endif
    enddo
  enddo

  ! now dump
  call write_nloper( iunout, nloper_atom )


  write(stdout,*) ' end shirley_nloper'
  
  end program shirley_nloper
