  module epsilon_module

! All the variables associated with epsilon and routines that generate
! them. Used mainly by shirley_epsilon.

#include "f_defs.h"

  use kinds, only : dp
  use kpt_module
  use corerepair_module
  use hamq_shirley

  implicit none

  integer,parameter :: stdin =5
  integer,parameter :: maxchar=255

  type(kpt_type) :: kpt, qpt
  type(corerepair_type) :: corerep
  type(matrix_list),allocatable :: vnl_atom(:)

  integer :: iunout
  character(maxchar) :: tmpfile
  real(dp) :: epsilon_cutoff
  real(dp) :: efermi

  contains

  subroutine epsilon_init

  use io_global, only : stdout, ionode, ionode_id
  use mp, only : mp_bcast
  use mp_global, only : mpime
  use pwmat_module

  integer :: iunhq, iunvnl, ierr
  character(maxchar) :: hamqfile
  character(maxchar) :: vnlfile
  character(maxchar) :: pwifile
  character(maxchar) :: pwmfile
  character(maxchar) :: outfile

  character(len=3) :: nodenumber
  character(maxchar) :: outfilenode
  integer :: iatom, it, np, icore

  integer,external :: freeunit

  namelist / input / hamqfile, vnlfile, outfile, pwifile, pwmfile, &
                     epsilon_cutoff, efermi, tmpfile


! Initialize for shirley_epsilon

  ! Initialize MPI
  call start_shirley( nodenumber )

  ! Modify stdout
  call init_stdout( stdout )

  ! read name list
  if( ionode ) then
    ! default value
    epsilon_cutoff = 0.d0
    tmpfile = 'tmpfile'
    read(stdin,nml=input,iostat=ierr)
    if( ierr /= 0 ) &
      call errore('epsilon_init','problem reading namelist &input',101)

    ! read k-points
    call kpt_read( stdin, kpt )
    write(stdout,*) 'done reading k-points'

    ! read q-points
    call kpt_read( stdin, qpt )
    write(stdout,*) 'done reading q-points'

    ! read atomic matrix elements for core-valence position
    call read_corerepair( stdin, corerep )
  endif
  ! broadcast some input
  call mp_bcast( pwmfile, ionode_id )
  call mp_bcast( outfile, ionode_id )
  call mp_bcast( tmpfile, ionode_id )
  call mp_bcast( epsilon_cutoff, ionode_id )
  call mp_bcast( efermi, ionode_id )
  call kpt_bcast( kpt, ionode_id )
  call kpt_bcast( qpt, ionode_id )
  call bcast_corerepair( corerep, ionode_id )


  if( ionode ) then
  ! open file for Hamiltonian input
    iunhq = freeunit()
    open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) &
      call errore('epsilon_init','problem opening file '//trim(hamqfile),102)

  ! open file for Non-local potential input
    iunvnl = freeunit()
    open(iunvnl,file=trim(vnlfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) &
      call errore('epsilon_init','problem opening file '//trim(vnlfile),103)
  endif

  ! output file for all nodes
  iunout = freeunit()
  outfilenode=trim(outfile)
  if( trim(nodenumber) /= '' ) &
    outfilenode = trim(outfilenode)//'.'//trim(nodenumber)
  open(iunout,file=trim(outfilenode),form='formatted',iostat=ierr)
  if( ierr /= 0 ) &
    call errore('epsilon_init','problem opening file '//trim(outfilenode),102)

  ! now start reading files
  write(stdout,*) ' reading hamiltonian ...'
  call read_hamq( iunhq )

  ! read atomic matrix elements for non-local potential
  allocate( vnl_atom(natom) )
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type(it)
    write(stdout,*) iatom, it, np
    allocate( vnl_atom(iatom)%matrix(np,np) )
  enddo
  call read_nloper( iunvnl, vnl_atom )

  ! checks on core-repair
  if( corerep%nspecies /= ntype ) &
    call errore('epsilon_init','core-repair has incorrect no. of species', &
                abs(corerep%nspecies))
  if( corerep%natom /= natom ) &
    call errore('epsilon_init','core-repair has incorrect no. of atoms', &
                abs(corerep%natom))
  do icore=1,corerep%ncore
    iatom = corerep%core(icore)%atom
    it = type_atom(iatom)
    np = nproj_type(corerep%core(icore)%species)
    if( it /= corerep%core(icore)%species ) then
      call errore('epsilon_init','species mismatch with core-repair',icore)
    endif
    if( np /= corerep%core(icore)%nproj1 .or. &
        np /= corerep%core(icore)%nproj2 ) then
      write(stdout,*) ' no. valence projections from core-repair = ', &
                      corerep%core(icore)%nproj1, &
                      corerep%core(icore)%nproj2
      write(stdout,*) ' no. valence projections for species ', &
                      corerep%core(icore)%species, ' = ', np
      call errore('epsilon_init','projector mismatch with core-repair',icore)
    endif
  enddo

  ! pwmatrix elements
  call pwi_read( pwifile, ionode_id, mpime )
  call pwm_open( pwmfile )

  return
  end subroutine epsilon_init


  end module epsilon_module
