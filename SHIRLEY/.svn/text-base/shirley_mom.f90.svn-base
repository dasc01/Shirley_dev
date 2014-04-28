  program shirley_mom

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and produce x-ray absorption spectra
  ! by averaging over a k-point grid

  ! David Prendergast, UCB, Jul 2007

#include "f_defs.h"
  use hamq_shirley
!  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime
  use mp, only : mp_bcast, mp_end, mp_barrier
  use kpt_module
  use corerepair_module

  implicit none

!  integer,parameter :: stdout=6
  integer,parameter :: stdin =5
  integer,parameter :: dp = kind(1.d0)
  integer,parameter :: maxchar=255
  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  type(kpt_type) :: kpt
  type(corerepair_type) :: corerep

  character(len=3) :: nodenumber
  character(maxchar) :: hamqfile
  character(maxchar) :: vnlfile
  character(maxchar) :: outfile
  character(maxchar) :: fmtstr

  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: zeigval(:)
  complex(dp),allocatable :: eigvec(:,:)
  complex(dp),allocatable :: betaq(:,:)
  integer :: iunhq, iunout, iunvnl, ierr
  integer :: iunoutgnu, iunoutstd
  character(maxchar) :: outfilenode

  integer :: nblock, lwork
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)

  type(matrix_list),allocatable :: vnl_atom(:)
  integer :: iatom, it, np, npm
  complex(dp),allocatable :: beta_block(:,:), vnl(:,:), ztmp(:,:)
  complex(dp),allocatable :: mom(:,:,:), mom_cr(:,:,:), mtmp(:,:)

  integer,external :: freeunit
  integer,external :: ilaenv

  integer :: ik, i, j
  integer :: icore, ibasis, ncm, ixyz

  complex(dp),allocatable :: betaq_tmp(:,:)
  complex(dp),allocatable :: ctmp(:,:)

  namelist / input / hamqfile, vnlfile, outfile


  ! initialize mpi
  CALL start_shirley (nodenumber) 

! debugging
!  if( .not. ionode ) stdout = 500+mpime
  call init_stdout( stdout )

  write(stdout,*) ' shirley_mom'
  write(stdout,*)

  ! read name list
  if( ionode ) then
    read(stdin,nml=input,iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_mom','problem reading namelist &input',101)
  endif

  if( ionode ) call kpt_read( stdin, kpt )
  write(stdout,*) 'done reading k-points'
  ! divide evenly among processes
  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

  ! broadcast namelist input first
  call mp_bcast( hamqfile, ionode_id )
  call mp_bcast( vnlfile, ionode_id )
  call mp_bcast( outfile, ionode_id )

  ! open files
  if( ionode ) then
    iunhq = freeunit()
    open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) call errore('shirley_mom','problem opening file '//trim(hamqfile),102)

    iunvnl = freeunit()
    open(iunvnl,file=trim(vnlfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) stop 104
  endif

  if( ionode .and. trim(kpt%param%grid_type) == 'bandstructure' ) then
    iunoutgnu = freeunit()
    open(iunoutgnu,file=trim(outfile)//'.gnu',form='formatted',iostat=ierr)
    if( ierr /= 0 ) call errore('shirley_mom','problem opening file '//trim(outfile)//'.gnu',103)
  endif

  iunout = freeunit()
  outfilenode=trim(outfile)
  if( trim(nodenumber) /= '' ) &
    outfilenode = trim(outfilenode)//'.'//trim(nodenumber)
  open(iunout,file=trim(outfilenode),form='unformatted',iostat=ierr)
  if( ierr /= 0 ) &
    call errore('shirley_mom','problem opening file '//trim(outfilenode),102)

!  iunoutstd = freeunit()
!  outfilenode=trim(outfile)
!  if( trim(nodenumber) /= '' ) &
!    outfilenode = trim(outfilenode)//'.std'//trim(nodenumber)
!  open(iunoutstd,file=trim(outfilenode),form='formatted',iostat=ierr)
!  if( ierr /= 0 ) &
!    call errore('shirley_mom','problem opening file '//trim(outfilenode),102)

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

  ! read atomic matrix elements for core-valence position
  if( ionode ) call read_corerepair( stdin, corerep )
  call bcast_corerepair( corerep, ionode_id )

  ! checks on core-repair
  if( corerep%nspecies /= ntype ) &
    call errore('shirley_mom','core-repair has incorrect no. of species', &
                abs(corerep%nspecies))
  if( corerep%natom /= natom ) &
    call errore('shirley_mom','core-repair has incorrect no. of atoms', &
                abs(corerep%natom))
  do icore=1,corerep%ncore
    iatom = corerep%core(icore)%atom
    it = type_atom(iatom)
    np = nproj_type(corerep%core(icore)%species)
    if( it /= corerep%core(icore)%species ) then
      call errore('shirley_mom','species mismatch with core-repair',icore)
    endif
    if( np /= corerep%core(icore)%nproj1 .or. &
        np /= corerep%core(icore)%nproj2 ) then
      write(stdout,*) ' no. valence projections from core-repair = ', &
                      corerep%core(icore)%nproj1, &
                      corerep%core(icore)%nproj2
      write(stdout,*) ' no. valence projections for species ', &
                      corerep%core(icore)%species, ' = ', np
      call errore('shirley_mom','projector mismatch with core-repair',icore)
    endif

  enddo


  allocate( eigval(nbasis), eigvec(nbasis,nbasis) )
  allocate( zeigval(nbasis) )
  allocate( betaq(nproj,nbasis) )
  allocate( betaq_tmp(nproj,nbasis) )
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


  write(fmtstr,'(a,i,a)') '(i,4e,',nbasis,'e)'

  npm = maxval( nproj_type(1:ntype) )
  allocate( beta_block(npm,nbasis), vnl(npm,npm), ztmp(npm,nbasis), ctmp(npm,npm) )
  ncm = maxval( corerep%core(1:corerep%ncore)%nproj2 )
  allocate( mom(nbasis,nbasis,3), &
            mom_cr(nbasis,nbasis,3), &
            mtmp(nbasis,nbasis) )

! ======================================================================
  do ik=1,kpt%list%nk
! ======================================================================

    write(stdout,'(a,i,a,i,a,a)') ' k-point ', ik, ' of ', kpt%list%nk, &
                                  ' on node ', trim(nodenumber)
    write(stdout,'(3f,4x,f)') kpt%list%kvec(1:3,ik), kpt%list%wk(ik)

    write(iunout) kpt%list%kvec(1:3,ik), kpt%list%wk(ik)
    write(iunout) 7 ! eigenvalues and 2 sets of 3 Cartesian transition matrix

!    write(iunoutstd,*) kpt%list%kvec(1:3,ik), kpt%list%wk(ik)
!    write(iunoutstd,*) 4 ! eigenvalues and 3 Cartesian transition matrix

    ! build the Hamiltonian for this q-point
    ! local contribution
    call build_hamq_local( kpt%list%kvec(1:3,ik), kpt%param%cartesian, eigvec )

    ! non-local contribution
    call build_hamq_projs( kpt%list%kvec(1:3,ik), kpt%param%cartesian, betaq )

      do iatom=1,natom
        it = type_atom(iatom)
        np = nproj_type_nl(it)

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

!    endif

    ! diagonalize matrix to find eigenvalues and vectors
    CALL ZHEEV( 'V', 'U', nbasis, eigvec, nbasis, eigval, work, lwork, rwork, ierr )
    if( ierr /= 0 ) stop 105

    ! write eigenvalues
    zeigval = cmplx( eigval )

!    write(stdout,*) ' eigenvalues '
!    write(stdout,'(2i,a)') nbasis, 1, '  NELEC'
!    do ibasis=1,nbasis
!      write(stdout,'(i,3e)') ibasis, eigval(ibasis), eigval(ibasis)*rytoev, kpt%list%wk(ik)
!    enddo

    write(iunout) nbasis, 1
    write(iunout) zeigval

    ! construct momentum matrix element
    call mom_hamq( kpt%list%kvec(1:3,ik), kpt%param%cartesian, mom )
    
    ! expand in the eigenbasis
    ! mom
    do ixyz=1,3
      CALL ZGEMM( 'N', 'N', nbasis, nbasis, nbasis, one, &
                  mom(1,1,ixyz), nbasis, &
                  eigvec, nbasis, &
                  zero, mtmp, nbasis )
      CALL ZGEMM( 'C', 'N', nbasis, nbasis, nbasis, one, &
                  eigvec, nbasis, &
                  mtmp, nbasis, &
                  zero, mom(1,1,ixyz), nbasis )
    enddo

    ! dump betaq - in eigenbasis
    betaq_tmp = betaq
    CALL ZGEMM( 'N', 'N', nproj, nbasis, nbasis, one, &
                betaq_tmp, nproj, &
                eigvec, nbasis, &
                zero, betaq, nproj )
!    write(stdout,*) ' |betaq|^2'
!    do ibasis=1,20
!    do i=1,nproj
!      write(stdout,'(2i,3f)') i, ibasis, real(conjg(betaq(i,ibasis))*betaq(i,ibasis)), betaq(i,ibasis)
!    enddo
!      write(stdout,*)
!    enddo

    ! construct corerepair
    mom_cr = 0.d0  ! beginning value
    do icore=1,corerep%ncore

      iatom = corerep%core(icore)%atom
      it = type_atom(iatom)
      np = nproj_type(corerep%core(icore)%species)

      write(stdout,*) ' iatom = ', iatom
      write(stdout,*) ' it = ', it
      write(stdout,*) ' np = ', np
      write(stdout,*) ' index_betaq = ', index_betaq(1:np,iatom)
      !write(stdout,*) ' corerep_x = ', corerep%core(icore)%matrix(1:np,1:np,1)

      forall( i=1:np, j=1:nbasis ) &
        beta_block(i,j) = &
          betaq(index_betaq(i,iatom),j)

      ! sum <i|beta>mom_cr<beta|j>
      do ixyz=1,3
        CALL ZGEMM( 'N', 'N', np, nbasis, np, one, &
                    corerep%core(icore)%matrix(1,1,ixyz), &
                    corerep%core(icore)%nproj1, &
                    beta_block, npm, &
                    zero, ztmp, npm )
        CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, one, &
                    beta_block, npm, &
                    ztmp, npm, &
                    one, mom_cr(1,1,ixyz), nbasis )
      enddo ! loop over ixyz

    enddo ! loop over icore

    ! add regular momentum
    mom_cr = mom_cr + mom

    ! write transition matrix elements
    do ixyz=1,3
      write(iunout) nbasis, nbasis
      write(iunout) mom(:,:,ixyz)
      write(iunout) nbasis, nbasis
      write(iunout) mom_cr(:,:,ixyz)
    enddo

!    ! write transition matrix elements
!    write(stdout,*) ' mom'
!    do j=1,nbasis
!    do i=1,32
!      write(stdout,'(2i,6e)') i, j, mom(i,j,1:3)
!    enddo
!      write(stdout,*)
!    enddo
!    write(stdout,*) ' mom_cr'
!    do j=1,nbasis
!    do i=1,32
!      write(stdout,'(2i,6e)') i, j, mom_cr(i,j,1:3)
!    enddo
!      write(stdout,*)
!    enddo

    ! make sure iunout is written
    call flush(iunout)

! ======================================================================
  enddo ! loop over k-points ik
! ======================================================================

  ! close binary dump
  close(iunout)

  ! deallocate
  deallocate( eigval, eigvec, betaq )
  deallocate( beta_block, vnl, ztmp )
  deallocate( betaq_tmp )

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_mom'
  stop
  
  end program shirley_mom
