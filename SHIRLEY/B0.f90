! ---------------------------------------------------------------------- 
  subroutine B0( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode
  USE cell_base
  USE gvect  
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk
  USE wvfct
  use control_flags, only : gamma_only
  USE gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
! davegp
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm

  implicit none

  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)

  integer :: ibnd
  complex(dp),allocatable :: B0vector(:)
  real(dp) :: G0vector

  integer :: ierr, iunb0

  logical :: exst
  integer,external :: freeunit

  WRITE( stdout, '(/5x,"Calling B0 .... ",/)')

  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
  ! ======================================================================
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  ! ======================================================================
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)
!  !call sum_band
  write(stdout,*) ' ecutwfc = ', ecutwfc
  write(stdout,*) ' tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nktot = ', nks, nkstot
  write(stdout,*) ' xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '     npw = ', ngk(1:nks)
  write(stdout,*) '    npwx = ', npwx


!  ! gamma-point only
!  qvec = 0.d0
!  current_k = 1
!  if( lsda ) current_spin = isk(1)

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, nwordwfc, iunwfc, 1, - 1 )

  write(stdout,*)
  write(stdout,*) ' construct B0:'
  write(stdout,*)

  ! ======================================================================
  ! B0 matrix elements
  ! ======================================================================
  write(stdout,*) ' B0 matrix elements'

  allocate( B0vector(nbnd) )
  B0vector = zero
  G0vector = 0.d0
  if( gstart == 2 ) then
    B0vector(1:nbnd) = evc(1,1:nbnd)
    G0vector = g2kin(1)
  endif
  call mp_sum( B0vector, intra_pool_comm )
  call mp_sum( G0vector, intra_pool_comm )

  ! dump to disk
  if( ionode ) then
    write(stdout,'(a,f12.5)') '|G0vector|^2 = ', G0vector
    iunb0 = freeunit()
    call seqopn(iunb0,'B0','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping B0 to file ', trim(prefix) // '.B0'

    write(iunb0) B0vector
    close(iunb0)
  endif
  
  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( B0vector )

  write(stdout,*) ' Done with B0'
  return

  end subroutine B0
! ---------------------------------------------------------------------- 

