  subroutine wfc_reorder( ik, ik0, npw, npw0, igwx0, igk_l2g0 )

  ! reorder the wave function coefficients so that k-point ik has the
  ! same ordering as k-point ik0

  ! David Prendergast, UCB, Jan 2007

  USE kinds, ONLY : dp
  USE mp_wave
  USE mp,        ONLY : mp_max, mp_barrier, mp_bcast
  USE mp_global, ONLY : npool, nproc
  USE io_global, ONLY : ionode_id, stdout
  USE wvfct,     ONLY : nbnd, igk_l2g, npwx
  USE klist, ONLY : ngk
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module, ONLY : noncolin
  USE parallel_include

  USE basis_shirley, only : evcb

  implicit none

  integer,intent(in) :: ik, ik0, npw, npw0, igwx0
  integer,intent(in) :: igk_l2g0(igwx0)

  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),parameter :: one  = cmplx(1.d0,0.d0)

  integer :: j, i, ierr
  complex(dp),allocatable :: wtmp(:)
  complex(dp) :: dotprod, dotprodw

  ! error checking
  if( npool > 1 ) &
    call errore('wfc_reorder','npool must be set to 1',1)
  if( noncolin ) &
    call errore('wfc_reorder','not implemented for noncolin',1)

  write(stdout,*) ' igwx0 = ', igwx0
  write(stdout,*) ' npwx = ', npwx
  write(stdout,*) ' npw  = ', npw
  write(stdout,*) ' npw0 = ', npw0

  allocate( wtmp(igwx0), stat=ierr )
  if( ierr /= 0 ) &
    call errore('wfc_reorder','unable to allocate temporary space',igwx0)

  do j=1,nbnd
    ! zero the temp array
    wtmp = zero

    ! merge wave function from evc into wtmp
    do i=1,npw
      wtmp(igk_l2g(i,ik)) = evc(i,j)
    enddo
    call reduce( 2*igwx0, wtmp )

    ! check
    dotprod = dot_product( evc(1:npw,j), evc(1:npw,j) )
    call reduce( 2, dotprod )
    dotprodw = dot_product(wtmp, wtmp)
    write(stdout,'(i,3e)') j, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
    call flush_unit(stdout)

    ! zero the former wave function
    evcb(:,j) = zero

    ! split wave function from wtmp into evcb
    call mp_bcast( wtmp, ionode_id )
    do i=1,npw0
      evcb(i,j) = wtmp(igk_l2g0(i))
    enddo

    dotprod = dot_product( evcb(1:npw0,j), evcb(1:npw0,j) )
    call reduce( 2, dotprod )
    write(stdout,'(i,3e)') j, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
    call flush_unit(stdout)

  enddo

  deallocate( wtmp )

  end subroutine wfc_reorder
