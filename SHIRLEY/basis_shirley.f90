  module basis_shirley

  use kinds, only : dp
  USE mp_global,    ONLY : intra_pool_comm
  use mp, only : mp_sum
  use io_global, only : stdout

  implicit none

  private
  public :: eigU, indx_basis, indx_kpbnd, nbasis, nbasis_trunc
  public :: bandstruc
  public :: basis_init, construct_overlap, diagonalize_overlap, &
            truncate_basis, expand_basis

  complex(dp),parameter :: zero=(0.d0,0.d0)
  complex(dp),parameter :: one =(1.d0,0.d0)

  integer :: nbasis, nbasis_trunc
  real(dp),allocatable :: eigU(:)
  complex(dp),allocatable :: S_l(:,:), B_l(:,:), Sij(:)

  type bandstruc
    integer :: ik
    integer :: ibnd
  end type bandstruc

  ! mapping from one set of indices to the other
  type(bandstruc),allocatable :: indx_basis(:)
  integer,allocatable :: indx_kpbnd(:,:)
  integer,allocatable :: irev(:)

  contains

    subroutine basis_init( nkpt, nbnd )

    use scalapack_module

    integer,intent(in) :: nkpt, nbnd

    integer :: ierr, i, ik, ibnd
    integer :: nbasis_lr, nbasis_lc

    nbasis = nkpt*nbnd

    ! initialize scalapack
    write(stdout,*) ' scalapack_init'
    call scalapack_init
    write(stdout,*) ' scalapack_init', nprocs, nprow, npcol, myrow, mycol

    ! distribute using scalapack
    write(stdout,*) ' scalapack_distrib'
    call scalapack_distrib( nbasis, nbasis, nbasis_lr, nbasis_lc )
    write(stdout,*) ' allocating local part of overlap matrix '
    write(stdout,*) nbasis_lr, ' x ', nbasis_lc
    allocate( Sij(nbasis), &
              S_l(nbasis_lr, nbasis_lc), &
              B_l(nbasis_lr, nbasis_lc), &
              eigU(nbasis), &
              stat=ierr )

    if( ierr /= 0 ) &
      call errore('basis_init','unable to allocate space for Shirley basis',1)

    ! create map from i -> n,k
    allocate( indx_basis(nbasis), &
              indx_kpbnd(nbnd,nkpt) )
    i=0
    do ik=1,nkpt
      do ibnd=1,nbnd
        i=i+1
        indx_basis(i)%ik = ik
        indx_basis(i)%ibnd = ibnd
        indx_kpbnd(ibnd,ik) = i
      enddo
    enddo

    end subroutine basis_init


    subroutine construct_overlap( nkpt, nbnd )

    use wfc_shirley, only : read_shirley_wfc_1, read_shirley_wfc_2, &
                            dotprod_shirley_wfc
    USE mp_global, ONLY : nproc
    use scalapack_module

    integer,intent(in) :: nkpt, nbnd
    integer :: i, j, ik, jk, ibnd, jbnd, i_l, j_l
    logical :: islocal

    write(stdout,*)
    write(stdout,*) ' constructing the overlap matrix'
    write(stdout,*)

    ! construct overlap
    i=0
    do ik=1,nkpt
      do ibnd=1,nbnd
        i=i+1
        indx_kpbnd(ibnd,ik) = i
      enddo
    enddo

    do jk=1,nkpt
      call read_shirley_wfc_1( jk )
      do ik=1,jk
        call read_shirley_wfc_2( ik )
        if( ik /= jk ) then
          do jbnd=1,nbnd
            j = indx_kpbnd(jbnd,jk)
            do ibnd=1,nbnd
              Sij(ibnd) = dotprod_shirley_wfc( ibnd, jbnd )
            enddo
    ! sum contributions from all processors
#ifdef __PARA
            call mp_sum( Sij, intra_pool_comm )
#endif
            do ibnd=1,nbnd
              i = indx_kpbnd(ibnd,ik)
              call scalapack_localindex( i, j, i_l, j_l, islocal )
              if( .not. islocal ) cycle
              S_l(i_l,j_l)=Sij(ibnd)
            enddo
            do ibnd=1,nbnd
              i = indx_kpbnd(ibnd,ik)
              if( i/=j ) then
                call scalapack_localindex( j, i, j_l, i_l, islocal )
                if( .not. islocal ) cycle
                S_l(j_l,i_l)=conjg(Sij(ibnd))
              endif
            enddo
          enddo
        else
          do jbnd=1,nbnd
            j = indx_kpbnd(jbnd,jk)
            do ibnd=1,jbnd
              Sij(ibnd) = dotprod_shirley_wfc( ibnd, jbnd )
            enddo
    ! sum contributions from all processors
#ifdef __PARA
            call mp_sum( Sij, intra_pool_comm )
#endif
            do ibnd=1,jbnd
              i = indx_kpbnd(ibnd,ik)
              call scalapack_localindex( i, j, i_l, j_l, islocal )
              if( .not. islocal ) cycle
              S_l(i_l,j_l) = Sij(ibnd)
            enddo
            do ibnd=1,jbnd
              i = indx_kpbnd(ibnd,ik)
              if( i/=j ) then
                call scalapack_localindex( j, i, j_l, i_l, islocal )
                if( .not. islocal ) cycle
                S_l(j_l,i_l) = conjg(Sij(ibnd))
              endif
            enddo
          enddo
        endif
      enddo  ! ik
    enddo  ! jk

    write(stdout,*) ' done with overlap'

    end subroutine construct_overlap


    subroutine diagonalize_overlap()

    use scalapack_module


    write(stdout,*) ' diagonalize overlap'

    write(stdout,*) ' scalapack_diag'
    call scalapack_diag( nbasis, S_l, eigU, B_l )

    end subroutine diagonalize_overlap


    subroutine truncate_basis( trace_tol )

    real(dp),intent(in) :: trace_tol
    integer :: i
    real(dp) :: trace, pctrc

    ! redefine the eigenvalues
    eigU = - eigU**2
    write(stdout,*) ' done diagonalize overlap'

    ! report and truncate
    ! note that diagonalizing S instead of U means that the eigenvalues
    ! are in increasing order, so we need to reverse them
    allocate( irev(nbasis) )
    do i=1,nbasis
      irev(i)=nbasis-i+1
    enddo
    trace = sum(eigU)
    write(stdout,*)
    write(stdout,*) ' -trace(U=-S^2) = ', -trace
    pctrc = 0.d0
    write(stdout,'(2x,a10,3a19)') 'i', 'eigU(i)', 'coverage', '%cover'
    do i=1,nbasis
      pctrc=pctrc+eigU(irev(i))
      write(stdout,'(2x,i10,2e19.10,3f19.10)') i, eigU(irev(i)), -pctrc, pctrc/trace*100.d0
    enddo
  
    ! truncation
    nbasis_trunc = nbasis
    if( trace_tol > 0 ) then
      pctrc=0.d0
      do i=1,nbasis
        pctrc=pctrc+eigU(irev(i))
        if( abs(pctrc/trace) > 1.d0-abs(trace_tol) ) exit
      enddo
      nbasis_trunc = min( i, nbasis )
      write(stdout,*) 'truncating to ', nbasis_trunc, ' basis functions'
      write(stdout,*)
    endif

    end subroutine truncate_basis

    subroutine expand_basis( npw )

    use wavefunctions_module, only : evc
    use wfc_shirley, only : read_shirley_wfc_1, evc_1
    use scalapack_module

    integer,intent(in) :: npw
    complex(dp),allocatable :: ztmp(:,:)
    integer :: i, ik, ibnd, j, k
    integer :: nk, nbnd, i_l, j_l
    logical :: islocal
    real(dp) :: invnorm

    write(stdout,*) ' expand new basis functions'

    nbnd = size(indx_kpbnd,1)
    nk = size(indx_kpbnd,2)

    evc = zero

    allocate( ztmp(nbnd,nbasis_trunc) )
    do ik=1,nk
      call read_shirley_wfc_1( ik )

      ztmp=zero
      do j=1,nbasis_trunc
      do ibnd=1,nbnd
!        i=ibnd+(ik-1)*nbnd
        i = indx_kpbnd(ibnd,ik)
        ! note that we want the truncated index in reverse order
        call scalapack_localindex( i, irev(j), i_l, j_l, islocal )
        if( .not. islocal ) cycle
 
        ztmp(ibnd,j) = B_l(i_l,j_l)
      enddo
      enddo
#ifdef __PARA
      call mp_sum( ztmp, intra_pool_comm )
#endif
 
      call ZGEMM( 'N', 'N', npw, nbasis_trunc, nbnd, one, &
                  evc_1, size(evc_1,1), ztmp, nbnd, one, evc, size(evc,1) )

    enddo
    deallocate( ztmp )

    deallocate( irev )

    ! normalize
    write(stdout,*) ' normalize new basis functions'
    do j=1,nbasis_trunc
      invnorm = real(dot_product(evc(:,j),evc(:,j)))
#ifdef __PARA
      call mp_sum( invnorm, intra_pool_comm )
#endif
      invnorm = 1.d0 / sqrt(invnorm)
      evc(:,j) = evc(:,j) * invnorm
    enddo

    end subroutine expand_basis

  end module basis_shirley
