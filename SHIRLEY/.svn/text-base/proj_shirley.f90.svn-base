  module proj_shirley

  USE kinds, ONLY : DP
  
  implicit none

  integer :: nbnd, natomwfc, nkb, nat, nsym
  integer,allocatable :: irt(:,:)

  real(dp) :: d1(3,3,48), d2(5,5,48), d3(7,7,48)

  TYPE wfc_label
     INTEGER na, n, l, m 
  END TYPE wfc_label 
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
  
  COMPLEX (DP), ALLOCATABLE :: projk (:,:)
  complex(dp),allocatable :: atomic_proj_matrix(:,:)

  contains


! ----------------------------------------------------------------------
  subroutine init_proj( natomwfc_, nbnd_, nkb_, nat_, nsym_ )
! ----------------------------------------------------------------------

  integer :: natomwfc_, nkb_, nbnd_, nat_, nsym_

  natomwfc = natomwfc_
  nkb = nkb_
  nbnd = nbnd_
  nat = nat_
  nsym = nsym_

  if( allocated( irt ) ) deallocate( irt )
  allocate( irt(48,nat) )

  if( allocated( nlmchi ) ) deallocate( nlmchi )
  allocate( nlmchi(natomwfc) )

  if( allocated( atomic_proj_matrix ) ) deallocate( atomic_proj_matrix )
  allocate( atomic_proj_matrix(natomwfc,nkb) )

  end subroutine init_proj


! ----------------------------------------------------------------------
  subroutine write_proj( iunprj )
! ----------------------------------------------------------------------

  integer :: iunprj

  write(iunprj) natomwfc, nkb, nbnd, nat, nsym
  write(iunprj) irt
  write(iunprj) d1, d2, d3
  write(iunprj) nlmchi
  write(iunprj) atomic_proj_matrix

  end subroutine write_proj


! ----------------------------------------------------------------------
  subroutine read_proj( iunprj )
! ----------------------------------------------------------------------

  integer :: iunprj

  read(iunprj) natomwfc, nkb, nbnd, nat, nsym
  call init_proj( natomwfc, nkb, nbnd, nat, nsym )
  read(iunprj) irt
  read(iunprj) d1, d2, d3
  read(iunprj) nlmchi
  read(iunprj) atomic_proj_matrix

  end subroutine read_proj


! ----------------------------------------------------------------------
  subroutine bcast_proj( mpime, root )
! ----------------------------------------------------------------------

  use mp, only : mp_bcast

  integer :: mpime, root
  integer :: i

  call mp_bcast( natomwfc, root )
  call mp_bcast( nkb, root )
  call mp_bcast( nbnd, root )
  call mp_bcast( nat, root )
  call mp_bcast( nsym, root )
  if( mpime /= root ) then
    call init_proj( natomwfc, nkb, nbnd, nat, nsym )
  endif
  call mp_bcast( irt, root )
  call mp_bcast( d1, root )
  call mp_bcast( d2, root )
  call mp_bcast( d3, root )
  do i=1,natomwfc
    call mp_bcast( nlmchi(i)%na, root )
    call mp_bcast( nlmchi(i)%n,  root )
    call mp_bcast( nlmchi(i)%l,  root )
    call mp_bcast( nlmchi(i)%m,  root )
  enddo
  call mp_bcast( atomic_proj_matrix, root )

  end subroutine bcast_proj


! ----------------------------------------------------------------------
  subroutine eigk_proj( beta )
! ----------------------------------------------------------------------

  complex(dp),intent(in) :: beta(nkb,nbnd)
  integer :: ibnd

  if( .not. allocated(projk) ) allocate( projk(natomwfc,nbnd) )

  call zgemm( 'N', 'N', natomwfc, nbnd, nkb, (1.d0, 0.d0), &
              atomic_proj_matrix, natomwfc, beta, nkb, (0.d0, 0.d0), &
              projk, natomwfc )

  do ibnd=1,nbnd
    write(*,*) ' projk ', ibnd
    write(*,*) projk(:,ibnd)
  enddo

  end subroutine eigk_proj


! ----------------------------------------------------------------------
  subroutine symm_projk( proj )
! ----------------------------------------------------------------------
  ! symmetrize the projections 
  ! make sure to have expanded the basis projections proj0 at some
  ! k-point to give non-zero projk

  real(dp),intent(out) :: proj(natomwfc,nbnd)

  INTEGER :: ibnd, na, nb, nt, isym, n,  m, m1, l, lm, nwfc, nwfc1
  complex(dp) :: work1(nbnd)

  if( .not. allocated( projk ) ) &
    call errore('symm_projk','projections not yet expanded with eigk_proj',1)

  proj=0.d0
  DO nwfc = 1, natomwfc
     ! 
     !  atomic wavefunction nwfc is on atom na 
     ! 
     na= nlmchi(nwfc)%na
     n = nlmchi(nwfc)%n
     l = nlmchi(nwfc)%l
     m = nlmchi(nwfc)%m
     ! 
     DO isym = 1, nsym
        nb = irt (isym, na)
        DO nwfc1 =1, natomwfc
           IF (nlmchi(nwfc1)%na == nb             .AND. &
                nlmchi(nwfc1)%n == nlmchi(nwfc)%n .AND. &
                nlmchi(nwfc1)%l == nlmchi(nwfc)%l .AND. &
                nlmchi(nwfc1)%m == 1 ) go to 10
        END DO
        CALL errore('proj_basis','cannot symmetrize',1)
 10     nwfc1=nwfc1-1
        ! 
        !  nwfc1 is the first rotated atomic wfc corresponding to nwfc 
        ! 
        IF (l == 0) THEN
           work1(:) = projk (nwfc1 + 1,:)
        ELSE IF (l == 1) THEN
           work1(:) = 0.d0
           DO m1 = 1, 3
              work1(:) = work1(:) + d1 (m1, m, isym) * projk (nwfc1 + m1,:) 
           ENDDO
        ELSE IF (l == 2) THEN
           work1(:) = 0.d0
           DO m1 = 1, 5
              work1(:) = work1(:) + d2 (m1, m, isym) * projk (nwfc1 + m1,:) 
           ENDDO
        ELSE IF (l == 3) THEN
           work1(:) = 0.d0
           DO m1 = 1, 7
              work1(:) = work1(:) + d3 (m1, m, isym) * projk (nwfc1 + m1,:) 
           ENDDO
        ENDIF
        DO ibnd = 1, nbnd
           proj (nwfc, ibnd) = proj (nwfc, ibnd) + &
              work1(ibnd) * CONJG (work1(ibnd)) / nsym
        ENDDO
     ENDDO
  ENDDO

  end subroutine symm_projk


  end module proj_shirley
