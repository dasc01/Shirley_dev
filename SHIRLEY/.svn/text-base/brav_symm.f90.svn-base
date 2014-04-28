  subroutine brav_symm( ibrav, at, noinv, nosym, nrot, tipo, s, sname, symm_type )

  use kinds, only : dp

  implicit none

  integer,intent(in) :: ibrav
  real(dp),intent(in) :: at(3,3)
  logical,intent(in) :: noinv, nosym

  integer,intent(out) :: nrot, tipo
  integer,intent(out) :: s(3,3,48)
  CHARACTER(LEN=45),intent(out) :: sname(48)
  CHARACTER(LEN=9),intent(out) :: symm_type

  integer :: irot, isym

  !
  !  ... generate transformation matrices for the crystal point group
  !  ... First we generate all the symmetry matrices of the Bravais lattice
  !
  IF ( ibrav == 4 .OR. ibrav == 5 ) THEN
     !
     ! ... here the hexagonal or trigonal bravais lattice
     !
     CALL hexsym( at, s, sname, nrot )
     !
     tipo = 2
     !
  ELSE IF ( ibrav >=1  .AND. ibrav <= 14 ) THEN
     !
     ! ... here for the cubic bravais lattice
     !
     CALL cubicsym( at, s, sname, nrot )
     !
     tipo = 1
     !
  ELSE IF ( ibrav == 0 ) THEN
     !
     IF ( symm_type == 'cubic' ) THEN
        !
        tipo = 1
        !
        CALL cubicsym( at, s, sname, nrot )
        !
     ELSE IF ( symm_type == 'hexagonal' ) THEN
        !
        tipo = 2
        !
        CALL hexsym( at, s, sname, nrot )
        !
     END IF
     !
  ELSE
     !
     CALL errore( 'setup', 'wrong ibrav', 1 )
     !
  END IF
  !
  ! ... if noinv is .TRUE. eliminate all symmetries which exchange z with -z
  !
  IF ( noinv ) THEN
     !
     irot = 0
     !
     DO isym = 1, nrot
        IF ( s(1,3,isym) == 0 .AND. s(3,1,isym) == 0 .AND. &
             s(2,3,isym) == 0 .AND. s(3,2,isym) == 0 .AND. &
             s(3,3,isym) == 1) THEN
           !
           irot = irot + 1
           !
           s(:,:,irot) = s(:,:,isym)
           !
           sname(irot) = sname(isym)
           !
        END IF
        !
     END DO
     !
     nrot = irot
     !
  END IF
  !
  ! ... If nosym is true do not use any point-group symmetry
  !
  IF ( nosym ) nrot = 1

  end subroutine brav_symm

