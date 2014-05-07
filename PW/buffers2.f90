!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE buffers2
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  PRIVATE
  PUBLIC :: open_buffer2, init_buffer2, get_buffer2, save_buffer2, close_buffer2
  !
  SAVE
  !
  ! ... global variables
  !
  COMPLEX(DP), ALLOCATABLE :: buffer1(:,:)
  INTEGER :: nword_
  CHARACTER(LEN=80) :: extension_
  !
  CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE open_buffer2 (unit, extension, nword, maxrec, exst)
  !-----------------------------------------------------------------------
  !
  !     unit > 6 : connect unit "unit" to a file "prefix"."extension" in
  !     tmp_dir for direct I/O access, record length nword complex numbers;
  !     maxrec is ignored, exst=T(F) if the file (does not) exists
  !
  !     unit =-1 : allocate a buffer for storing up to maxrec records
  !     of length nword complex numbers; extension is saved but ignored
  !     exst=T(F) if the buffer is already allocated
  !
  USE io_files,  ONLY : diropn
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: extension
  INTEGER, INTENT(IN) :: unit, nword, maxrec
  LOGICAL, INTENT(OUT) :: exst
  !
  INTEGER :: ierr
#ifdef __GFORTRAN
  complex(dp),allocatable :: buf(:)
#endif
  !
  IF ( unit == -1 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        CALL errore ('open_buffer', 'buffer already allocated',1)
        !
     ELSE
        !
        nword_ = nword
        extension_ = extension
        ALLOCATE ( buffer1 ( nword, maxrec ) )
        call init_buffer2 (unit, exst, ierr)
        !
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     CALL diropn (unit, extension, 2*nword, exst)

#ifdef __GFORTRAN
! davegp
! Hack - there seems to be a bug in gfortran - if you save buffers larger
!        4kB out of order with direct access then the file size is wrong
!        so, not all records are written correctly - this becomes a problem
!        when nspin>1, since evc are written out of order.
!        Writing the last (zeroed) record seems to fix the file size.
!
     if( .not. exst ) then
       allocate( buf(nword) )
       buf=(0.d0,0.d0)
       call save_buffer( buf, nword, unit, maxrec )
       deallocate( buf )
     endif
! davegp
#endif
     !
  ELSE
     !
     CALL errore ('open_buffer', 'incorrect unit specified', ABS(unit))
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE open_buffer2
!
!----------------------------------------------------------------------------
SUBROUTINE save_buffer2( vect, nword, unit, nrec )
  !----------------------------------------------------------------------------
  !
  ! ... copy vect(1:nword) into the "nrec"-th record of
  ! ... - a previously allocated buffer, if unit = -1
  ! ... - a previously opened direct-access file with unit > 6
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec
  COMPLEX(DP), INTENT(IN) :: vect(nword)
  !
  IF ( unit == -1 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        IF ( nrec > SIZE ( buffer1, 2) )  &
           CALL errore ('save_buffer', 'too many records', ABS(nrec))
        !
        IF ( nword /= SIZE ( buffer1, 1) )  &
           CALL errore ('save_buffer', 'record length mismatch', ABS(nword))
        !
        buffer1(:,nrec) = vect(:)
        !
     ELSE
        !
        CALL errore ('save_buffer', 'buffer not allocated', ABS(unit))
        !
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     CALL davcio ( vect, 2*nword, unit, nrec, +1 )
     !
  ELSE
     !
     CALL errore ('save_buffer', 'incorrect unit specified', ABS(unit))
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE save_buffer2
!
!----------------------------------------------------------------------------
SUBROUTINE get_buffer2( vect, nword, unit, nrec )
  !----------------------------------------------------------------------------
  !
  ! ... copy vect(1:nword) from the "nrec"-th record of
  ! ... - a previously allocated buffer, if unit = -1
  ! ... - a previously opened direct-access file with unit > 6
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec
  COMPLEX(DP), INTENT(OUT) :: vect(nword)
  !
  IF ( unit == -1 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        IF ( nrec > SIZE ( buffer1, 2) )  &
           CALL errore ('get_buffer', 'no such record', ABS(nrec))
        !
        IF ( nword /= SIZE ( buffer1, 1) )  &
           CALL errore ('get_buffer', 'record length mismatch', ABS(nword))
        !
        vect(:) = buffer1(:,nrec)
        !
     ELSE
        !
        CALL errore ('get_buffer', 'buffer not allocated', ABS(unit))
        !
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     CALL davcio ( vect, 2*nword, unit, nrec, -1 )
     !
  ELSE
     !
     CALL errore ('get_buffer', 'incorrect unit specified', ABS(unit))
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE get_buffer2
!
SUBROUTINE close_buffer2 ( unit, status )
  !
  !     unit > 6 : close unit with status "status" ('keep' or 'delete')
  !     unit =-1 : deallocate buffer; if "status='keep'" save to file
  !                (using saved value of extension)
  !
  USE io_files,       ONLY : find_free_unit, diropn
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: unit
  CHARACTER(LEN=*), INTENT(IN) :: status
  !
  INTEGER :: unit_, i
  LOGICAL :: exst, opnd
  !
  IF ( unit == -1 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        IF ( TRIM(status) == 'KEEP' .OR. TRIM(status) == 'keep') THEN
           !
           unit_ = find_free_unit () 
           CALL diropn (unit_, extension_, 2*nword_, exst)
           DO i = 1, SIZE (buffer1, 2)
              CALL davcio ( buffer1(1,i), 2*nword_, unit_, i, +1 )
           END DO
           CLOSE( UNIT = unit_, STATUS = status )
           !
        END IF
        !
        DEALLOCATE (buffer1)
        !
     ELSE
        !
        CALL infomsg ('close_buffer', 'buffer not allocated')
        !
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     INQUIRE( UNIT = unit, OPENED = opnd )
     !
     IF ( opnd ) CLOSE( UNIT = unit, STATUS = status )
     !
  ELSE
     !
     CALL infomsg ('get_buffer', 'incorrect unit specified')
     !
  END IF
  !
END SUBROUTINE close_buffer2
!
SUBROUTINE init_buffer2 ( unit, exst, ierr )
  !
  !     unit > 6 : ignored
  !     unit =-1 : read into buffer the array previously saved to file
  !                when the buffer was closed (used in NEB calculations)
  !     exst     : T if the file where to read from is present
  !     ierr     : 0 if everything ok, 1 otherwise
  !
  USE io_files,       ONLY : find_free_unit, diropn
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: unit
  INTEGER, INTENT(OUT) :: ierr
  LOGICAL, INTENT(OUT) :: exst
  !
  INTEGER :: unit_, i
  !
  ierr = 1 
  !
  IF ( unit == -1 ) THEN
     !
     IF ( .NOT. ALLOCATED ( buffer1 ) ) THEN
        CALL infomsg ('init_buffer2', 'buffer not allocated')
        RETURN
     END IF
     !
     unit_ = find_free_unit () 
     CALL diropn (unit_, extension_, 2*nword_, exst)
     IF ( .NOT. exst ) THEN
        CLOSE (UNIT = unit_ , STATUS = 'delete')
        RETURN
     END IF
     !
     DO i = 1, SIZE (buffer1, 2)
        CALL davcio ( buffer1(1,i), 2*nword_, unit_, i, -1 )
     END DO
     CLOSE( UNIT = unit_, STATUS = 'keep' )
     ierr = 0
     !
  ELSE
     !
     CALL infomsg ('init_buffer2', 'incorrect unit specified')
     !
  END IF
  !
END SUBROUTINE init_buffer2
!
END MODULE buffers2
