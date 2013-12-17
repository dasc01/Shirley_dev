!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE vxcinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes the exchange correlation potential in the array
  ! ... vr
  !
  ! ...    read rho from the file produced by the scf calculation
  ! 
  ! ... the scf potential is recalculated and saved in vr
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE cell_base,        ONLY : alat, omega
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE basis,            ONLY : starting_pot
  USE klist,            ONLY : nelec
  USE lsda_mod,         ONLY : lsda, nspin
  USE gvect,            ONLY : ngm, gstart, nl, g, gg
  use grid_dimensions,  only : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE gvecs,            ONLY : doublegrid
  USE control_flags,    ONLY : lscf
  USE scf,              ONLY : rho, rho_core, vltot, vr, vrs
  USE ener,             ONLY : ehart, etxc, vtxc
  USE ldaU,             ONLY : niter_with_fixed_ns
  USE ldaU,             ONLY : lda_plus_u, Hubbard_lmax, ns, nsnew
  USE noncollin_module, ONLY : noncolin, factlist, pointlist, pointnum, &
                               mcons, i_cons, lambda, vtcon, report
  USE io_files,         ONLY : prefix, iunocc, input_drho
  USE spin_orb,         ONLY : domag
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : intra_image_comm
  USE io_global,        ONLY : ionode, ionode_id
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL (DP) :: charge           ! the starting charge
  REAL (DP) :: etotefield       ! 
  INTEGER        :: ios
  INTEGER        :: ldim             ! integer variable for I/O control
  LOGICAL        :: exst 
  !
  !
  IF ( ionode ) THEN
     !
     CALL seqopn( 4, 'rho', 'UNFORMATTED', exst )
     !
     IF ( exst ) THEN
        !
        CLOSE( UNIT = 4, STATUS = 'KEEP' )
        !
     ELSE
        !
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
     END IF
     !
  END IF
  !
  CALL mp_bcast( exst, ionode_id, intra_image_comm )
  !
  IF ( starting_pot == 'file' .AND. exst ) THEN
     ! 
     ! ... Cases a) and b): the charge density is read from file
     !
     CALL io_pot( -1, 'rho', rho, nspin )
     !       
     IF ( lscf ) THEN
        WRITE( stdout, '(/5X,"The initial density is read from file ", A20)' )&
           TRIM( prefix ) // '.rho'
     ELSE
        WRITE( stdout, '(/5X,"The potential is recalculated from file ",A20)')&
           TRIM( prefix ) // '.rho'
     END IF
     !
     ! ... The occupations ns also need to be read in order to build up 
     ! ... the potential
     !
     IF ( lda_plus_u ) THEN  
        !
        ldim = 2 * Hubbard_lmax + 1
        !
        IF ( ionode ) THEN
           !
           CALL seqopn( iunocc, 'occup', 'FORMATTED', exst )
           READ( UNIT = iunocc, FMT = * ) ns
           CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
           !
        ELSE
           !  
           ns(:,:,:,:) = 0.D0
           !
        END IF
        !
        CALL reduce( ( ldim * ldim * nspin * nat ), ns )  
        CALL poolreduce( ( ldim * ldim * nspin * nat ), ns )  
        !
        nsnew = ns
        !
     END IF
     !
  ELSE
     !   
     call errore( 'vxcinit','density not found in .rho file',1)
     !
  END IF
  !
  ! ... check the integral of the starting charge
  !
  IF ( nspin == 2 ) THEN
     !
     charge = SUM ( rho (:, 1:nspin) ) * omega / ( nr1 * nr2 * nr3 )
     !
  ELSE
     !
     charge = SUM ( rho (:, 1) ) * omega / ( nr1 * nr2 * nr3 )
     !
  END IF
  !
  call reduce (1, charge)
  !
  IF ( lscf .AND. ABS( charge - nelec ) / charge > 1.D-6 ) THEN
     !
     WRITE( stdout, &
          '(/,5X,"starting charge ",F10.5,", renormalised to ",F10.5)') &
          charge, nelec
     !
     rho = rho / charge * nelec
     !
  ELSE IF ( .NOT. lscf .AND. ABS( charge - nelec ) / charge > 1.D-6 ) THEN
     !
     CALL errore ( 'potinit', 'starting and expected charges differ', 1 )
     !
  END IF
  !
  ! ... calculate exchange-correlation potential - store in vr
  !
  CALL v_xc( rho, rho_core, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
             nrxx, nl, ngm, g, nspin, alat, omega, etxc, vtxc, vr )
  !
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
  !
  ! ... write on output the parameters used in the lda+U calculation
  !
  IF ( lda_plus_u ) THEN
     !
     WRITE( stdout, '(/5X,"Parameters of the lda+U calculation:")')
     WRITE( stdout, '(5X,"Number of iteration with fixed ns =",I3)') &
         niter_with_fixed_ns
     WRITE( stdout, '(5X,"Starting ns and Hubbard U :")')
     !
     CALL write_ns()
     !
  END IF
  !
  IF ( report /= 0 .AND. noncolin .AND. domag .AND. lscf ) CALL report_mag()
  !
  RETURN
  !
END SUBROUTINE vxcinit
