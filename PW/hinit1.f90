!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hinit1()
  !----------------------------------------------------------------------------
  !
  ! ... Atomic configuration dependent hamiltonian initialization
  !
  USE ions_base,     ONLY : nat, nsp, ityp, tau
  USE cell_base,     ONLY : at, bg, omega, tpiba2
  USE grid_dimensions,ONLY: nrxx
  USE gvect,         ONLY : ngm, g
  USE gvecs,       ONLY : doublegrid
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : vrs, vltot, v, kedtau
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : qpointlist
  USE wannier_new,   ONLY : use_wannier
  USE martyna_tuckerman, ONLY : tag_wg_corr_as_obsolete
  USE dfunct,                 only : newd
  !
  IMPLICIT NONE
  !
  !
  ! ... update the wavefunctions, charge density, potential
  ! ... update_pot initializes structure factor array as well
  !
  CALL update_pot()
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid )
  !
  IF ( tqr ) CALL qpointlist()
  !
  ! ... update the D matrix
  !
  CALL newd()
  !
  ! ... and recalculate the products of the S with the atomic wfcs used 
  ! ... in LDA+U calculations
  !
  IF ( lda_plus_u .OR. use_wannier ) CALL orthoatwfc()
  !
  call tag_wg_corr_as_obsolete
  !
  RETURN
  !
END SUBROUTINE hinit1

