!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the  acfdt program
!
!
MODULE acfdtest
 USE kinds
 SAVE
 LOGICAL :: acfdt_is_active=.FALSE.
 LOGICAL :: acfdt_num_der=.FALSE.
 LOGICAL :: acfdt_term1 = .FALSE.
 LOGICAL :: acfdt_term2 = .FALSE.
 LOGICAL :: acfdt_term3 = .FALSE.
 LOGICAL :: test_oep=.FALSE. 
 INTEGER :: ir_point=0
 REAL(DP):: delta_vrs=0.0_DP
 REAL(DP):: f1=1.0_DP
 REAL(DP):: f2=1.0_DP
 REAL(DP):: f3=1.0_DP
END MODULE acfdtest

