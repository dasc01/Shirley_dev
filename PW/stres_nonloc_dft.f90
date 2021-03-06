!
! Copyright (C) 2010- Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
subroutine stres_nonloc_dft( rho, rho_core, nspin, sigma_nonloc_dft )

  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  use funct,            ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
                               dft_is_gradient, get_igcc, dft_is_vdW
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE grid_dimensions,  ONLY : nrxx
  USE vdW_DF,           ONLY : stress_vdW_DF, print_sigma 
  !
  IMPLICIT NONE
  !
  real(DP), intent(in)     :: rho (nrxx, nspin), rho_core (nrxx)
  real(DP), intent(inout)  :: sigma_nonloc_dft (3, 3)
  integer ::nspin

  integer :: l, m


  sigma_nonloc_dft(:,:) = 0.d0
  
  if (dft_is_vdW()) then
     if (nspin>1) call errore('stres_vdW_DF', &
                  'vdW+DF non implemented in spin polarized calculations',1)
     CALL stress_vdW_DF(rho, rho_core, sigma_nonloc_dft)

  end if

  return

end subroutine stres_nonloc_dft

