!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine eh_weights (nks, wk, nbnd, nelec, nholes, eg_min, degauss, ngauss, &
     et, ef, demet, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss, nholes, eg_min
  real(DP), intent(out) :: wg (nbnd, nks), ef, demet
  integer, intent(in) :: is, isk(nks)
  !
  integer :: kpoint, ibnd
  real(DP) , external :: wgauss, w1gauss, efermig
  real(DP):: nvbel,ncbel,efvb,efcb
  real(DP):: ecb(nbnd,nks),vbocc(nbnd,nks),cbocc(nbnd,nks)
  integer :: homo(nks), lumo(nks)
  ! Calculate the Fermi energy ef
  ef = efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)

  efvb=efermig (et, nbnd, nks, nelec-nholes, wk, degauss, ngauss, is, isk)
  efcb=efermig (et, nbnd, nks, nelec+nholes, wk, degauss, ngauss, is, isk)

  write(*,'(A,3f12.6)') "efvb,ef,efcb=",efvb*13.6056925,ef*13.6056925,efcb*13.6056925

  do  kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        if(et(ibnd,kpoint) .lt. ef+eg_min) then
           homo(kpoint) = ibnd
        else
           exit
        endif
     enddo
     do ibnd = nbnd, 1, -1
        if(et(ibnd,kpoint) .gt. ef+eg_min) then
           lumo(kpoint) = ibnd
        else
           exit
        endif
     enddo
     write(*,'(A,3i5)') "ik,homo,lumo=",kpoint,homo(kpoint),lumo(kpoint)
  enddo

  demet = 0.d0
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, homo(kpoint)
        ! Calculate the gaussian weights
        wg (ibnd, kpoint) = wk (kpoint) * &
                            wgauss ( (efvb-et(ibnd,kpoint)) / degauss, ngauss)
        !
        ! The correct (i.e. variational) form of the band energy is 
        !    Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        ! This differs by the term "demet" from the sum of KS eigenvalues:
        !    Eks = \sum wg(n,k) et(n,k)
        ! which is non variational. When a Fermi-Dirac function is used
        ! for a given T, the variational energy is really the free energy F,
        ! and F = E - TS , with E = non variational energy, -TS = demet
        !
        demet = demet + wk (kpoint) * &
                 degauss * w1gauss ( (efvb-et(ibnd,kpoint)) / degauss, ngauss)
     enddo
     do ibnd = lumo(kpoint), nbnd
        ! Calculate the gaussian weights
        wg (ibnd, kpoint) = wk (kpoint) * &
                            wgauss ( (efcb-et(ibnd,kpoint)) / degauss, ngauss)
        !
        ! The correct (i.e. variational) form of the band energy is 
        !    Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        ! This differs by the term "demet" from the sum of KS eigenvalues:
        !    Eks = \sum wg(n,k) et(n,k)
        ! which is non variational. When a Fermi-Dirac function is used
        ! for a given T, the variational energy is really the free energy F,
        ! and F = E - TS , with E = non variational energy, -TS = demet
        !
        demet = demet + wk (kpoint) * &
                 degauss * w1gauss ( (efcb-et(ibnd,kpoint)) / degauss, ngauss)
     enddo

  enddo
  return
end subroutine eh_weights
