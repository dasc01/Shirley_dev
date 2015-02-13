!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine def_weights (nks, wk, nbnd, nelec, nspin, nholes, ncorex, eg_min, e_excit, degauss, ngauss, &
     et, ef, demet, netot, wg, elaste, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss, nspin
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss, nholes, ncorex, eg_min, e_excit
  real(DP), intent(out) :: wg (nbnd, nks), ef, demet, netot, elaste
  integer, intent(in) :: is, isk(nks)
  !
  integer :: kpoint, ibnd, numk
  real(DP) , external :: wgauss, w1gauss
  real(DP):: nvbel,ncbel,efvb,efcb, efgs,dwg(nbnd,nks),scal,temp_tot,kwgt,gaus, wkk(nks), efermig2,demet_k
  integer :: homo(nks), lumo(nks)
  logical:: excit
  real(DP)::efk(nks),efgsk(nks),efvbk(nks),efcbk(nks),temp
  real(DP), parameter :: eps= 1.0d-8
  ! Calculate the Fermi energy ef
  wkk(:)=2.0d0/nspin  !currently works only for non-spin polarized systems because only nelec is specified
  demet=0.0d0
  netot=0.0d0
  ef=0.0d0
  numk=0
  do  kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     efk(kpoint) = efermig2 (et(1,kpoint), nbnd, 1, nelec, wkk(kpoint), degauss, ngauss, is, isk(kpoint))
     efgsk(kpoint) = efermig2 (et(1,kpoint), nbnd, 1, nelec-ncorex, wkk(kpoint), degauss, ngauss, is, isk(kpoint))
     efvbk(kpoint) = efermig2 (et(1,kpoint), nbnd, 1, nelec-ncorex-nholes, wkk(kpoint), degauss, ngauss, is, isk(kpoint))
     efcbk(kpoint) = efermig2 (et(1,kpoint), nbnd, 1, nelec+nholes, wkk(kpoint), degauss, ngauss, is, isk(kpoint))
     nvbel=0.0d0
     demet_k=0.0d0
     do ibnd = 1, nbnd
        ! Calculate the gaussian weights
        wg (ibnd, kpoint) = wkk (kpoint) * &
             wgauss ( (efvbk(kpoint)-et(ibnd,kpoint)) / degauss, ngauss)
        !
        ! The correct (i.e. variational) form of the band energy is 
        !    Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        ! This differs by the term "demet" from the sum of KS eigenvalues:
        !    Eks = \sum wg(n,k) et(n,k)
        ! which is non variational. When a Fermi-Dirac function is used
        ! for a given T, the variational energy is really the free energy F,
        ! and F = E - TS , with E = non variational energy, -TS = demet
        !
        demet_k = demet_k + wkk (kpoint) * &
             degauss * w1gauss ( (efvbk(kpoint)-et(ibnd,kpoint)) / degauss, ngauss)
        nvbel = nvbel + wg (ibnd, kpoint)
     enddo
     ncbel=0.0d0
     do ibnd=nbnd, 1, -1
        temp=wkk (kpoint) * wgauss ( (efcbk(kpoint)-et(ibnd,kpoint)) / degauss, ngauss)
        if((nvbel+ncbel+temp).le. nelec+eps) then
           wg (ibnd, kpoint) = temp
           demet_k = demet_k + wkk (kpoint) * &
                degauss * w1gauss ( (efcbk(kpoint)-et(ibnd,kpoint)) / degauss, ngauss)
           ncbel = ncbel + wg (ibnd, kpoint)
        else
           exit
        endif
     enddo
     netot=netot+(nvbel+ncbel)*wk(kpoint)/wkk(kpoint)
     wg(:, kpoint)=wg(:, kpoint)*wk(kpoint)/wkk(kpoint)
     demet_k=demet_k*wk(kpoint)/wkk(kpoint)
     demet=demet+demet_k
     numk=numk+1
     ef=ef+efk(kpoint)
!     write(*,'(A,I3,4F14.6)') "ik, netot, nvbel, ncbel, demet_k=",kpoint, netot, nvbel, ncbel, demet_k
  enddo
  if(numk .gt. 0) ef=ef/numk   !average Fermi energy over all k-points

!Calculate average energy of last electron
  elaste=0.0d0
  do  kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     temp_tot=0.0d0
     do ibnd=nbnd, 1, -1
        temp=wkk(kpoint) * wgauss ( (efcbk(kpoint)-et(ibnd,kpoint)) / degauss, ngauss)
        if((temp_tot+temp).le. 1.0d0+eps) then
           temp_tot = temp_tot + temp
           elaste=elaste+temp*wk(kpoint)/wkk(kpoint)*et(ibnd,kpoint)
        else
           temp=1.0d0-temp_tot
           elaste=elaste+temp*wk(kpoint)/wkk(kpoint)*et(ibnd,kpoint)
           exit
        endif
     enddo
  enddo

  return
end subroutine def_weights

!--------------------------------------------------------------------
FUNCTION efermig2 (et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - Gaussian Broadening
  !     (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !
  USE io_global, ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY: rytoev
  USE mp, ONLY : mp_max, mp_min
  USE mp_global, ONLY : inter_pool_comm
  implicit none
  !  I/O variables
  integer, intent(in) :: nks, nbnd, Ngauss, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), Degauss, nelec
  real(DP) :: efermig2
  !
  real(DP), parameter :: eps= 1.0d-10
  integer, parameter :: maxiter = 300
  ! internal variables
  real(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  real(DP)::  sumkg2
  integer :: i, kpoint
  !
  !      find bounds for the Fermi energy. Very safe choice!
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  do kpoint = 2, nks
     Elw = min (Elw, et (1, kpoint) )
     Eup = max (Eup, et (nbnd, kpoint) )
  enddo
  Eup = Eup + 2 * Degauss
  Elw = Elw - 2 * Degauss
#ifdef __PARA
  !
  ! find min and max across pools
  !
  call mp_max( eup, inter_pool_comm )
  call mp_min( elw, inter_pool_comm )
#endif
  !
  !      Bisection method
  !
  sumkup = sumkg2 (et, nbnd, nks, wk, Degauss, Ngauss, Eup, is, isk)
  sumklw = sumkg2 (et, nbnd, nks, wk, Degauss, Ngauss, Elw, is, isk)
  if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       call errore ('efermig2', 'internal error, cannot bracket Ef', 1)
  do i = 1, maxiter
     Ef = (Eup + Elw) / 2.d0
     sumkmid = sumkg2 (et, nbnd, nks, wk, Degauss, Ngauss, Ef, is, isk)
     if (abs (sumkmid-nelec) < eps) then
        efermig2 = Ef
        return
     elseif ( (sumkmid-nelec) < -eps) then
        Elw = Ef
     else
        Eup = Ef
     endif
  enddo
  if (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig2 = Ef
  return
end FUNCTION efermig2

!-----------------------------------------------------------------------
function sumkg2 (et, nbnd, nks, wk, degauss, ngauss, e, is, isk)
  !-----------------------------------------------------------------------
  !
  !     This function computes the number of states under a given energy e
  !
  !
  USE kinds
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_sum
  implicit none
  ! Output variable
  real(DP) :: sumkg2
  ! Input variables
  integer, intent(in) :: nks, nbnd, ngauss
  ! input: the total number of K points
  ! input: the number of bands
  ! input: the type of smearing
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), degauss, e
  ! input: the weight of the k points
  ! input: the energy eigenvalues
  ! input: gaussian broadening
  ! input: the energy to check
  integer, intent(in) :: is, isk(nks)
  !
  ! local variables
  !
  real(DP), external :: wgauss
  ! function which compute the smearing
  real(DP) ::sum1
  integer :: ik, ibnd
  ! counter on k points
  ! counter on the band energy
  !
  sumkg2 = 0.d0
  do ik = 1, nks
     sum1 = 0.d0
     if (is /= 0) then
        if (isk(ik).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        sum1 = sum1 + wgauss ( (e-et (ibnd, ik) ) / degauss, ngauss)
     enddo
     sumkg2 = sumkg2 + wk (ik) * sum1
  enddo
!!$#ifdef __PARA
!!$  call mp_sum ( sumkg, inter_pool_comm )
!!$#endif
  return
end function sumkg2
