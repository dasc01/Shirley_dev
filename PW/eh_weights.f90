!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine eh_weights (nks, wk, nbnd, nelec, nholes, ncorex, eg_min, e_excit, degauss, ngauss, &
     et, ef, demet, netot, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss, nholes, ncorex, eg_min, e_excit
  real(DP), intent(out) :: wg (nbnd, nks), ef, demet, netot
  integer, intent(in) :: is, isk(nks)
  !
  integer :: kpoint, ibnd
  real(DP) , external :: wgauss, w1gauss, efermig
  real(DP):: nvbel,ncbel,efvb,efcb, eflu,efgs,dwg(nbnd,nks),scal,temp_tot,kwgt,gaus
  integer :: homo(nks), lumo(nks)
  logical:: excit
  ! Calculate the Fermi energy ef
  ef = efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)

  efgs=efermig (et, nbnd, nks, nelec-ncorex, wk, degauss, ngauss, is, isk)
  efvb=efermig (et, nbnd, nks, nelec-ncorex-nholes, wk, degauss, ngauss, is, isk)
  eflu=efermig (et, nbnd, nks, nelec-ncorex+nholes, wk, degauss, ngauss, is, isk)
  efcb=efermig (et, nbnd, nks, nelec+nholes, wk, degauss, ngauss, is, isk)
  excit=.false.
  if(e_excit .gt. efcb-efvb) then
     excit=.true.
     efcb=efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)
  endif

  write(*,'(A,5f12.6)') "efvb,efgs,ef,eflu,efcb=",efvb*13.6056925,efgs*13.6056925,ef*13.6056925,eflu*13.6056925,efcb*13.6056925

  do  kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        if(et(ibnd,kpoint) .lt. efgs+eg_min) then
           homo(kpoint) = ibnd
        else
           exit
        endif
     enddo
     do ibnd = nbnd, 1, -1
        if(et(ibnd,kpoint) .gt. efgs+eg_min) then
           lumo(kpoint) = ibnd
        else
           exit
        endif
     enddo
     write(*,'(A,3i5)') "ik,homo,lumo=",kpoint,homo(kpoint),lumo(kpoint)
  enddo
  if(excit) dwg(:,:)=0.0d0
  netot = 0.0d0
  demet = 0.d0
  temp_tot=0.0d0
  kwgt=0.0d0
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
        netot = netot + wg (ibnd, kpoint)
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
        netot = netot + wg (ibnd, kpoint)
     enddo
     if(excit) then
        kwgt=kwgt+wk (kpoint)
        do ibnd = lumo(kpoint), nbnd
!!$           write(*,'(A,3F12.7)') 'et(ibnd,kpoint),(efvb+e_excit),diff=',&
!!$                &et(ibnd,kpoint)*13.6056925,(efvb+e_excit)*13.6056925,(et(ibnd,kpoint)-(efvb+e_excit))*13.6056925
           dwg(ibnd, kpoint)= wk (kpoint) * gaus(et(ibnd,kpoint)-(efvb+e_excit), max(degauss*2.0d0,0.073498648))
           temp_tot=temp_tot+dwg(ibnd, kpoint)
        enddo
     endif
  enddo
  if(excit) then
     scal=(nholes*kwgt*0.5d0)/(temp_tot+1.0d-16)
     write(*,'(A,3F12.7)') "temp_tot,scal,kwgt=",temp_tot,scal,kwgt
     do kpoint = 1, nks
        if (is /= 0) then
           if (isk(kpoint).ne.is) cycle
        end if
        do ibnd = lumo(kpoint), nbnd
           wg(ibnd, kpoint)=wg(ibnd, kpoint)+dwg(ibnd, kpoint)*scal
           netot = netot + dwg(ibnd, kpoint)*scal
        enddo
     enddo
  endif
  return
end subroutine eh_weights

function gaus(de,dg)
  USE kinds
  implicit none
  real(DP)::de,dg,gaus,pi

  pi=4.0d0*atan(1.0d0)

  gaus=(1/(dg*sqrt(2*pi)))*exp((-1.0d0/2)*((de/dg)**2))

  return
end function gaus


!--------------------------------------------------------------------
subroutine e_last_elec(nk, nbnd, nelec, nspin, eig, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  USE constants, ONLY: rytoev
  implicit none
  !
  integer, intent(in) :: nk, nbnd, nspin
  real(DP), intent(in) :: eig (nbnd, nk), nelec
  real(DP), intent(out) :: wg (nbnd, nk) 
  integer, intent(in) :: is, isk(nk)
  !
  integer :: kpoint, ibnd, ispin, pk, qk, nbk, numk
  real(DP):: etmp, sumw, avge, dw
  real(DP),allocatable::eig_srtd(:,:),wg_srtd(:,:),ik(:)
  logical:: break

  numk=nk/nspin
  ! Calculate the Fermi energy ef
  allocate(eig_srtd(nbnd*numk,nspin),wg_srtd(nbnd*numk,nspin),ik(nspin))
  nbk=nbnd*numk
  ik(:)=1
  do ibnd=nbnd,1,-1
     do  kpoint = 1, nk
        if(isk(kpoint) .eq. 0 .or. isk(kpoint) .eq. 1) then
           eig_srtd(ik(1),1)=eig(ibnd,kpoint)
           wg_srtd(ik(1),1)=wg(ibnd,kpoint)
           ik(1)=ik(1)+1
        else if(isk(kpoint) .eq. 2) then
           eig_srtd(ik(2),2)=eig(ibnd,kpoint)
           wg_srtd(ik(2),2)=wg(ibnd,kpoint)
           ik(2)=ik(2)+1
        endif
     enddo
  enddo
  
  do ispin=1, nspin
     do pk=nbk,1,-1
        do qk=1, pk-1
           if(eig_srtd(qk,ispin) .lt. eig_srtd(qk+1,ispin)) then
              etmp=eig_srtd(qk+1,ispin)
              eig_srtd(qk+1,ispin)=eig_srtd(qk,ispin)
              eig_srtd(qk,ispin)=etmp
              
              etmp=wg_srtd(qk+1,ispin)
              wg_srtd(qk+1,ispin)=wg_srtd(qk,ispin)
              wg_srtd(qk,ispin)=etmp
           endif
        enddo
     enddo
  enddo

  do ispin=1, nspin
     sumw=0.0d0
     avge=0.0d0
     pk=1
     break=.false.
     do while(sumw <= 1.0d0)
        if(sumw+wg_srtd(pk,ispin) .gt. 1.0d0) then
           dw=(1.0d0-sumw)
           break=.true.
        else
           dw=wg_srtd(pk,ispin)
        endif
        avge=avge+eig_srtd(pk,ispin)*dw
        sumw=sumw+dw
        if(break) exit
        pk=pk+1
     enddo
     avge=avge/sumw
     write(*,'(A,I4,F16.10,A)') "     Energy of last electron with spin", ispin,avge*rytoev," eV"  
  enddo
  
  return
end subroutine e_last_elec
