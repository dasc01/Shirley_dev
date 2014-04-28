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
SUBROUTINE rhoG()
  !----------------------------------------------------------------------------
  !
  ! ... This routine dumps the charge density fourier coefficients
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE cell_base,        ONLY : alat, omega, tpiba, tpiba2, at, bg
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE basis,            ONLY : startingpot
  USE klist,            ONLY : nelec
  USE klist,            ONLY : xk, ngk
  USE lsda_mod,         ONLY : lsda, nspin
  USE gvect,            ONLY : ngm, gstart, nl, g, gg, ig_l2g
  use grid_dimensions,  only : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE gvecs,            ONLY : doublegrid
  use wvfct,     only : ectuwfc, npw, igk, g2kin, igk_l2g, npwx
  USE control_flags,    ONLY : lscf
  USE scf,              ONLY : rho, rho_core, vltot, vr, vrs
  USE ener,             ONLY : ehart, etxc, vtxc
  USE ldaU,             ONLY : niter_with_fixed_ns
  USE ldaU,             ONLY : lda_plus_u, Hubbard_lmax, ns, nsnew
  USE noncollin_module, ONLY : noncolin, factlist, pointlist, pointnum, &
                               mcons, i_cons, lambda, vtcon, report
  USE io_files,         ONLY : prefix, iunocc, input_drho
  USE spin_orb,         ONLY : domag
  USE mp,               ONLY : mp_bcast, mp_sum, mp_max
  USE mp_global,        ONLY : intra_image_comm, me_pool, nproc_pool, &
                               root_pool, intra_pool_comm
  USE mp_wave,   ONLY : mergewf
  USE io_global,        ONLY : ionode, ionode_id
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  REAL (DP) :: charge           ! the starting charge
  REAL (DP) :: etotefield       ! 
  INTEGER        :: ios
  INTEGER        :: ldim             ! integer variable for I/O control
  LOGICAL        :: exst 
  REAL (DP), ALLOCATABLE :: aux(:,:)
  INTEGER                     :: ig, ixyz
  integer :: ngk_g, igwx
  real(dp), allocatable :: g_g(:,:), glatt_g(:,:)
  complex(dp),allocatable :: rho_g(:), jtmp(:)
  integer :: iunrhoG

  integer,external :: freeunit
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
  IF ( startingpot == 'file' .AND. exst ) THEN
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
  ELSE
     !   
     call errore( 'rhoG','density not found in .rho file',1)
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
  ! ... convert charge density to fourier space and merge components
  !
  !
  ALLOCATE( aux( 2, nrxx ) )
  !
  ! ... copy total rho in aux
  !
  aux(2,:) = 0.D0
  aux(1,:) = rho(:,1)
  !
  IF ( nspin == 2 ) aux(1,:) = aux(1,:) + rho(:,2)
  !
  ! ... bring rho (aux) to G space
  !
  CALL cft3( aux, nr1, nr2, nr3, nr1x, nr2x, nr3x, -1 )
  !
  ! ----------------------------------------------------------------------
  ! now merge according to the wave function ordering
  ! ----------------------------------------------------------------------
  !
  ! determine the local number of plane-waves less then pwmtxel_cutoff
  call n_plane_waves (ecutwfc, tpiba2, 1, xk, g, ngm, npwx, ngk)

  ! find the global number
  ngk_g = ngk(1)
  call mp_sum( ngk_g )

  write(stdout,*) ' nkg_g = ', ngk_g
  allocate( g_g(3,ngk_g), rho_g(ngk_g) )

  ! generate sorted local list
  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, ngk(1), igk, g2kin )
  call gk_l2gmap (ngm, ig_l2g(1), ngk(1), igk, igk_l2g(1,1))

  ! merge vectors
  igwx = maxval( igk_l2g(1:ngk(1),1) )
  call mp_max( igwx )
  allocate( jtmp(npwx), rho_g(igwx) )
  do ixyz=1,3
    jtmp(1:ngk(1)) = cmplx(g(ixyz,igk(1:ngk(1))),0.d0)
    rho_g = zero
    call mergewf(jtmp, rho_g, ngk(1), igk_l2g(:,1), &
                 me_pool, nproc_pool, root_pool, intra_pool_comm)
    g_g(ixyz,1:ngk_g) = rho_g(1:ngk_g)
  enddo

  ! convert G-vectors back to lattice coords
  glatt_g = matmul( transpose(at), g_g )

  ! merge rho magn
  jtmp(1:ngk(1)) = cmplx(aux(1,nl(1:ngk(1))),aux(2,nl(1:igk(1))))
  rho_g = zero
  call mergewf(jtmp, rho_g, ngk(1), igk_l2g(:,1), &
               me_pool, nproc_pool, root_pool, intra_pool_comm)

  if( ionode ) then
    iunrhoG=freeunit()
    open(iunrhoG,file=trim(prefix)//'.rhoG',form='unformatted')
    write(stdout,*) ' density coefficients saved to file: ', &
                    trim(prefix)//'.rhoG'
    write(iunrhoG) ecutwfc  ! cut-off
    write(iunrhoG) ngk_g    ! number of G vectors
    write(iunrhoG) at       ! transformation matrix of bravais lattice
    write(iunrhoG) bg       ! transformation matrix of reciprocal lattice
    write(iunrhoG) omega    ! cell volume
    write(iunrhoG) tpiba    ! 2 * pi / alat

    write(iunrhoG) g_g      ! the G-vectors in Cartesian coords
    write(iunrhoG) glatt_g  ! the G-vectors in Cartesian coords
    write(iunrhoG) rho_g(1:ngk_g) ! the density coefficients
    close(iunrhoG)
  endif
  !
  DEALLOCATE( aux )
  deallocate( jtmp, rho_g )
  !
  RETURN
  !
END SUBROUTINE rhoG
