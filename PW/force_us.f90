!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE force_us( forcenl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  ! ... wrapper
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq, deeq, qq_so, deeq_nc
  USE uspp_param,           ONLY : upf, nh, newpseudo, nhm
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : bec_type, becp, allocate_bec_type, deallocate_bec_type, becunit !DAS
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm, mpime, root !DAS
  USE mp,                   ONLY : mp_sum
  USE io_files, ONLY:find_free_unit !DAS
  !
  IMPLICIT NONE
  !
  ! ... the dummy variable
  !
  REAL(DP) :: forcenl(3,nat)
  ! output: the nonlocal contribution
  !
  CALL allocate_bec_type ( nkb, nbnd, becp )   
  !
  !DASb
  if(mpime==root) then
     becunit=find_free_unit()
     open(becunit,file='becp.dmp',form='formatted',status='replace')
     write(becunit,'(1I12)') nks
  endif
  !DASe
  !
  IF ( gamma_only ) THEN
     !
     CALL force_us_gamma( forcenl )
     !
  ELSE
     !
     CALL force_us_k( forcenl )
     !
  END IF  
  !
  CALL deallocate_bec_type ( becp )   
  !
  !DASb
  if(mpime==root) then
     close(becunit)     
  endif
  !DASe
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma( forcenl )
       !-----------------------------------------------------------------------
       !
       ! ... calculation at gamma
       !
       USE becmod, ONLY : calbec, dump_becp
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       REAL(DP), ALLOCATABLE    :: rdbecp (:,:,:)
       ! auxiliary variable, contains <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
       ! auxiliary variable contains g*|beta>
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       !
       ALLOCATE( rdbecp( nkb, nbnd, 3 ) )    
       ALLOCATE( vkb1(  npwx, nkb ) ) 
       !   
       IF ( nks > 1 ) REWIND iunigk
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk (ik)
          IF ( nks > 1 ) THEN
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec ( npw, vkb, evc, becp )
          !
          DO ipol = 1, 3
             DO jkb = 1, nkb
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk(ig))
                END DO
             END DO
             !
             CALL calbec ( npw, vkb1, evc, rdbecp(:,:,ipol) )
             !
          END DO
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         DO ipol = 1, 3
                            forcenl(ipol,na) = forcenl(ipol,na) - &
                                       ps * wg(ibnd,ik) * 2.D0 * tpiba * &
                                       rdbecp(ikb,ibnd,ipol) *becp%r(ikb,ibnd)
                         END DO
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq(ih,jh,nt)
                               DO ipol = 1, 3
                                  forcenl(ipol,na) = forcenl(ipol,na) - &
                                     ps * wg(ibnd,ik) * 2.d0 * tpiba * &
                                     (rdbecp(ikb,ibnd,ipol) *becp%r(jkb,ibnd) + &
                                      rdbecp(jkb,ibnd,ipol) *becp%r(ikb,ibnd) )
                               END DO
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
!
!DASb
          call dump_becp(becp,ik,current_spin)
!DASe
       END DO
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl )
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces
       !
       CALL symvector ( nat, forcenl )
       !
       DEALLOCATE( vkb1 )
       DEALLOCATE(rdbecp ) 
       !
       RETURN
       !
     END SUBROUTINE force_us_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k( forcenl )
       !-----------------------------------------------------------------------
       !  
       USE becmod, ONLY : calbec, dump_becp
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       COMPLEX(DP), ALLOCATABLE :: dbecp(:,:,:), dbecp_nc(:,:,:,:)
       ! auxiliary variable contains <beta|psi> and <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
       ! auxiliary variable contains g*|beta>
       COMPLEX(DP) :: psc(2,2), fac
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0, &
                        is, js, ijs
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       !
       IF (noncolin) then
          ALLOCATE( dbecp_nc(nkb,npol,nbnd,3) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       ELSE
          ALLOCATE( dbecp( nkb, nbnd, 3 ) )    
          ALLOCATE( deff(nhm,nhm,nat) )
       ENDIF
       ALLOCATE( vkb1( npwx, nkb ) )   
       ! 
       IF ( nks > 1 ) REWIND iunigk
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          IF ( nks > 1 ) THEN
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec ( npw, vkb, evc, becp)
          !
          DO ipol = 1, 3
             DO jkb = 1, nkb
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk(ig))
                END DO
             END DO
             !
             IF (noncolin) THEN
                IF ( nkb > 0 ) &
                   CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw, ( 1.D0, 0.D0 ),&
                            vkb1, npwx, evc, npwx, ( 0.D0, 0.D0 ),    &
                            dbecp_nc(1,1,1,ipol), nkb )
             ELSE
                IF ( nkb > 0 ) &
                   CALL ZGEMM( 'C', 'N', nkb, nbnd, npw, ( 1.D0, 0.D0 ),   &
                            vkb1, npwx, evc, npwx, ( 0.D0, 0.D0 ),      &
                            dbecp(1,1,ipol), nkb )
             END IF
          END DO
          !
          DO ibnd = 1, nbnd
             IF (noncolin) THEN
                CALL compute_deff_nc(deff_nc,et(ibnd,ik))
             ELSE
                CALL compute_deff(deff,et(ibnd,ik))
             ENDIF
             fac=wg(ibnd,ik)*tpiba
             ijkb0 = 0
             DO nt = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      DO ih = 1, nh(nt)
                         ikb = ijkb0 + ih
                         IF (noncolin) THEN
                            DO ipol=1,3
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     forcenl(ipol,na) = forcenl(ipol,na)- &
                                         deff_nc(ih,ih,na,ijs)*fac*( &
                                         CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                         becp%nc(ikb,js,ibnd)+ &
                                         CONJG(becp%nc(ikb,is,ibnd))* &
                                         dbecp_nc(ikb,js,ibnd,ipol) )
                                  END DO
                               END DO
                            END DO
                         ELSE
                            DO ipol=1,3
                               forcenl(ipol,na) = forcenl(ipol,na) - &
                                  2.D0 * fac * deff(ih,ih,na)*&
                                      DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                            becp%k(ikb,ibnd) )
                            END DO
                         END IF
                         !
                         IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                            DO jh = ( ih + 1 ), nh(nt)
                               jkb = ijkb0 + jh
                               IF (noncolin) THEN
                                  DO ipol=1,3
                                     ijs=0
                                     DO is=1,npol
                                        DO js=1,npol
                                           ijs=ijs+1
                                           forcenl(ipol,na)=forcenl(ipol,na)- &
                                           deff_nc(ih,jh,na,ijs)*fac*( &
                                          CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                                 becp%nc(jkb,js,ibnd)+ &
                                          CONJG(becp%nc(ikb,is,ibnd))* &
                                                dbecp_nc(jkb,js,ibnd,ipol))- &
                                           deff_nc(jh,ih,na,ijs)*fac*( &
                                          CONJG(dbecp_nc(jkb,is,ibnd,ipol))* &
                                                becp%nc(ikb,js,ibnd)+ &
                                          CONJG(becp%nc(jkb,is,ibnd))* &
                                                dbecp_nc(ikb,js,ibnd,ipol) )
                                        END DO
                                     END DO
                                  END DO
                               ELSE
                                  DO ipol = 1, 3
                                     forcenl(ipol,na) = forcenl (ipol,na) - &
                                          2.D0 * fac * deff(ih,jh,na)* &
                                       DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                             becp%k(jkb,ibnd) +       &
                                             dbecp(jkb,ibnd,ipol) * &
                                             CONJG( becp%k(ikb,ibnd) ) )
                                  END DO
                               END IF
                            END DO !jh
                         END IF ! tvanp
                      END DO ! ih = 1, nh(nt)
                      ijkb0 = ijkb0 + nh(nt)
                   END IF ! ityp(na) == nt
                END DO ! nat
             END DO ! ntyp
          END DO ! nbnd
!
!DASb
          call dump_becp(becp,ik,current_spin)
!DASe
       END DO ! nks
       !
#ifdef __PARA
       CALL mp_sum(  forcenl , intra_pool_comm )
#endif
       !
       DEALLOCATE( vkb1 )
       IF (noncolin) THEN
          DEALLOCATE( dbecp_nc )
          DEALLOCATE( deff_nc )
       ELSE
          DEALLOCATE( dbecp )
          DEALLOCATE( deff )
       ENDIF
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl )
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces.
       !
       CALL symvector ( nat, forcenl )
       !
       RETURN
       !
     END SUBROUTINE force_us_k
     !     
END SUBROUTINE force_us
