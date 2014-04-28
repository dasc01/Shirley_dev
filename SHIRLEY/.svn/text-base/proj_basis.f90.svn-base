! ----------------------------------------------------------------------
  subroutine proj_basis
! ----------------------------------------------------------------------

! This is a blatant rip-off of proj_basis in ../PP/projwfc.f90
! I don't want the fancy formatting output, I just want the matrix proj
! Note that proj has no k-point index since the basis is entirely at the
! Gamma-point

  USE io_global,  ONLY : stdout, ionode
  USE atom
  USE ions_base, ONLY : nat, ityp, atm, ntyp => nsp
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU, only : lda_plus_u, swfcatom
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symm_base, ONLY: nsym, irt
  USE wvfct
  use control_flags, only : gamma_only
  USE uspp, ONLY: nkb, vkb
  USE becmod,   ONLY: becp, rbecp
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc

  !use proj_shirley, only : init_proj, write_proj, nlmchi, irt, &
  !                         d1, d2, d3, atomic_proj_matrix
  use shirley_ham_input, only : nkgrid, ikgrid, ksplord
  use splines_module


  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, lm, nwfc,&
       nwfc1, lmax_wfc, is
  integer :: ikb
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ... 
  REAL   (DP), ALLOCATABLE ::roverlap(:,:),          rwork1(:),rproj0(:,:)
  ! ... or for gamma-point. 

  complex(dp),allocatable :: qbp(:,:)
  real(dp),allocatable :: rqbp(:,:)
  ! 
  integer :: iunproj
  logical :: exst
  integer,external :: freeunit

  TYPE wfc_label
     INTEGER na, n, l, m
  END TYPE wfc_label
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)

  real(dp) :: d1(3,3,48), d2(5,5,48), d3(7,7,48)

  type(spline_dataset_type) :: proj_dataset
  type(spline_fit_type) :: proj_fit

  integer :: nkr
  real(dp) :: xk_cart(3)
  

  if( ionode ) then
    iunproj = freeunit()
    call seqopn(iunproj,'proj','unformatted',exst)
    write(stdout,*)
    write(stdout,*) ' dumping atomic angular momentum projections to file ', trim(prefix) // '.proj'
  endif
  ! 
  ! 
  WRITE( stdout, '(/5x,"Calling proj_basis .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  END IF
  ! 
  !! store dimensions in proj_shirley
  !call init_proj( natomwfc, nkb, nbnd, nat, nsym )
  !irt=irt_
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1 
  ! 
  CALL d_matrix (d1, d2, d3)
  ! 
  ! fill structure nlmchi 
  ! 
  ALLOCATE (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, nchi (nt)
        IF (oc (n, nt) >= 0.d0) THEN
           l = lchi (n, nt)
           lmax_wfc = MAX (lmax_wfc, l )
           DO m = 1, 2 * l + 1
              nwfc=nwfc+1
              nlmchi(nwfc)%na = na
              nlmchi(nwfc)%n  =  n
              nlmchi(nwfc)%l  =  l
              nlmchi(nwfc)%m  =  m
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  ! 
  IF (lmax_wfc > 3) CALL errore ('proj_basis', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('proj_basis', 'wrong # of atomic wfcs?', 1)
  ! 
  IF (.NOT. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) )
  ALLOCATE(wfcatom (npwx, natomwfc) )
  ALLOCATE(overlap (natomwfc, natomwfc) )
  overlap= (0.d0,0.d0)
  IF ( gamma_only ) THEN
     ALLOCATE(roverlap (natomwfc, natomwfc) )
     roverlap= 0.d0
     ALLOCATE (rbecp (nkb,natomwfc))
  ELSE
     ALLOCATE ( becp (nkb,natomwfc))
  END IF
  ALLOCATE(e (natomwfc) )

  ! Only one k-point - Gamma
  CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  ! load eigenfunctions
  CALL davcio (evc, 2*nwordwfc, iunwfc, 1, - 1)

  ! allocate space for spline data
  call create_spline_dataset( natomwfc, nbnd, nkgrid, ikgrid, &
                              bg*tpiba, transpose(at)/tpiba, &
                              proj_dataset )
  ! allocate space for the spline fit
  call create_spline_fit( natomwfc, nbnd, ksplord, proj_dataset, proj_fit )

  ! k-point grid - uniform only
  write(stdout,*) ' generate atomic projectors on the following k-point grid:'
  write(stdout,'(2(2x,a,3i6))') 'nkgrid =', nkgrid(1:3), 'ikgrid =', ikgrid(1:3)

  nkr = numgridcoord_spline_dataset( proj_dataset )

  do ik=1,nkr

     call getcartcoord_spline_dataset( ik, xk_cart, proj_dataset )
     xk_cart = xk_cart / tpiba
 
     write(stdout,*) ' load atomic projections for ik = ', ik, ' of ', nkr
     write(stdout,'(x,a,3f12.5)') ' xk_cart = ', xk_cart(1:3)
     call flush_unit(stdout)
 
     CALL atomic_wfc_shirley (npw, igk, xk_cart, wfcatom)
 
     CALL init_us_2_shirley (npw, igk, xk_cart, vkb)
  
     ! construct <beta_i|phi_j>
     IF ( gamma_only ) THEN
        CALL pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, wfcatom, npwx, rbecp, nkb)
     ELSE
        CALL ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom)
     END IF

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j> 
     ! 
      IF ( gamma_only ) THEN
        CALL pw_gemm ('Y', natomwfc, natomwfc, npw, wfcatom, npwx, swfcatom, &
             npwx, roverlap, natomwfc)
        overlap(:,:)=CMPLX(roverlap(:,:))
        ! TEMP: diagonalization routine for real matrix should be used instead 
     ELSE
        CALL ccalbec (natomwfc, npwx, npw, natomwfc, overlap, wfcatom, swfcatom) 
     END IF
     ! 
     ! calculate O^{-1/2} 
     ! 
     ALLOCATE(work (natomwfc, natomwfc) )
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO
     DO i = 1, natomwfc
        DO j = i, natomwfc
           overlap (i, j) = (0.d0, 0.d0)
           DO k = 1, natomwfc
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * CONJG (work (i, k) )
           ENDDO
           IF (j /= i) overlap (j, i) = CONJG (overlap (i, j))
        ENDDO
     ENDDO
     DEALLOCATE (work)
     ! 
     ! calculate wfcatom = O^{-1/2} \hat S | phi> 
     ! 
     IF ( gamma_only ) THEN
        roverlap(:,:)=REAL(overlap(:,:),DP)
        ! TEMP: diagonalization routine for real matrix should be used instead 
        CALL DGEMM ('n', 't', 2*npw, natomwfc, natomwfc, 1.d0 , &
             swfcatom, 2*npwx,  roverlap, natomwfc, 0.d0, wfcatom, 2*npwx)
     ELSE
        CALL ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , &
             swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx)
     END IF

     ! 
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j> 
     ! 
     IF ( gamma_only ) THEN
        ALLOCATE(rproj0(natomwfc,nbnd), rwork1 (nbnd) )
        CALL pw_gemm ('Y', natomwfc, nbnd, npw, wfcatom, npwx, evc, npwx, rproj0, natomwfc)
        ALLOCATE(proj0(natomwfc,nbnd) )
        proj0 = cmplx( rproj0 )
     ELSE
        ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) )
        CALL ccalbec (natomwfc, npwx, npw, nbnd, proj0, wfcatom, evc)
     END IF

     ! now dump proj into spline_dataset
     do ibnd=1,nbnd
     do nwfc=1,natomwfc
       call putelement_spline_dataset( nwfc, ibnd, ik, proj0(nwfc,ibnd), proj_dataset )
     enddo
     enddo

     IF ( gamma_only ) THEN
        DEALLOCATE (rwork1)
        DEALLOCATE (rproj0)
        DEALLOCATE (proj0)
     ELSE
        DEALLOCATE (work1)
        DEALLOCATE (proj0)
     END IF

     ! on k-points 
  ENDDO
  !
  DEALLOCATE (e)
  IF ( gamma_only ) THEN
     DEALLOCATE (roverlap)
     DEALLOCATE (rbecp)
  ELSE
     DEALLOCATE ( becp)
  END IF
  DEALLOCATE (overlap)
  DEALLOCATE (wfcatom)
  IF (.NOT. lda_plus_u) DEALLOCATE (swfcatom)

  ! fit the spline
  call fit_spline_to_dataset( proj_dataset, proj_fit )

  ! dump to disk
  if( ionode ) then
    call dump_spline_dataset( proj_dataset )

    call write_spline_fit( iunproj, proj_fit )
  endif
  !
  end subroutine proj_basis


     !-----------------------------------------------------------------------
     SUBROUTINE qbetaphi_gamma( m, ps )
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : vkb, nkb, qq, okvan
  USE uspp_param, ONLY : nh, tvanp
  USE basis,     ONLY : natomwfc
  USE wvfct,      ONLY : igk, g2kin
  use gvecs,      only : nls
  USE smooth_grid_dimensions,    ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE ldaU,       ONLY : lda_plus_u
  USE ions_base,  ONLY : nat, nsp, ityp
       USE becmod, ONLY : rbecp
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       integer,intent(in) :: m
       REAL(DP),intent(out) :: ps(nkb,m)
         ! the product vkb and phi
       !
       !
       ! ... The product with the beta functions
       !
       IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
       !
       ps(:,:) = 0.D0
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( tvanp (nt) ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd = 1, m
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                           qq(ih,jh,nt) * rbecp(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          END IF
       END DO
       !
       RETURN
       !
     END SUBROUTINE qbetaphi_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE qbetaphi_k( m, ps )
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : vkb, nkb, qq, okvan
  USE uspp_param, ONLY : nh, tvanp
  USE basis,     ONLY : natomwfc
  USE wvfct,      ONLY : igk, g2kin
  use gvecs,      only : nls
  USE smooth_grid_dimensions,    ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE ldaU,       ONLY : lda_plus_u
  USE ions_base,  ONLY : nat, nsp, ityp
       USE becmod,  ONLY : becp
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       integer,intent(in) :: m
       COMPLEX(DP),intent(out) :: ps(nkb,m)
         ! the product vkb and psi
       !
       !
       ! ... The product with the beta functions
       !
       IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
       !
       ps(:,:) = ( 0.D0, 0.D0 )
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( tvanp(nt) ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd = 1, m
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                           qq(ih,jh,nt) * becp(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          END IF
       END DO
       !
       RETURN
       !
     END SUBROUTINE qbetaphi_k     
