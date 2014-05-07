!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
module proj_weights
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data produced by PWscf
  !
  ! ... written by Carlo Sbraccia (2005)
  !
  USE iotk_module
  USE xml_io_base
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, delete_if_present, &
                        qexml_version, qexml_version_init,find_free_unit
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global,            ONLY : kunit, nproc, nproc_pool, my_pool_id, me_pool, &
       root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm 
  USE mp,        ONLY : mp_bcast, mp_sum, mp_max
  USE parser,    ONLY : version_compare

  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: set_proj_weights,set_ext_weights
  !
  INTEGER, PRIVATE :: iunout
  COMPLEX(DP), ALLOCATABLE, TARGET :: auxevc(:,:),aux2(:,:)
!!$  USE basis,    ONLY : natomwfc
!!$  USE lsda_mod, ONLY : lsda, isk
!!$  USE klist,    ONLY : nkstot, wk, nelec
!!$  USE wvfct,    ONLY : et, wg, nbnd
!!$  USE ener,     ONLY : ef, ef_up, ef_dw
  REAL(DP)::ef_aux,ef_up_aux,ef_dw_aux,nelec_aux
  REAL(DP), ALLOCATABLE, SAVE:: et_aux(:,:),wg_aux(:,:)
  INTEGER::natomwfc_aux,nbnd_aux
  INTEGER, ALLOCATABLE,SAVE::isk_aux(:)
  
  !
contains

subroutine set_proj_weights
  USE kinds,             ONLY : DP
  USE control_flags,     ONLY : io_level, twfcollect
  USE noncollin_module,  ONLY : npol
  USE wvfct,             ONLY : nbnd, npwx
  USE klist,             ONLY : nks, nkstot, wk, xk, ngk
  USE io_files,          ONLY : nwordwfc, iunwfc, iunigk
  USE buffers2,          ONLY : open_buffer2, get_buffer2
  USE buffers,           ONLY : get_buffer
  USE input_parameters,  ONLY : read_extwfc, auxprefix
  USE fft_base,          ONLY : dffts
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, btype
  USE uspp, ONLY              : okvan
  USE funct,                ONLY : dft_is_meta
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  !
  IMPLICIT NONE
  !
  LOGICAL            :: exst
  LOGICAL, save      :: first=.true.
  INTEGER            :: ierr, incr, ik, ibnd, jk, jspin, jbnd
  TYPE (bec_type)    ::psipsi
  REAL(DP)           :: totel
  REAL(DP),allocatable::auxwg(:,:)

  if(.not. allocated(auxevc)) allocate(auxevc(SIZE(evc,1),SIZE(evc,2)),aux2(SIZE(evc,1),SIZE(evc,2)))
  if(.not. allocated(auxwg)) allocate(auxwg(nbnd,nkstot))

  if(first) then
     CALL open_buffer2( iunwfc, 'wfc2', nwordwfc, nks, exst )
     CALL read_wavefunctions2(trim(auxprefix), ierr)
     CALL  read_band_structure2( trim(auxprefix), ierr )
     first=.false.
  endif
  auxwg(:,:)=0

  !  RETURN
  !
  !
  !IF ( nks > 1 ) 
  REWIND( iunigk )
  !
  !
  incr = 1
  !
  if(dffts%have_task_groups) stop "Task groups not supported in proj_weights!!"
  if(okvan) stop "Wavefunction projection not supported with ultrasoft pseudo!"
  if(dft_is_meta()) stop "Meta ggas not supported on proj_weights!"
  totel=0.0d0
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     npw = ngk(ik)
     !
     !
     if(nks > 1) READ( iunigk ) igk
     CALL get_buffer2 ( auxevc, nwordwfc, iunwfc, ik )
     !
     CALL get_buffer ( aux2, nwordwfc, iunwfc, ik )

     CALL allocate_bec_type(nbnd, nbnd, psipsi)

     CALL calbec(npw, auxevc, aux2, psipsi)
     do ibnd=1, nbnd  !loop over target states

        do jbnd=1, nbnd !loop over source states
           auxwg(ibnd,ik)=auxwg(ibnd,ik)+wg(jbnd,ik)*psipsi%k(jbnd,ibnd)*CONJG(psipsi%k(jbnd,ibnd))
        enddo
!!$        write(*,'(A,2I6,F16.12)') "ik, ibnd, psipsi=", ik, ibnd, REAL(psipsi%k(ibnd,ibnd))
     enddo
     CALL deallocate_bec_type ( psipsi )

!check completeness of projection
     do ibnd=1, nbnd
       totel=totel+auxwg(ibnd,ik)
       write(*,'(A,2I6,F16.12)') "ik, ibnd, auxwg =", ik, ibnd, auxwg(ibnd, ik)
     enddo
     
!!$     do jk=1, nks
!!$
!!$        jspin = 1
!!$
!!$        IF ( lsda ) jspin = isk(jk)
!!$        
!!$        if(jspin .eq. current_spin) then
!!$           CALL get_buffer2 ( aux2, nwordwfc, iunwfc, jk )
!!$        endif
!!$        
!!$        CALL allocate_bec_type(nbnd, nbnd, psipsi)
!!$        CALL calbec(npw, auxevc, aux2, psipsi)
!!$        do ibnd=1, nbnd
!!$           write(*,'(A,3I6,F16.12)') "ik, jk, ibnd, psipsi=", ik, jk, ibnd, REAL(psipsi%k(ibnd,ibnd))
!!$        enddo
!!$        CALL deallocate_bec_type ( psipsi )
!!$        
!!$     enddo
     !
     DO ibnd = 1, nbnd, incr
        !
        IF( dffts%have_task_groups ) THEN
        ELSE
        END IF
        !
        IF (dft_is_meta()) THEN
        END IF
        !
        !
     END DO
     !
     IF( dffts%have_task_groups ) THEN
     END IF
     !
     ! ... If we have a US pseudopotential we compute here the becsum term
     !
     IF ( .NOT. okvan ) CYCLE k_loop
     !
     !
  END DO k_loop
  write(*,'(A,F12.8)') "totel = ", totel
  !
  
  !
end subroutine set_proj_weights

!------------------------------------------------------------------------
SUBROUTINE read_wavefunctions2( dirname, ierr )
  !------------------------------------------------------------------------
  !
  ! ... This routines reads wavefunctions from the new file format and
  ! ... writes them into the old format
  !
  USE control_flags,        ONLY : twfcollect, lkpoint_dir
  USE cell_base,            ONLY : tpiba2
  USE lsda_mod,             ONLY : nspin, isk
  USE klist,                ONLY : nkstot, wk, nelec, nks, xk, ngk
  USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, nbnd, ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE gvect,   ONLY : ig_l2g
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers2,             ONLY : save_buffer2
  USE gvect,                ONLY : ngm, ngm_g, g
  USE noncollin_module,     ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  INTEGER,          INTENT(OUT) :: ierr
  !
  CHARACTER(LEN=256)   :: filename
  INTEGER              :: ik, ipol, ik_eff, num_k_points
  INTEGER, ALLOCATABLE :: kisort(:)
  INTEGER              :: npool, nkbl, nkl, nkr, npwx_g
  INTEGER              :: ike, iks, npw_g, ispin
  INTEGER, ALLOCATABLE :: ngk_g(:)
  INTEGER, ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:)
  LOGICAL              :: opnd
  REAL(DP)             :: scalef
  !
  !
  IF ( iunwfc > 0 ) THEN
     !
     INQUIRE( UNIT = iunwfc, OPENED = opnd )
     !
     IF ( .NOT. opnd ) CALL errore( 'read_wavefunctions', &
          & 'wavefunctions unit (iunwfc) is not opened', 1 )
  END IF
  !
  IF ( ionode ) THEN
     !
     CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
          & TRIM( xmlpun ), IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, ionode_id, intra_image_comm )
  !
  IF ( ierr > 0 ) RETURN
  !
  IF ( nkstot > 0 ) THEN
     !
     ! ... find out the number of pools
     !
     npool = nproc / nproc_pool
     !
     ! ... find out number of k points blocks
     !
     nkbl = nkstot / kunit
     !
     !  k points per pool
     !
     nkl = kunit * ( nkbl / npool )
     !
     ! ... find out the reminder
     !
     nkr = ( nkstot - nkl * npool ) / kunit
     !
     ! ... Assign the reminder to the first nkr pools
     !
     IF ( my_pool_id < nkr ) nkl = nkl + kunit
     !
     ! ... find out the index of the first k point in this pool
     !
     iks = nkl * my_pool_id + 1
     !
     IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
     !
     ! ... find out the index of the last k point in this pool
     !
     ike = iks + nkl - 1
     !
  END IF
  !
  ! ... find out the global number of G vectors: ngm_g  
  !
  ngm_g = ngm
  !
  CALL mp_sum( ngm_g, intra_pool_comm )
  !
  ! ... build the igk_l2g array, yielding the correspondence between
  ! ... the local k+G index and the global G index - see also ig_l2g
  !
  ALLOCATE ( igk_l2g( npwx, nks ) )
  igk_l2g = 0
  !
  ALLOCATE( kisort( npwx ) )
  !
  DO ik = 1, nks
     !
     kisort = 0
     npw    = npwx
!!$     !
     CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
          ecutwfc/tpiba2, npw, kisort(1), g2kin )
     !
     CALL gk_l2gmap( ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
!!$     !
     ngk(ik) = npw
     !
  END DO
  !
  DEALLOCATE( kisort )
  !
  ! ... compute the global number of G+k vectors for each k point
  !
  ALLOCATE( ngk_g( nkstot ) )
  !
  ngk_g = 0
  ngk_g(iks:ike) = ngk(1:nks)
  !
  CALL mp_sum( ngk_g, intra_image_comm )
  !
  ! ... compute the Maximum G vector index among all G+k an processors
  !
  npw_g = MAXVAL( igk_l2g(:,:) )
  !
  CALL mp_max( npw_g, intra_image_comm )
  !
  ! ... compute the Maximum number of G vector among all k points
  !
  npwx_g = MAXVAL( ngk_g(1:nkstot) )
  !
  ! 
  ! ... define a further l2g map to read gkvectors and wfc coherently 
  ! 
  ALLOCATE( igk_l2g_kdip( npwx_g, nks ) )
  igk_l2g_kdip = 0
  !
  DO ik = iks, ike
     !
     CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
          igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
  END DO
  !
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin( iunpun, "EIGENVECTORS" )
     !
  END IF
  !
  num_k_points = nkstot
  !
  IF ( nspin == 2 ) num_k_points = nkstot / 2
  !
  k_points_loop: DO ik = 1, num_k_points
     !
     IF ( ionode ) THEN
        !
        CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
        !
        IF ( nspin == 2 .OR. noncolin ) THEN
           !
           CALL iotk_scan_begin( iunpun, "WFC.1", FOUND = twfcollect  )
           IF ( twfcollect ) CALL iotk_scan_end( iunpun, "WFC.1" )
           !
        ELSE
           !
           CALL iotk_scan_begin( iunpun, "WFC", FOUND = twfcollect  )
           IF ( twfcollect ) CALL iotk_scan_end( iunpun, "WFC" )
           !
        ENDIF
        !
     END IF
     !
     CALL mp_bcast( twfcollect, ionode_id, intra_image_comm )
     !
     IF ( .NOT. twfcollect ) THEN
        !
        IF ( ionode ) THEN
           !
           CALL iotk_scan_end( iunpun, &
                "K-POINT" // TRIM( iotk_index( ik ) ) )
           !
        END IF
        !
        EXIT k_points_loop
        !
     END IF
     !
     IF ( nspin == 2 ) THEN
        !
        ispin = 1 
        auxevc=(0.0_DP, 0.0_DP)
        !
        ! ... no need to read isk here: they are read from band structure
        ! ... and correctly distributed across pools in read_file
!!! isk(ik) = 1
        !
        IF ( ionode ) THEN
           !
           filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin, &
                DIR=lkpoint_dir ) )
           !
        END IF
        !
        CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,      &
             auxevc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
             ngk(ik-iks+1), filename, scalef, &
             ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
        !
        IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
           !
           CALL save_buffer2 ( auxevc, nwordwfc, iunwfc, (ik-iks+1) )
           !
        END IF
        !
        ispin = 2
        ik_eff = ik + num_k_points
        auxevc=(0.0_DP, 0.0_DP)
        !
        ! ... no need to read isk here (see above why)
        !isk(ik_eff) = 2
        !
        IF ( ionode ) THEN
           !
           filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin, &
                DIR=lkpoint_dir ) )
           !
        END IF
        !
        CALL read_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin,      &
             auxevc, npw_g, nbnd, igk_l2g_kdip(:,ik_eff-iks+1),   &
             ngk(ik_eff-iks+1), filename, scalef, &
             ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
        !
        IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
           !
           CALL save_buffer2 ( auxevc, nwordwfc, iunwfc, (ik_eff-iks+1) )
           !
        END IF
        !
     ELSE
!!$        !
!!$        ! ... no need to read isk here (see above why)
!!$        !isk(ik) = 1
!!$        !
        auxevc=(0.0_DP, 0.0_DP)
        IF ( noncolin ) THEN
           !
           DO ipol = 1, npol
              !
              IF ( ionode ) THEN
                 !
                 filename = TRIM( wfc_filename( dirname, 'evc', ik, ipol, &
                      DIR=lkpoint_dir ) )
                 !
              END IF
              !
!!! TEMP
              nkl=(ipol-1)*npwx+1
              nkr= ipol   *npwx
              CALL read_wfc( iunout, ik, nkstot, kunit, ispin,          &
                   npol, auxevc(nkl:nkr,:), npw_g, nbnd,         &
                   igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1),   &
                   filename, scalef, & 
                   ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
              !
           END DO
           !
        ELSE
!!$           !
           IF ( ionode ) THEN
              !
              filename = TRIM( wfc_filename( dirname, 'evc', ik, &
                   DIR=lkpoint_dir ) )
              !
           END IF
           !
           CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,         &
                auxevc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),      &
                ngk(ik-iks+1), filename, scalef, &
                ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
           !
        END IF
        !
        IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
!!$           !
           CALL save_buffer2 ( auxevc, nwordwfc, iunwfc, (ik-iks+1) )
!!$           !
!!$           ! the following two line can be used to debug read_wfc
!!$           ! WRITE(200+10*ik+me_pool,fmt="(2D18.10)") evc
!!$           ! CLOSE(200+10*ik+me_pool )
!!$           !
        END IF
        !
     END IF
!!$     !
     IF ( ionode ) THEN
        !
        CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
        !
     END IF
     !
  END DO k_points_loop
  !
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( igk_l2g_kdip )
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_end( iunpun, "EIGENVECTORS" )
     !
     CALL iotk_close_read( iunpun )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE read_wavefunctions2
!
!
!----------------------------------------------------------------------------
SUBROUTINE gk_l2gmap( ngm, ig_l2g, ngk, igk, igk_l2g )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine maps local G+k index to the global G vector index
  ! ... the mapping is used to collect wavefunctions subsets distributed
  ! ... across processors.
  ! ... Written by Carlo Cavazzoni
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER, INTENT(IN)  :: ngm, ngk, igk(ngk), ig_l2g(ngm)
  INTEGER, INTENT(OUT) :: igk_l2g(ngk)
  INTEGER              :: ig
  !
  ! ... input: mapping between local and global G vector index
  !
  DO ig = 1, ngk
     !
     igk_l2g(ig) = ig_l2g(igk(ig))
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE gk_l2gmap
!
!-----------------------------------------------------------------------
SUBROUTINE gk_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip, igwk )
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine maps local G+k index to the global G vector index
  ! ... the mapping is used to collect wavefunctions subsets distributed
  ! ... across processors.
  ! ... This map is used to obtained the G+k grids related to each kpt
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER,           INTENT(IN)  :: npw_g, ngk_g, ngk
  INTEGER,           INTENT(IN)  :: igk_l2g(ngk)
  INTEGER, OPTIONAL, INTENT(OUT) :: igwk(ngk_g), igk_l2g_kdip(ngk)
  !
  INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
  INTEGER              :: ig, ig_, ngg
  !
  !
  ALLOCATE( itmp( npw_g ) )
  ALLOCATE( igwk_( ngk_g ) )
  !
  itmp(:)  = 0
  igwk_(:) = 0
  !
  !
  DO ig = 1, ngk
     !
     itmp(igk_l2g(ig)) = igk_l2g(ig)
     !
  END DO
  !
  CALL mp_sum( itmp, intra_pool_comm )
  !
  ngg = 0
  DO ig = 1, npw_g
     !
     IF ( itmp(ig) == ig ) THEN
        !
        ngg = ngg + 1
        !
        igwk_(ngg) = ig
        !
     END IF
     !
  END DO
  !
  IF ( ngg /= ngk_g ) &
       CALL errore( 'igk_l2g_kdip', 'unexpected dimension in ngg', 1 )
  !
  IF ( PRESENT( igwk ) ) THEN
     !
     igwk(1:ngk_g) = igwk_(1:ngk_g)
     !
  END IF
  !
  IF ( PRESENT( igk_l2g_kdip ) ) THEN
     !
     ALLOCATE( igwk_lup( npw_g ) )
     !
     !$omp parallel private(ig_, ig)
     !$omp workshare
     igwk_lup = 0
     !$omp end workshare
     !$omp do
     do ig_ = 1, ngk_g
        igwk_lup(igwk_(ig_)) = ig_
     end do
     !$omp end do
     !$omp do
     do ig = 1, ngk
        igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
     end do
     !$omp end do
     !$omp end parallel
     !
     DEALLOCATE( igwk_lup )

  END IF
  !
  DEALLOCATE( itmp, igwk_ )
  !
  RETURN
  !
END SUBROUTINE gk_l2gmap_kdip
!
!------------------------------------------------------------------------
SUBROUTINE read_band_structure2( dirname, ierr )
  !------------------------------------------------------------------------
  !
  USE control_flags, ONLY : lkpoint_dir
!  USE basis,    ONLY : natomwfc
  USE lsda_mod, ONLY : lsda, isk
  USE klist,    ONLY : nkstot, wk, nks !, nelec
  USE wvfct,    ONLY : et, wg , nbnd
  USE constants, ONLY :e2
!  USE ener,     ONLY : ef, ef_up, ef_dw
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  INTEGER,          INTENT(OUT) :: ierr
  !
  INTEGER :: ik, ik_eff, num_k_points, iunout
  LOGICAL :: found, two_fermi_energies_
  LOGICAL , SAVE :: lbs_read = .false.
  !
  if(.not. allocated(isk_aux)) allocate(isk_aux(SIZE(isk,1)))
  if(.not. allocated(et_aux)) allocate(et_aux(nbnd,nkstot))
  if(.not. allocated(wg_aux)) allocate(wg_aux(nbnd,nkstot))
  ierr = 0
  IF ( lbs_read ) RETURN
  !
!!$  IF ( .NOT. lspin_read ) &
!!$       CALL errore( 'read_band_structure', 'read spin first', 1 )
!!$  IF ( .NOT. lbz_read ) &
!!$       CALL errore( 'read_band_structure', 'read band_structure first', 1 )
  !
  IF ( ionode ) &
       CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
       & TRIM( xmlpun ), IERR = ierr )
  !
  CALL mp_bcast( ierr, ionode_id, intra_image_comm )
  !
  IF ( ierr > 0 ) RETURN
  !
  IF (.NOT.lkpoint_dir) THEN
     !
     IF ( ionode ) &
          CALL iotk_open_read( iunout, FILE = TRIM( dirname ) // '/' // &
          & TRIM( xmlpun )//'.eig', IERR = ierr )
     !
     CALL mp_bcast( ierr, ionode_id, intra_image_comm )
     !
     IF ( ierr > 0 ) RETURN
     !
  END IF
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
     !
     CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_aux )
     !
     CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc_aux, &
          FOUND = found )
     IF ( .NOT. found ) natomwfc_aux = 0
     !
     CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_aux )
     !
     CALL iotk_scan_dat( iunpun, "FERMI_ENERGY", ef_aux, FOUND = found )
     !
     IF ( found ) THEN
        ef_aux = ef_aux * e2
     ELSE
        ef_aux = 0.d0
     END IF
     !
     CALL iotk_scan_dat( iunpun, "TWO_FERMI_ENERGIES", two_fermi_energies_ )
     !
     IF ( two_fermi_energies_ ) THEN
        !
        CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_UP", ef_up_aux )
        CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw_aux )
        !
        ef_up_aux = ef_up_aux * e2
        ef_dw_aux = ef_dw_aux * e2
        !
     ENDIF
     !
     CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
     !
  END IF
  !
  num_k_points = nkstot
  !
  IF ( lsda ) num_k_points = nkstot / 2
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin( iunpun, "EIGENVALUES" )
     !
     k_points_loop: DO ik = 1, num_k_points
        !
        CALL iotk_scan_begin( iunpun, &
             "K-POINT" // TRIM( iotk_index( ik ) ) )
        !
        IF ( lsda ) THEN
           !
           isk_aux(ik) = 1
           !
           IF (lkpoint_dir) THEN
              CALL iotk_scan_begin(iunpun, "DATAFILE"//TRIM(iotk_index(1)))
              CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et_aux(:,ik)  )
              CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg_aux(:,ik) )
              CALL iotk_scan_end(iunpun, "DATAFILE"//TRIM(iotk_index(1)) )
           ELSE
              CALL iotk_scan_begin( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP")
              CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et_aux(:,ik)  )
              CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg_aux(:,ik) )
              CALL iotk_scan_end( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP")
           ENDIF
           !
           ik_eff = ik + num_k_points
           !
           isk_aux(ik_eff) = 2
           !
           IF (lkpoint_dir) THEN
              CALL iotk_scan_begin(iunpun,"DATAFILE"//TRIM(iotk_index(2)) )
              CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et_aux(:,ik_eff) )
              CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg_aux(:,ik_eff) )
              CALL iotk_scan_end( iunpun, "DATAFILE"//TRIM(iotk_index(2)) )
           ELSE
              CALL iotk_scan_begin( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW")
              CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et_aux(:,ik_eff) )
              CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg_aux(:,ik_eff) )
              CALL iotk_scan_end( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW")
           ENDIF
           !
           !
        ELSE
           !
           isk_aux(ik) = 1
           !
           IF (lkpoint_dir) THEN
              CALL iotk_scan_begin( iunpun, "DATAFILE" )
              CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et_aux(:,ik) )
              CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg_aux(:,ik) )
              CALL iotk_scan_end  ( iunpun, "DATAFILE" )
           ELSE
              CALL iotk_scan_begin( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) ))
              CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et_aux(:,ik) )
              CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg_aux(:,ik) )
              CALL iotk_scan_end( iunout, &
                   "DATA_EIG"//TRIM( iotk_index( ik ) ))
           ENDIF
           !
           !
        END IF
        !
        CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
        !
     END DO k_points_loop
     !
     et_aux(:,:) = et_aux(:,:) * e2
     !
     FORALL( ik = 1:nkstot ) wg_aux(:,ik) = wg_aux(:,ik)*wk(ik)
     !
     CALL iotk_scan_end( iunpun, "EIGENVALUES" )
     !
     CALL iotk_close_read( iunpun )
     !
     IF (.NOT.lkpoint_dir) CALL iotk_close_read( iunout )
     !
  END IF
  !
  CALL mp_bcast( nelec_aux,    ionode_id, intra_image_comm )
  CALL mp_bcast( natomwfc_aux, ionode_id, intra_image_comm )
  CALL mp_bcast( nbnd_aux,     ionode_id, intra_image_comm )
  CALL mp_bcast( isk_aux,      ionode_id, intra_image_comm )
  CALL mp_bcast( et_aux,       ionode_id, intra_image_comm )
  CALL mp_bcast( wg_aux,       ionode_id, intra_image_comm )
  CALL mp_bcast( ef_aux,       ionode_id, intra_image_comm )

  CALL poolscatter( nbnd, nkstot, et_aux, nks, et_aux )
  CALL poolscatter( nbnd, nkstot, wg_aux, nks, wg_aux )

  !
  lbs_read = .TRUE.
  !
  RETURN
  !
END SUBROUTINE read_band_structure2
!
!!$!--------------------------------------------------------------------
subroutine set_ext_weights (nks, nkstot, wk, nbnd, nelec, degauss, ngauss, &
     et, ef, demet, netot, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  USE lsda_mod,       ONLY : lsda, nspin, current_spin
  implicit none
  !
  integer, intent(in) :: nks, nkstot, nbnd, ngauss
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss
  integer, intent(in) :: is, isk(nks)
  real(DP), intent(out) :: wg (nbnd, nks), ef, demet, netot
!  real(DP) :: wg (nbnd, nks), ef, demet, netot
  !
  integer :: kpoint, ibnd, unit, ikd, ibd
  real(DP) , external :: wgauss, w1gauss, efermig
  real(DP):: nvbel,ncbel,efvb,efcb, efgs,dwg(nbnd,nks),scal,temp_tot,kwgt,&
       &gaus,dne,maxocc,dnefrac
  integer :: homo(nks), lumo(nks)
  logical:: excit
  logical, save:: first=.true.
  character(64)::filename

  if(.not. allocated(et_aux)) allocate(et_aux(nbnd,nkstot))
  if(.not. allocated(wg_aux)) allocate(wg_aux(nbnd,nkstot))

  if(first) then
     if(ionode) then
        unit = find_free_unit() 
        filename='weights.dat'
        OPEN( UNIT = unit,    FILE = filename,    STATUS ='OLD' , FORM = 'FORMATTED')
        read(unit,*)
        read(unit,*)
        do kpoint = 1, nkstot
           do ibnd=1, nbnd
              read(unit,*) ikd, ibd, et_aux(ibnd,kpoint), wg_aux(ibnd,kpoint)
           enddo
        enddo
        CLOSE(unit)
     endif
     CALL mp_bcast( et_aux,       ionode_id, intra_image_comm )
     CALL mp_bcast( wg_aux,       ionode_id, intra_image_comm )
     !
     CALL poolscatter( nbnd, nkstot, et_aux, nks, et_aux )
     CALL poolscatter( nbnd, nkstot, wg_aux, nks, wg_aux )
     !
     first=.false.
  endif
  wg(1:nbnd,1:nks)=wg_aux(1:nbnd,1:nks)
  ef = efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)
  demet = 0.d0
  netot = 0.0d0
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        netot = netot + wg (ibnd, kpoint)
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
             degauss * w1gauss ( (ef-et(ibnd,kpoint)) / degauss, ngauss)
     enddo
     
  enddo
  !
  CALL mp_sum( demet, inter_pool_comm )
  !
  !
  CALL mp_sum( netot, inter_pool_comm )
  !
  dne=(nelec-netot)
  write(*,'(A,2F18.12)') "netot=",netot,dne
  if(dne/nelec > 1.0d-7) then
     if(me_pool .eq. root_pool) then
        maxocc=(2.0d0/nspin)
        do kpoint = 1, nks
           if (is /= 0) then
              if (isk(kpoint).ne.is) cycle
           end if
           dnefrac=dne*wk(kpoint)*(nspin/2.0d0)
           bnd_loop: do ibnd=1, nbnd
              if((wg(ibnd, kpoint)+dnefrac) .lt. maxocc*wk(kpoint)) then
                 wg(ibnd, kpoint)=wg(ibnd, kpoint)+dnefrac
                 exit bnd_loop
              endif
           enddo bnd_loop
        enddo
     endif
     CALL mp_bcast( wg,  root_pool, intra_pool_comm )
     netot=0.0d0
     do kpoint = 1, nks
        if (is /= 0) then
           if (isk(kpoint).ne.is) cycle
        end if
        do ibnd = 1, nbnd
           netot = netot + wg (ibnd, kpoint)
        enddo
     enddo
     CALL mp_sum( netot, inter_pool_comm )
     !
     write(*,'(A,F18.12)') "netot=",netot     
  endif
  
  return
end subroutine set_ext_weights


end module proj_weights

