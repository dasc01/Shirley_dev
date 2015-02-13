!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
! ----------------------------------------------------------------------
  module corerepair_mod
! ----------------------------------------------------------------------
  use kinds, only : dp

  implicit none
  private
  public :: core_type, corerepair_type, corerep
  public :: read_corerepair, bcast_corerepair

  type core_type
    integer :: species
    integer :: atom
    character(255) :: filename
    character(255) :: label
    integer :: nwfc1, nwfc2
    integer,pointer :: lwfc1(:), lwfc2(:)
    integer :: nproj1, nproj2
    integer :: nnonzero
    complex(dp),pointer :: matrix(:,:,:)
  end type core_type

  type corerepair_type
    integer :: nspecies, natom
    integer :: ncore
    type(core_type),pointer :: core(:)    
  end type corerepair_type

  type(corerepair_type) :: corerep

  contains

! ----------------------------------------------------------------------
  subroutine read_corerepair( iunit, corerep )
! ----------------------------------------------------------------------
  ! for XAS
  USE io_files,     ONLY : find_free_unit
  !
  !
  integer,intent(in) :: iunit
  type(corerepair_type) :: corerep

  integer :: i, iuncore, ierr
!  integer,external :: freeunit

  read(iunit,*) ! heading
  read(iunit,*) corerep%nspecies, corerep%natom
  read(iunit,*) corerep%ncore
  allocate( corerep%core(corerep%ncore) )
  do i=1,corerep%ncore
    read(iunit,*) corerep%core(i)%species, corerep%core(i)%atom, &
                  corerep%core(i)%filename

    iuncore=find_free_unit()
    open(iuncore,file=trim(corerep%core(i)%filename),form='formatted', &
         iostat=ierr)
    if( ierr/=0 ) call errore( 'read_corerepair','unable to open file '// &
                               trim(corerep%core(i)%filename), 1 )

    call read_core( iuncore, corerep%core(i) )
    close(iuncore)
  enddo

  end subroutine read_corerepair

! ----------------------------------------------------------------------
  subroutine read_core( iunit, core )
! ----------------------------------------------------------------------
  integer,intent(in) :: iunit
  type(core_type) :: core

  integer :: i
  integer :: ip1, ip2, ixyz
  real(dp) :: cR, cI

  read(iunit,*) core%label
  read(iunit,*) core%nwfc1, core%nwfc2

  allocate( core%lwfc1(core%nwfc1) )
  allocate( core%lwfc2(core%nwfc2) )
  read(iunit,*) core%lwfc1(1:core%nwfc1)
  read(iunit,*) core%lwfc2(1:core%nwfc2)

  core%nproj1 = sum(2*core%lwfc1(1:core%nwfc1)+1)
  core%nproj2 = sum(2*core%lwfc2(1:core%nwfc2)+1)

  allocate( core%matrix(core%nproj1,core%nproj2,3) )
  core%matrix = 0.d0
  read(iunit,*) core%nnonzero
  do i=1,core%nnonzero
    read(iunit,*) ip1, ip2, ixyz, cR, cI
    core%matrix(ip1,ip2,ixyz) = cmplx(cR,cI,kind=dp)
  enddo
  
  end subroutine read_core

! ----------------------------------------------------------------------
  subroutine bcast_corerepair( corerep, root )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast
  use mp_global, only : mpime
  type(corerepair_type) :: corerep
  integer,intent(in) :: root

  integer :: i

  call mp_bcast( corerep%nspecies, root )
  call mp_bcast( corerep%natom,    root )
  call mp_bcast( corerep%ncore,    root )

  if( mpime/=root ) allocate( corerep%core(corerep%ncore) )
  do i=1,corerep%ncore
    call bcast_core( corerep%core(i), root )
  enddo

  end subroutine bcast_corerepair

! ----------------------------------------------------------------------
  subroutine bcast_core( core, root )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast
  use mp_global, only : mpime
  type(core_type) :: core
  integer,intent(in) :: root

  integer :: i

  call mp_bcast( core%species, root )
  call mp_bcast( core%atom, root )
  call mp_bcast( core%filename, root )

  call mp_bcast( core%label, root )

  call mp_bcast( core%nwfc1, root )
  call mp_bcast( core%nwfc2, root )
  if( mpime/=root ) allocate(core%lwfc1(core%nwfc1))
  if( mpime/=root ) allocate(core%lwfc2(core%nwfc2))
  call mp_bcast( core%lwfc1, root )
  call mp_bcast( core%lwfc2, root )

  call mp_bcast( core%nproj1, root )
  call mp_bcast( core%nproj2, root )
  call mp_bcast( core%nnonzero, root )
  if( mpime/=root ) allocate(core%matrix(core%nproj1,core%nproj2,3) )
  call mp_bcast( core%matrix, root )

  end subroutine bcast_core


  end module corerepair_mod
!---------------------------------------------------------------------
!
  module broaden

    implicit none

    public::add_gauss

  contains
    
    subroutine add_gauss( nener, ener, e, sigma, weight, spec )
      use kinds,only: dp
      implicit none
      
      integer :: nener
      real(dp) :: ener(nener), e, sigma, weight(:), spec(:,:)
      
      real(dp) :: its2, pre(size(weight))
      real(dp) :: arg(nener)
      integer :: i, j
      
      its2=2.d0*sigma*sigma
      pre=1.d0/sqrt(acos(-1.d0)*its2) * weight
      its2 = 1.d0/its2
      arg = ener - e
      arg = arg ** 2.d0
      arg = - arg * its2
      forall( i=1:nener, j=1:size(weight) ) &
           spec(i,j) = spec(i,j) + pre(j) * exp( arg(i) )
      
    end subroutine add_gauss
    
  end module broaden
!-------------------------------------------------------------------
MODULE projections
  USE kinds, ONLY : DP

  TYPE wfc_label
     INTEGER na, n, l, m
  END TYPE wfc_label
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:,:)

END MODULE projections
!
!
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
  COMPLEX(DP), ALLOCATABLE, TARGET :: auxevc(:,:),aux2(:,:),projDM(:,:,:)
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
  COMPLEX(DP),allocatable::wgpsps(:,:)

  if(.not. allocated(auxevc)) allocate(auxevc(SIZE(evc,1),SIZE(evc,2)),aux2(SIZE(evc,1),SIZE(evc,2)))
  if(.not. allocated(auxwg)) allocate(auxwg(nbnd,nks))
  if(.not. allocated(projDM)) allocate(projDM(nbnd,nbnd,nks))
  if(first) then
     CALL open_buffer2( iunwfc, 'wfc2', nwordwfc, nks, exst )
     CALL read_wavefunctions2(trim(auxprefix), ierr)
     CALL  read_band_structure2( trim(auxprefix), ierr )
     first=.false.
  endif

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
  auxwg(:,:)=0.0d0
  CALL allocate_bec_type(nbnd, nbnd, psipsi)
  allocate(wgpsps(nbnd, nbnd))

  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     npw = ngk(ik)
     !
     !
     if(nks > 1) READ( iunigk ) igk
     CALL get_buffer2 ( auxevc, nwordwfc, iunwfc, ik )  !these are the physical states of interest
     !
     CALL get_buffer ( aux2, nwordwfc, iunwfc, ik )    !these are the basis states
          
     CALL calbec(npw, auxevc, aux2, psipsi)   !get projection of physical states (auxevc) onto basis states

!!$     do ibnd=1, nbnd  !loop over target states
!!$
!!$        do jbnd=1, nbnd !loop over source states
!!$           auxwg(ibnd,ik)=auxwg(ibnd,ik)+wg_aux(jbnd,ik)*psipsi%k(jbnd,ibnd)*CONJG(psipsi%k(jbnd,ibnd))
!!$        enddo
!!$        write(*,'(A,2I6,F16.12)') "ik, ibnd, psipsi=", ik, ibnd, REAL(psipsi%k(ibnd,ibnd))
!!$     enddo

     !multiply the rows by apporpriate band occupancy weights 
     do ibnd=1, nbnd
        wgpsps(ibnd,:)=wg_aux(ibnd,ik)*psipsi%k(ibnd,:)
     enddo
     
     CALL ZGEMM( 'C', 'N', nbnd, nbnd, &
          nbnd, (1.0d0,0.0d0), psipsi%k(1,1), nbnd, &
          wgpsps(1,1), nbnd, &
          (0.0d0,0.0d0), projDM(1,1,ik), nbnd )

     
     do ibnd=1, nbnd  !loop over target states
        totel=totel+projDM(ibnd,ibnd,ik)
     enddo

!check completeness of projection
!!$     do ibnd=1, nbnd
!!$       totel=totel+auxwg(ibnd,ik)
!!$       write(*,'(A,2I6,F16.12)') "ik, ibnd, auxwg =", ik, ibnd, auxwg(ibnd, ik)
!!$     enddo
     
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
  CALL deallocate_bec_type ( psipsi )
  deallocate(wgpsps)
  CALL mp_sum( totel, inter_pool_comm )
  write(*,'(A,F12.8)') "totel = ", totel
  !
  call setup_xas()
  !
end subroutine set_proj_weights
!-----------------------------------------------------------------------
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine setup_xas()
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : degauss, ngauss, lgauss
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir, trimcheck
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup, nproc_ortho
  USE environment,      ONLY : environment_start
  !
  ! for GWW
  USE io_files,     ONLY : find_free_unit
  !
  !
  ! for XAS
  use corerepair_mod
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos, filproj, io_choice, outdir, chapprox, spectype, filename
  REAL (DP)      :: Emin, Emax, DeltaE, degauss1, smoothing, broadening
  INTEGER :: ngauss1, ios, Nener
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes, fixocc
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)

  !
  ! for GWW
  INTEGER :: iun, idum, unit
  REAL(DP) :: rdum1,rdum2,rdum3
  LOGICAL :: lex, lgww
  !
  !
  !for XAS
  LOGICAL :: readcorerep
  !
  NAMELIST / inputpp / outdir, prefix, ngauss1, degauss1, lsym, &
             Emin, Emax, Nener, DeltaE, io_choice, smoothing, filpdos, filproj, &
             lgww, & !if .true. use GW QP energies from file bands.dat
             kresolveddos,  irmin, irmax,  &
             readcorerep, chapprox, broadening, spectype, fixocc
  !
  ! initialise environment
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.0d0  !Right Shift of spectrum in eV
  Nener  = 1000
  ngauss1 = 0
  lsym   = .true.
  degauss1= 0.d0
  lgww   = .false.
  kresolveddos = .false.
  irmin(:,:)  = 1
  irmax(:,:)  = 0
  readcorerep = .false.
  chapprox = 'XCH'
  broadening = 0.0
  spectype='RAW'
  fixocc=.false.
  !
  ios = 0
  !
  IF ( ionode )  THEN
     unit = find_free_unit() 
     filename='pw2xas.in'
     OPEN( UNIT = unit,    FILE = filename,    STATUS ='OLD' , FORM = 'FORMATTED')
     !    
     READ (unit, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id )
  !
  IF (ios /= 0) CALL errore ('setup_xas', 'reading inputpp namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix,  ionode_id )
  CALL mp_bcast( filproj,  ionode_id )
  CALL mp_bcast( ngauss1, ionode_id )
  CALL mp_bcast( degauss1,ionode_id )
  CALL mp_bcast( DeltaE,  ionode_id )
  CALL mp_bcast( Nener,  ionode_id )
  CALL mp_bcast( lsym,  ionode_id )
  CALL mp_bcast( Emin, ionode_id )
  CALL mp_bcast( Emax, ionode_id )
  ! for GWW
  CALL mp_bcast( lgww, ionode_id )
  ! for projection on boxes
  CALL mp_bcast( irmin, ionode_id )
  CALL mp_bcast( irmax, ionode_id )
  CALL mp_bcast( readcorerep, ionode_id )
  CALL mp_bcast( chapprox, ionode_id )
  CALL mp_bcast( broadening, ionode_id )
  CALL mp_bcast( spectype, ionode_id )
  CALL mp_bcast( fixocc, ionode_id )
!!$  !
!!$  !   Now allocate space for pwscf variables, read and check them.
!!$  !
!!$  CALL read_file ( )
!!$  !
!!$  CALL openfil_pp ( )
!!$  !
  if( readcorerep ) then
     ! read atomic matrix elements for core-valence position
     if(ionode) then
        call read_corerepair( unit, corerep )
        write(stdout,*) 'done reading corerepair'
     endif
     call bcast_corerepair( corerep, ionode_id )
  endif

  !
  !   decide Gaussian broadening
  !
!!$  lgauss=.true.
  IF (degauss1/=0.d0) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss1,degauss1
!!$     lgauss=.true.
  ELSEIF (lgauss) THEN
     degauss1=degauss
     ngauss1 = ngauss
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss1,degauss1
  ELSE
     degauss1=0.1d0/rytoev
     ngauss1 =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss1,degauss1
!!$     lgauss=.true.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !
  CALL projwave (filproj, lsym, lgww)

  CALL partialdos (Emin, Emax, DeltaE, Nener, broadening, chapprox, spectype, fixocc, kresolveddos, filpdos)

END subroutine setup_xas
!

!-----------------------------------------------------------------------
SUBROUTINE projwave( filproj, lsym, lgww )
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE grid_dimensions, ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symm_base, ONLY: nsym, irt, d1, d2, d3
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf, nh
  USE becmod,   ONLY: bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, find_free_unit
  USE spin_orb, ONLY: lspinorb
  USE wavefunctions_module, ONLY: evc
  USE buffers,           ONLY : get_buffer
  !
  USE projections
  !
  !
  ! for XAS
  use corerepair_mod
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*) :: filproj
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, nwfc,&
       nwfc1, lmax_wfc, is, ios, iunproj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE ::roverlap(:,:), rwork1(:),rproj0(:,:)
  ! ... or for gamma-point.
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2)
  INTEGER  :: nksinit, nkslast
  CHARACTER(len=256) :: filename
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  INTEGER, ALLOCATABLE :: idx(:)
  LOGICAL :: lsym
  !
  !
  ! for GWW
  INTEGER :: iun, idum
  REAL(DP) :: rdum1,rdum2,rdum3
  LOGICAL :: lex, lgww
  !
  ! for XAS
  INTEGER :: b_ptr, nbeta, nbnt , ixyz, ncp
  COMPLEX(DP),ALLOCATABLE::psibeta(:,:),temp(:,:),tcr(:,:)
  TYPE (bec_type)    ::becp
  
  !
  WRITE( stdout, '(/5x,"Calling projwave .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  !
  if(.not. allocated(auxevc)) allocate(auxevc(SIZE(evc,1),SIZE(evc,2)))
  ! fill structure nlmchi
  !
  if(.not. allocated(nlmchi)) ALLOCATE (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           lmax_wfc = max (lmax_wfc, l )
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
  IF (lmax_wfc > 3) CALL errore ('projwave', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projwave', 'wrong # of atomic wfcs?', 1)
  !
  ncp = corerep%core(1)%nproj2
  if(.not. allocated(proj_aux))  ALLOCATE( proj_aux (ncp, nbnd, 3, nks) )
!!$  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj_aux      = (0.d0, 0.0d0)
  CALL allocate_bec_type (nkb, nbnd, becp )
!!$  ALLOCATE(e (natomwfc) )
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  DO ik = 1, nks
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

     CALL get_buffer ( auxevc, nwordwfc, iunwfc, ik )    !these are the basis states
!!$     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)

!!$     CALL atomic_wfc (ik, wfcatom)

     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL calbec ( npw, vkb, auxevc, becp)

!Pick out the <psi|beta> on the core-excited atom
     nbeta=corerep%core(1)%nproj1
     if(.not. ALLOCATED(psibeta)) ALLOCATE(psibeta(nbnd,nbeta))
     psibeta(:,:)=CMPLX(0.0d0)
!loop in the right order over becp
     b_ptr=0
     DO nt=1, ntyp !loop over types       
        DO na = 1, nat !loop over atoms
           IF(ityp(na) .eq. nt ) THEN  !atom is of current type
              IF(corerep%core(1)%atom .eq. na) THEN !is the core-excited atom
                 IF(nbeta .ne. nh(nt)) then
                    write(*,*) "atom=",na,"type=",ityp(na),"spec=",corerep%core(1)%species
                    write(*,*) 'error! wrong number of nbeta, nh(nt)',nbeta, nh(nt)
                    STOP
                 endif
                 IF(gamma_only) THEN
                    DO nb=1, nbeta
                       DO ibnd=1, nbnd
                          psibeta(ibnd,nb) = CMPLX(becp%r(b_ptr+nb,ibnd))
                       ENDDO
                    ENDDO
                 ELSE
                    DO nb=1, nbeta
                       DO ibnd=1, nbnd
                          psibeta(ibnd,nb) = becp%k(b_ptr+nb,ibnd)
                       ENDDO
                    ENDDO                    
                 ENDIF
              ENDIF
              b_ptr=b_ptr+nh(nt)
           ENDIF
        ENDDO
     ENDDO
    ! sum <phi|r|psi><beta|nk>
     DO ixyz=1, 3
        CALL ZGEMM( 'T', 'C', ncp, nbnd, &
                    nbeta, (1.0d0,0.0d0), corerep%core(1)%matrix(1,1,ixyz), nbeta, &
                    psibeta, nbnd, &
                    (0.0d0,0.0d0), proj_aux(1,1,ixyz,ik), ncp )
     ENDDO
     ! on k-points
  ENDDO
  call deallocate_bec_type(becp)
!!$  CALL poolrecover (proj_aux, ncp*nbnd*3, nkstot, nks)
  !
  !
  RETURN
  !
END SUBROUTINE projwave
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE  partialdos (Emin, Emax, DeltaE, Nener, broadening, chapprox, spectype, fixocc, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nks, nkstot, degauss, ngauss, lgauss, nelec
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd, wg
  USE constants, ONLY: rytoev
  !
  USE projections
  ! for XAS
  USE kinds, only : DP
  USE broaden, only : add_gauss
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos, chapprox, spectype
  REAL(DP) :: Emin, Emax, DeltaE, broadening
  LOGICAL :: kresolveddos, fixocc
  INTEGER :: Nener
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  INTEGER :: ik, ibnd, jbnd, m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff
  REAL(DP) :: etev, delta, Elw, Eup, wkeff
  REAL(DP), ALLOCATABLE :: dostot(:,:,:), pdos(:,:,:,:), pdostot(:,:,:), &
       ldos(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !for XAS
  REAL(DP), ALLOCATABLE::occ(:,:),xas_xyz_sp(:,:),spec(:,:,:),tmp_spec(:,:)
  REAL(DP)::Ener(Nener), dE, ef, ehomo, elumo
  INTEGER::ncp,j,iunout,ispin
  real(DP) , external :: efermig, wgauss


!!$  if(trim(chapprox) .eq. 'XCH') then
!!$     ef = efermig (et, nbnd, nkstot, nelec-1, wk, degauss, ngauss, 0, isk)
!!$  else
!!$     ef = efermig (et, nbnd, nkstot, nelec, wk, degauss, ngauss, 0, isk)
!!$  endif
!!$
!!$  ALLOCATE(occ(nbnd,nkstot))
!!$  occ=0.0d0
!!$  ehomo=et(1,1)
!!$  elumo=et(nbnd,1)
!!$  DO ik=1, nkstot
!!$     DO ibnd=1,nbnd
!!$  !      write(*,*) "occ: ik, ibnd=", ik, ibnd
!!$        if(et(ibnd,ik) .gt. ehomo .and. et(ibnd,ik) .le. ef ) ehomo=et(ibnd,ik)
!!$        if(et(ibnd,ik) .lt. elumo .and. et(ibnd,ik) .ge. ef ) elumo=et(ibnd,ik)
!!$        if(trim(spectype) .eq. 'RAW') then
!!$           occ(ibnd,ik)=1.0d0
!!$        else if (trim(spectype) .eq. 'XAS') then
!!$           if(fixocc) then
!!$              occ(ibnd,ik)= 1 - (wg(ibnd,ik)/wk(ik))
!!$           else
!!$              occ(ibnd,ik)=1-wgauss((ef-et(ibnd,ik)) / degauss, ngauss)
!!$           endif
!!$        else if (trim(spectype) .eq. 'XES') then
!!$           if(fixocc) then
!!$              occ(ibnd,ik)= wg(ibnd,ik)/wk(ik)
!!$           else
!!$              occ(ibnd,ik)=wgauss((ef-et(ibnd,ik)) / degauss, ngauss)
!!$           endif
!!$        endif
!!$     ENDDO
!!$  ENDDO
!!$  write(*,*) "done occ..."
  DeltaE =0.0d0
  DeltaE = DeltaE/rytoev   !convert to Ryd
  broadening = broadening/rytoev !convert to Ryd
!!$  if(trim(spectype) .eq. 'XAS' .or. trim(spectype) .eq. 'XES') then
!!$     DeltaE=DeltaE-elumo
!!$  endif
  !
  !
  ! find band extrema
  !
  if(ionode) then
     Elw = et (1, 1)
     Eup = et (nbnd, 1)
     DO ik = 2, nkstot
        Elw = min (Elw, et (1, ik) )
        Eup = max (Eup, et (nbnd, ik) )
     ENDDO
     IF (broadening/=0.d0) THEN
        Eup = Eup + 5d0 * broadening
        Elw = Elw - 5d0 * broadening
     ENDIF
     Emin=Elw
     Emax=Eup
  endif
  CALL mp_bcast(Emin, ionode_id)
  CALL mp_bcast(Emax, ionode_id)
!!$  Emin=Emin/rytoev
!!$  Emax=Emax/rytoev
!!$  if(Emin .le. -1000.d0) Emin = max (Emin, Elw)
!!$  if(Emax .ge. +1000.d0) Emax = min (Emax, Eup)
!!$  if(trim(spectype) .eq. 'RAW') then
!!$     Emin = Elw
!!$     Emax = Eup
!!$  endif
  dE = (Emax - Emin)/dble(Nener)
  do ie=1,nener
    ener(ie) = Emin + dble(ie-1)*dE
  enddo
  write(*,*) "done erange..."

!
  ALLOCATE(xas_xyz_sp(3,nspin),spec(Nener,4,nspin),tmp_spec(Nener,4))
  spec(1:nener,1:4,1:nspin)=0.0d0
  tmp_spec(1:nener,1:4)=0.0d0
  xas_xyz_sp = 0.0d0
  ncp=size(proj_aux,1)
  current_spin = 1
!!$  ie_delta = 5 * degauss / DeltaE + 1

  DO ik = 1,nks
     ! use true weights
     wkeff=wk(ik)
     ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
     ikeff=1
     !
     IF ( nspin == 2 ) current_spin = isk ( ik )
     DO ibnd = 1, nbnd
        etev = et(ibnd,ik)
        forall(j=1:3) xas_xyz_sp(j,current_spin) = sum(conjg(proj_aux(1:ncp,ibnd,j,ik))*proj_aux(1:ncp,ibnd,j,ik))*wk(ik)
        do jbnd=1, nbnd
           forall(j=1:3) xas_xyz_sp(j,current_spin) = xas_xyz_sp(j,current_spin)- sum(conjg(proj_aux(1:ncp,ibnd,j,ik))*proj_aux(1:ncp,jbnd,j,ik))*projDM(ibnd,jbnd,ik)
        enddo
        tmp_spec(:,:)=0.0d0
        call add_gauss( Nener, ener, etev+DeltaE, &
             broadening, xas_xyz_sp(:,current_spin), tmp_spec(1:nener,2:4) )
        spec(1:nener,2:4,current_spin)=spec(1:nener,2:4,current_spin)+tmp_spec(1:nener,2:4)
     ENDDO
  ENDDO
  if(allocated(xas_xyz_sp)) deallocate(xas_xyz_sp)
  if(allocated(tmp_spec)) deallocate(tmp_spec)
  write(*,*) "done spec..."

  CALL mp_sum( spec, inter_pool_comm )
  

  iunout=4
!!$  fileout = trim(filpdos)//".pw.raw"
!!$  if(trim(spectype) .eq. 'XAS') then 
!!$     fileout = trim(filpdos)//".pw.xas"
!!$  else if(trim(spectype) .eq. 'XES') then
!!$     fileout = trim(filpdos)//".pw.xes"
!!$  endif
  fileout = trim(filpdos)//".XAS"
  if(ionode) then
     OPEN (iunout,file=fileout,form='formatted', status='unknown')
  endif

  do ispin=1, nspin
     forall( ie=1:nener) spec(ie,1,ispin) = sum(spec(ie,2:4,ispin))/3.d0
  enddo
  if(ionode) then
     write(iunout,'(a,f12.5,a)') '# Applied a delta shift of ', DeltaE*rytoev, ' eV'
     do ispin=1, nspin
        write(iunout,*) '#SPIN=',ispin
        do ie=1,nener
           write(iunout,'(5e14.5e3)') ener(ie)*rytoev,spec(ie,1:4,ispin)
        enddo
        write(iunout,*)
     enddo
     close(iunout)
  endif
  DEALLOCATE (proj_aux, spec)
  !
  RETURN
END SUBROUTINE partialdos
!
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
  INTEGER :: ik, ik_eff, num_k_points, iunout, unit, ibnd
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
  IF ( ionode )  THEN
     unit = find_free_unit() 
     OPEN( UNIT = unit,    FILE = 'OccVsE.dat',    STATUS ='UNKNOWN' , FORM = 'FORMATTED')
     !    
     DO ik=1,nkstot
        DO ibnd=1,nbnd
           write(unit,'(2F16.8)') et_aux(ibnd,ik)*13.60569, wg_aux(ibnd,ik)/wk(ik)
        ENDDO
     ENDDO
     CLOSE(unit)
  ENDIF

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

