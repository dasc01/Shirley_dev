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
!
!-----------------------------------------------------------------------
PROGRAM pw2xas 
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! calculates Lowdin charges, spilling parameter, projected DOS
  ! or computes the LDOS in a volume given in input as function of energy
  !
  !
  ! Input (namelist &inputpp ... / ):                         Default value
  !
  !    prefix        prefix of input file produced by pw.x    'pwscf'
  !                    (wavefunctions are needed)
  !    outdir        directory containing the input file       ./
  !    ngauss        type of gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            =-99  Fermi-Dirac function
  !    degauss       gaussian broadening, Ry (not eV!)          0.0
  !    Emin, Emax    min, max energy (eV) for DOS plot          band extrema
  !    DeltaE        energy grid step (eV)                      none
  !    lsym          if true the projections are symmetrized    .true.
  !    filproj       file containing the projections            none
  !    filpdos       prefix for output files containing PDOS(E) prefix
  !    lgww          if .true. take energies from previous GWW calculation
  !                  (file bands.dat)
  !    kresolveddos  if .true. the DOS is written as function   .false.
  !                  of the k-point (not summed over all of them)
  !                  all k-points but results
  !    tdosinboxes   if .true., the local DOS in specified      .false.
  !                  volumes (boxes) is computed
  !    n_proj_boxes  number of volumes for the local DOS        0
  !    irmin, irmax  first and last point of the FFT grid       1, 0
  !                  included in the volume
  !    plotboxes     if .true., the volumes are given in output .false.
  !
  !
  ! Output:
  !
  !   Projections are written to standard output,
  !   and also to file filproj if given as input.
  !   The total DOS and the sum of projected DOS are written to file
  !   "filpdos".pdos_tot.
  !   The format for the collinear, spin-unpolarized case and the
  !   non-collinear, spin-orbit case is
  !        E DOS(E) PDOS(E)
  !   The format for the collinear, spin-polarized case is
  !        E DOSup(E) DOSdw(E)  PDOSup(E) PDOSdw(E)
  !   The format for the non-collinear, non spin-orbit case is
  !        E DOS(E) PDOSup(E) PDOSdw(E)
  !
  !   In the collinear case and the non-collinear, non spin-orbit case
  !   projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   (one file per atomic wavefunction found in the pseudopotential file)
  !   - The format for the collinear, spin-unpolarized case is
  !        E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
  !     where LDOS = \sum m=1,2l+1 PDOS_m(E)
  !     and PDOS_m(E) = projected DOS on atomic wfc with component m
  !   - The format for the collinear, spin-polarized case and the
  !     non-collinear, non spin-orbit case is as above with
  !     two components for both  LDOS(E) and PDOS_m(E)
  !
  !   In the non-collinear, spin-orbit case (i.e. if there is at least one
  !   fully relativistic pseudopotential) wavefunctions are projected
  !   onto eigenstates of the total angular-momentum.
  !   Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l_j),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   and j is the value of the total angular momentum.
  !   In this case the format is
  !      E LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
  !
  !   All DOS(E) are in states/eV plotted vs E in eV
  !
  !   If the kresolveddos option is used, the k-point index is prepended
  !   to the formats above, e.g. (collinear, spin-unpolarized case)
  !        ik E DOS(E) PDOS(E)
  !
  !   If the local DOS(E) is computed (tdosinboxes),
  !   projections are written to file "filproj" if given as input.
  !   Volumes are written as xsf files with 3D datagrids, valued 1.0
  !   inside the box volume and 0 outside, named filpdos.box#n.xsf
  !   The local DOS(E) is written to file "filpdos".ldos_boxes, with format
  !      E totDOS(E) SumLDOS(E) LDOS_1(E) ... LDOS_n(E)
  !
  !  Order of m-components for each l in the output:
  !
  !  1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)
  !
  !  where phi is the polar angle:x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
  !  This is determined in file flib/ylmr2.f90 that calculates spherical harm.
  !      L=1 :
  !  1 pz     (m=0)
  !  2 px     (real combination of m=+/-1 with cosine)
  !  3 py     (real combination of m=+/-1 with sine)
  !      L=2 :
  !  1 dz2    (m=0)
  !  2 dzx    (real combination of m=+/-1 with cosine)
  !  3 dzy    (real combination of m=+/-1 with sine)
  !  4 dx2-y2 (real combination of m=+/-2 with cosine)
  !  5 dxy    (real combination of m=+/-1 with sine)
  !
  ! Important notice:
  !
  !    The tetrahedron method is presently not implemented.
  !    Gaussian broadening is used in all cases:
  !    - if degauss is set to some value in namelist &inputpp, that value
  !      (and the optional value for ngauss) is used
  !    - if degauss is NOT set to any value in namelist &inputpp, the
  !      value of degauss and of ngauss are read from the input data
  !      file (they will be the same used in the pw.x calculations)
  !    - if degauss is NOT set to any value in namelist &inputpp, AND
  !      there is no value of degauss and of ngauss in the input data
  !      file, degauss=DeltaE (in Ry) and ngauss=0 will be used
  ! Obsolete variables, ignored:
  !   io_choice
  !   smoothing
  !
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
  CHARACTER (len=256) :: filpdos, filproj, io_choice, outdir, chapprox, spectype
  REAL (DP)      :: Emin, Emax, DeltaE, degauss1, smoothing, broadening
  INTEGER :: ngauss1, ios, Nener
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes, fixocc
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)

  !
  ! for GWW
  INTEGER :: iun, idum
  REAL(DP) :: rdum1,rdum2,rdum3
  LOGICAL :: lex, lgww
  !
  !
  !for XAS
  LOGICAL :: readcorerep
  !
  NAMELIST / inputpp / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, Nener, DeltaE, io_choice, smoothing, filpdos, filproj, &
             lgww, & !if .true. use GW QP energies from file bands.dat
             kresolveddos, tdosinboxes, n_proj_boxes, irmin, irmax, plotboxes, &
             readcorerep, chapprox, broadening, spectype, fixocc
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PROJWFC' )
  !
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
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  lgww   = .false.
  kresolveddos = .false.
  tdosinboxes = .false.
  plotboxes   = .false.
  n_proj_boxes= 1
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
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id )
  !
  IF (ios /= 0) CALL errore ('projwfc', 'reading inputpp namelist', abs (ios) )
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
  CALL mp_bcast( tdosinboxes, ionode_id )
  CALL mp_bcast( n_proj_boxes, ionode_id )
  CALL mp_bcast( irmin, ionode_id )
  CALL mp_bcast( irmax, ionode_id )
  CALL mp_bcast( readcorerep, ionode_id )
  CALL mp_bcast( chapprox, ionode_id )
  CALL mp_bcast( broadening, ionode_id )
  CALL mp_bcast( spectype, ionode_id )
  CALL mp_bcast( fixocc, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  CALL openfil_pp ( )
  !
  if( readcorerep ) then
     ! read atomic matrix elements for core-valence position
     if(ionode) then
        call read_corerepair( 5, corerep )
        write(stdout,*) 'done reading corerepair'
     endif
     call bcast_corerepair( corerep, ionode_id )
  endif

  !
  !   decide Gaussian broadening
  !
  lgauss=.true.
  IF (degauss1/=0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ELSEIF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  ELSE
     degauss=0.1d0/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !
!!$  IF ( tdosinboxes ) THEN
!!$     IF( nproc_ortho > 1 ) THEN
!!$        CALL errore ('projwfc', 'nproc_ortho > 1 not yet implemented', 1)
!!$     ELSE
!!$        CALL projwave_boxes (filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes)
!!$     ENDIF
!!$  ELSE
!!$     IF (noncolin) THEN
!!$        CALL projwave_nc(filproj, lsym )
!!$     ELSE
!!$        IF( nproc_ortho > 1 ) THEN
!!$           CALL pprojwave (filproj, lsym)
!!$        ELSE
  CALL projwave (filproj, lsym, lgww)
!!$        ENDIF
!!$     ENDIF
!!$  ENDIF
  !
!!$  IF ( ionode ) THEN
!!$     IF ( tdosinboxes ) THEN
!!$        CALL partialdos_boxes (Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
!!$     ELSE
!!$        IF ( lsym ) THEN
!!$           !
!!$           IF (noncolin) THEN
!!$              CALL partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
!!$           ELSE
     write(*,*) "calling partialdos..."
     CALL partialdos (Emin, Emax, DeltaE, Nener, broadening, chapprox, spectype, fixocc, kresolveddos, filpdos)
!!$           ENDIF
!!$           !
!!$        ENDIF
!!$     ENDIF
!!$  ENDIF
  !
  CALL stop_pp
  !
END PROGRAM pw2xas 
!
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
MODULE projections_nc
  USE kinds, ONLY : DP

  TYPE wfc_label_nc
     INTEGER na, n, l, m, ind
     REAL (DP) jj
  END TYPE wfc_label_nc
  TYPE(wfc_label_nc), ALLOCATABLE :: nlmchi(:)

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)

END MODULE projections_nc
!
MODULE projections_ldos
  USE kinds, ONLY : DP
  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
END MODULE projections_ldos
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
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, find_free_unit
  USE spin_orb, ONLY: lspinorb
  USE wavefunctions_module, ONLY: evc
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
  !
  WRITE( stdout, '(/5x,"Calling projwave .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
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
  ALLOCATE( proj_aux (ncp, nbnd, 3, nkstot) )
!!$  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj_aux      = (0.d0, 0.0d0)
!!$  proj_aux  = (0.d0, 0.d0)
  !
!!$  IF (.not. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) )
!!$  ALLOCATE(wfcatom (npwx, natomwfc) )
!!$  ALLOCATE(overlap (natomwfc, natomwfc) )
!!$  overlap= (0.d0,0.d0)
!!$  !
!!$  IF ( gamma_only ) THEN
!!$     ALLOCATE(roverlap (natomwfc, natomwfc) )
!!$     roverlap= 0.d0
!!$  ENDIF
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
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)

!!$     CALL atomic_wfc (ik, wfcatom)

     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL calbec ( npw, vkb, evc, becp)

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
!!$    allocate(temp(nbnd,ncp),tcr(nbeta,ncp))
     DO ixyz=1, 3
        CALL ZGEMM( 'T', 'C', ncp, nbnd, &
                    nbeta, (1.0d0,0.0d0), corerep%core(1)%matrix(1,1,ixyz), nbeta, &
                    psibeta, nbnd, &
                    (0.0d0,0.0d0), proj_aux(1,1,ixyz,ik), ncp )
!!$        tcr(:,:)=corerep%core(1)%matrix(:,:,ixyz)
!!$        proj_aux(:,:,ixyz,ik) = MATMUL(psibeta, tcr)

!!$        CALL ZGEMM( 'C', 'N', ncp, nbnd, nbeta, &
!!$             CMPLX(1.0d0), corerep%core(1)%matrix(1,1,ixyz), &
!!$             corerep%core(1)%nproj1, psibeta, nbeta,  &
!!$             CMPLX(0.0d0), proj_aux (1, 1, ixyz, ik), ncp )
     ENDDO
!!$     deallocate(temp,tcr)
     ! on k-points
  ENDDO
!!$  CALL poolrecover (proj_aux, ncp*nbnd*3, nkstot, nks)
  !
  !
  RETURN
  !
END SUBROUTINE projwave
!-----------------------------------------------------------------------
SUBROUTINE  partialdos (Emin, Emax, DeltaE, Nener, broadening, chapprox, spectype, fixocc, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, nks, degauss, ngauss, lgauss, nelec
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd, wg
  USE constants, ONLY: rytoev
  USE mp_global, ONLY:mp_bcast, mp_sum, inter_pool_comm
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
  INTEGER :: ik, ibnd,  m, &
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
!!$  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
!!$  !
!!$  IF (kresolveddos) THEN
!!$     IF ( nspin==2 ) THEN
!!$        nkseff=nkstot/2
!!$     ELSE
!!$        nkseff=nkstot
!!$     ENDIF
!!$  ELSE
!!$     nkseff=1
!!$  ENDIF


  if(trim(chapprox) .eq. 'XCH') then
     ef = efermig (et, nbnd, nkstot, nelec-1, wk, degauss, ngauss, 0, isk)
  else
     ef = efermig (et, nbnd, nkstot, nelec, wk, degauss, ngauss, 0, isk)
  endif
  CALL mp_bcast(ef, ionode_id)
  ALLOCATE(occ(nbnd,nkstot))
  occ=0.0d0
  ehomo=et(1,1)
  elumo=et(nbnd,1)
  if(ionode) then
     DO ik=1, nkstot
        DO ibnd=1,nbnd
!!$        write(*,*) "occ: ik, ibnd=", ik, ibnd
           if(et(ibnd,ik) .gt. ehomo .and. et(ibnd,ik) .le. ef ) ehomo=et(ibnd,ik)
           if(et(ibnd,ik) .lt. elumo .and. et(ibnd,ik) .ge. ef ) elumo=et(ibnd,ik)
           if(trim(spectype) .eq. 'RAW') then
              occ(ibnd,ik)=1.0d0
           else if (trim(spectype) .eq. 'XAS') then
              if(fixocc) then
                 occ(ibnd,ik)= 1 - (wg(ibnd,ik)/wk(ik))
              else
                 occ(ibnd,ik)=1-wgauss((ef-et(ibnd,ik)) / degauss, ngauss)
              endif
           else if (trim(spectype) .eq. 'XES') then
              if(fixocc) then
                 occ(ibnd,ik)= wg(ibnd,ik)/wk(ik)
              else
                 occ(ibnd,ik)=wgauss((ef-et(ibnd,ik)) / degauss, ngauss)
              endif
           endif
        ENDDO
     ENDDO
  write(*,*) "done occ..."
  endif
  CALL poolscatter(nbnd, nkstot, occ, nks, occ)
  CALL mp_bcast(ehomo, ionode_id)
  CALL mp_bcast(elumo, ionode_id)

  DeltaE = DeltaE/rytoev   !convert to Ryd
  broadening = broadening/rytoev !convert to Ryd
!  if(trim(spectype) .eq. 'XAS' .or. trim(spectype) .eq. 'XES') then
!     DeltaE=DeltaE-elumo
!  endif
  !
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  if(ionode) then
     DO ik = 2, nkstot
        Elw = min (Elw, et (1, ik) )
        Eup = max (Eup, et (nbnd, ik) )
     ENDDO
     IF (broadening/=0.d0) THEN
        Eup = Eup + 5d0 * broadening
        Elw = Elw - 5d0 * broadening
     ENDIF
     Emin=Emin/rytoev
     Emax=Emax/rytoev
     if(Emin .le. -1000.d0) Emin = max (Emin, Elw)
     if(Emax .ge. +1000.d0) Emax = min (Emax, Eup)
     if(trim(spectype) .eq. 'RAW') then
        Emin = Elw
        Emax = Eup
     endif
  endif
  CALL mp_bcast(Emin, ionode_id)
  CALL mp_bcast(Emax, ionode_id)
  dE = (Emax - Emin)/dble(Nener)
  do ie=1,nener
    ener(ie) = Emin + dble(ie-1)*dE
  enddo
  write(*,*) "done erange..."


  !
!!$  ALLOCATE (pdos(0:ne,natomwfc,nspin,nkseff))
!!$  ALLOCATE (dostot(0:ne,nspin,nkseff), pdostot(0:ne,nspin,nkseff), ldos(0:ne,nspin,nkseff) )
!!$  pdos(:,:,:,:) = 0.d0
!!$  dostot(:,:,:) = 0.d0
!!$  pdostot(:,:,:)= 0.d0
!!$  !

  ALLOCATE(xas_xyz_sp(3,nspin),spec(Nener,4,nspin),tmp_spec(Nener,4))
  spec(1:nener,1:4,1:nspin)=0.0d0
  tmp_spec(1:nener,1:4)=0.0d0
  xas_xyz_sp = 0.0d0
  ncp=size(proj_aux,1)
  current_spin = 1
!!$  ie_delta = 5 * degauss / DeltaE + 1

  DO ik = 1,nks
!!$     write(*,*) "ik=",ik

     ! use true weights
     wkeff=wk(ik)
     ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
     ikeff=1
     !
     IF ( nspin == 2 ) current_spin = isk ( ik )
     DO ibnd = 1, nbnd
!!$     write(*,*) "ibnd=",ibnd
        etev = et(ibnd,ik)
        forall(j=1:3) xas_xyz_sp(j,current_spin) = sum(conjg(proj_aux(1:ncp,ibnd,j,ik))*proj_aux(1:ncp,ibnd,j,ik))*wk(ik)*occ(ibnd,ik)
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
  fileout = trim(filpdos)//".pw.raw"
  if(trim(spectype) .eq. 'XAS') then 
     fileout = trim(filpdos)//".pw.xas"
  else if(trim(spectype) .eq. 'XES') then
     fileout = trim(filpdos)//".pw.xes"
  endif
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
     enddo
     close(iunout)
  endif
  DEALLOCATE (proj_aux, spec)
  !
  RETURN
END SUBROUTINE partialdos
!
!!$!-----------------------------------------------------------------------
!!$SUBROUTINE write_io_header(filplot, iunplot, title, nr1x, nr2x, nr3x, &
!!$           nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
!!$           nkstot,nbnd,natomwfc)
!!$   !-----------------------------------------------------------------------
!!$
!!$   USE kinds, ONLY: DP
!!$   USE ions_base, ONLY : zv, atm, tau, ityp
!!$   USE noncollin_module, ONLY: noncolin
!!$   USE spin_orb, ONLY: lspinorb
!!$
!!$   IMPLICIT NONE
!!$   CHARACTER (len=*) :: filplot
!!$   CHARACTER (len=*) :: title
!!$   INTEGER :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav
!!$   REAL(DP) :: celldm (6), gcutm, dual, ecutwfc, at(3,3)
!!$   INTEGER :: iunplot, ios, na, nt, i
!!$   INTEGER :: nkstot,nbnd,natomwfc
!!$   !
!!$   IF (filplot == ' ') CALL errore ('write_io_h', 'filename missing', 1)
!!$
!!$   OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
!!$         STATUS = 'unknown', ERR = 101, IOSTAT = ios)
!!$   101     CALL errore ('write_io_h', 'opening file '//trim(filplot), abs (ios) )
!!$   WRITE (iunplot, '(a)') title
!!$   WRITE (iunplot, '(8i8)') nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
!!$   WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
!!$   IF  (ibrav == 0) THEN
!!$       WRITE ( iunplot, * ) at(:,1)
!!$       WRITE ( iunplot, * ) at(:,2)
!!$       WRITE ( iunplot, * ) at(:,3)
!!$   ENDIF
!!$   WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, 9
!!$   WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
!!$         (nt, atm (nt), zv (nt), nt=1, ntyp)
!!$   WRITE (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
!!$         (tau (i, na), i = 1, 3), ityp (na), na = 1, nat)
!!$   WRITE (iunplot, '(3i8)') natomwfc, nkstot, nbnd
!!$   WRITE (iunplot, '(2l5)') noncolin, lspinorb
!!$
!!$   RETURN
!!$END SUBROUTINE write_io_header
!!$!
!!$!-----------------------------------------------------------------------
!!$FUNCTION compute_mj(j,l,m)
!!$   !-----------------------------------------------------------------------
!!$   USE kinds, ONLY: DP
!!$   IMPLICIT NONE
!!$   !
!!$   REAL(DP) :: compute_mj, j
!!$   INTEGER  :: l, m
!!$
!!$   IF (abs(j-l-0.5d0)<1.d-4) THEN
!!$       compute_mj=m+0.5d0
!!$   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
!!$      compute_mj=m-0.5d0
!!$   ELSE
!!$      CALL errore('compute_mj','l and j not compatible',1)
!!$   ENDIF
!!$
!!$   RETURN
!!$END FUNCTION compute_mj
!!$!
!!$!-----------------------------------------------------------------------
!!$SUBROUTINE  write_proj (filename, projs)
!!$  !-----------------------------------------------------------------------
!!$  !
!!$  USE kinds
!!$  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir
!!$  USE basis,            ONLY : natomwfc
!!$  USE cell_base
!!$  USE klist,            ONLY : wk, xk, nkstot, nelec
!!$  USE noncollin_module, ONLY : noncolin
!!$  USE lsda_mod,         ONLY : nspin, isk
!!$  USE ener,             ONLY : ef
!!$  USE wvfct,            ONLY : et, nbnd
!!$  USE iotk_module
!!$  IMPLICIT NONE
!!$
!!$  CHARACTER(*),  INTENT(in) :: filename
!!$  COMPLEX(DP),   INTENT(in) :: projs(natomwfc,nbnd,nkstot)
!!$  !
!!$  CHARACTER(256)          :: tmp
!!$  CHARACTER(iotk_attlenx) :: attr
!!$  INTEGER :: ik, ik_eff, ia, ierr, num_k_points
!!$
!!$!
!!$! subroutine body
!!$!
!!$
!!$  tmp = trim( tmp_dir ) // trim( prefix ) // '.save/' //trim(filename)
!!$  !
!!$  CALL iotk_open_write(iun, FILE=trim(tmp), ROOT="ATOMIC_PROJECTIONS", IERR=ierr )
!!$  IF ( ierr /= 0 ) RETURN
!!$  !
!!$  !
!!$  num_k_points = nkstot
!!$  IF ( nspin == 2 ) num_k_points = nkstot / 2
!!$  !
!!$  CALL iotk_write_begin(iun, "HEADER")
!!$  !
!!$  CALL iotk_write_dat(iun, "NUMBER_OF_BANDS", nbnd)
!!$  CALL iotk_write_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
!!$  CALL iotk_write_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin)
!!$  CALL iotk_write_dat(iun, "NON-COLINEAR_CALCULATION",noncolin)
!!$  CALL iotk_write_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
!!$  CALL iotk_write_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
!!$  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.true.  )
!!$  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
!!$  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.true.  )
!!$  CALL iotk_write_empty (iun,  "UNITS_FOR_ENERGY", ATTR=attr)
!!$  CALL iotk_write_dat(iun, "FERMI_ENERGY", ef )
!!$  !
!!$  CALL iotk_write_end(iun, "HEADER")
!!$  !
!!$  !
!!$  CALL iotk_write_dat(iun, "K-POINTS", xk(:,1:num_k_points) , COLUMNS=3 )
!!$  CALL iotk_write_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points), COLUMNS=8 )
!!$  !
!!$  CALL iotk_write_begin(iun, "EIGENVALUES")
!!$  !
!!$  DO ik=1,num_k_points
!!$     !
!!$     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
!!$     !
!!$     IF ( nspin == 2 ) THEN
!!$        !
!!$        ik_eff = ik + num_k_points
!!$        !
!!$        CALL iotk_write_dat( iun, "EIG.1", et(:,ik) )
!!$        CALL iotk_write_dat( iun, "EIG.2", et(:,ik_eff) )
!!$        !
!!$     ELSE
!!$        !
!!$        CALL iotk_write_dat( iun, "EIG", et(:,ik) )
!!$        !
!!$     ENDIF
!!$     !
!!$     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
!!$     !
!!$  ENDDO
!!$  !
!!$  CALL iotk_write_end(iun, "EIGENVALUES")
!!$
!!$  !
!!$  ! main loop atomic wfc
!!$  !
!!$  CALL iotk_write_begin(iun, "PROJECTIONS")
!!$  !
!!$  DO ik=1,num_k_points
!!$     !
!!$     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
!!$     !
!!$     IF ( nspin == 2 ) THEN
!!$        !
!!$        CALL iotk_write_begin ( iun, "SPIN.1" )
!!$           !
!!$           DO ia = 1, natomwfc
!!$               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
!!$           ENDDO
!!$           !
!!$        CALL iotk_write_end ( iun, "SPIN.1" )
!!$        !
!!$        ik_eff = ik + num_k_points
!!$        !
!!$        CALL iotk_write_begin ( iun, "SPIN.2" )
!!$           !
!!$           DO ia = 1, natomwfc
!!$               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
!!$           ENDDO
!!$           !
!!$        CALL iotk_write_end ( iun, "SPIN.2" )
!!$        !
!!$     ELSE
!!$        !
!!$        DO ia = 1,natomwfc
!!$            CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
!!$        ENDDO
!!$        !
!!$     ENDIF
!!$     !
!!$     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
!!$     !
!!$  ENDDO
!!$  !
!!$  CALL iotk_write_end(iun, "PROJECTIONS")
!!$
!!$  !
!!$  ! closing the file
!!$  !
!!$  CALL iotk_close_write(iun)
!!$
!!$END SUBROUTINE write_proj
!!$
!!$!
