  module wfc_shirley

  ! this module contains details of the Fourier space wave functions
  ! used in expanding the optimal Shirley basis set

  ! David Prendergast

  use kinds, only : dp


  implicit none

  public

  ! two arrays used to store wave functions from different k-points
  complex(dp),allocatable :: evc_1(:,:)
  complex(dp),allocatable :: evc_2(:,:)
  integer :: nwordwfc_gamma, iunwfcr

  ! G-space details at the Gamma-point
  integer :: npw_gamma, nbnd_gamma
  integer :: igwx_gamma
  real(dp) :: ecutwfc_gamma
  integer,allocatable :: igk_l2g(:,:)
  integer,allocatable :: igk_l2g_gamma(:)

  ! inverse mapped to unmapped
  type kmap
    integer :: nmap
    real(dp),pointer :: xmap(:,:)
  end type kmap
  type(kmap),allocatable :: xk_map(:)


  contains


! ----------------------------------------------------------------------
  subroutine expand_Gspace( ecut_in )
! ----------------------------------------------------------------------

  ! determines the largest G-vector index over all k-points and expands
  ! the G-space cut-off to include this G-vector at the Gamma-point

  use gvect,     only : ngm, gcutm, g, ig_l2g
  use cell_base, only : bg, tpiba2
  use klist,     only : nkstot, nks, xk, wk, ngk
  use wvfct,     only : ecutwfc, npw, igk, g2kin, npwx
  use mp,        only : mp_max
  use io_global, only : stdout
  use io_files,  only : iunigk

  real(dp),intent(in) :: ecut_in

  integer :: ik, igwx, i, npwx_store
  integer :: ikm, imap
  real(dp) :: ecut
  real(dp) :: xk_gamma(3)
  real(dp),allocatable :: xk_exp(:,:)
  integer :: ierr


  ! summary of input details
  call summary

  ! old cut-off in Gspace units
  ecut = ecutwfc  / tpiba2

  if( .NOT. allocated(igk_l2g) ) call errore('expand_Gspace','igk_l2g not allocated',1)

  ! max G index
  igwx = max( 1, maxval( igk_l2g ) )
  call mp_max( igwx )

  ! Gamma-point vector
  xk_gamma = 0.d0

  ! allocate these arrays larger than necessary
  if( allocated( igk ) ) deallocate( igk )
  if( allocated( g2kin ) ) deallocate( g2kin )
  allocate( igk(ngm), g2kin(ngm) )

  ! ---------------------------------------------
  ! increase npwx to the max possible
  npwx = ngm
  ! ---------------------------------------------

  write(stdout,*)
  write(stdout,*) ' expand_Gspace -'
  write(stdout,*) '   increasing cut-off to include G-vectors from all k-points'
  write(stdout,*) '     at the Gamma-point'
  write(stdout,*)
  write(stdout,'(2x,a6,2a12,2a10,a10)') ' iter', ' max-G-index', ' max-k-index', ' new Ecut', ' old Ecut', ' npwx'
  call flush_unit(stdout)

  if( allocated(igk_l2g_gamma) ) deallocate( igk_l2g_gamma )
  allocate( igk_l2g_gamma(ngm) )

  i = 0
  igwx_gamma = 0
  do while( igwx_gamma < igwx .and. ecut <= gcutm )
    CALL gk_sort (xk_gamma, ngm, g, ecut, npw_gamma, igk, g2kin)
    call gk_l2gmap (ngm, ig_l2g(1), npw_gamma, igk, igk_l2g_gamma(1))

    ! for the Gamma-point the ordering of ig_l2g is correct
    ! or is it?
    igwx_gamma = maxval( igk_l2g_gamma(1:npw_gamma) )
    call mp_max( igwx_gamma )

    i=i+1
    ! new cut-off
    ecutwfc_gamma = ecut * tpiba2

    ! report
    write(stdout,'(2x,i6,2i12,2f10.5,i10)') &
      i, igwx_gamma, igwx, ecutwfc_gamma, ecutwfc, npw_gamma
  call flush_unit(stdout)

    ecut = ecut + 1.d0
  enddo
  ecut = ecut - 1.d0 ! correct for additional 1.d0
  write(stdout,*)
  write(stdout,*) ' all-encompassing cut-off      = ', ecutwfc_gamma, npw_gamma
  write(stdout,*)
  call flush_unit(stdout)

  if( ecut > gcutm ) call errore('expand_Gspace','expanded cut-off too large',1)

  ! from this point on we assume only one k-point k=0
  ! and for the purposes of generating a new FFT grid 
  ! we set nkstot and nks to 1
  nkstot = 1
  nks = 1
  xk(1:3,1) = 0.d0  ! N.B. Important to set this to Gamma point
                    !      It is stored in the save file
  wk(1) = 2.d0      ! the weight is two for non-spin-polarized
  ! spin should have no meaning for the optimal basis - it should work for both
  ! spin-channels, so really wk is just a dummy variable now
  ecutwfc = ecutwfc_gamma

  !
  ! if the input cut-off is non-zero
  if( ecut_in > 0.d0 .and. ecut_in < ecutwfc ) then
    ecutwfc = min( ecut_in, ecutwfc )
    ecut = ecutwfc / tpiba2
    write(stdout,*) ' using user input cut-off instead: ecutwfc = ', &
      ecut_in
  endif
  !
  ! update FFT grids based on the new cut-off
  write(stdout,*) ' update FFT grids '
  call flush_unit(stdout)
  call update_shirley_fft

  write(stdout,*) 'local number of G-vectors before FFT update: ', igwx_gamma

    igwx_gamma = maxval( igk_l2g_gamma(1:npw_gamma) )
    call mp_max( igwx_gamma )

  write(stdout,*) 'local number of G-vectors after FFT update: ', igwx_gamma

  CALL gk_sort (xk_gamma, ngm, g, ecut, npw_gamma, igk, g2kin)

  write(stdout,*) ' igk, g2kin, igk_l2g_gamma ', size(igk), size(g2kin), size(igk_l2g_gamma)

  ! is this even necessary?
  ! allocate the correct amount of space and reinitialize indexing arrays
  if( allocated( igk ) ) deallocate( igk )
  if( allocated( g2kin ) ) deallocate( g2kin )
  if( allocated(igk_l2g_gamma) ) deallocate( igk_l2g_gamma )
  allocate( igk(npw_gamma), g2kin(npw_gamma), igk_l2g_gamma(npw_gamma) )
  ! ------------------------------
  ! set npwx - very important
  ! it seems that gk_sort uses npwx to size incoming arrays, so it should match the sizes of igk, g2kin
  npwx = npw_gamma
  ! ------------------------------
  write(stdout,*) ' igk, g2kin, igk_l2g_gamma ', size(igk), size(g2kin), size(igk_l2g_gamma)
  CALL gk_sort (xk_gamma, ngm, g, ecut, npw_gamma, igk, g2kin)
  call gk_l2gmap (ngm, ig_l2g(1), npw_gamma, igk, igk_l2g_gamma(1))

  write(stdout,*) ' all-encompassing cut-off = ', ecutwfc_gamma, npwx
  write(stdout,*)

  write(stdout,*) 'leaving expand_Gspace'
  return
  end subroutine expand_Gspace


! ----------------------------------------------------------------------
  subroutine update_shirley_fft( )
! ----------------------------------------------------------------------
  !
  ! update all fft information based on the new cut-off ecutwfc_gamma
  !
  use wvfct,                  only : ecutwfc
  use gvect
  use grid_dimensions,        only : nr1, nr2, nr3
  use grid_subroutines,       only : realspace_grids_init
  use gvecs,                  only : dual, doublegrid, gcutms
  use smooth_grid_dimensions, only : nr1s, nr2s, nr3s
  use recvec_subs,            only : ggen
  use cell_base,              only : tpiba2, at, bg
  use control_flags,          only : gamma_only


  !
  ! deallocate space used for FFT's and g-vectors
  !
  call deallocate_fft
  !
  ! ... Compute the cut-off of the G vectors
  !
  gcutm = dual * ecutwfc / tpiba2
  !
  doublegrid = ( dual > 4.D0 )
  !
  IF ( doublegrid ) THEN
     !
     gcutms = 4.D0 * ecutwfc / tpiba2
     !
  ELSE
     !
     gcutms = gcutm
     !
  END IF
  !
  ! initialize FFT grid dimensions
  nr1=0 ; nr2=0; nr3=0
  nr1s=0 ; nr2s=0; nr3s=0
  !
  ! ... calculate dimensions of the FFT grid
  !
  !call set_fft_dim
  call realspace_grids_init( at, bg, gcutm, gcutms )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  call allocate_fft
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  call ggen ( gamma_only, at, bg )
  !
  call summary
  !
  return
  end subroutine update_shirley_fft

! ----------------------------------------------------------------------
  subroutine init_shirley_wfc_1( )
! ----------------------------------------------------------------------

  ! Initialize the memory and disk for shirley wave functions

  use wvfct,            only : nbnd
  use noncollin_module, only : npol
  use io_files,         only : diropn
  use io_global, only : stdout

  logical :: exst, opend
  integer :: ierr
  integer,external :: freeunit

  ! copy number of bands
  nbnd_gamma = nbnd

  allocate( evc_1(npw_gamma,nbnd_gamma), stat=ierr )
  if( ierr/=0 ) call errore('init_shirley_wfc','unable to allocate space for shirley wave function #1',1)

  ! open a new set of files for reordered wave functions
  nwordwfc_gamma = 2 * nbnd_gamma * npw_gamma * npol
  iunwfcr = freeunit()
  call diropn( iunwfcr, 'wfcr', nwordwfc_gamma, exst )

  return
  end subroutine init_shirley_wfc_1


! ----------------------------------------------------------------------
  subroutine init_shirley_wfc_2( )
! ----------------------------------------------------------------------

  ! Initialize the memory and disk for shirley wave functions

  integer :: ierr

  allocate( evc_2(npw_gamma,nbnd_gamma), stat=ierr )
  if( ierr/=0 ) call errore('init_shirley_wfc','unable to allocate space for shirley wave function #2',1)

  return
  end subroutine init_shirley_wfc_2

! ----------------------------------------------------------------------
  subroutine reorder_shirley_wfc( ik )
! ----------------------------------------------------------------------

  ! reorder the wave function coefficients so that k-point ik has the
  ! same ordering as the Gamma-point

  ! David Prendergast, UCB, Jan 2007

  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_global, ONLY : npool, mpime
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  USE mp_wave,   ONLY : mergewf, splitwf
  USE io_global, ONLY : ionode_id, stdout
  USE wvfct,     ONLY : nbnd
  USE klist, ONLY : ngk
  USE gvect, ONLY : ig_l2g
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module, ONLY : noncolin
  USE parallel_include

  implicit none

  integer,intent(in) :: ik

  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),parameter :: one  = cmplx(1.d0,0.d0)

  integer :: ikp
  integer :: j, i, ierr
  complex(dp),allocatable :: wtmp(:)
  complex(dp) :: dotprod, dotprodw

  write(stdout,*) 'reorder_shirley_wfc', ik

  ! error checking
  if( npool > 1 ) &
    call errore('reorder_shirley_wfc','npool must be set to 1',1)
  if( noncolin ) &
    call errore('reorder_shirley_wfc','not implemented for noncolin',1)

  if( allocated(wtmp) ) deallocate(wtmp)
  allocate( wtmp(igwx_gamma), stat=ierr )
  if( ierr /= 0 ) &
    call errore('reorder_shirley_wfc','unable to allocate temporary space',igwx_gamma)

! check for igk_l2g here
  if( .NOT. allocated(igk_l2g) ) call errore('reorder_shirley_wfc','igk_l2g not allocated',1)

  write(stdout,*) ' size ngk ', size(ngk), ngk
  do j=1,nbnd
    ! zero the temp array
    wtmp = zero

    ! merge wave function from evc into wtmp
    call mergewf(evc(:,j), wtmp, ngk(ik), igk_l2g(:,ik), me_pool, nproc_pool, root_pool, intra_pool_comm)

    ! check
    dotprod = dot_product( evc(1:ngk(ik),j), evc(1:ngk(ik),j) )
#ifdef __PARA
    call mp_sum( dotprod, intra_pool_comm )
#endif
    dotprodw = dot_product(wtmp, wtmp)
    if( abs(real(dotprodw-dotprod)) > 1.d-12 ) then
      write(stdout,'(a,i6,3e12.5)') 'warning: bad merged norm for band ', &
        j, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
    endif

    ! zero the former wave function
    evc_1(:,j) = zero

    ! split wave function from wtmp into evc_1
    call splitwf(evc_1(:,j), wtmp, npw_gamma, igk_l2g_gamma(:), me_pool, nproc_pool, root_pool, intra_pool_comm)

    dotprod = dot_product( evc_1(1:npw_gamma,j), evc_1(1:npw_gamma,j) )
#ifdef __PARA
    call mp_sum( dotprod, intra_pool_comm )
#endif
    if( abs(real(dotprodw-dotprod)) > 1.d-12 ) then
      write(stdout,'(a,i6,3e12.5)') 'warning: bad split norm for band ', &
        j, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
    endif
  enddo

  deallocate( wtmp )

  return
  end subroutine reorder_shirley_wfc


! ----------------------------------------------------------------------
  subroutine write_shirley_wfc( ik )
! ----------------------------------------------------------------------

  integer,intent(in) :: ik

  call davcio( evc_1, nwordwfc_gamma, iunwfcr, ik, +1 )

  return
  end subroutine write_shirley_wfc


! ---------------------------------------------------------------------- 
  subroutine read_shirley_wfc_1( ik )
! ---------------------------------------------------------------------- 
  integer,intent(in) :: ik
  integer,save :: ik_current=0

  if( ik/=ik_current ) then
    call davcio( evc_1, nwordwfc_gamma, iunwfcr, ik, -1 )
    ik_current = ik
  endif

  return
  end subroutine read_shirley_wfc_1


! ---------------------------------------------------------------------- 
  subroutine read_shirley_wfc_2( jk )
! ---------------------------------------------------------------------- 
  integer,intent(in) :: jk
  integer,save :: jk_current=0

  if( jk/=jk_current ) then
    call davcio( evc_2, nwordwfc_gamma, iunwfcr, jk, -1 )
    jk_current = jk
  endif

  return
  end subroutine read_shirley_wfc_2


! ---------------------------------------------------------------------- 
  subroutine delete_shirley_wfc
! ---------------------------------------------------------------------- 

  close(iunwfcr,status='delete')

  return
  end subroutine delete_shirley_wfc


! ---------------------------------------------------------------------- 
  function dotprod_shirley_wfc( i, j )
! ---------------------------------------------------------------------- 
  complex(dp) :: dotprod_shirley_wfc
  integer,intent(in) :: i, j
  complex(dp) :: ZDOTC

  !dotprod_shirley_wfc = dot_product( evc_2(:,i), evc_1(:,j) )
  dotprod_shirley_wfc = ZDOTC( size(evc_2,1), evc_2(1,i), 1, evc_1(1,j), 1 )

  return
  end function dotprod_shirley_wfc


! ---------------------------------------------------------------------- 
  subroutine fix_shirley_planewaves
! ---------------------------------------------------------------------- 

  ! note that nkstot should be 1 now and ecutwfc should be larger

  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! and that igk_l2g is set up for indexing wave function components

  use io_global, only : stdout
  use cell_base, only : tpiba2
  use klist,     only : xk, nks, nkstot, ngk
  USE gvect,     ONLY : ngm, g, ig_l2g
  USE wvfct,     ONLY : ecutwfc, npw, npwx, igk, g2kin

  integer :: ik

  write(stdout,*) ' establishing correct dimensions for Gamma-point functions'
  write(stdout,*) ' nkstot  = ', nkstot
  write(stdout,*) ' nks     = ', nks   
  write(stdout,*) ' ecutwfc = ', ecutwfc   
  write(stdout,*) ' ngm = ', ngm   
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)
  if( allocated( igk_l2g ) ) deallocate( igk_l2g )
  allocate( igk_l2g(npwx,nkstot) )
  do ik=1,nkstot
    call gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,ik))
  enddo

  return
  end subroutine fix_shirley_planewaves


! ---------------------------------------------------------------------- 
  subroutine map_shirley_bz( ndim, band_subset )
! ---------------------------------------------------------------------- 

  ! This is an important mapping of the available wave functions from
  ! their original k-points to k in [0,1]
  ! Note that we include the zone boundary at 1 aswell
  ! This means for example that the Gamma-point will be mapped
  ! to itself (000) and also (001), (010), (011), (100), (101), (110), (111)
  ! If we have a system of reduced dimensions (e.g. 1D-periodic) then the
  ! possible mappings is reduced, e.g. (000) and (001) for 1D along z

  ! I need to add inversion symmetry to this routine

  use constants, only : tpi
  use klist,     only : xk, nks, nkstot, ngk, wk
  use cell_base, only : ibrav, at, bg, symm_type, tpiba2, tpiba, alat
  use symm_base, only : s, sname, ftau, nsym, t_rev, irt, ftau, invsym
  USE gvect, ONLY : ngm, ngm_g, gcutm, g, ig_l2g, nl, gg
  use grid_dimensions, only : nrxx
  use gvecs, only : nls, dual, doublegrid, gcutms
  use smooth_grid_dimensions, only : nrxxs
  USE wvfct,         ONLY : ecutwfc, nbnd, npwx, npw, igk, et, wg, nbndx, g2kin
  USE wavefunctions_module, ONLY : evc, psic
  use noncollin_module, only : npol
  USE io_files,      ONLY : iunwfc, nwordwfc, iunigk, prefix, diropn
  use io_global, only : stdout
  use control_flags, only: twfcollect
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  USE mp_wave,   ONLY : mergewf, splitwf
  USE mp,        ONLY : mp_max, mp_barrier
  use fft_base,  only : dffts, dfftp
  use fft_interfaces, only : fwfft, invfft
  USE mp_global, ONLY : mpime


  real(dp),parameter :: eps=1.d-12
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)

  integer,intent(in) :: ndim(3)
  integer,intent(in) :: band_subset(2)

  integer :: ik, isym, i, j
  integer :: imap
  integer :: ixyz

  integer :: nkstot_orig, npw_orig, npwx_orig
  real(dp),allocatable :: xk_orig(:,:)
  integer,allocatable :: igk_orig(:)
  integer,allocatable :: igk_l2g_orig(:,:)
!  integer,allocatable :: ig_l2g_orig(:)
  integer,allocatable :: ngk_orig(:)

  integer :: igwx, igwx_orig
  integer :: nwordwfcr, iunwfcr

  real(dp),allocatable :: r(:,:)
  complex(dp),allocatable :: expikr(:), wtmp(:)
  integer :: ibnd, jbnd, ikmapped
  complex(dp),allocatable :: evc_tmp(:,:), evc_tmpm(:,:)
  integer :: nbndm, nwordwfcm, iunwfcm

  real(dp) :: x(3)

  logical :: exst, opend
  integer,external :: freeunit


  write(stdout,*)
  write(stdout,*) ' map_shirley_bz '
  write(stdout,*) ' Map wave functions on input k-point set to all equivalent'
  write(stdout,*) ' points in the positive BZ [0,1]^3'
  write(stdout,*)

  ! expand kpoints using symmetry, inversion and translation
  call kpoint_expand( ndim )

! this is new
  ! I no longer have access to igk_l2g - need to regenerate
  if( allocated(igk_l2g) ) deallocate( igk_l2g )
  allocate( igk_l2g( npwx, nks ) )
  igk_l2g = 0
  if( nks > 1 ) rewind(iunigk)
  do ik=1,nks
    npw=ngk(ik)
    if( nks > 1 ) read(iunigk) igk
    call gk_l2gmap( ngm, ig_l2g(1), npw, igk(1), igk_l2g(1,ik) )
  enddo
! this is new

  write(stdout,*) '    nks = ', nks
  write(stdout,*) ' nkstot = ', nkstot, size(xk,2), size(ngk), size(igk_l2g,2)
  write(stdout,*) '   npwx = ', npwx, size(ig_l2g), size(igk_l2g,1)

  ! store original input k-point info
  nkstot_orig = nkstot
  npwx_orig = npwx
  allocate( xk_orig(3,nkstot_orig) )
  xk_orig = xk(1:3,1:nkstot_orig)
!  allocate( igk_orig(npwx_orig) )
  allocate( ngk_orig(nkstot_orig) )
  ngk_orig = ngk(1:nkstot_orig)
!  allocate( ig_l2g_orig(npwx_orig) )
!  ig_l2g_orig = ig_l2g
  allocate( igk_l2g_orig(npwx_orig,nkstot_orig) )
  igk_l2g_orig = igk_l2g
  igwx_orig = maxval( igk_l2g_orig )

  ! update k-points based on expansion from kpoint_expand
  nkstot = sum(xk_map(:)%nmap)
  nks = nkstot
  write(stdout,*) ' new number of mapped k-points = ', nkstot
  ikmapped=0
  do ik=1,nkstot_orig
    do imap=1,xk_map(ik)%nmap
      ikmapped=ikmapped+1
      xk(:,ikmapped) = matmul( bg, xk_map(ik)%xmap(:,imap) ) + xk_orig(:,ik)
      wk(ikmapped) = 0.d0 ! I'm zeroing this just in case
    enddo
  enddo

  !! update FFT grids based on the new set of k-points
  !write(stdout,*) ' update FFT grids '
  !call update_shirley_fft


  ! find the new npwx
  call n_plane_waves (ecutwfc, tpiba2, nkstot_orig, xk_orig, g, ngm, npwx, ngk)
  if( npwx > npwx_orig ) then
    ! these two allocations are important - gk_sort assumes that the dimension
    ! of igk and g2kin is npwx as stored in module wvfct
    if( allocated( igk ) ) deallocate( igk )
    allocate( igk(npwx) )
    if( allocated( g2kin ) ) deallocate( g2kin )
    allocate( g2kin(npwx) )
    if( allocated( igk_l2g ) ) deallocate( igk_l2g )
    allocate( igk_l2g(npwx,nkstot_orig) )
    if( allocated(evc) ) deallocate(evc)
    allocate( evc(npwx,nbnd) )
  endif

  igk_l2g = 0  ! important initialization
  do ik=1,nkstot_orig
    CALL gk_sort (xk_orig(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,ik))

    write(stdout,*) ' npwx = ', npwx, ' npw = ', npw, ' ngk(ik) = ', ngk(ik)
  enddo
  igwx = maxval( igk_l2g )

  ! new output file for reordered wave functions
  write(stdout,*) ' wfcr '
  call flush_unit(stdout)
  iunwfcr = freeunit()
  nwordwfcr = 2 * nbnd * npwx * npol
  call diropn( iunwfcr, 'wfcr', nwordwfcr, exst )

  ! reorder *original* k-points with respect to new FFT grids
  igwx = max( igwx, igwx_orig )
  call mp_max( igwx )
  write(stdout,*) ' igwx = ', igwx
  allocate( wtmp(igwx) )

  ! create pointer to evc subarray
  allocate( evc_tmp(1:npwx_orig,1:nbnd), &
            evc_tmpm(1:npwx,1:nbnd)      )
  write(stdout,*) ' sizes'
  write(stdout,*) ' evc      = ', size(evc,1), size(evc,2)
  write(stdout,*) ' evc_tmp  = ', size(evc_tmp,1), size(evc_tmp,2)
  write(stdout,*) ' evc_tmpm = ', size(evc_tmpm,1), size(evc_tmpm,2)

  do ik=1,nkstot_orig

    write(stdout,*) ' reordering original k-point ', ik

    ! note that nwordwfc needs a factor of 2 in v4.3
    write(stdout,*) size(evc_tmp)*2*npol, 2*nwordwfc

    ! load wave functions for ik
    call davcio( evc_tmp, 2*nwordwfc, iunwfc, ik, -1 )
    evc(1:npwx_orig,1:nbnd) = evc_tmp

    do ibnd=1,nbnd
      wtmp = 0.d0
      call mergewf(evc(:,ibnd), wtmp, ngk_orig(ik), igk_l2g_orig(:,ik), me_pool, nproc_pool, root_pool, intra_pool_comm)
      call splitwf(evc(:,ibnd), wtmp, ngk(ik), igk_l2g(:,ik), me_pool, nproc_pool, root_pool, intra_pool_comm)
    enddo

    ! save re-ordered wave functions for ik
    evc_tmpm = evc(1:npwx,1:nbnd)
    call davcio( evc_tmpm, nwordwfcr, iunwfcr, ik, +1 )

  enddo
  if( allocated(wtmp) ) deallocate( wtmp )
  if( allocated(evc_tmp) ) deallocate( evc_tmp )
  if( allocated(evc_tmpm) ) deallocate( evc_tmpm )

  ! adjust dimensions for new k-point set
  !if( npwx > npwx_orig ) then
  !  if( allocated( igk_orig ) ) deallocate( igk_orig )
  !  allocate( igk_orig(npwx) )
  !endif
  npwx_orig = npwx

  ! I need to reallocate ngk here to reflect the (potentially) larger number of k-points
  if( allocated( ngk ) ) deallocate( ngk )
  allocate( ngk(nkstot) )

  ! find the new npwx
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)

  write(stdout,*) ' npwx_orig = ', npwx_orig, ' npwx = ', npwx
  write(stdout,*) ' ngk = ', ngk

  if( npwx > npwx_orig ) then
    write(stdout,*) npwx
    write(stdout,*) ngk(1:nkstot)
    if( allocated(igk) ) deallocate(igk)
    allocate( igk(npwx) )
    if( allocated(g2kin) ) deallocate(g2kin)
    allocate( g2kin(npwx) )
    if( allocated(evc) ) deallocate(evc)
    allocate( evc(npwx,nbnd) )
    ! hope against hope that this is correct
  endif
  allocate( igk_orig(npwx) )

  ! we will redefine this below
  if( allocated(igk_l2g) ) deallocate( igk_l2g )
  allocate( igk_l2g(npwx,nkstot) )
  igk_l2g = 0  ! important initialization
 

  ! new output file for mapped wave functions
  write(stdout,*) ' wfcm '
  call flush_unit(stdout)
  iunwfcm = freeunit()
  nbndm = band_subset(2)-band_subset(1)+1
  nwordwfcm = 2 * nbndm * npwx * npol
  call diropn( iunwfcm, 'wfcm', nwordwfcm, exst )

  ! make real-space grid for r
  write(stdout,*) ' rgrid '
  call flush_unit(stdout)
  allocate( r(nrxx,3), expikr(nrxx) )
  write(stdout,*) ' nrxx = ', nrxx
  write(stdout,*) ' dfftp%nnr = ', dfftp%nnr
  call rgrid( r )
!  if( doublegrid ) then
!    do ixyz=1,3
!      call interpolate( r(1,ixyz), r(1,ixyz), -1 )
!    enddo
!  endif
  allocate( wtmp(nrxx) )

  ! begin mapping wave functions
  ikmapped = 0
  do ik=1,nkstot_orig
    
    ! determine igk
    ! why is this done here? are we wanting to keep _orig variables pristine
    ! or only their dimensions - because they will be overwritten here...
    CALL gk_sort (xk_orig(1,ik), ngm, g, ecutwfc/tpiba2, npw_orig, igk_orig, g2kin)
    ! Gah!
    if( npw_orig > npwx ) stop
 
    ! pointers to array sections of evc
    allocate( evc_tmp(1:npwx_orig,1:nbnd), &
              evc_tmpm(1:npwx,1:nbndm)     )

    write(stdout,*) ' mapping k-point ', ik

    ! make phases
    do imap=1,xk_map(ik)%nmap
      
      ikmapped=ikmapped+1

      ! determine igk and igk_l2g
      CALL gk_sort (xk(1,ikmapped), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,ikmapped))

      write(stdout,*) ikmapped, npw, maxval(igk_l2g(:,ikmapped))
 
      ! load wave functions for ik
      call davcio( evc_tmp, nwordwfcr, iunwfcr, ik, -1 )
      evc(1:npwx_orig,1:nbnd) = evc_tmp

      ! e^(-i 2 pi xk_map.r )
      expikr(:) = matmul( r, xk_map(ik)%xmap(:,imap) )
      expikr(:)= exp( ( - iota * tpi ) * expikr(:) )

      x = matmul( bg, xk_map(ik)%xmap(:,imap) )

      do ibnd=1,nbndm
        !
        jbnd = ibnd-1+band_subset(1)
        !
        CALL start_clock( 'firstfft' )
        !
        psic(1:nrxxs) = ( 0.D0, 0.D0 )
        psic(nls(igk_orig(1:npw_orig))) = evc(1:npw_orig,jbnd)
        !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
        CALL invfft ('Wave', psic, dffts)
        !
        CALL stop_clock( 'firstfft' )
        !
        ! ... product with the phase factor expikr =  exp(i xmap.r) 
        !
        !!!!psic(1:nrxxs) = psic(1:nrxxs) * expikr(1:nrxxs)
        !
        ! ... make the realspace product at higher resolution
        call cinterpolate( wtmp, psic, +1 )
        wtmp(1:nrxx) = wtmp(1:nrxx) * expikr(1:nrxx)
        call cinterpolate( wtmp, psic, -1 )
        !
        ! ... back to reciprocal space
        !
        CALL start_clock( 'secondfft' )
        !
        !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
        call fwfft ('Wave', psic, dffts)
        !
        ! ... store with correct ordering
        !
        !evc_map(1:npw,ibnd) = psic(nls(igk(1:npw)))
        evc(1:npw,jbnd) = psic(nls(igk(1:npw)))
        if( npw < npwx ) evc(npw+1:npwx,jbnd) = 0.d0
        !
        CALL stop_clock( 'secondfft' )
        !
      enddo

      ! write to file
      write(stdout,*) '     dumping mapped k-point ', ikmapped
      evc_tmpm = evc(1:npwx,band_subset(1):band_subset(2))
      call davcio( evc_tmpm, nwordwfcm, iunwfcm, ikmapped, +1 )

    enddo
    
    if( allocated(evc_tmp) ) deallocate( evc_tmp )
    if( allocated(evc_tmpm) ) deallocate( evc_tmpm )

  enddo

  if( allocated(wtmp) ) deallocate( wtmp )
  if( allocated(r) ) deallocate( r )
  if( allocated(expikr) ) deallocate( expikr )
  if( allocated(xk_orig) ) deallocate( xk_orig )
  deallocate( igk_orig )

  ! close first reordered wave functions
  !call delete_shirley_wfc()
  close( iunwfcr, status='delete' )

  ! modify the wave function variables for future reference
  inquire(unit=iunwfc,opened=opend)
  if( opend ) then
    if( twfcollect ) then
      close( iunwfc, status='delete' )
    else
      close( iunwfc )
    endif
  endif
  ! switch wave function unit to wfcm's
  iunwfc = iunwfcm
  ! new convention for nwordwfc in v4.3
  nwordwfc = nwordwfcm/2
  nbnd = nbndm
  nbndx = nbnd
  if( allocated(evc) ) deallocate( evc )
  allocate( evc(npwx,nbnd) )


  end subroutine map_shirley_bz


  end module wfc_shirley
