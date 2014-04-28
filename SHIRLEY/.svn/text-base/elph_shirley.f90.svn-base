  module elph_shirley

#include "f_defs.h"

  USE io_global,  ONLY : stdout, ionode, ionode_id

  use kinds, only : dp
  implicit none

  integer :: nbnddv, nspindv, npedv
  real(dp),allocatable :: w2(:)
  complex(dp),allocatable :: dyn(:,:)
  complex(dp),allocatable :: elph(:,:,:)

  contains

! ---------------------------------------------------------------------- 
  subroutine elphmat( )
! ---------------------------------------------------------------------- 

  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base
  USE gvect  
  use grid_dimensions, only : nr1, nr2, nr3
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  USE control_flags, only : gamma_only
  use gvecs, only : nls, doublegrid
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, iunigk
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin
  use mp,        only : mp_sum, mp_bcast
  use mp_global, only : intra_pool_comm
  use fft_base, only : dffts
  use fft_interfaces, only : fwfft, invfft

  use hamq_shirley
  use bspline, only : dbsnak, dbs3in, dbs3gd
  use shirley_ham_input, only : debug, band_subset, &
                                nrdv1, nrdv2, nrdv3, &
                                nspindv_in=>nspindv, npedv_in=>npedv, &
                                dvfile, elphfile, dynfile

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)

  integer :: ios, iundv, nrecdv, unf_nrecdv
  integer :: ipedv
  integer :: iundyn

  ! spline data
  integer :: kxord, kyord, kzord
  integer :: nxcoef, nycoef, nzcoef
  real(dp),allocatable :: rx(:), ry(:), rz(:)
  real(dp),allocatable :: rxp(:), ryp(:), rzp(:)
  real(dp),allocatable :: xknot(:), yknot(:), zknot(:)

  complex(dp),allocatable :: dv(:,:,:,:)
  real(dp),allocatable :: dvgrid(:,:,:)

  real(dp),allocatable :: dvsplR(:), dvsplI(:)
  real(dp),allocatable :: dvintR(:), dvintI(:)
  complex(dp),allocatable :: dvint(:)

  ! band data
  integer :: ibnd, jbnd
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)

  ! tmp data
  integer :: ij, i, j, k, n, ispin, i1,i2,i3
  complex(dp),allocatable :: jtmp(:)
  integer :: iunepm
  integer,external :: freeunit

#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8 
#endif

  nspindv = nspindv_in
  npedv = npedv_in

  WRITE( stdout, '(/5x,"Calling elphmat .... ",/)')
  write(stdout,*) ' npwx  = ', npwx
  write(stdout,*) ' npw   = ', npw
  write(stdout,*) ' nbnd  = ', nbnd
  write(stdout,*) ' nbndx = ', nbndx

  ! sort out which band subset we will work with
  if( band_subset(1) > band_subset(2) ) then
    i=band_subset(2)
    band_subset(2) = band_subset(1)
    band_subset(1) = i
  endif
  if( band_subset(1)>=nbnd .or. band_subset(1)<=0 ) band_subset(1) = 1
  if( band_subset(2)>=nbnd .or. band_subset(2)<=0 ) band_subset(2) = nbnd

  if( band_subset(2)-band_subset(1)+1 < nbnd ) then
    write(stdout,*) ' Requested band subset:', band_subset(1), &
                    ' ... ', band_subset(2)
    write(stdout,*) ' Reducing total number of bands from ', nbnd, &
                    ' to ', band_subset(2)-band_subset(1)+1
  endif

  nbnd = band_subset(2)-band_subset(1)+1
  allocate( ibnd_indx(nbnd) )
  do i=1,nbnd
    ibnd_indx(i) = band_subset(1)+i-1
  enddo

  nbnddv=nbnd

  ! check dvfile dimensions
  if( nrdv1==0 ) nrdv1 = nr1
  if( nrdv2==0 ) nrdv2 = nr2
  if( nrdv3==0 ) nrdv3 = nr3
  if( nspindv==0 ) nspindv=nspin
  if( nspindv/=nspin ) then
    call errore('elphmat','number of perturbation spin channels, nspindv, is different',abs(nspindv))
  endif
  if( npedv==0 ) then
    call errore('elphmat','number of perturbations, npedv, is zero',1)
  endif
  ! report
  write(stdout,*) ' dimensions of input dvfile: ', trim(dvfile)
  write(stdout,*) ' nrdv1 = ', nrdv1, ' nr1 = ', nr1
  write(stdout,*) ' nrdv2 = ', nrdv2, ' nr2 = ', nr2
  write(stdout,*) ' nrdv3 = ', nrdv3, ' nr3 = ', nr3
  write(stdout,*) ' number of perturbations: ', npedv

  if( ionode ) then
    write(stdout,*) ' generating electron-phonon matrix elements in this basis'
    write(stdout,*) '   < B_i | dV_SCF | B_j > '

    ! try to open input file containing dvscf
    nrecdv = 2 * nrdv1 * nrdv2 * nrdv3 * nspindv
    unf_nrecdv = DIRECT_IO_FACTOR * nrecdv
    iundv = freeunit()
    open( unit=iundv, file=trim(dvfile), iostat=ios, form='unformatted', &
          status='unknown', access='direct', recl=unf_nrecdv )
    if (ios /= 0) call errore ('elphmat', 'error opening '//trim(dvfile), iundv)

    ! output file
    iunepm = freeunit()
    open(iunepm,file=trim(elphfile),form='unformatted')
    write(stdout,*) ' will be saved to file: ', trim(elphfile)
    !
  endif

  
  if( ionode ) then
    ! try to open the dynamical matrix
    iundyn = freeunit()
    open( unit=iundyn, file=trim(dynfile), iostat=ios, form='formatted' )
    if (ios /= 0) call errore ('elphmat', 'error opening '//trim(dynfile), iundv)
    write(stdout,*) ' reading dynamical matrix', trim(dynfile)

    ! diagonalize dynamical matrix and get the eigenvalues
    call readmat( iundyn )
    !
    write(stdout,*) ' done reading dynamical matrix', trim(dynfile)
  endif
  call mp_bcast( npedv, ionode_id )
  if( .not. ionode ) then
    allocate( w2(npedv), dyn(npedv,npedv) )
  endif
  call mp_bcast( w2, ionode_id )
  call mp_bcast( dyn, ionode_id )
  write(stdout,*) ' dynamical matrix broadcasted', trim(dynfile)


  call flush_unit( stdout )

  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)

  write(stdout,*) ' ecutwfc = ', ecutwfc
  write(stdout,*) ' tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nktot = ', nks, nkstot
  write(stdout,*) ' xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '     npw = ', ngk(1:nks)
  write(stdout,*) '    npwx = ', npwx


!  ! gamma-point only
!  qvec = 0.d0
!  current_k = 1
!  if( lsda ) current_spin = isk(1)

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )

! davegp - looks like igk_l2g is not needed so skip this line
!  call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
!
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! report norms
  allocate( norm(nbnd) )
  do ibnd=1,nbnd
    norm(ibnd) = dot_product( evc(:,ibnd_indx(ibnd)), evc(:,ibnd_indx(ibnd)) )
  enddo
  call mp_sum( norm, intra_pool_comm )
  do ibnd=1,nbnd
    if( abs(norm(ibnd)-1.d0) > eps ) then
      write(stdout,'(a,i6,a,f14.10)') ' band ', ibnd, ' norm = ', norm(ibnd)
      call errore('hamq','wave function norm is not 1',1)
    endif
  enddo
  deallocate( norm )
 
  ! tmp space for wave functions
  allocate( jtmp(npw) )


  ! ======================================================================
  ! spline stuff
  ! ======================================================================
  allocate( dv(nrdv1,nrdv2,nrdv3,nspin), dvint(nr1*nr2*nr3), &
            elph(nbnd,nbnd,nspin), &
            stat=ios )
  if( ios/= 0 ) call errore('elphmat','unable to allocate space for dv',abs(ios))

  ! spline orders - fixed here at 3
  kxord=3
  kyord=3
  kzord=3

  nxcoef=nrdv1+1
  nycoef=nrdv2+1
  nzcoef=nrdv3+1

  ! dv grid coordinates
  allocate( rx(nxcoef), xknot(nxcoef+kxord) )
  allocate( ry(nycoef), yknot(nycoef+kyord) )
  allocate( rz(nzcoef), zknot(nzcoef+kzord) )

  forall( i=1:nxcoef ) rx(i)=dble(i-1)/dble(nxcoef-1)
  forall( i=1:nycoef ) ry(i)=dble(i-1)/dble(nycoef-1)
  forall( i=1:nzcoef ) rz(i)=dble(i-1)/dble(nzcoef-1)

  ! make knots
  call dbsnak( nxcoef, rx, kxord, xknot )
  call dbsnak( nycoef, ry, kyord, yknot )
  call dbsnak( nzcoef, rz, kzord, zknot )

  write(stdout,*) ' made knots on dv grid '

  ! wave function grid coordinates
  allocate( rxp(nr1) )
  allocate( ryp(nr2) )
  allocate( rzp(nr3) )

  forall( i=1:nr1 ) rxp(i)=dble(i-1)/dble(nr1)
  forall( i=1:nr2 ) ryp(i)=dble(i-1)/dble(nr2)
  forall( i=1:nr3 ) rzp(i)=dble(i-1)/dble(nr3)

  ! space for interpolation
  allocate( dvgrid(nxcoef,nycoef,nzcoef), &
            dvsplR(nxcoef*nycoef*nzcoef), &
            dvsplI(nxcoef*nycoef*nzcoef), &
            dvintR(nr1*nr2*nr3), dvintI(nr1*nr2*nr3) )

  write(stdout,*) ' data grid allocated'


  ! ======================================================================
  ! electron-phonon matrix elements
  ! ======================================================================
  write(stdout,*) ' electron-phonon matrix elements'
  pert_loop: do ipedv=1,npedv

    ! read this dv
    if( ionode ) then
    read(iundv,rec=ipedv,iostat=ios) dv
    if( ios/=0 ) call errore('elphmat','problem reading record from '//trim(dvfile),ipedv)
    endif
    call mp_bcast( dv, ionode_id )

    write(stdout,*) ' perturbation : ', ipedv, ' sq. freq. = ', w2(ipedv)

    spin_loop: do ispin=1,nspin

      ! put on grid separated in real ...
      forall( i1=1:nxcoef, i2=1:nycoef, i3=1:nzcoef ) &
        dvgrid(i1,i2,i3) = &
          real ( dv( mod(i1-1,nxcoef)+1, &
                     mod(i2-1,nxcoef)+1, &
                     mod(i3-1,nzcoef)+1, ispin ) &
               )

      ! B-spline interpolate
      call dbs3in( nxcoef, rx, nycoef, ry, nzcoef, rz, &
                   dvgrid, nxcoef, nycoef, &
                   kxord, kyord, kzord, &
                   xknot, yknot, zknot, dvsplR )
    
      !                           ... and imaginary
      forall( i1=1:nxcoef, i2=1:nycoef, i3=1:nzcoef ) &
        dvgrid(i1,i2,i3) = &
          aimag( dv( mod(i1-1,nxcoef)+1, &
                     mod(i2-1,nxcoef)+1, &
                     mod(i3-1,nzcoef)+1, ispin ) &
               )

      ! B-spline interpolate
      call dbs3in( nxcoef, rx, nycoef, ry, nzcoef, rz, &
                   dvgrid, nxcoef, nycoef, &
                   kxord, kyord, kzord, &
                   xknot, yknot, zknot, dvsplI )

      ! interpolate values on the current grid
      ! zeros correspond to derivative values
      call dbs3gd( 0, 0, 0, nr1, rxp, nr2, ryp, nr3, rzp,          &
                   kxord, kyord, kzord, xknot, yknot, zknot,       &
                   nxcoef, nycoef, nzcoef, dvsplR, dvintR, nr1, nr2)
      call dbs3gd( 0, 0, 0, nr1, rxp, nr2, ryp, nr3, rzp,          &
                   kxord, kyord, kzord, xknot, yknot, zknot,       &
                   nxcoef, nycoef, nzcoef, dvsplI, dvintI, nr1, nr2)
      forall(i=1:nr1*nr2*nr3) dvint(i) = cmplx( dvintR(i), dvintI(i) )

      if( doublegrid ) call cinterpolate( dvint, dvint, -1 )

      ! now generate matrix elements
      do ibnd=1,nbnd

        ! transform to real space and complex conjugate
        psic(1:nrxxs) = ( 0.D0, 0.D0 )
        psic(nls(igk(1:npw))) = conjg( evc(1:npw,ibnd_indx(ibnd)) )
        !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
        CALL invfft ('Wave', psic, dffts)

        ! multiply by the dv
        psic(1:nrxxs) = psic(1:nrxxs) * dvint(1:nrxxs)

        ! convert back to Fourier space
        !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, -2 )
        CALL fwfft ('Wave', psic, dffts)
        jtmp(1:npw) = psic(nls(igk(1:npw)))

        do jbnd=1,nbnd

          ! dot-product with band j
          elph(ibnd,jbnd,ispin) = &
            sum( jtmp(1:npw) * evc(1:npw,ibnd_indx(jbnd)) )

        enddo ! jbnd

      enddo ! ibnd
         
    enddo spin_loop

    ! sum over all processes
    call mp_sum( elph )

    if( ionode ) call write_elph( iunepm, ipedv )

  enddo pert_loop

  if( ionode ) close(iunepm)

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( dvgrid, dvsplR, dvsplI, dvintR, dvintI )
  deallocate( rx, ry, rz )
  deallocate( rxp, ryp, rzp )
  deallocate( dv, dvint, elph )
  deallocate( ibnd_indx, jtmp )
  deallocate( w2, dyn )

  return

  contains

!-----------------------------------------------------------------------
SUBROUTINE readmat ( iudyn )
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  ! Input
  INTEGER :: iudyn 
  ! local
  integer :: ibrav, nat, ntyp
  REAL(DP) :: celldm (6), q (3)
  integer,allocatable :: ityp(:)
  real(dp),allocatable :: amass(:), tau(:,:)
  REAL(DP),allocatable :: dynr (:, :, :, :, :)
  CHARACTER(len=80) :: line
  CHARACTER(len=3)  :: atm
  INTEGER :: nt, na, nb, naa, nbb, nu, mu, i, j
  !
  ! tmp space
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer :: nblock, lwork, ierr, ii, jj
  integer,external :: ilaenv
  !
  !
  REWIND (iudyn)
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, * ) ntyp, nat, ibrav, celldm

  ! allocate space
  allocate( w2(3*nat), dyn(3*nat,3*nat) )
  allocate( dynr (2, 3, nat, 3, nat), amass(ntyp), tau(3,nat), ityp(nat) )

  DO nt = 1, ntyp
     READ (iudyn, * ) i, atm, amass
  ENDDO
  DO na = 1, nat
     READ (iudyn, * ) i, ityp(na), tau(:,na)
  ENDDO
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (line (11:80), * ) (q (i), i = 1, 3)
  READ (iudyn, '(a)') line
  DO na = 1, nat
     DO nb = 1, nat
        READ (iudyn, * ) naa, nbb
        IF (na.NE.naa.OR.nb.NE.nbb) CALL errore ('readmat', 'error reading &
             &file', nb)
        READ (iudyn, * ) ( (dynr (1, i, na, j, nb), dynr (2, i, na, j, nb) &
             , j = 1, 3), i = 1, 3)
     ENDDO
  ENDDO
  !
  ! divide the dynamical matrix by the masses
  !
  DO nb = 1, nat
     DO j = 1, 3
        DO na = 1, nat
           DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! solve the eigenvalue problem.
  ! NOTA BENE: eigenvectors are overwritten on dyn
  !
  ! tmp space for diagonalization
  nblock = ILAENV( 1, 'ZHETRD', 'U', 3*nat, - 1, - 1, - 1 )
  IF ( nblock < 1 ) nblock = MAX( 1, 3*nat )
  IF ( nblock == 1 .OR. nblock >= 3*nat ) THEN
     lwork = 2 * 3*nat - 1
  ELSE
     lwork = ( nblock + 1 ) * 3*nat
  END IF
  allocate( work(lwork), rwork(3*3*nat-2) )
  
  jj=0
  DO nb = 1, nat
     DO j = 1, 3
        jj=jj+1
        ii=0
        DO na = 1, nat
           DO i = 1, 3
              ii=ii+1
              dyn (ii,jj) = cmplx( dynr (1, i, na, j, nb),  &
                                   dynr (2, i, na, j, nb) )
           enddo
        enddo
     enddo
  enddo

  CALL ZHEEV( 'V', 'U', 3*nat, dyn, 3*nat, w2, work, lwork, rwork, ierr )
  if( ierr /= 0 ) call errore('readmat','error diagonalizing dynamical matrix',abs(ierr))

  deallocate( work, rwork )

!!  CALL cdiagh (3 * nat, dynr, 3 * nat, w2, dyn)
  !
  ! divide by sqrt(mass) to get displacements
  !
  DO nu = 1, 3 * nat
     DO mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu) = dyn (mu, nu) / SQRT (amass (ityp (na) ) )
     ENDDO
  ENDDO
  !
  !
  RETURN
END SUBROUTINE readmat

  end subroutine elphmat
! ---------------------------------------------------------------------- 


! ---------------------------------------------------------------------- 
  subroutine write_elph( iunepm, ipedv )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunepm, ipedv

  ! save to file
  if( ipedv==1 ) then
    rewind(iunepm)
    ! some dimensions
    write(iunepm) nbnddv, nspindv, npedv
    write(stdout,*) ' writing the first elph record'
    write(stdout,*) ' dimensions: ', nbnddv, nspindv, npedv
  endif
  ! now everything else
  write(iunepm) w2(ipedv)
  write(iunepm) dyn(1:npedv,ipedv)
  write(iunepm) elph

  end subroutine write_elph


! ---------------------------------------------------------------------- 
  subroutine read_elph( iunepm, ipedv )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunepm, ipedv

  ! read from file
  if( ipedv==1 ) then
    rewind(iunepm)
    ! some dimensions
    read(iunepm) nbnddv, nspindv, npedv
    write(stdout,*) 'elph dimensions: ', nbnddv, nspindv, npedv
    ! allocate
    if( allocated(w2) ) deallocate(w2)
    if( allocated(dyn) ) deallocate(dyn)
    if( allocated(elph) ) deallocate(elph)
    allocate( w2(npedv), dyn(npedv,npedv), elph(nbnddv,nbnddv,nspindv) )
  endif
  ! now everything else
  read(iunepm) w2(ipedv)
  read(iunepm) dyn(1:npedv,ipedv)
  read(iunepm) elph

  end subroutine read_elph
  

  end module elph_shirley
