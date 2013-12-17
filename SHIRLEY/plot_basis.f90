! ---------------------------------------------------------------------- 
  subroutine plot_basis( )
! ---------------------------------------------------------------------- 

  use constants, only : tpi
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base
  USE ions_base, ONLY : ntyp=>nsp, nat, ityp, tau, atm
  USE gvect  
  USE grid_dimensions,  ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  use control_flags, only : gamma_only
  use gvecs, only : nls
  USE smooth_grid_dimensions,  ONLY : nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nrxxs
  USE io_files, ONLY: prefix, tmp_dir, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc, psic
  USE lsda_mod, ONLY : current_spin, lsda, isk, nspin

  use shirley_ham_input, only : debug, band_subset, &
                                nplot, plot_prefix, plot_style, &
                                plot_center_atom
  use plot_module, only : plot_basis_write_header, &
                          plot_basis_save_header
  use mp, only : mp_bcast, mp_barrier
  use fft_base, only : dffts
  use fft_interfaces, only : invfft

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)
  integer,parameter :: maxchar=255

  integer :: ibnd, jbnd, ik
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)

  integer :: i, j, ierr
  integer :: iunplt
  integer :: nwordpwm
  logical :: exst

  !integer :: nbasis
  real(dp),allocatable :: allplot(:), psisq(:), allplot_sort(:)
  integer,allocatable :: indplot_sort(:)
  complex(dp),allocatable :: plotevc(:,:)
  real(dp),allocatable :: ploteigval(:), plotde(:)
  complex(dp),allocatable :: ploteigvec(:,:)
  integer,allocatable :: iunout(:)
  character(maxchar) :: fmtstr, isofile
  real(dp),allocatable :: rvec(:,:)
  real(dp) :: integral, tot_integral, fac

  real(dp) :: tauc(3,nat)
  real(dp),allocatable :: grid(:,:,:)
  real(dp) :: shift(3)

  integer,external :: freeunit


  WRITE( stdout, '(/5x,"Calling plot_basis .... ",/)')
  write(stdout,*) ' npwx  = ', npwx
  write(stdout,*) ' npw   = ', npw
  write(stdout,*) ' nbnd  = ', nbnd
  write(stdout,*) ' nbndx = ', nbndx
  write(stdout,*) '  nr1x = ', nr1x
  write(stdout,*) '  nr2x = ', nr2x
  write(stdout,*) '  nr3x = ', nr3x

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

  if( ionode ) then
    write(stdout,*) ' plotting the optimal basis on a real-space grid'
    write(stdout,*) '   < r | B_i > '
  endif

  call flush_unit( stdout )

  ! open plt file
  if( ionode ) then
    iunplt = freeunit()
    open(iunplt,file=trim(plot_prefix)//'.plt',form='formatted',iostat=ierr)
    write(stdout,*) ' reading from file: '
    write(stdout,*) trim(plot_prefix)//'.plt'
    call errore( 'plot_basis', &
      'unable to open file '//trim(plot_prefix)//'.plt',abs(ierr) )

    !read(iunplt,*) nplot, nbasis, plot_style
    !if( nbasis /= nbnd ) call errore('plot_basis', &
    !  'number of basis functions does not match',1)

    write(stdout,*) ' style of plot is ', trim(plot_style)
    if( trim(plot_style) /= 'cube' .and. trim(plot_style) /= 'xsf' ) then
      call errore('plot_basis','unrecognized plot_style',1)
    endif

    allocate( plotevc(npwx,nplot), ploteigval(nplot), plotde(nplot), &
              ploteigvec(nbnd,nplot), iunout(nplot) )

    allocate( rvec(1:2,nbnd) )
    do i=1,nplot
      read(iunplt,*) ibnd, ik
      read(iunplt,*) ploteigval(i), plotde(i)

      write(stdout,*) ' eigen-energy band= ', ibnd, ' k= ', ik, &
                      ' e= ', ploteigval(i), ' de= ', plotde(i)

      write(fmtstr,'(a,i6,a,i6,a)') '(a,i', int(log10(dble(ik)))+1, &
                         ',a,i', int(log10(dble(ibnd)))+1, ',a)'

      if( trim(plot_style) == 'cube' ) then
        write(isofile,trim(fmtstr)) trim(plot_prefix)//'.', ik, '.', ibnd, '.cube'
      else if( trim(plot_style) == 'xsf' ) then
        write(isofile,trim(fmtstr)) trim(plot_prefix)//'.', ik, '.', ibnd, '.xsf'
      else
        call errore('plot_basis','unrecognized plot_style',1)
      endif

      iunout(i)=freeunit()
      open(iunout(i),file=trim(isofile),form='formatted',iostat=ierr)

      read(iunplt,*) rvec
      ploteigvec(:,i) = cmplx( rvec(1,:), rvec(2,:), kind=dp )
    enddo
    deallocate( rvec )
    close(iunplt)
  endif

  call mp_bcast( nplot, ionode_id )
  if( .not. ionode ) then
    allocate( plotevc(npwx,nplot), ploteigval(nplot), &
              ploteigvec(nbnd,nplot) )
  endif
  call mp_bcast( ploteigvec, ionode_id )
  call mp_bcast( ploteigval, ionode_id )


  ! allocate space
  allocate( psisq( size(psic) ), &
            allplot( nr1x*nr2x*nr3x ), &
            allplot_sort( nr1x*nr2x*nr3x ), &
            indplot_sort( nr1x*nr2x*nr3x ), &
            grid(nr1x,nr2x,nr3x), stat=ierr )
  call errore('plot_basis','unable to allocate space for output',abs(ierr))

  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
!  ! gamma-point only
!  qvec = 0.d0
!  current_k = 1
!  if( lsda ) current_spin = isk(1)

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
!  call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
!  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! ======================================================================
  ! plane-wave matrix elements
  ! ======================================================================
  write(stdout,*) ' plotting'


!  call plot_basis_save_header( alat, at, nat, ntyp, &
!               tau, atm, ityp, nr1, nr2, nr3, nbnd )
!
!  if( ionode ) then
!    call  plot_basis_write_header( iunplt )
!  endif

  ! expand the eigenvector
  write(stdout,*) ' matmul ...'
  plotevc = matmul( evc(:,ibnd_indx(1:nbnd)), ploteigvec(1:nbnd,:) )
!  call ZGEMM( 'N', 'N', npw, nplot, nbnd, one, &
!              evc(:,ibnd_indx(1:nbnd)), npwx, &
!              ploteigvec, nbnd, &
!              zero, plotevc, npwx )
  write(stdout,*) ' ... done'

  do i=1,nplot
    !
    call mp_barrier
    !
    psic(:) = zero
    !
    psic(nls(igk(1:npw))) = plotevc(1:npw,i)
    !
    !CALL cft3s( psic, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, 2 )
    CALL invfft ('Wave', psic, dffts)
    !
    psisq = dble(conjg(psic)*psic)
    call gather( psisq, allplot )
    !
    if( ionode ) then

      ! isosurfaces
      tot_integral = sum(allplot)
      allplot_sort = allplot
      indplot_sort = 0
      call hpsort_eps( nr1x*nr2x*nr3x, allplot_sort, indplot_sort, 1.d-10 )
      fac = 10.d0
      integral = 0.d0
      write(stdout,'(a,i6)') '# Isosurface dump ', i
      do j=1,nr1x*nr2x*nr3x
        integral = integral + allplot_sort(j)
        if( integral/tot_integral*100.d0 > fac ) then
          write(stdout,'(a,f7.2,a,f12.5)') '#', 100.d0-fac, &
                                          '% isosurface at ', &
                                          allplot_sort(j)
          fac = fac + 10.d0
          if( fac > 90.d0 ) exit
        endif
      enddo

      ! mapping of atoms and grid such that it is centered
      ! convert to crystal coords
      tauc = tau
      call cryst_to_cart( nat, tauc,  bg, -1 )
      ! center the system
      if( plot_center_atom < 1 .or. plot_center_atom > nat) then
        call all_center( tauc, shift )
      else
        call atom_center( tauc, plot_center_atom, shift )
      endif
      call grid_shift( tauc, shift )
      ! convert back to Cartesian coords / alat
      call cryst_to_cart( nat, tauc,  at, +1 )


      if( trim(plot_style) == 'cube' ) then
        call write_cubefile ( alat, at, bg, nat, tauc, atm, ityp, allplot, &
                              nr1, nr2, nr3, nr1x, nr2x, nr3x, iunout(i) )
      else if( trim(plot_style) == 'xsf' ) then
        call xsf_struct (alat, at, nat, tauc, atm, ityp, iunout(i))
        call xsf_fast_datagrid_3d &
          (allplot, nr1, nr2, nr3, nr1x, nr2x, nr3x, at, alat, iunout(i))
      else
        call errore('plot_basis','unrecognized plot_style',1)
      endif

      close(iunout(i))
    endif
    !
  enddo

  deallocate( psisq, psic )

  return

  contains


  SUBROUTINE gather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers dfftp%nproc distributed data on the first processor of that subset
  !
  ! ... real(dp)  f_in  = distributed variable (nxx)
  ! ... real(dp)  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  USE fft_base,  ONLY : dfftp
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  real (DP) :: f_in(:), f_out(:)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfftp%root,    &
                    dfftp%comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
#endif
  !
  RETURN
  !
  END SUBROUTINE gather


  SUBROUTINE cgather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers dfftp%nproc distributed data on the first processor of that subset
  !
  ! ... complex(dp)  f_in  = distributed variable (nxx)
  ! ... complex(dp)  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  USE fft_base,  ONLY : dfftp
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  complex (DP) :: f_in(:), f_out(:)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfftp%root,    &
                    dfftp%comm, info )
  !
  CALL errore( 'cgather', 'info<>0', info )
  !
#endif
  !
  RETURN
  !
  END SUBROUTINE cgather


  subroutine all_center( tauc, cvec )

  real(dp) :: tauc(3,nat)
  real(dp) :: avec(3), cvec(3)
  integer :: n

  write(stdout,*) ' cell_center '

  avec=0.d0
  do n=1,nat
    avec = avec + tauc(:,n)
  enddo
  avec = avec / dble(nat)
  cvec = 0.5d0 - avec

  end subroutine all_center


  subroutine atom_center( tauc, iatom, cvec )

  real(dp) :: tauc(3,nat)
  real(dp) :: cvec(3)
  integer :: iatom

  cvec = 0.5d0 - tauc(:,iatom)

  end subroutine atom_center


  subroutine grid_shift( tauc, shift )

  real(dp) :: tauc(3,nat)
  real(dp) :: shift(3), cvec(3)
  integer :: n, nr1c, nr2c, nr3c, i1, i2, i3, j1, j2, j3, j123

  write(stdout,*) ' grid_shift '

  cvec = shift

  nr1c = cvec(1) * dble(nr1)
  nr2c = cvec(2) * dble(nr2)
  nr3c = cvec(3) * dble(nr3)

  ! make sure these shifts are positive in the cell
  nr1c = nr1c - floor( dble(nr1c)/dble(nr1) )*nr1
  nr2c = nr2c - floor( dble(nr2c)/dble(nr2) )*nr2
  nr3c = nr3c - floor( dble(nr3c)/dble(nr3) )*nr3

  cvec(1) = dble(nr1c)/dble(nr1)
  cvec(2) = dble(nr2c)/dble(nr2)
  cvec(3) = dble(nr3c)/dble(nr3)

  do n=1,nat
    tauc(:,n) = tauc(:,n) + cvec
  enddo
  
  grid = reshape( allplot, (/ nr1x, nr2x, nr3x /) )
  allplot = 0.d0

  do i3=1,nr3
    j3=mod(i3+nr3c-1,nr3)+1
    do i2=1,nr2
      j2=mod(i2+nr2c-1,nr2)+1
      do i1=1,nr1
        j1=mod(i1+nr1c-1,nr1)+1
        j123=j1+(j2-1)*nr1+(j3-1)*nr1*nr2

        allplot(j123) = grid(i1,i2,i3)
      enddo
    enddo
  enddo

  end subroutine grid_shift

  end subroutine plot_basis
! ---------------------------------------------------------------------- 

