! ---------------------------------------------------------------------- 
  module pot_rho_G
! ---------------------------------------------------------------------- 

! This module is useful for keeping a Fourier space copy of the potential
! and charge density along with the original mapping from local to global
! G-vectors, just in case you should wish to change the FFT grid and thus
! introduce a new mapping ig_l2g

  use kinds, only : dp
  use mp,        only : mp_max, mp_sum
  use mp_global, only : mpime
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft

  complex(dp),allocatable :: rhoG(:,:)
  complex(dp),allocatable :: vrG(:,:)
  real(dp) :: ref_charge
  integer,allocatable :: ig_l2g_ref(:)

  real(dp),allocatable :: ns(:,:,:,:)

  contains


! ---------------------------------------------------------------------- 
  subroutine read_pot_rho_G ( )
! ---------------------------------------------------------------------- 

  use lsda_mod,  only : nspin
  use scf,       only : rho, v
  use gvect,     only : ngm, nl, ig_l2g
  use cell_base, only : omega
  USE mp_global, ONLY : intra_pool_comm
  use io_rho_xml, only : read_rho
  use io_global, only : stdout
  use ldaU,      only : lda_plus_u


  implicit none

  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),allocatable :: aux(:)
  integer :: ispin, ig

  ! copy map from local to global G-vectors
  allocate( ig_l2g_ref(size(ig_l2g)) )
  ig_l2g_ref = ig_l2g

  ! read potential and charge from existing files

  call read_rho( rho, nspin )

  ! store total charge for future reference

  IF ( nspin == 2 ) THEN
     ref_charge = SUM ( rho%of_r (:, 1:nspin) ) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  ELSE
     ref_charge = SUM ( rho%of_r (:, 1) ) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  END IF
#ifdef __PARA
  call mp_sum( ref_charge, intra_pool_comm )
#endif

  WRITE( stdout, &
          '(/,5X,"starting charge ",F10.5)') &
          ref_charge

  if( lda_plus_u ) then
    write(stdout,*) 'LDA+U - reading occupations'
    write(stdout,*) rho%ns

    ! copy occupancy
    if( allocated(ns) ) deallocate(ns)
    allocate( ns(size(rho%ns,1),size(rho%ns,2),size(rho%ns,3),size(rho%ns,4)) )
    ns = rho%ns
  endif

  ! transform charge and density to G-space

  allocate( rhoG(ngm,2), vrG(ngm,2) )
  allocate( aux(dfftp%nnr) )

  ! density
  rhoG = zero
  do ispin=1,nspin
    aux(:) = cmplx(rho%of_r(:,ispin))
    !CALL cft3( aux, nr1, nr2, nr3, nr1x, nr2x, nr3x, -1 )
    CALL fwfft ('Dense', aux, dfftp)
    do ig=1,ngm
      rhoG(ig,ispin) = aux(nl(ig))
    enddo
  enddo

  ! potential
  vrG = zero
  do ispin=1,nspin
    aux(:) = cmplx(v%of_r(:,ispin))
    !CALL cft3( aux, nr1, nr2, nr3, nr1x, nr2x, nr3x, -1 )
    CALL fwfft ('Dense', aux, dfftp)
    do ig=1,ngm
      vrG(ig,ispin) = aux(nl(ig))
    enddo
  enddo

  deallocate( aux )

  return
  end subroutine read_pot_rho_G


! ---------------------------------------------------------------------- 
  subroutine write_pot_rho_G ( )
! ---------------------------------------------------------------------- 

  use lsda_mod,  only : nspin
  use ldaU,      only : lda_plus_u
  use scf,       only : rho, v, create_scf_type, destroy_scf_type
  use gvect,     only : ngm, nl, ig_l2g
  use cell_base, only : omega
  USE mp_global, ONLY : intra_pool_comm
  use io_rho_xml, only : write_rho
  use io_global, only : stdout
  use ldaU,      only : lda_plus_u

  implicit none

  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),allocatable :: aux(:), aux1(:)
  integer :: ispin, ig, ierr, igwx
  real(dp) :: charge

  ! reallocate space
  ! charge
  call destroy_scf_type( rho )
  call create_scf_type( rho )

  if( lda_plus_u ) then
    if( .not. allocated(ns) ) call errore('write_pot_rho_G','occupancies not saved',1)
    ! restore ns
    rho%ns = ns
    deallocate( ns )
  endif

  ! potential
  call destroy_scf_type( v )
  call create_scf_type( v )

  ! reorder G-vectors globally
  ! what is the largest index?
  igwx = max( maxval( ig_l2g_ref ), maxval( ig_l2g ) )
  call mp_max( igwx )

  ! transform charge and density to R-space

  allocate( aux(dfftp%nnr) )
  allocate( aux1(igwx) )

  ! density
  do ispin=1,nspin
    ! copy rhoG to original global ordering provided by ig_l2g_ref
    aux1 = zero
    do ig=1,min(size(rhoG,1),ngm)
      aux1(ig_l2g_ref(ig)) = rhoG(ig,ispin)
    enddo
    call mp_sum( aux1 )
    ! distribute rhoG using new global ordering ig_l2g
    aux = zero
    do ig=1,ngm
      aux(nl(ig)) = aux1(ig_l2g(ig))
    enddo
    !CALL cft3( aux, nr1, nr2, nr3, nr1x, nr2x, nr3x, +1 )
    CALL invfft ('Dense', aux, dfftp)
    rho%of_r(:,ispin) = real(aux(:))
  enddo

  ! potential
  do ispin=1,nspin
    ! copy vrG to original global ordering provided by ig_l2g_ref
    aux1 = zero
    do ig=1,min(size(vrG,1),ngm)
      aux1(ig_l2g_ref(ig)) = vrG(ig,ispin)
    enddo
    call mp_sum( aux1 )
    ! distribute rhoG using new global ordering ig_l2g
    aux = zero
    do ig=1,ngm
      aux(nl(ig)) = aux1(ig_l2g(ig))
    enddo
    !CALL cft3( aux, nr1, nr2, nr3, nr1x, nr2x, nr3x, +1 )
    CALL invfft ('Dense', aux, dfftp)
    v%of_r(:,ispin) = real(aux(:))
  enddo

  deallocate( aux, aux1 )
  deallocate( rhoG, vrG )

  ! renormalize total charge?

  IF ( nspin == 2 ) THEN
     charge = SUM ( rho%of_r (:, 1:nspin) ) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  ELSE
     charge = SUM ( rho%of_r (:, 1) ) * omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  END IF
#ifdef __PARA
  call mp_sum( charge, intra_pool_comm )
#endif

  WRITE( stdout, &
          '(/,5X,"starting charge ",F10.5,", resized charge ",F10.5)') &
          charge, ref_charge
  write(stdout,*) ' ave abs diff = ', abs(charge-ref_charge)/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  ! warning if different
  if( abs(charge - ref_charge)/(dfftp%nr1*dfftp%nr2*dfftp%nr3) > 1.d-10 ) &
    call errore('write_pot_rho_G','resizing of charge density has led to modified total charge',-1)

  ! write potential and charge to output files

  call write_rho( rho, nspin )

  if( lda_plus_u ) then
    write(stdout,*) 'LDA+U - writing occupations'
    write(stdout,*) rho%ns
  endif

  return
  end subroutine write_pot_rho_G


  end module pot_rho_G
