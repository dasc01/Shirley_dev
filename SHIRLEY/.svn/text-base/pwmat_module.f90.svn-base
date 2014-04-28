! ----------------------------------------------------------------------
  module pwmat_module
! ----------------------------------------------------------------------

  use kinds, only : dp

  implicit none
  public

  real(dp) :: pwmtxel_cutoff
  integer :: ngk_pwg
  integer :: nbnd
  real(dp) :: at(3,3), bg(3,3)
  real(dp) :: omega
  real(dp) :: tpiba, tpiba2
  real(dp),allocatable :: g_pw(:,:)
  
  integer :: nqpg
  real(dp),allocatable :: qpgvec(:,:), qpgvec_len(:)

  integer :: iunpwm, pwmrecl

  contains


! ----------------------------------------------------------------------
  subroutine pwi_read( pwifile, root, mpime )
! ----------------------------------------------------------------------

  character(*),intent(in) :: pwifile
  integer,intent(in) :: root, mpime
  integer,external :: freeunit

  integer :: iunpwi

  if( mpime == root ) then

  iunpwi = freeunit()
  open(iunpwi,file=trim(pwifile),form='unformatted')
  read(iunpwi) pwmtxel_cutoff ! cut-off
  read(iunpwi) ngk_pwg  ! number of G vectors
  read(iunpwi) nbnd     ! number of bands
  read(iunpwi) at       ! transformation matrix of bravais lattice
  read(iunpwi) bg       ! transformation matrix of reciprocal lattice
  read(iunpwi) omega    ! cell volume
  read(iunpwi) tpiba    ! 2 * pi / alat

  tpiba2 = tpiba*tpiba

  allocate( g_pw(3,ngk_pwg) )
  read(iunpwi) g_pw  ! the G-vectors in Cartesian coords
  close(iunpwi)

  endif
  call pwi_bcast( root, mpime )

  end subroutine pwi_read


! ----------------------------------------------------------------------
  subroutine pwi_bcast( root, mpime )
! ----------------------------------------------------------------------

  use mp, only : mp_bcast

  integer,intent(in) :: root, mpime

  call mp_bcast( pwmtxel_cutoff, root )
  call mp_bcast( ngk_pwg, root )
  call mp_bcast( nbnd, root )
  call mp_bcast( at, root )
  call mp_bcast( bg, root )
  call mp_bcast( omega, root )
  call mp_bcast( tpiba, root )

  if( mpime /= root ) tpiba2 = tpiba*tpiba

  if( mpime /= root ) allocate( g_pw(3,ngk_pwg) )
  call mp_bcast( g_pw, root )

  end subroutine pwi_bcast


! ----------------------------------------------------------------------
  subroutine pw_qpgcut( qvec, ecut )
! ----------------------------------------------------------------------
  real(dp) :: qvec(3)
  real(dp) :: ecut

  real(dp) :: ecutg
  real(dp) :: qpg(3), qpg2
  integer :: i

  ! renormalize ecut to G-vector units
  ecutg = ecut / tpiba2
  ! cut off the list of G-vectors + qvec
  nqpg=0
  do i=1,ngk_pwg
    qpg = qvec + g_pw(:,i)
    qpg2 = dot_product( qpg, qpg )
    if( qpg2 <= ecutg ) nqpg = nqpg + 1
  enddo

  allocate( qpgvec(3,nqpg), qpgvec_len(nqpg) )
  nqpg=0
  do i=1,ngk_pwg
    qpg = qvec + g_pw(:,i)
    qpg2 = dot_product( qpg, qpg )
    if( qpg2 <= ecutg ) then
      nqpg = nqpg + 1
      qpgvec(:,nqpg) = qpg
      qpgvec_len(nqpg) = sqrt( qpg2 )
    endif
  enddo
    
  end subroutine pw_qpgcut


! ----------------------------------------------------------------------
  subroutine pwm_open( pwmfile )
! ----------------------------------------------------------------------
  character(*),intent(in) :: pwmfile
  integer,external :: freeunit

  iunpwm=freeunit()
  pwmrecl = nbnd*nbnd*2*8  ! 8 bytes per word
  open(unit=iunpwm,file=trim(pwmfile),form='unformatted',status='old', &
       access='direct',action='read',recl=pwmrecl)

  end subroutine pwm_open


! ----------------------------------------------------------------------
  subroutine pwmtxel_int( ig, eigk1, eigk2, iblk1, iblk2, pwm1, pwm2 )
! ----------------------------------------------------------------------
! Assumes that we only want two blocks from the entire matrix
  use io_global, only : stdout
  integer,intent(in) :: ig, iblk1, iblk2
  complex(dp),intent(in) :: eigk1(nbnd,nbnd), eigk2(nbnd,nbnd)
  complex(dp),intent(out) :: pwm1(nbnd-iblk1,iblk2), pwm2(iblk1,nbnd-iblk2)

  complex(dp),allocatable :: pwm_in(:,:)
  complex(dp),allocatable :: pwm_tmp(:,:)
  complex(dp),parameter :: ONE =cmplx(1.d0,0.d0)
  complex(dp),parameter :: ZERO=cmplx(0.d0,0.d0)

  allocate( pwm_in(nbnd,nbnd) )
  read(unit=iunpwm,rec=ig) pwm_in

  ! block 1 (lower left)
  allocate( pwm_tmp(nbnd,iblk2) )
  call ZGEMM('N','N',nbnd,iblk2,nbnd,ONE,pwm_in,nbnd, &
             eigk2,nbnd,ZERO,pwm_tmp,nbnd)
  call ZGEMM('C','N',nbnd-iblk1,iblk2,nbnd,ONE,eigk1(1,iblk1+1),nbnd, &
             pwm_tmp,nbnd,ZERO,pwm1,nbnd-iblk1)

  ! block 2 (upper right)
  deallocate( pwm_tmp )
  allocate( pwm_tmp(nbnd,nbnd-iblk2) )
  call ZGEMM('N','N',nbnd,nbnd-iblk2,nbnd,ONE,pwm_in,nbnd, &
             eigk2(1,iblk2+1),nbnd,ZERO,pwm_tmp,nbnd)
  call ZGEMM('C','N',iblk1,nbnd-iblk2,nbnd,ONE,eigk1,nbnd, &
             pwm_tmp,nbnd,ZERO,pwm2,iblk1)

  deallocate( pwm_in, pwm_tmp )

  end subroutine pwmtxel_int


  end module pwmat_module

