  program lifetime

  ! compute lifetime/broadening due to el-ph coupling

  ! David Prendergast, UCB, Feb 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use hamq_shirley, only : init_stdout, read_hamq, nbasis
  use diag_shirley
  use elph_shirley
  USE io_global,  ONLY : stdout
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_bcast, mp_barrier, mp_end, mp_sum

  implicit none

!  integer,parameter :: stdout=6
  integer,parameter :: stdin =5
  integer,parameter :: maxchar=255
  REAL(DP), PARAMETER :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  character(len=3) :: nodenumber
  character(maxchar) :: hamqfile
  character(maxchar) :: elphfile
  character(maxchar) :: outfile
  character(maxchar) :: kpt_type
  character(maxchar) :: fmtstr
  integer :: nq
  integer :: nqgrid(3), iqgrid(3)
  real(dp),allocatable :: qvec(:,:), wq(:)
  integer :: iqs
  integer :: nqtot
  integer,allocatable :: nq_local(:)
  integer,allocatable :: iqs_local(:)
  logical :: cartesian

  real(dp),allocatable :: eigval_q(:,:)
  complex(dp),allocatable :: eigvec_q(:,:,:)

  integer :: iunhq, iunepm, iunout, ierr

  real(dp) :: ef, nelecf, def, kT
  real(dp),allocatable :: wg(:,:)

  real(dp) :: efermi=0.d0 ! eV
  real(dp) :: eta = 0.1d0 ! eV
  real(dp) :: temp = 100.d0 ! = 100 K
  logical :: debug=.false.

  real(dp) :: de1, de2, phocc1, phocc2, bosefac, fermifac, omega, wm
  complex(dp),allocatable :: elph_dv(:,:,:)
  complex(dp),allocatable :: dv1(:,:)
  complex(dp),allocatable :: dmat(:,:)
  real(dp),allocatable :: gamma_q(:,:,:)

  real(dp) :: pi, invpi

  integer,external :: freeunit

  integer :: i,j,k
  integer :: iq, ik, ipedv


  namelist / input / hamqfile, elphfile, &
                     efermi, eta, temp, &
                     outfile, debug 


  ! constant pi
  pi = acos(-1.d0)
  invpi = 1.d0 / pi

  ! initialize mpi
  CALL start_shirley (nodenumber) 

  write(stdout,*) ' shirley_efermi'
  write(stdout,*)

  if( mpime==root ) then

    read(stdin,nml=input,iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_qdiagp','problem reading namelist &input',101)

    ! conversion to internal units
    kT = temp * kelvin2rydberg
    efermi = efermi /rytoev
    eta = eta / rytoev

    read(stdin,*) kpt_type ! type of coordinates
    kpt_type = adjustl(kpt_type)

    cartesian = .false.
    if( trim(kpt_type) == 'tpiba' .or. trim(kpt_type) == 'crystal' ) then
      read(stdin,*) nq
      allocate( qvec(3,nq), wq(nq) )
      do iq=1,nq
        read(stdin,*) qvec(1:3,iq), wq(iq)
      enddo
      if( trim(kpt_type) == 'tpiba' ) cartesian = .true.
    else if( trim(kpt_type) == 'automatic' ) then
      read(stdin,*) nqgrid, iqgrid
      if( any(iqgrid > 1) .or. any(iqgrid < 0) ) &
        call errore('shirley_qdiagp','grid shift values should be 0 or 1',201)
      if( any(nqgrid < 1) ) &
        call errore('shirley_qdiagp','grid values should be positive',202)
      nq = product(nqgrid)
      allocate( qvec(3,nq), wq(nq) )
      nq=0
      do i=1,nqgrid(1)
        do j=1,nqgrid(2)
          do k=1,nqgrid(3)
            nq=nq+1
            !  xk are the components of the complete grid in crystal axis
            qvec(1,nq) = DBLE(i-1)/nqgrid(1) - DBLE(iqgrid(1))/2/nqgrid(1)
            qvec(2,nq) = DBLE(j-1)/nqgrid(2) - DBLE(iqgrid(2))/2/nqgrid(2)
            qvec(3,nq) = DBLE(k-1)/nqgrid(3) - DBLE(iqgrid(3))/2/nqgrid(3)
          enddo
        enddo
      enddo
      wq = 2.d0/dble(nq)
    else
      write(stdout,*) ' kpoints flag unrecognized'
      write(stdout,*) ' should be: tpiba, crystal, or automatic'
      call errore('shirley_qdiagp','stopping',1)
    endif
 
    ! now divide k-point set over processors
    nqtot = nq
    allocate( nq_local(nproc), iqs_local(nproc) )
    nq_local = nqtot / nproc
    do i=1,nqtot-nq_local(nproc)*nproc
      nq_local(i)=nq_local(i)+1
    enddo

    iqs_local(1)=0
    do i=2,nproc
      iqs_local(i)=iqs_local(i-1)+nq_local(i-1)
    enddo

  else  ! mpime==root

    allocate( nq_local(nproc), iqs_local(nproc) )

  endif ! mpime==root

  ! broadcast namelist input first
  call mp_bcast( hamqfile, root )
  call mp_bcast( outfile, root )
  call mp_bcast( debug, root )
  call mp_bcast( kT, root )
  call mp_bcast( efermi, root )
  call mp_bcast( eta, root )

  ! broadcast k-point stuff next
  call mp_bcast( nqtot, root )
  call mp_bcast( cartesian, root )
  if( mpime/=root ) allocate( qvec(3,nqtot), wq(nqtot) )
  call mp_bcast( nq_local, root )
  call mp_bcast( iqs_local, root )
  call mp_bcast( qvec, root )
  call mp_bcast( wq, root )
  nq = nq_local(mpime+1)
  iqs = iqs_local(mpime+1)
  
! debugging
  if( debug .and. mpime/=root ) stdout = 500+mpime

  write(stdout,*) ' # q-points         = ', nqtot
  write(stdout,*) ' # q-points locally = ', nq, ':', iqs+1, ' to ', iqs+nq


  if( mpime==root ) then
    iunhq = freeunit()
    open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) call errore('shirley_qdiagp','problem opening file '//trim(hamqfile),102)

    iunepm = freeunit()
    open(iunepm,file=trim(elphfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) call errore('shirley_qdiagp','problem opening file '//trim(elphfile),103)
  endif

  iunout = freeunit()
  open(iunout,file=trim(outfile)//'.'//trim(nodenumber),form='formatted',iostat=ierr)
  if( ierr /= 0 ) call errore('shirley_qdiagp','problem opening file '//trim(outfile)//'.'//trim(nodenumber),102)

! initialize stdout in hamq_shirley
  call init_stdout( stdout )

  write(stdout,*) ' reading hamiltonian ...'
  call read_hamq( iunhq )
  write(stdout,*) ' ... done'

  if( mpime==root ) then
    write(stdout,*) ' reading electron-phonon matrix element dimensions ...'
    call read_elph( iunepm, 1 )
    write(stdout,*) ' nbnddv = ', nbnddv, ' nbasis = ', nbasis
    write(stdout,*) ' nspindv = ', nspindv
    write(stdout,*) ' npedv = ', npedv
    if( nspindv/= 1 ) goto 999
    allocate( elph_dv(nbasis,nbasis,npedv) )
    do ipedv=1,npedv
      call read_elph( iunepm, ipedv )
      ! I consider spin unpolarized
      elph_dv(:,:,ipedv) = elph(:,:,1)
    enddo
  endif
  call mp_barrier
  call mp_bcast( npedv, root )
  if( mpime/=root ) then
    allocate( elph_dv(nbasis,nbasis,npedv) )
    allocate( w2(npedv), dyn(npedv,npedv) )
  endif
  call mp_bcast( elph_dv, root )
  call mp_bcast( w2, root )
  call mp_bcast( dyn, root )


  write(stdout,*) ' about to allocate eigenspace: ', nbasis
  call flush_unit(stdout)
  allocate( eigvec_q(nbasis,nbasis,nqtot), eigval_q(nbasis,nqtot), stat=ierr )
  allocate( dv1(nbasis,nbasis), dmat(nbasis,nbasis) )
  allocate( gamma_q(nbasis,npedv,nqtot) )
  write(fmtstr,'(a,i,a)') '(i,3e,',nbasis,'e)'

  ! find Hamiltonian solutions first
  eigval_q = 0.d0
  eigvec_q = zero
  do iq=iqs+1,iqs+nq
    write(stdout,*) ' diag_hamq for q-point ', iq-iqs, ' of ', nq
    call diag_hamq( qvec(:,iq), eigval_q(:,iq), eigvec_q(:,:,iq), cartesian )
  enddo
  ! reduce over all processors
  do iq=1,nqtot
    write(stdout,*) 'iq = ', iq
  call flush_unit(stdout)
    call mp_sum( eigval_q(:,iq) )
    call mp_sum( eigvec_q(:,:,iq) )
  enddo


  gamma_q = 0.d0

  do iq=iqs+1,iqs+nq
    write(stdout,*) ' q-point ', iq-iqs, ' of ', nq
  call flush_unit(stdout)

    ! loop over perturbations
    do ipedv=1,npedv
!!!      write(stdout,*) ' perturbation ', ipedv, ' of ', npedv

      ! if we have a translational acoustic mode
      if( w2(ipedv) < 1.d-10 ) cycle

      ! dyn is actually a displacement with 1/sqrt(M) factored in
      ! this wm looks like 1/(M w(q))
      omega = sqrt(w2(ipedv))
      wm = (1.d0/omega) * sum( conjg(dyn(:,ipedv)) * dyn(:,ipedv) )

      phocc1 = bosefunc( omega, kT )  ! absorption
      phocc2 = phocc1 + 1.d0               ! emission

      call ZGEMM( 'N', 'N', nbasis, nbasis, nbasis, one, &
                  elph_dv(1,1,ipedv), nbasis, eigvec_q(1,1,iq), nbasis, &
                  zero, dv1, nbasis )

      ! now a second loop over k
      ! check previous results - there was a bug here
      !    i had used nq instead of nqtot
      do ik=1,nqtot
!!!        write(stdout,*) ' k-point ', ik, ' of ', nq

        call ZGEMM( 'C', 'N', nbasis, nbasis, nbasis, one, &
                    eigvec_q(1,1,ik), nbasis, dv1, nbasis, &
                    zero, dmat, nbasis )

        do i=1,nbasis ! loop over q bands
          do j=1,nbasis ! loop over k bands
            ! excited quasiparticle
            fermifac = 1.d0
            if( eigval_q(i,iq) <= efermi ) fermifac = 0.d0
            fermifac = fermifunc( eigval_q(j,ik), efermi, kT ) &
                     - fermifac

            ! if there is no allowed transition
            if( abs(fermifac) < 1.d-10 ) cycle

            de1 = eigval_q(j,ik) - eigval_q(i,iq)
            de2 = de1 + omega
            de1 = de1 - omega
            bosefac = phocc1 * deltafunc( de1, eta ) &
                    + phocc2 * deltafunc( de2, eta )

            ! if no allowed phonon scattering
            if( abs(bosefac) < 1.d-10 ) cycle

            gamma_q(i,ipedv,iq) = gamma_q(i,ipedv,iq) &
                                + real( conjg(dmat(j,i))*dmat(j,i) ) & 
                                * bosefac * fermifac * wm * wq(ik)
          enddo ! k bands
        enddo ! q bands
        
      enddo ! k-point
      
    enddo ! perturbation

    write(fmtstr,*) '(f,', npedv, 'f)'
    do i=1,nbasis
      write(iunout,trim(fmtstr)) (eigval_q(i,iq)-efermi)*rytoev, &
                                  (gamma_q(i,ipedv,iq), ipedv=1,npedv)
    enddo

  enddo ! q-point

  call mp_barrier

  ! deallocate
  deallocate( qvec, wq )
  deallocate( nq_local, iqs_local )
  deallocate( eigvec_q, eigval_q )
  deallocate( gamma_q, dv1, dmat )


  999 write(stdout,*) ' end shirley_qdiag'
  call mp_barrier
  call mp_end
  stop
  
  contains


    pure function fermifunc( e, ef, kT )
    real(dp),intent(in) :: e, ef, kT
    real(dp) :: fermifunc

    fermifunc = 1.d0 / ( exp( (e-ef)/kT ) + 1 )
    end function fermifunc

    pure function bosefunc( e, kT )
    real(dp),intent(in) :: e, kT
    real(dp) :: bosefunc

    bosefunc = 1.d0 / ( exp( e/kT ) - 1 )
    end function bosefunc

    pure function deltafunc( e, eta )
    ! actually a lorentzian
    real(dp) deltafunc
    real(dp),intent(in) :: e, eta
    real(dp) :: df

    df = eta*0.5d0
    df = df / ( df*df + e*e )
    deltafunc = invpi * df
    end function deltafunc

  end program lifetime
