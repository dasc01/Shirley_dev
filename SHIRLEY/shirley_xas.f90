  program shirley_xas

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and produce x-ray absorption spectra
  ! by averaging over a k-point grid

  ! David Prendergast, UCB, Jul 2007

#include "f_defs.h"
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
!  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_scatter_size, mp_scatter, mp_sum
  use mpio
  use scalapack_module
  use kpt_module
  use corerepair_module
  use shirley_input_module, nelec_in=>nelec
  use diag_module
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, desc_striped, context_cyclic, &
                        local_striped_dim, cyclic_localindex, &
                        create_striped_desc2, local_striped_dim2, &
                        context_global

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)


  character(255) :: eigval_file, info_file, xmat_file,eigval_file2

  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: fheigval, fheigvec, fhxmat
  integer :: reqeigval, reqreigvec, reqxmat
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  integer :: group_ranks(1), group_size, ionode_group, ionode_comm
  integer :: wcomm, wgroup
  integer :: nbasis_stripe
  integer :: nbasis_stripe2

  complex(dp),allocatable :: ztmp(:,:), beta_block(:,:)
  real(dp) :: fermi_energy

  complex(dp),allocatable :: posn(:,:,:), posn_l(:,:)

  integer :: iuninf

  integer,external :: freeunit
  integer,external :: ilaenv

  logical :: islocal

  integer :: ik, i, j, nbnd, nk, ispin,iband
  integer :: icore, ibasis, ncp, ncm, ixyz, core_root
  integer :: iatom, it, np, npm

  complex(dp),allocatable :: beta_local(:,:), eigvec_s(:,:)
  complex(dp),allocatable :: beta_local2(:,:)
  integer :: ind, ims
  logical :: first_diag

  integer :: desc1(DLEN_), desc2(DLEN_)

  real(dp) :: nelec, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin


!DASb
  integer :: becunit, b_nk, b_ik, b_iik, b_is, b_nkb, b_nbnd, b_ib, b_ikb
  integer :: b_iat, b_it, b_npra, b_ptr, ityp, b_ip
  integer :: gw_nqp, gw_nspin, iqp, xasunit
  complex(dp),allocatable ::becp(:,:,:),b_ztmp(:,:)
  real (dp), allocatable ::gwqp(:,:,:)
  integer, allocatable::bkbs(:),ksindx(:,:)
!DASe

  namelist /info/ nk, nbnd, ncp, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u, gwrun


  call shirley_input

  call dump_system( nelec, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*)
  write(stdout,*) ' System details: '
  write(stdout,*) '       nelec = ', nelec
  write(stdout,*) '        alat = ', alat
  write(stdout,*) ' cell volume = ', volume
  write(stdout,*) '       nspin = ', nspin
  write(stdout,*)

  !DASb
  if(gwrun) goto 111
  !DASe

  ! fix band_subset if undefined
  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. &
      band_subset(2) > nbasis ) band_subset(2)=nbasis


  if( any( band_subset > 0 ) ) then
    write(stdout,*) ' band_subset = ', band_subset

    call diagx_init( band_subset(1), band_subset(2) )
  else
    call diag_init
  endif

111 continue
  
  write(stdout,*) ' shirley_xas'
  write(stdout,*)

!  call kpt_scatter( kpt%list, kpt%param, ionode_id ) 

  if( corerep%ncore /= 1 ) then
    write(*,*) ' the number of cores for an xas calculation should be one'
    goto 999
  endif

!  npm = maxval( nproj_type_nl(1:ntype) )
  npm = maxval( nproj_type(1:ntype) )
  ncm = maxval( corerep%core(1:corerep%ncore)%nproj2 )

  ! store total number of k-points
  nk=kpt%list%nk

  if( mpime==root ) then
    info_file=trim(outfile)//'.info'
    if( any( band_subset > 0 ) ) then
      nbnd=band_subset(2)-band_subset(1)+1
    else
      nbnd=nbasis
    endif
    ncp = corerep%core(1)%nproj2

    fermi_energy = efermi

    write(stdout,*) ' Fermi energy = ', fermi_energy

    iuninf=freeunit()
    open(iuninf,file=info_file,form='formatted')
    write(iuninf,nml=info)
    write(iuninf,*) kpt%param%cartesian
    write(iuninf,*) trim(kpt%param%grid_type)
    write(iuninf,*) kpt%list%wk
    write(iuninf,*) kpt%list%kvec
    if( trim(kpt%param%grid_type) == 'tetrahedra' ) then
      write(iuninf,*) kpt%tetra%ntetra
      write(iuninf,*) kpt%tetra%tetra
    endif
    close(iuninf)
  endif
  call mp_bcast( nbnd, root )

  !DASb
  if(gwrun) then
     if( mpime==root ) then
        becunit=freeunit()
        open(becunit,file='becp.dmp',form='formatted',status='old')
        read(becunit,'(1I12)') b_nk
        do b_ik=1, b_nk
           read(becunit,'(4I12)') b_iik, b_is, b_nkb, b_nbnd
           if(b_ik .eq. 1) then
              allocate(becp(b_nkb,b_nbnd,b_nk))
              allocate(bkbs(b_nk))
           endif
           bkbs(b_ik)=b_is
           do b_ib=1,b_nbnd
              read(becunit,'(2F18.14)') (becp(b_ikb,b_ib,b_ik),b_ikb=1,b_nkb)
           enddo
        enddo
        close(becunit)
        becunit=freeunit()
        open(becunit,file='bands.dat',form='formatted',status='old')
        read(becunit,'(i8)') gw_nqp
        read(becunit,'(i8)') gw_nspin
        allocate(gwqp(gw_nqp,4,gw_nspin),ksindx(gw_nqp,gw_nspin))
        do ispin=1,gw_nspin
           do iqp=1,gw_nqp
              read(becunit,'(i5,4f10.5)') ksindx(iqp,ispin),gwqp(iqp,1,ispin),gwqp(iqp,2,ispin), &
                   gwqp(iqp,3,ispin),gwqp(iqp,4,ispin)
           enddo
        enddo
        gwqp(:,:,:)=gwqp(:,:,:)/RYTOEV
        close(becunit)
        neig_found=gw_nqp
     endif
     goto 112
  endif
  !DASe

  ! make another descriptor for matrix multiplies determined by splitting
  ! data across those bands we want
  call create_striped_desc2( nbnd )

!DASb
112 continue
     if(  mypoolid==mypoolroot ) then
        xasunit=freeunit()
        write(*,*) "xasunit=",xasunit
     endif
!DASe

  ! MPI-IO
  if( mypoolid==mypoolroot ) then
    eigval_file=trim(outfile)//'.eigval'
    xmat_file=trim(outfile)//'.xmat'
    call mp_file_open_dp( eigval_file, fheigval, rootpool, cross_pool_comm )
    call mp_file_open_dp( xmat_file, fhxmat, rootpool, cross_pool_comm )
  endif

!!$  if( mpime==root ) then
!!$     fheigval2=freeunit()
!!$     eigval_file2=trim(outfile)//'.eigval2'
!!$     open( fheigval2, file=eigval_file2,form='formatted',status='REPLACE')
!!$  endif


  first_diag=.true.
! ======================================================================
  do ispin=1,nspin
  do ik=1,kpt%list%nk
! ======================================================================

    if( mod(ik-1,npool)/=mypool ) cycle

    write(stdout,'(a,3f12.5,a,f12.5,i6,a,i6,a,a)') ' k-point ', &
      kpt%list%kvec(1:3,ik), ' weight ', kpt%list%wk(ik), ik,   &
      ' of ', kpt%list%nk, ' on node ', trim(nodenumber)

    if(gwrun) goto 113 !DAS

    ! build the Hamiltonian for this q-point
    ! local contribution
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )

    call mp_barrier( intra_pool_comm )

    ! diagonalize matrix to find eigenvalues and vectors
    if( any( band_subset > 0 ) ) then
      call diagx_ham
    else
      call diag_ham
    endif

    call mp_bcast( neig_found, mypoolroot, intra_pool_comm )
    if ( neig_found /= nbnd ) then
      write(stdout,*) ' only found ', neig_found, ' eigenvalues of a total ', nbnd
      call errore('shirley_xas','problem diagonalizing hamiltonian',1)
    endif

113 continue

    if( mypoolid==mypoolroot ) then
!DASb
       if(gwrun) then
          offset = ((ispin-1)*kpt%list%nk + ik-1)*neig_found
          call mpi_file_write_at( fheigval, offset, &
               gwqp(1,qpwrite,ispin), neig_found, &
               MPI_DOUBLE_PRECISION, status, ierr )
          write(*,*) "offset=",offset, neig_found
          write(*,'(F12.8)') (gwqp(b_ib,qpwrite,ispin),b_ib=1,neig_found)
       else
!DASe
    ! write eigenvalues
          offset = ((ispin-1)*kpt%list%nk + ik-1)*neig_found
          call mpi_file_write_at( fheigval, offset, &
               eigval, neig_found, &
               MPI_DOUBLE_PRECISION, status, ierr )
          write(*,*) "offset, neig_found, ik=",offset, neig_found, ik
          write(xasunit,'(F12.8)') (eigval(b_ib)*rytoev,b_ib=1,neig_found)
       endif
    endif


    if( first_diag ) then
      ! allocations
      icore=1
      iatom = corerep%core(icore)%atom
      core_root=mod(iatom-1,nproc_per_pool)

      ! This process will have access to the corresponding beta
      if( core_root==mypoolid ) then
        ! convert iatom to local index
        iatom=(iatom-1)/nproc_per_pool+1
        it = type_atom(iatom)
        np = nproj_type(it)
        ncp = corerep%core(icore)%nproj2
      endif

      if( nproc_per_pool > 1 ) call mp_bcast( np, core_root, intra_pool_comm )
      if( nproc_per_pool > 1 ) call mp_bcast( ncp, core_root, intra_pool_comm )

!DASb
      if(gwrun)then
         allocate( b_ztmp(npm,b_nbnd) )
         b_ptr=0   !initialize pointer
         do ityp=1, ntype   !loop over types
            do b_iat=1, natom !loop over atoms
               b_it=type_atom(b_iat) 
               if(b_it .eq. ityp) then  !if atom if of given type
                  b_npra=nproj_type(b_it)
                  if(iatom .eq. b_iat) then
                     do b_ip=1,b_npra       !loop over projs on atom
                        do b_ik=1, b_nk
                           if(b_ik .eq. ik .and. bkbs(b_ik) .eq. ispin) then
                              do b_ib=1,b_nbnd
                                 b_ztmp(b_ip,b_ib)=becp(b_ptr+b_ip,b_ib,b_ik)
                              enddo              !band
                           endif
                        enddo           !k-points
                     enddo              !projs
                  endif
                  b_ptr=b_ptr+b_npra    !increment pointer
               endif                    !check atom type
            enddo                       !atoms
         enddo                          !types
         allocate( ztmp(np,neig_found) )
         allocate( posn(neig_found,ncp,3) )
         do b_ib=1, b_nbnd
            ztmp(1:np,b_ib)=b_ztmp(1:np,b_ib)
         enddo
         deallocate(b_ztmp)
      else
!DASe
         allocate( beta_block(np,nbasis) )
         allocate( ztmp(np,neig_found) )
         allocate( posn(neig_found,ncp,3) )
         
         if( nproc_per_pool > 1 ) then
            call local_striped_dim( nbasis_stripe )
            
            call local_striped_dim2( nbasis_stripe2 )
            
            write(stdout,*) 'nbasis_stripe = ', nbasis_stripe, nbasis_stripe2
            
            allocate( beta_local(np,nbasis_stripe), &
                 eigvec_s(nbasis_stripe,nbasis), & 
                 beta_local2(np,nbasis_stripe2), &
                 posn_l(nbasis_stripe2,ncp) )
         endif
!DASb
      endif
!DASe
      first_diag=.false.  
    endif

    if(gwrun) goto 114  !DAS

    if( core_root==mypoolid ) then
      forall( i=1:np, j=1:nbasis ) &
        beta_block(i,j) = &
          betaq(index_betaq(i,iatom),j)
    endif

    ! construct <nk|beta>
    if( nproc_per_pool > 1 ) then

      call mp_scatter( beta_block, beta_local, core_root, intra_pool_comm )

      ! redistribute eigenbasis
      call PZGEMR2D(nbasis,nbasis, &
                    eigvec,1,1,desc_cyclic, &
                    eigvec_s,1,1,desc_striped,context_cyclic)
                   
      if( nbasis_stripe > 0 ) then

      ! expand projectors in the eigenbasis
      CALL ZGEMM( 'N', 'N', np, neig_found, nbasis_stripe, one, &
                  beta_local, np, &
                  eigvec_s, nbasis_stripe, &
                  zero, ztmp, np )

      else

      ztmp = zero

      endif

      ! sum over pool
      call mp_sum( ztmp, intra_pool_comm )

     ! now split ztmp and reuse beta_local - the choice of root?
      call mp_scatter( ztmp, beta_local2, mypoolroot, intra_pool_comm )

    else

      CALL ZGEMM( 'N', 'N', np, neig_found, nbasis, one, &
                  beta_block, np, &
                  eigvec, nbasis, &
                  zero, ztmp, np )

    endif

114 continue

    ! sum <nk|beta><psi|r|phi>
    if( nproc_per_pool > 1 ) posn=zero

    do ixyz=1,3
      if( nproc_per_pool > 1 ) then

        if( nbasis_stripe2 > 0 ) then

        CALL ZGEMM( 'C', 'N', nbasis_stripe2, ncp, &
                    np, one, beta_local2, np, &
                    corerep%core(icore)%matrix(1,1,ixyz), &
                    corerep%core(icore)%nproj1, &
                    zero, posn_l, nbasis_stripe2 )

        ind=mypoolid*nbasis_stripe2
        ims=neig_found-nbasis_stripe2*nproc_per_pool
        if( ims > 0 ) ind=ind+ims
        forall(i=1:nbasis_stripe2, j=1:ncp) &
          posn(ind+i,j,ixyz) = posn_l(i,j)

        else

        posn = zero

        endif

      else

        CALL ZGEMM( 'C', 'N', neig_found, ncp, &
                    np, one, ztmp, np, &
                    corerep%core(icore)%matrix(1,1,ixyz), &
                    corerep%core(icore)%nproj1, &
                    zero, posn(1,1,ixyz), neig_found )

      endif
    enddo

    ! accumulate posn
    if( nproc_per_pool > 1 ) call mp_sum( posn, intra_pool_comm )

    if( mypoolid==mypoolroot ) then
      ! write xas
      offset = ((ispin-1)*kpt%list%nk + ik-1)*neig_found*ncp*3*2
        call mpi_file_write_at( fhxmat, offset, &
                                posn, neig_found*ncp*3*2, &
                                MPI_DOUBLE_PRECISION, status, ierr )
    endif

! ======================================================================
  enddo ! loop over k-points ik
  enddo ! loop over spin ispin
! ======================================================================

  ! deallocate space
  if( nproc_per_pool > 1 ) deallocate( beta_local, eigvec_s, beta_local2, posn_l )
  if(allocated(posn)) deallocate( posn )
  if(allocated(beta_block)) deallocate( beta_block, ztmp )


  if( mypoolid==mypoolroot ) then
    ! close binary dump
    call mpi_file_close( fheigval, ierr )
    if( ierr/=0 ) &
      call errore('shirley_xas','problem closing eigval file',abs(ierr))

    ! close xmat file
    call mpi_file_close( fhxmat, ierr )
    if( ierr/=0 ) &
      call errore('shirley_xas','problem closing xmat file',abs(ierr))
  endif

!  close(iunout)

  ! end
  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier
  call mp_end
  write(stdout,*) ' end shirley_xas'
  stop
  
  end program shirley_xas
