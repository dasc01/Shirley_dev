!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
program addlocalproj_upf
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !  Add a local projector to the UPF file
  !  Dump to a new file in UPF again
  !
  use upf

  implicit none
  !
  real(8),parameter :: eps12=1.d-12
  real(8),parameter :: lambda=6.d0 ! exponent for local projector
  !
  integer :: ios, iunps = 4
  character (len=256) :: filein
  integer :: nb, lloc, ichi_loc, i, ikk_loc
  real(8) :: rcut_loc, epsloc, dion_loc
  !
  ! pp_beta
  real(8), allocatable :: betar_old(:,:)
  integer, allocatable:: lll_old(:), ikk2_old(:)
  ! pp_dij
  real(8), allocatable :: dion_old(:,:)
  ! pp_pswfc
  real(8), allocatable :: chi_loc(:)
  real(8), allocatable :: betar_loc(:)
  !
  real(8),external :: intradialprod
  !
  !---------------------------------------------------------------------- 
  ! read
  print '(''  Input PP file in UPF format > '',$)'
  read (5, '(a)', end = 20, err = 20) filein
  open(unit=iunps,file=filein,status='old',form='formatted',iostat=ios)
  if (ios.ne.0) stop
  call read_upf(iunps)
  close (unit=iunps)

  !---------------------------------------------------------------------- 
  ! do something to pseudo
  if( pseudotype /= 'NC' ) then
    write(*,*) ' error: this code is not yet ready to handle non-normconserving pseudopotentials'
    goto 20
  endif

  write(*,*) ' number of projectors = ', nbeta
  do nb=1,nbeta
    write(*,*) ' proj ', nb, ' l=', lll(nb)
  enddo
  write(*,*) ' number of atomic wave functions = ', ntwfc
  do nb=1,ntwfc
    write(*,*) '  chi ', nb, ' l=', lchiw(nb)
  enddo
  print '(''  Angular momentum of local channel > '',$)'
  read (*,'(i)',end=20,err=20) lloc

  ! check on existing wave functions
  if( ntwfc == 0 ) then
    write(*,*) ' There are no atomic wave functions for generating this projector'
    if( lloc > 0 ) then
      write(*,*) ' hydrogenic wave implemented only for s-waves'
      goto 20
    endif
    write(*,*) ' Using a hydrogenic wave instead'
    ichi_loc=0
  else
    do nb=1,ntwfc
      if( lloc==lchiw(nb) ) write(*,*) ' chi ', nb, ' has the correct ang mom'
    enddo
    print '(''  Pick which chi to use in constructing the local projector > '',$)'
    read(*,'(i)',end=20,err=20) ichi_loc
    if( ichi_loc < 1 .or. ichi_loc > ntwfc ) then
      write(*,*) ' error: chosen chi is out of bounds' ; goto 20
    else if( lchiw(ichi_loc)/=lloc ) then
      write(*,*) ' error: chosen chi has the wrong ang momentum'; goto 20
    endif
  endif

  ! cut-off
  ! recompute ikk2
  do nb=1,nbeta
    i=size(betar,1)
    do while( i > 1 .and. abs(betar(i,nb)) < eps12 )
      i=i-1
    enddo
    ikk2(nb) = i
  enddo
  if( nbeta > 0 ) then
    ikk_loc = maxval(ikk2(1:nbeta)) 
    rcut_loc = r(ikk_loc)
  else
    rcut_loc = 0.5d0
    do i=1,size(r)
      if( r(i) > rcut_loc ) then
        ikk_loc = i-1
        exit
      endif
    enddo
  endif
  write(*,*) ' chosen cut-off radius = ', rcut_loc

  ! copy old projectors and resize
  ! pp_beta
  allocate( betar_old(size(betar,1),size(betar,2)), &
            lll_old(size(lll,1)), &
            ikk2_old(size(ikk2,1)) )
  betar_old = betar
  lll_old = lll
  ikk2_old = ikk2
  deallocate( betar ) ; allocate( betar(size(betar_old,1),size(betar_old,2)+1) )
  deallocate( lll ) ; allocate( lll(size(lll_old,1)+1) )
  deallocate( ikk2 ) ; allocate( ikk2(size(ikk2_old,1)+1) )
  betar(:,1:size(betar,2)-1) = betar_old
  lll(1:size(lll)-1) = lll_old
  ikk2(1:size(ikk2)-1) = ikk2_old
  ! pp_dij
  allocate( dion_old(size(dion,1),size(dion,2)) )
  dion_old = dion
  deallocate( dion ) ; allocate( dion(size(dion_old,1)+1,size(dion_old,2)+1) )
  dion=0.d0
  dion(1:size(dion,1)-1,1:size(dion,2)-1) = dion_old

  ! generate chi_loc and betar_loc
  allocate( chi_loc(size(r)), betar_loc(size(r)) )
  if( ichi_loc<1 ) then
  ! compute chi_loc
    epsloc = (-log(eps12))
    ! compute ikk_loc
    i=size(r)
    do while( i > 1 .and. r(i) > epsloc )
      i=i-1
    enddo
    ikk_loc = i
    chi_loc = 0.d0
    do i=1,ikk_loc
      ! hydrogenic s-wave
      chi_loc(i) = exp( -r(i) )
    enddo
    ! why the factor of 2 here - surely it's irrelevant
    chi_loc = 2.d0 * r * chi_loc  ! dont forget that we store r*chi_loc

  ! compute betar_loc
    epsloc = ((-log(eps12))**(1.d0/lambda))*rcut_loc
    ! compute ikk_loc
    i=size(r)
    do while( i > 1 .and. r(i) > epsloc )
      i=i-1
    enddo
    ikk_loc = i
    betar_loc = 0.d0
    do i=1,ikk_loc
      ! hydrogenic s-wave
      betar_loc(i) = exp( -(r(i)/rcut_loc)**lambda )*chi_loc(i)
    enddo
  else
    ! modify rcut_loc to give beta_loc=0 outside the cut-off radius
    rcut_loc = rcut_loc / ((-log(eps12) )**(1.d0/lambda))
    chi_loc = chi(:,ichi_loc)
    betar_loc = 0.d0
    do i=1,ikk_loc
      betar_loc(i) = exp(-(r(i)/rcut_loc)**lambda)*chi_loc(i)
    enddo
    
  endif

  ! dion = 1 / integral beta * chi
  dion_loc = 1.d0 / intradialprod( ikk_loc, r, betar_loc, chi_loc )

  ! update pseudo-variables
  nbeta = size(betar,2)
  betar(:,size(betar,2)) = betar_loc
  lll(size(lll)) = lloc
  ikk2(size(ikk2)) = ikk_loc
  dion(size(dion,1),size(dion,2)) = dion_loc
  
  write(*,*) ' new local projector has been made'

  !---------------------------------------------------------------------- 
  ! write
  print '(''  Output PP file in UPF format > '',$)'
  read (5, '(a)', end = 20, err = 20) filein
  open(unit=iunps,file=filein,form='formatted',iostat=ios)
  if (ios.ne.0) stop
  call write_upf(iunps)
  close (unit=iunps)

20 stop
end program addlocalproj_upf
!
