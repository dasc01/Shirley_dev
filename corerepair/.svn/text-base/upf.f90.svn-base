!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module upf
  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format)
  !
  integer,parameter :: lqmax=5
  ! pp_info
  integer :: rel
  real(8) :: rcloc
  integer :: nwfs
  real(8), allocatable :: oc(:), rcut(:), rcutus(:), epseu(:)
  character(len=2), allocatable :: els(:)
  integer, allocatable:: lchi (:), nns (:)
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd, pseudotype
  integer :: nv = 0
  integer :: iexch, icorr, igcx, igcc
  integer :: lmax, mesh, nbeta, ntwfc
  logical :: nlcc, isus
  character (len=20):: dft
  real(8) :: zp, ecutrho, ecutwfc, etotps
  real(8), allocatable :: ocw(:)
  character(len=2), allocatable :: elsw(:)
  integer, allocatable:: lchiw(:)
  !
  ! pp_mesh
  real(8), allocatable :: r(:), rab(:)
  !
  ! pp_nlcc
  real(8), allocatable :: rho_atc(:)
  !
  ! pp_local
  real(8), allocatable ::  vloc0(:)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8), allocatable :: betar(:,:)
  integer, allocatable:: lll(:), ikk2(:)  
  ! pp_dij
  real(8), allocatable :: dion(:,:)
  ! pp_qij
  integer ::  nqf, nqlc
  real(8), allocatable :: rinner(:), qqq(:,:), qfunc(:,:,:)
  ! pp_qfcoef
  real(8), allocatable :: qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(8), allocatable :: chi(:,:)
  !
  ! pp_rhoatom
  real(8), allocatable :: rho_at(:)
end module upf
!
!---------------------------------------------------------------------
subroutine read_upf (iunps)  
  !---------------------------------------------------------------------
  !
  !  Read pseudopotential in the Unified Pseudopotential Format (UPF)
  !
  use upf
  implicit none
  !
  integer :: iunps  
  ! iunps: unit connected with pseudopotential file
  !
  write ( *, * ) " Reading pseudopotential file in UPF format..."  
  !------->Search for Header
  write(*,*) 'info'
  call scan_begin (iunps, "INFO", .true.)  
  call read_pseudo_comment (iunps)  
  call scan_end (iunps, "INFO")  
  write(*,*) 'info'

  !------->Search for Header
  write(*,*) 'header'
  call scan_begin (iunps, "HEADER", .true.)  
  call read_pseudo_header (iunps)  
  call scan_end (iunps, "HEADER")  
  write(*,*) 'header'

  !-------->Search for mesh information
  write(*,*) 'mesh'
  call scan_begin (iunps, "MESH", .true.)  
  call read_pseudo_mesh (iunps)  
  call scan_end (iunps, "MESH")  
  write(*,*) 'mesh'
  !-------->If  present, search for nlcc
  if (nlcc  ) then  
  write(*,*) 'nlcc'
     call scan_begin (iunps, "NLCC", .true.)  
     call read_pseudo_nlcc (iunps)  
     call scan_end (iunps, "NLCC")  
  write(*,*) 'nlcc'
  endif
  !-------->Search for Local potential
  write(*,*) 'local'
  call scan_begin (iunps, "LOCAL", .true.)  
  call read_pseudo_local (iunps)  
  call scan_end (iunps, "LOCAL")  
  write(*,*) 'local'
  !-------->Search for Nonlocal potential
  write(*,*) 'nonlocal'
  call scan_begin (iunps, "NONLOCAL", .true.)  
  call read_pseudo_nl (iunps)  
  call scan_end (iunps, "NONLOCAL")  
  write(*,*) 'nonlocal'
  !-------->Search for atomic wavefunctions
  write(*,*) 'pswcf'
  call scan_begin (iunps, "PSWFC", .true.)  
  call read_pseudo_pswfc (iunps)  
  call scan_end (iunps, "PSWFC")  
  write(*,*) 'pswcf'
  !-------->Search for atomic charge
  write(*,*) 'rhoatom'
  call scan_begin (iunps, "RHOATOM", .true.)  
  call read_pseudo_rhoatom (iunps)  
  call scan_end (iunps, "RHOATOM")  
  write(*,*) 'rhoatom'
  !
  write ( *, * ) " ...done"
  return
end subroutine read_upf
!---------------------------------------------------------------------

subroutine scan_begin (iunps, string, rew)  
  !---------------------------------------------------------------------
  !
  implicit none
  ! Unit of the input file
  integer :: iunps  
  ! Label to be matched
  character (len=*) :: string  
  logical :: rew  
  ! Flag: if .true. rewind the file
  character (len=80) :: rstring  
  ! String read from file
  integer :: ios
  logical, external :: matches 

  ios = 0
  if (rew) rewind (iunps)  
  do while (ios.eq.0)  
     read (iunps, *, iostat = ios, err = 300) rstring  
     if (matches ("<PP_"//string//">", rstring) ) return  
  enddo
300 call errore ('scan_begin', 'No '//string//' block', abs (ios) )  

end subroutine scan_begin
!---------------------------------------------------------------------

subroutine scan_end (iunps, string)  
  !---------------------------------------------------------------------
  implicit none
  ! Unit of the input file
  integer :: iunps
  ! Label to be matched
  character (len=*) :: string  
  ! String read from file
  character (len=80) :: rstring
  integer :: ios
  logical, external :: matches 

  read (iunps, '(a)', iostat = ios, err = 300) rstring  
  if (matches ("</PP_"//string//">", rstring) ) return  
300 call errore ('scan_end', &
       'No '//string//' block end statement, possibly corrupted file',  - 1)
end subroutine scan_end
!
  !---------------------------------------------------------------------
  subroutine read_pseudo_comment (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine reads the comments of the UPF file
    !
    use upf
    implicit none
    integer :: ounps  

    integer :: nb, ios  
    character (len=75) :: dummy  

    read (ounps, '(a)', err = 100, iostat = ios) generated
    read (ounps, '(a)', err = 100, iostat = ios) date_author
    read (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       read (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &dummy
    else if (rel==1) then  
       read (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &dummy
    else if (rel==0) then 
       read (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & dummy
    endif

    if (rcloc > 0.d0) &
       read (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
              rcloc, dummy

    if (nwfs>0) &
       read (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) dummy
    do nb = 1, nwfs  
       read (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    enddo

    return
100 call errore ('read_pseudo_comment', 'Reading pseudo file', abs ( &
         ios))   
  end subroutine read_pseudo_comment

!---------------------------------------------------------------------

subroutine read_pseudo_header (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps  
  !
  integer :: ios, nw  
  character (len=75) :: dummy  
  logical, external :: matches 

  read (iunps, *, err = 100, iostat = ios) nv, dummy  
  read (iunps, *, err = 100, iostat = ios) psd , dummy  
  read (iunps, *, err = 100, iostat = ios) pseudotype
  if (matches (pseudotype, "US") ) isus = .true.  
  read (iunps, *, err = 100, iostat = ios) nlcc , dummy  
  read (iunps, '(a20,t24,a)', err = 100, iostat = ios) dft, dummy
  read (iunps, * ) zp , dummy  
  read (iunps, * ) etotps, dummy  
  read (iunps, * ) ecutwfc, ecutrho, dummy
  read (iunps, * ) lmax , dummy  
  read (iunps, *, err = 100, iostat = ios) mesh , dummy  
  read (iunps, *, err = 100, iostat = ios) ntwfc, nbeta , dummy
  read (iunps, '(a)', err = 100, iostat = ios) dummy
  allocate( elsw(ntwfc), lchiw(ntwfc), ocw(ntwfc) )
  do nw = 1, ntwfc
     read (iunps, * ) elsw (nw), lchiw (nw), ocw (nw)  
  enddo
  return  
100 call errore ('read_pseudo_header', 'Reading pseudo file', abs (ios))
end subroutine read_pseudo_header
!
!---------------------------------------------------------------------
subroutine read_pseudo_local (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps  
  !
  integer :: ir, ios  
  !
  allocate( vloc0(mesh) )
  read (iunps, *, err=100, iostat=ios) (vloc0(ir) , ir=1,mesh)

100 call errore ('read_pseudo_local','Reading pseudo file', abs(ios) )

  return  
end subroutine read_pseudo_local
!
!---------------------------------------------------------------------

subroutine read_pseudo_mesh (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps  
  !
  integer :: ir, ios
  !
  call scan_begin (iunps, "R", .false.)  
  allocate( r(mesh), rab(mesh) )
  read (iunps, *, err = 100, iostat = ios) (r(ir), ir=1,mesh )
  call scan_end (iunps, "R")  
  call scan_begin (iunps, "RAB", .false.)  
  read (iunps, *, err = 100, iostat = ios) (rab(ir), ir=1,mesh )
  call scan_end (iunps, "RAB")  

  return  

100 call errore ('read_pseudo_mesh', 'Reading pseudo file', abs (ios) )  
end subroutine read_pseudo_mesh
!
!---------------------------------------------------------------------

subroutine read_pseudo_nl (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps
  !
  integer :: nb, mb, n, ir, nd, ios, idum, ldum, icon, lp, i
  ! counters
  character (len=75) :: dummy  
  !
  allocate( lll(nbeta), ikk2(nbeta) )
  allocate( betar(mesh,nbeta) )
  allocate( dion(nbeta,nbeta) )
  do nb = 1, nbeta   
     call scan_begin (iunps, "BETA", .false.)  
     read (iunps, *, err = 100, iostat = ios) idum, lll(nb), dummy
     read (iunps, '(i6)', err = 100, iostat = ios) ikk2(nb)  
     read (iunps, *, err = 100, iostat = ios) &
          (betar(ir,nb), ir=1,ikk2(nb))
     do ir = ikk2(nb) + 1, mesh   
        betar (ir, nb) = 0.d0  
     enddo
     call scan_end (iunps, "BETA")  
  enddo

  call scan_begin (iunps, "DIJ", .false.)  
  read (iunps, *, err = 100, iostat = ios) nd, dummy  
  dion (:,:) = 0.d0
  do icon = 1, nd  
     read (iunps, *, err = 100, iostat = ios) nb, mb, dion(nb,mb)
     dion (mb,nb) = dion (nb,mb)  
  enddo
  call scan_end (iunps, "DIJ")  

  if (isus  ) then  
     call scan_begin (iunps, "QIJ", .false.)  
     read (iunps, *, err = 100, iostat = ios) nqf
     nqlc = 2 * lmax  + 1
     if (nqlc.gt.lqmax .or. nqlc.lt.0) &
          call errore (' read_pseudo_nl', 'Wrong  nqlc', nqlc  )
     if (nqf.ne.0) then
        allocate( rinner(nqlc) )
        call scan_begin (iunps, "RINNER", .false.)  
        read (iunps,*,err=100,iostat=ios) &
             (idum,rinner(i),i=1,nqlc)
        call scan_end (iunps, "RINNER")  
     end if
     allocate( qqq(nbeta,nbeta), qfunc(mesh,nbeta,nbeta), &
               qfcoef(nqf,nqlc,nbeta,nbeta) )
     do nb = 1, nbeta  
        do mb = nb, nbeta

           read (iunps,*,err=100,iostat=ios) idum, idum, ldum, dummy
           !"  i    j   (l)"
           if (ldum.ne.lll(mb) ) call errore ('read_pseudo_nl', &
                'inconsistent angular momentum for Q_ij', 1)

           read (iunps,*,err=100,iostat=ios) qqq(nb,mb), dummy
           ! "Q_int"
           qqq(mb,nb) = qqq(nb,mb)  

           read (iunps,*,err=100,iostat=ios) &
                        (qfunc(n,nb,mb), n=1,mesh)
           do n = 0, mesh   
              qfunc(n,mb,nb) = qfunc(n,nb,mb)  
           enddo

           if (nqf.gt.0) then
              call scan_begin (iunps, "QFCOEF", .false.)  
              read (iunps,*,err=100,iostat=ios) &
                        ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
              call scan_end (iunps, "QFCOEF")  
           end if

        enddo
     enddo
     call scan_end (iunps, "QIJ")  
  else  
     qqq (:,:) = 0.d0
     qfunc(:,:,:) =0.d0
  endif

100 call errore ('read_pseudo_nl', 'Reading pseudo file', abs (ios) )  
  return  
end subroutine read_pseudo_nl
!
!---------------------------------------------------------------------
subroutine read_pseudo_nlcc (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps  
  !
  integer :: ir, ios  

  allocate( rho_atc(mesh) )
  read (iunps, *, err = 100, iostat = ios) (rho_atc(ir), ir=1,mesh )
  !
100 call errore ('read_pseudo_nlcc', 'Reading pseudo file', abs (ios) )  
  return  
end subroutine read_pseudo_nlcc
!
!---------------------------------------------------------------------
subroutine read_pseudo_pswfc (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps
  !
  character (len=75) :: dummy  
  integer :: nb, ir, idum, ios  
  !
  allocate( chi(mesh,ntwfc) )
  do nb = 1, ntwfc  
     read (iunps,*,err=100,iostat=ios) dummy  !Wavefunction labels
     read (iunps,*,err=100,iostat=ios) (chi(ir,nb), ir=1,mesh)
  enddo
100 call errore ('read_pseudo_pswfc', 'Reading pseudo file', abs(ios))
  return  

end subroutine read_pseudo_pswfc
!
!---------------------------------------------------------------------
subroutine read_pseudo_rhoatom (iunps)  
  !---------------------------------------------------------------------
  !
  use upf
  implicit none
  !
  integer :: iunps
  !
  integer :: ir, ios  

  allocate( rho_at(mesh) )
  read (iunps,*,err=100,iostat=ios) (rho_at(ir), ir=1,mesh)
  return  

100 call errore ('read_pseudo_rhoatom','Reading pseudo file',abs(ios))

end subroutine read_pseudo_rhoatom

subroutine write_upf(ounps)

  use upf, only: nlcc

  integer :: ounps

  call write_pseudo_comment(ounps)  
  call write_pseudo_header(ounps)  
  call write_pseudo_mesh(ounps)
  if (nlcc)  call write_pseudo_nlcc(ounps)  
  call write_pseudo_local(ounps)  
  call write_pseudo_nl(ounps)  
  call write_pseudo_pswfc(ounps)
  call write_pseudo_rhoatom(ounps)  
  !
  print '("*** PLEASE TEST BEFORE USING!!! ***")'
  print '("review the content of the PP_INFO fields")'
  !
end subroutine write_upf

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_comment (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the comments of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  

    integer :: nb, ios  

    write (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"  

    write (ounps, '(a)', err = 100, iostat = ios) generated
    write (ounps, '(a)', err = 100, iostat = ios) date_author
    write (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Full-Relativistic Calculation"
    else if (rel==1) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    else if (rel==0) then 
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    endif

    if (rcloc > 0.d0) &
       write (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
              rcloc, "Local Potential cutoff radius"

    if (nwfs>0) &
       write (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
              &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    do nb = 1, nwfs  
       write (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    enddo

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"  
    return
100 call errore ('write_pseudo_comment', 'Writing pseudo file', abs ( &
         ios))   
  end subroutine write_pseudo_comment

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_header (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    character (len=4) :: shortname
    integer :: nb, ios  
    !
    !
    write (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"  

    write (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    write (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) psd , &
         "Element"
    if (pseudotype == 'NC') then  
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    else if (pseudotype == 'US') then
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    else
       call errore ('write_pseudo_header',&
            'Unknown PP type: '//pseudotype, 1)
    endif
    write (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    !call dftname (iexch, icorr, igcx, igcc, dft, shortname)
    !write (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
    !     dft, shortname," Exchange-Correlation functional"
    write (ounps, '(a,t24,a)', err = 100, iostat = ios) &
         dft," Exchange-Correlation functional"
    write (ounps, '(f17.11,t24,a)') zp , "Z valence"  
    write (ounps, '(f17.11,t24,a)') etotps, "Total energy"  
    write (ounps, '(2f11.7,t24,a)') ecutrho, ecutwfc, &
         "Suggested cutoff for wfc and rho"  

    write (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"  
    write (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    write (ounps, '(2i5,t24,a)', err = 100, iostat = ios) ntwfc, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    write (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    do nb = 1, ntwfc
       write (ounps, '(t24,a2,i3,f6.2)') elsw(nb), lchiw(nb), ocw(nb)
    enddo
    !---> End header writing

    write (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    return   
100 call errore ('write_pseudo_header','Writing pseudo file', abs(ios) )

  end subroutine write_pseudo_header

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_mesh (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  
    !
    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"  

    write (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    write (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"  
    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"  

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"  

    return

100 call errore ('write_pseudo_rhoatom','Writing pseudo file',abs(ios))

  end subroutine write_pseudo_mesh

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nlcc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the core charge for the nonlinear core
    !     correction of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"  

    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rho_atc(ir), ir = 1, mesh )
    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"  
    return

100 call errore ('write_pseudo_nlcc', 'Writing pseudo file', abs (ios))

  end subroutine write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_local (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vloc0(ir), ir = 1, mesh )
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"  
    return
100 call errore ('write_pseudo_local', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_local

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nl (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, mb, n, ir, nd, i, lp, ios  

    write (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"  
    do nb = 1, nbeta  
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"  
       write (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                    nb, lll(nb), "Beta    L"
       write (ounps, '(i6)', err=100, iostat=ios) ikk2 (nb)  
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betar(ir,nb), ir=1,ikk2(nb) )
       write (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"  
    enddo

    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"  
    nd = 0  
    do nb = 1, nbeta  
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 )  nd = nd + 1 
       enddo
    enddo
    write (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    do nb = 1, nbeta
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 ) &
             write(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, dion(nb,mb)
       enddo
    enddo
    write (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"  

    if (pseudotype == 'US') then  
       write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"  
       write (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       if (nqf.gt.0) then
          write (ounps, '(t5,a11)', err=100, iostat=ios) "<PP_RINNER>"  
          write (ounps,'(i5,1pe19.11)', err=100, iostat=ios) &
                                        (i, rinner(i), i = 1, nqlc)
          write (ounps, '(t5,a12)', err=100, iostat=ios) "</PP_RINNER>"  
       end if
       do nb = 1, nbeta 
          do mb = nb, nbeta
             write (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lll(mb) , "i  j  (l(j))"
             write (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qqq(nb,mb), "Q_int"
             write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qfunc (n,nb,mb), n=1,mesh )
             if (nqf.gt.0) then
                write (ounps, '(t5,a11)', err=100, iostat=ios) &
                                          "<PP_QFCOEF>"  
                write(ounps,'(1p4e19.11)', err=100, iostat=ios) &
                                 ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
                write (ounps, '(t5,a12)', err=100, iostat=ios) &
                                          "</PP_QFCOEF>"
             end if
          enddo
       enddo
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"  

    endif
    write (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"  
    return

100 call errore ('write_pseudo_nl', 'Writing pseudo file', abs (ios) )  

  end subroutine write_pseudo_nl

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_pswfc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the pseudo atomic functions
    !     of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"  
    do nb = 1, ntwfc
       write (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elsw(nb), lchiw(nb), ocw(nb), "Wavefunction"
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( chi(ir,nb), ir=1,mesh )
    enddo
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"  
    return

100 call errore ('write_pseudo_pswfc', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_rhoatom (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"  
    write (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rho_at(ir), ir=1,mesh )
    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"  
    return

100 call errore('write_pseudo_rhoatom','Writing pseudo file',abs(ios))
  end subroutine write_pseudo_rhoatom

  !---------------------------------------------------------------------
  subroutine dftname(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  implicit none
  integer iexch, icorr, igcx, igcc
  character (len=4) :: shortname
  character (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer, parameter :: nxc = 1, ncc = 9, ngcx = 4, ngcc = 5
  character (len=20) :: exc, corr, gradx, gradc  
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc)
  data exc / 'NOX ', 'SLA ' /  
  data corr / 'NOC ', 'PZ  ', 'VWN ', 'LYP ', 'PW  ', 'WIG ', 'HL  ',&
              'OBZ ', 'OBW ', 'GL  ' /
  data gradx / 'NOGX', 'B88 ', 'GGX ', 'PBE ', 'TPSS' /  
  data gradc / 'NOGC', 'P86 ', 'GGC ', 'BLYP', 'PBE ', 'TPSS' /  

  if (iexch==1.and.igcx==0.and.igcc==0) then
     shortname = corr(icorr)
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==5) then
     shortname = 'TPSS'
  else
     shortname = ' '
  end if
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  return
end subroutine dftname
