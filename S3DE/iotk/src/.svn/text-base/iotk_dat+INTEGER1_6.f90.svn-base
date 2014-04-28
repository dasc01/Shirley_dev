# 48 "iotk_dat.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 65 "iotk_dat.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_dat.spp"

# 78 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 6 <= __IOTK_MAXRANK
# 82 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                    intent(in)  :: unit
  character(len=*),                           intent(in)  :: name
  INTEGER (kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:) 
  type(iotk_dummytype),             optional              :: dummy
  character(len=*),                 optional, intent(in)  :: attr
  integer,                          optional, intent(in)  :: columns
  character(len=*),                 optional, intent(in)  :: sep
  character(len=*),                 optional, intent(in)  :: fmt
  integer,                          optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw,stream
  integer :: lcolumns
!-<
  integer :: dat_rank,i
  integer, dimension(:), allocatable :: dat_shape
  logical :: qe_syntax
!->
  integer(iotk_header_kind), parameter :: idummy=0
  character(100) :: lsep
  character(300) :: usefmt
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
!-<
  character(len=10) :: tmpstr
!->
  type (iotk_unit), pointer :: this
# 123 "iotk_dat.spp"
  INTEGER (kind=this_kind),allocatable :: dattmp(:)
# 125 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lcolumns = 1
  lsep(1:2) = " "//iotk_eos
  if(present(columns)) lcolumns = columns
  if(present(sep)) then
    call iotk_strcpy(lsep,sep,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
!-<  
  qe_syntax = .false.
  if (associated(this)) then
     qe_syntax = this%qe_syntax
  end if
!->
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 152 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 157 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 162 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 166 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 166 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  call iotk_write_attr(lattr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 174 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  end if
!->
# 198 "iotk_dat.spp"
!-<
  if (qe_syntax) then
    dat_rank = size(shape(dat))
    if (dat_rank>0) then
       allocate(dat_shape(dat_rank))
       dat_shape = shape(dat)
       call iotk_write_attr(lattr,"rank",dat_rank,ierr=ierrl,first=.true.)
       if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 206 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
          goto 1
       end if
       do i=1,dat_rank
     	  write(tmpstr,'(i20)') i
	  tmpstr = adjustl(tmpstr)
          call iotk_write_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),dat_shape(i),ierr=ierrl)
          if(ierrl/=0) then
             call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
             goto 1
          end if
       end do
       deallocate(dat_shape)
    end if
  end if
!->
  if(binary) then
    call iotk_write_attr(lattr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 225 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
# 230 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(lcolumns/=1) call iotk_write_attr(lattr,"columns",lcolumns,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 237 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 252 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 257 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 262 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 267 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"columns",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 272 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 277 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 283 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 288 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  else
    if(present(attr)) then
      call iotk_write_begin(unit,name,attr,ierr=ierrl)
    else
      call iotk_write_begin(unit,name,ierr=ierrl)
    end if
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 299 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
# 305 "iotk_dat.spp"
    call iotk_write_begin(unit,iotk_toupper("INTEGER"),lattr,ierr=ierrl) 
# 307 "iotk_dat.spp"
  end if
!->

  allocate(dattmp(size(dat)))
# 314 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 318 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 320 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 324 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 329 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 335 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 344 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 346 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 347 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 353 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 367 "iotk_dat.spp"
     write(lunit,fmt=trim(iotk_wfmt("INTEGER",kind(dattmp),lcolumns,-1,lsep)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 369 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
# 373 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 376 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  end if
!-<
  if (qe_syntax) then
# 386 "iotk_dat.spp"
    call iotk_write_end(unit,iotk_toupper("INTEGER"),ierr=ierrl) 
# 388 "iotk_dat.spp"
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 389 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
!->
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 396 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_6


# 753 "iotk_dat.spp"

# 755 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
!-<
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
!->
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
# 787 "iotk_dat.spp"
  INTEGER (kind=this_kind),              allocatable :: tmpdat(:)
# 789 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
!-<
! ... necessary because i need the syntax to use
  type (iotk_unit), pointer :: this
! ... necessary to read the tag describing the type
  character(iotk_taglenx) :: ltag
  character(iotk_namlenx) :: r_name
  character(iotk_attlenx) :: lattr2
  character(len=20) :: tmpstr
! ... necessary for scan_tag
  logical :: binary,stream,qe_syntax
  integer :: r_control,rrank
  integer :: rshape,i
!->
  integer :: columns
  logical :: inside,foundl
!-<
  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)
  qe_syntax = .false.
  if (associated(this)) THEN
     qe_syntax = this%qe_syntax
  end if    
!->
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 826 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!-<
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,columns,ierrl)
! Note that columns is not effectively used
  if(ierrl/=0) goto 1
!-<
  else
     call iotk_inquire(iotk_phys_unit(unit),binary,stream,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 839 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
     end if
     do
        call iotk_scan_tag(unit,+1,r_control,ltag,binary,stream,ierrl)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 845 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           goto 1
        end if
       if (r_control==4) then
	   cycle
	else
	   exit
	end if
     end do

     call iotk_tag_parse(ltag,r_name,lattr2,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 857 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
	goto 1
     end if
     call iotk_strcpy(rtype,r_name,ierrl)
     rtype = iotk_toupper(rtype)

      rlen = -1

     if (rtype(1:iotk_strlen(rtype))/="STRING") then
        call iotk_scan_attr(lattr2,"rank",rrank,ierr=ierrl,default=0)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 868 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           return
        end if
        rsize = 1
        if (rrank>0) then
           do i=1,rrank
     	      write(tmpstr,'(i20)') i
	      tmpstr = adjustl(tmpstr)
	      call iotk_scan_attr(lattr2,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
              if(ierrl/=0) then
                 call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 878 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
                 return
              end if
	      rsize = rsize*rshape
	   end do
	end if
     else
        call iotk_strcpy(rtype,"CHARACTER",ierrl)
        rsize = -1        
     end if

     call iotk_scan_attr(lattr2,"kind",rkind,ierr=ierrl,default=-1)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 891 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
     if(ierrl/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
# 896 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"columns",columns,ierr=ierrl,default=1)
     if(ierrl/=0) then
       call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
# 901 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
       return
     end if
  end if
!->
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 907 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 907 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 907 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 907 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 911 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 918 "iotk_dat.spp"
  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 921 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 921 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Error reading data')
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",name)
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"rkind",rkind)
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if
# 927 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 931 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dat,tmpdat,size(tmpdat),1)
# 933 "iotk_dat.spp"
#else
     dat = reshape(tmpdat,shape(dat))
#endif
# 937 "iotk_dat.spp"
1 continue
!-<
  if ( allocated(tmpdat) ) deallocate(tmpdat)
!->
  if(inside) then
!-<
    if (qe_syntax) then 
# 947 "iotk_dat.spp"
      call iotk_scan_end(unit,iotk_toupper("INTEGER"),ierr=ierrl2) 
# 949 "iotk_dat.spp"
       if(ierrl2/=0) then
          call iotk_error_clear(ierrl)
          ierrl=ierrl2
       end if
    end if
!->
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 964 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 964 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 964 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_INTEGER1_6

!-<
! ... new procedure: scan dat when you are already inside the data tag
# 980 "iotk_dat.spp"
subroutine iotk_scan_dat_inside_INTEGER1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif

  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
# 1005 "iotk_dat.spp"
  INTEGER (kind=this_kind),              allocatable :: tmpdat(:)
  character(len=20) :: tmpstr
  integer :: rrank
  integer :: rshape,i
# 1010 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  character(iotk_namlenx) :: rname
  type (iotk_unit), pointer :: this
  integer :: columns
  logical :: inside,foundl,qe_syntax


  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)

  qe_syntax = .true.
  IF (associated(this)) then
     qe_syntax = this%qe_syntax
  END IF

  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.

  if (.not.qe_syntax) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1034 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if

# 1042 "iotk_dat.spp"
  rname = iotk_tolower("INTEGER")
  call iotk_scan_begin(unit,iotk_toupper("INTEGER"),lattr,found=foundl,ierr=ierrl)
# 1045 "iotk_dat.spp"
  if(.not. foundl) goto 1


  foundl = .true.
  inside = .true.
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1051 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  
  rlen = -1


# 1062 "iotk_dat.spp"
     call iotk_scan_attr(lattr,"rank",rrank,ierr=ierrl,default=0)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1064 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     rsize = 1
     if (rrank>0) then
        do i=1,rrank
	   write(tmpstr,'(i20)') i
	   tmpstr = adjustl(tmpstr)
           call iotk_scan_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
           if(ierrl/=0) then
              call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1074 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
              return
           end if
           rsize = rsize*rshape
        end do
     end if
# 1081 "iotk_dat.spp"


  call iotk_scan_attr(lattr,"kind",rkind,ierr=ierrl,default=-1)
  if(ierrl/=0) then
     call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1085 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
  if(ierrl/=0) then
     call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1090 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"columns",columns,ierr=ierrl,default=1)
  if(ierrl/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1095 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if

  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1100 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 1107 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 1111 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Error reading data')
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rname",rname)
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rkind",rkind)
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if

# 1118 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 1122 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dat,tmpdat,size(tmpdat),1)
# 1124 "iotk_dat.spp"
#else
     dat = reshape(tmpdat,shape(dat))
#endif
# 1128 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
# 1134 "iotk_dat.spp"
     call iotk_scan_end(unit,iotk_toupper("INTEGER"),ierr=ierrl2) 
# 1136 "iotk_dat.spp"
     if(ierrl2/=0) then
        call iotk_error_clear(ierrl)
        ierrl=ierrl2
     end if
  end if

  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 1145 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 1145 "iotk_dat.spp"
call iotk_error_write(ierrl,"rname",rname)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_scan_dat_inside_INTEGER1_6
!->

#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_6
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_6


# 45 "iotk_dat.spp"

# 65 "iotk_dat.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_dat.spp"

# 78 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 7 <= __IOTK_MAXRANK
# 82 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                    intent(in)  :: unit
  character(len=*),                           intent(in)  :: name
  INTEGER (kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:) 
  type(iotk_dummytype),             optional              :: dummy
  character(len=*),                 optional, intent(in)  :: attr
  integer,                          optional, intent(in)  :: columns
  character(len=*),                 optional, intent(in)  :: sep
  character(len=*),                 optional, intent(in)  :: fmt
  integer,                          optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw,stream
  integer :: lcolumns
!-<
  integer :: dat_rank,i
  integer, dimension(:), allocatable :: dat_shape
  logical :: qe_syntax
!->
  integer(iotk_header_kind), parameter :: idummy=0
  character(100) :: lsep
  character(300) :: usefmt
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
!-<
  character(len=10) :: tmpstr
!->
  type (iotk_unit), pointer :: this
# 123 "iotk_dat.spp"
  INTEGER (kind=this_kind),allocatable :: dattmp(:)
# 125 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lcolumns = 1
  lsep(1:2) = " "//iotk_eos
  if(present(columns)) lcolumns = columns
  if(present(sep)) then
    call iotk_strcpy(lsep,sep,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
!-<  
  qe_syntax = .false.
  if (associated(this)) then
     qe_syntax = this%qe_syntax
  end if
!->
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 152 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 157 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 162 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 166 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 166 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 166 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  call iotk_write_attr(lattr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 174 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  end if
!->
# 198 "iotk_dat.spp"
!-<
  if (qe_syntax) then
    dat_rank = size(shape(dat))
    if (dat_rank>0) then
       allocate(dat_shape(dat_rank))
       dat_shape = shape(dat)
       call iotk_write_attr(lattr,"rank",dat_rank,ierr=ierrl,first=.true.)
       if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 206 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
          goto 1
       end if
       do i=1,dat_rank
     	  write(tmpstr,'(i20)') i
	  tmpstr = adjustl(tmpstr)
          call iotk_write_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),dat_shape(i),ierr=ierrl)
          if(ierrl/=0) then
             call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
             goto 1
          end if
       end do
       deallocate(dat_shape)
    end if
  end if
!->
  if(binary) then
    call iotk_write_attr(lattr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 225 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
# 230 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(lcolumns/=1) call iotk_write_attr(lattr,"columns",lcolumns,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 237 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 252 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 257 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 262 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 267 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"columns",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 272 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 277 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 283 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 288 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  else
    if(present(attr)) then
      call iotk_write_begin(unit,name,attr,ierr=ierrl)
    else
      call iotk_write_begin(unit,name,ierr=ierrl)
    end if
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 299 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
# 305 "iotk_dat.spp"
    call iotk_write_begin(unit,iotk_toupper("INTEGER"),lattr,ierr=ierrl) 
# 307 "iotk_dat.spp"
  end if
!->

  allocate(dattmp(size(dat)))
# 314 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 318 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 320 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 324 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 329 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 335 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 344 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 346 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 347 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 353 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 367 "iotk_dat.spp"
     write(lunit,fmt=trim(iotk_wfmt("INTEGER",kind(dattmp),lcolumns,-1,lsep)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 369 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
# 373 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 376 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  end if
!-<
  if (qe_syntax) then
# 386 "iotk_dat.spp"
    call iotk_write_end(unit,iotk_toupper("INTEGER"),ierr=ierrl) 
# 388 "iotk_dat.spp"
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 389 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
!->
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 396 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_7


# 753 "iotk_dat.spp"

# 755 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
!-<
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
!->
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
# 787 "iotk_dat.spp"
  INTEGER (kind=this_kind),              allocatable :: tmpdat(:)
# 789 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
!-<
! ... necessary because i need the syntax to use
  type (iotk_unit), pointer :: this
! ... necessary to read the tag describing the type
  character(iotk_taglenx) :: ltag
  character(iotk_namlenx) :: r_name
  character(iotk_attlenx) :: lattr2
  character(len=20) :: tmpstr
! ... necessary for scan_tag
  logical :: binary,stream,qe_syntax
  integer :: r_control,rrank
  integer :: rshape,i
!->
  integer :: columns
  logical :: inside,foundl
!-<
  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)
  qe_syntax = .false.
  if (associated(this)) THEN
     qe_syntax = this%qe_syntax
  end if    
!->
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 826 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!-<
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,columns,ierrl)
! Note that columns is not effectively used
  if(ierrl/=0) goto 1
!-<
  else
     call iotk_inquire(iotk_phys_unit(unit),binary,stream,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 839 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
     end if
     do
        call iotk_scan_tag(unit,+1,r_control,ltag,binary,stream,ierrl)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 845 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           goto 1
        end if
       if (r_control==4) then
	   cycle
	else
	   exit
	end if
     end do

     call iotk_tag_parse(ltag,r_name,lattr2,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 857 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
	goto 1
     end if
     call iotk_strcpy(rtype,r_name,ierrl)
     rtype = iotk_toupper(rtype)

      rlen = -1

     if (rtype(1:iotk_strlen(rtype))/="STRING") then
        call iotk_scan_attr(lattr2,"rank",rrank,ierr=ierrl,default=0)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 868 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           return
        end if
        rsize = 1
        if (rrank>0) then
           do i=1,rrank
     	      write(tmpstr,'(i20)') i
	      tmpstr = adjustl(tmpstr)
	      call iotk_scan_attr(lattr2,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
              if(ierrl/=0) then
                 call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 878 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
                 return
              end if
	      rsize = rsize*rshape
	   end do
	end if
     else
        call iotk_strcpy(rtype,"CHARACTER",ierrl)
        rsize = -1        
     end if

     call iotk_scan_attr(lattr2,"kind",rkind,ierr=ierrl,default=-1)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 891 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
     if(ierrl/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
# 896 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"columns",columns,ierr=ierrl,default=1)
     if(ierrl/=0) then
       call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
# 901 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
       return
     end if
  end if
!->
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 907 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 907 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 907 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 907 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 911 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 918 "iotk_dat.spp"
  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 921 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 921 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Error reading data')
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",name)
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"rkind",rkind)
# 921 "iotk_dat.spp"
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if
# 927 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 931 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dat,tmpdat,size(tmpdat),1)
# 933 "iotk_dat.spp"
#else
     dat = reshape(tmpdat,shape(dat))
#endif
# 937 "iotk_dat.spp"
1 continue
!-<
  if ( allocated(tmpdat) ) deallocate(tmpdat)
!->
  if(inside) then
!-<
    if (qe_syntax) then 
# 947 "iotk_dat.spp"
      call iotk_scan_end(unit,iotk_toupper("INTEGER"),ierr=ierrl2) 
# 949 "iotk_dat.spp"
       if(ierrl2/=0) then
          call iotk_error_clear(ierrl)
          ierrl=ierrl2
       end if
    end if
!->
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 964 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 964 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 964 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_INTEGER1_7

!-<
! ... new procedure: scan dat when you are already inside the data tag
# 980 "iotk_dat.spp"
subroutine iotk_scan_dat_inside_INTEGER1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif

  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
# 1005 "iotk_dat.spp"
  INTEGER (kind=this_kind),              allocatable :: tmpdat(:)
  character(len=20) :: tmpstr
  integer :: rrank
  integer :: rshape,i
# 1010 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  character(iotk_namlenx) :: rname
  type (iotk_unit), pointer :: this
  integer :: columns
  logical :: inside,foundl,qe_syntax


  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)

  qe_syntax = .true.
  IF (associated(this)) then
     qe_syntax = this%qe_syntax
  END IF

  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.

  if (.not.qe_syntax) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1034 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if

# 1042 "iotk_dat.spp"
  rname = iotk_tolower("INTEGER")
  call iotk_scan_begin(unit,iotk_toupper("INTEGER"),lattr,found=foundl,ierr=ierrl)
# 1045 "iotk_dat.spp"
  if(.not. foundl) goto 1


  foundl = .true.
  inside = .true.
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1051 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  
  rlen = -1


# 1062 "iotk_dat.spp"
     call iotk_scan_attr(lattr,"rank",rrank,ierr=ierrl,default=0)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1064 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     rsize = 1
     if (rrank>0) then
        do i=1,rrank
	   write(tmpstr,'(i20)') i
	   tmpstr = adjustl(tmpstr)
           call iotk_scan_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
           if(ierrl/=0) then
              call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1074 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
              return
           end if
           rsize = rsize*rshape
        end do
     end if
# 1081 "iotk_dat.spp"


  call iotk_scan_attr(lattr,"kind",rkind,ierr=ierrl,default=-1)
  if(ierrl/=0) then
     call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1085 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
  if(ierrl/=0) then
     call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1090 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"columns",columns,ierr=ierrl,default=1)
  if(ierrl/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1095 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if

  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1100 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 1107 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 1111 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Error reading data')
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rname",rname)
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rkind",rkind)
# 1111 "iotk_dat.spp"
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if

# 1118 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 1122 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dat,tmpdat,size(tmpdat),1)
# 1124 "iotk_dat.spp"
#else
     dat = reshape(tmpdat,shape(dat))
#endif
# 1128 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
# 1134 "iotk_dat.spp"
     call iotk_scan_end(unit,iotk_toupper("INTEGER"),ierr=ierrl2) 
# 1136 "iotk_dat.spp"
     if(ierrl2/=0) then
        call iotk_error_clear(ierrl)
        ierrl=ierrl2
     end if
  end if

  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
# 1145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
# 1145 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 1145 "iotk_dat.spp"
call iotk_error_write(ierrl,"rname",rname)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_scan_dat_inside_INTEGER1_7
!->

#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_7
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_7


# 45 "iotk_dat.spp"

