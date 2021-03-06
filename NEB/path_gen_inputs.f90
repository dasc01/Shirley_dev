subroutine path_gen_inputs(parse_file_name,engine_prefix,nimage,root,comm)
!
USE io_files, only : find_free_unit
USE mp_global, only : mp_rank

implicit none
!
character(len=*), intent(in) :: parse_file_name
character(len=*), intent(in) :: engine_prefix
integer, intent(out) :: nimage
integer, intent(in) :: root
integer, intent(in) :: comm
!
character(len=256) :: dummy, dummy2
integer :: i, j
integer :: parse_unit, neb_unit
integer, allocatable :: unit_tmp(:)
integer :: unit_tmp_i
character(len=10) :: a_tmp

integer :: myrank


myrank =  mp_rank(comm)

parse_unit = find_free_unit()
open(unit=parse_unit,file=trim(parse_file_name),status="old")


! ---------------------------------------------------
! NEB INPUT PART
! ---------------------------------------------------

i=0
nimage = 0

neb_unit = find_free_unit()

open(unit=neb_unit,file='neb.dat',status="unknown")
dummy=""
do while (LEN_TRIM(dummy)<1)
read(parse_unit,fmt='(A256)',END=10) dummy
enddo

if(trim(dummy)=="BEGIN") then
  do while (trim(dummy)/="END")
    read(parse_unit,*) dummy
    if(trim(dummy)=="BEGIN_PATH_INPUT") then

        read(parse_unit,'(A256)') dummy

      do while (trim(dummy)/="END_PATH_INPUT")
        if(myrank==root) write(neb_unit,*) trim(dummy)
        read(parse_unit,'(A256)') dummy
      enddo
    endif
    if(trim(dummy)=="FIRST_IMAGE") then
      nimage = nimage + 1
    endif
    if(trim(dummy)=="INTERMEDIATE_IMAGE") then
      nimage = nimage + 1
    endif
    if(trim(dummy)=="LAST_IMAGE") then
      nimage=nimage+1
    endif
  enddo
else
  write(0,*) "key word BEGIN missing"
endif
close(neb_unit)
!------------------------------------------------
!
!
! ------------------------------------------------
! ENGINE INPUT PART
! ------------------------------------------------

allocate(unit_tmp(1:nimage))
unit_tmp(:) = 0

do i=1,nimage
unit_tmp(i) = find_free_unit()
enddo

do i=1,nimage
if(i>=1.and.i<10) then
write(a_tmp,'(i1)') i
elseif(i>10.and.i<100) then
write(a_tmp,'(i2)') i
elseif(i>100.and.i<1000) then
write(a_tmp,'(i3)')
endif

unit_tmp_i = unit_tmp(i)

open(unit=unit_tmp_i,file=trim(engine_prefix)//trim(a_tmp)//".in")

REWIND(parse_unit)

dummy=""
do while (LEN_TRIM(dummy)<1)
read(parse_unit,fmt='(A256)',END=10) dummy
enddo

if(trim(dummy)=="BEGIN") then
  do while (trim(dummy)/="END")
    dummy=""
    do while (LEN_TRIM(dummy)<1)
    read(parse_unit,fmt='(A256)',END=10) dummy
    enddo

    if(trim(dummy)=="BEGIN_ENGINE_INPUT") then
       dummy=""
       do while (LEN_TRIM(dummy)<1)
       read(parse_unit,fmt='(A256)',END=10) dummy
       enddo
        
       do while (trim(dummy)/="BEGIN_POSITIONS")
          if(myrank==root) write(unit_tmp_i,*) trim(dummy)
          read(parse_unit,'(A256)') dummy
        enddo
        if(i==1) then
        do while (trim(dummy)/="FIRST_IMAGE")
          read(parse_unit,'(A256)') dummy
        enddo
        if(trim(dummy)=="FIRST_IMAGE") then
            read(parse_unit,'(A256)') dummy
          do while (trim(dummy)/="INTERMEDIATE_IMAGE".and.(trim(dummy)/="LAST_IMAGE"))
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
          do while (trim(dummy)/="END_POSITIONS")
            read(parse_unit,'(A256)') dummy
          enddo
          read(parse_unit,'(A256)') dummy
          do while (trim(dummy)/="END_ENGINE_INPUT")
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
        endif
        endif
        !
        if(i==nimage) then
        do while (trim(dummy)/="LAST_IMAGE")
          read(parse_unit,'(A256)') dummy
        enddo
        if(trim(dummy)=="LAST_IMAGE") then
            read(parse_unit,'(A256)') dummy
          do while (trim(dummy)/="END_POSITIONS")
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
          read(parse_unit,'(A256)') dummy
          do while (trim(dummy)/="END_ENGINE_INPUT")
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
        endif 
        endif
        !
        if(i/=nimage.and.i/=1) then
        do j=2,i
        do while (trim(dummy)/="INTERMEDIATE_IMAGE")
          read(parse_unit,'(A256)') dummy
        enddo
        enddo
        if(trim(dummy)=="INTERMEDIATE_IMAGE") then
          read(parse_unit,'(A256)') dummy
          do while ((trim(dummy)/="LAST_IMAGE").and.trim(dummy)/="INTERMEDIATE_IMAGE")
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
          do while (trim(dummy)/="END_POSITIONS")
            read(parse_unit,'(A256)') dummy
          enddo
          read(parse_unit,'(A256)') dummy
          do while (trim(dummy)/="END_ENGINE_INPUT")
            if(myrank==root) write(unit_tmp_i,*) trim(dummy)
            read(parse_unit,'(A256)') dummy
          enddo
        endif
        endif
        !
    endif
  enddo
endif

close(unit_tmp_i)
enddo

deallocate(unit_tmp)

close(parse_unit)
!
10 CONTINUE
!
end subroutine path_gen_inputs
