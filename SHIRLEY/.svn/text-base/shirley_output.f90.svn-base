  subroutine shirley_output_bandstructure( iunit, bandfile, nband, bandstructure, outfmt )

  use kpt_module, only : bandstructure_type

  implicit none
  integer,intent(in) :: iunit
  character(*),intent(in) :: bandfile
  integer,intent(inout) :: nband
  type(bandstructure_type) :: bandstructure
  character(*) :: outfmt

  character(255) :: ic
  integer :: i

  write(iunit,*) "set nokey"
  write(iunit,*) "set style data l"
  write(iunit,*) "bandcolor=-1"
  write(iunit,*) "efermi=0.0"

  if( trim(outfmt)=='eps' ) then
    write(iunit,*) "set term postscript eps enhanced"
    write(iunit,*) "set output 'qdiag.eps'"
  endif
  
  ! make axes and special k-points labels/lines
  write(iunit,*) "set xtics ( \"
  do i=1,bandstructure%nksp-1
    write(iunit,'(a,f,a)') &
      '"'//trim(bandstructure%labelsp(i))//'"', &
      bandstructure%kpathlen(i), ", \"
  enddo
  write(iunit,'(a,f,a)') &
    '"'//trim(bandstructure%labelsp(bandstructure%nksp))//'"', &
    bandstructure%kpathlen(bandstructure%nksp), ")"
  do i=2,bandstructure%nksp-1
    write(iunit,'(a,i,a,f,a,f,a)') "set arrow ", i, &
      " from first ", bandstructure%kpathlen(i), ", graph 0 to ", &
      bandstructure%kpathlen(i), ", graph 1 lt bandcolor lw 1 nohead" 
  enddo

  ! plot data
  write(iunit,'(a, a, a)') "plot '", trim(bandfile), &
    "' u 5:($6-efermi) lt bandcolor \"
  do i=2,nband-1
    write(ic,'(i)') i+5 ; ic=adjustl(ic)
    write(iunit,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor \"
  enddo
  write(ic,'(i)') nband+5 ; ic=adjustl(ic)
  write(iunit,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor"

  end subroutine shirley_output_bandstructure
