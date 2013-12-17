  program spline_plot

  use kinds,only : dp
  USE io_global,  ONLY : stdout
  use splines_module

  implicit none

  integer,parameter :: maxchar=255

  type(spline_fit_type) :: spl_fit
  complex(dp),allocatable :: spl_val(:,:)
  integer :: n,m
  character(maxchar) :: spl_file
  character(maxchar) :: cq1, cq2, cq3

  integer :: iunspl, ios, i, j
  real(dp) :: qvec(3)

  integer,external :: iargc
  integer,external :: freeunit

  call getarg( 1, cq1 )
  call getarg( 2, cq2 )
  call getarg( 3, cq3 )
  call getarg( 4, spl_file )

  read(cq1,*) qvec(1)
  read(cq2,*) qvec(2)
  read(cq3,*) qvec(3)
  write(stdout,'(a,3e)') '# evaluating spline for q=', qvec

  iunspl=freeunit()
  open(iunspl,file=trim(spl_file),form='unformatted',status='old',iostat=ios)
  if( ios/=0 ) then
    call errore('spline_plot','unable to open spline file',abs(ios))
  endif
  write(stdout,'(a)') '# reading spline'
  call read_spline_fit( iunspl, spl_fit )
  write(stdout,'(a)') '# dimension spline'
  call getdims_spline_fit( n, m, spl_fit )
  write(stdout,'(a)') '# allocate spline'
  allocate( spl_val(n,m) )

  write(stdout,'(a)') '# evaluate spline'
  call evaluate_spline_fit( qvec, spl_fit, spl_val )

  write(stdout,'(a)') '# dumping spline'
  do j=1,m
    do i=1,n
      write(stdout,'(2e)') spl_val(i,j)
    enddo
    write(stdout,*)
  enddo
  write(stdout,'(a)') '# done with spline'
  
  end program spline_plot
