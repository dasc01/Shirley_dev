  subroutine trnvect_wrapper( x, t, iflag )

  use kinds, only: dp
  use cell_base, only: at, bg, tpiba

  implicit none

  real(dp),intent(in) :: x(3)
  real(dp),intent(out) :: t(3)
  integer,intent(in) :: iflag

  t = x
  ! call trnvect(t, at, bg, iflag)
  t = t * tpiba
  call cryst_to_cart(1, t, bg, iflag)
  
  return
  end subroutine trnvect_wrapper
