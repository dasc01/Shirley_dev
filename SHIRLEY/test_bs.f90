program test_bs
  use bspline
!!$  integer:: nx, kxord, n2x
!!$  double precision,allocatable::y1(:),y2(:),x1(:),x2(:),x1k(:),bx(:)
!!$  double precision::dx1, dx2
!!$
!!$  integer::i,j
!!$  nx=100
!!$  kxord=4
!!$  allocate(x1(nx),y1(nx),x1k(nx+kxord))
!!$  dx1=10.0/nx
!!$  do i=1,100
!!$     x1(i)=(i-1)*dx1
!!$     y1(i)=sin(x1(i))
!!$  enddo
!!$
!!$  call dbsnak(nx,x1,kxord,x1k)
!!$  allocate(bx(nx))
!!$  call dbsint(nx,x1,y1,kxord,x1k,bx)
!!$
!!$  do i=1,nx
!!$     write(*,'(2F12.8)') x1(i), y1(i)
!!$  enddo
!!$
!!$  n2x=122
!!$  dx2=10.0/n2x
!!$  allocate(x2(n2x),y2(n2x))
!!$  do i=1,n2x
!!$     x2(i)=(i-1)*dx2
!!$    if(x2(i).lt.x1(nx)) y2(i)=dbsval(x2(i),kxord,x1k,nx,bx)
!!$  enddo
!!$  write(*,'(A)') '#'
!!$  do i=1,n2x
!!$     write(*,'(2F12.8)') x2(i), y2(i)
!!$  enddo

  integer::nx,ny,kx,ky,ldz,mx,my
  double precision, allocatable::xv(:),yv(:),zv(:,:),xk(:),yk(:),bxy(:),pv(:),qv(:),g(:,:)
  double precision::dx,dy,hx,hy,Lx,Ly

  integer::i,j

  Lx=4
  Ly=7
  nx=32
  ny=48
  kx=4
  ky=4
  dx=Lx/nx
  dy=Ly/ny

  allocate(xv(nx),yv(ny),zv(nx,ny),xk(nx+kx),yk(ny+ky),bxy(nx*ny))
  do i=1,nx
     xv(i)=(i-1)*dx
  enddo
  do j=1,ny
     yv(j)=(j-1)*dy
  enddo
  do j=1,ny
     do i=1, nx
        zv(i,j)=sin(xv(i))*cos(yv(j))
     enddo
  enddo

  do i=1,nx
     do j=1,ny
        write(66,'(3F12.6)') xv(i),yv(j),zv(i,j)
     enddo
     write(*,'(A)') ' '
  enddo
  call dbsnak(nx,xv,kx,xk)
  call dbsnak(ny,yv,ky,yk)

  call dbs2in(nx,xv,ny,yv,zv,nx,kx,ky,xk,yk,bxy)
  
  hx=min(dx,dy)
  hy=hx
  mx=ceiling(Lx/hx)
  my=ceiling(Ly/hy)
  allocate(pv(mx),qv(my),g(mx,my))
  do i=1,mx
     pv(i)=(i-1)*hx
  enddo
  do j=1,my
     qv(j)=(j-1)*hy
  enddo
  do i=1,mx
     do j=1,my
        g(i,j)=dbs2vl(min(pv(i),xv(nx)),min(qv(j),yv(ny)),kx,ky,xk,yk,nx,ny,bxy)
        write(77,'(3F12.6)') pv(i),qv(j),g(i,j)
     enddo
     write(*,'(A)') ' '
  enddo
  

end program test_bs
