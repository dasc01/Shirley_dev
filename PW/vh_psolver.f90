subroutine eh_psolver(etot,ehart)
  USE kinds,            ONLY : DP
  USE grid_dimensions,  ONLY : nrxx
  USE lsda_mod,         ONLY : nspin
  USE ions_base,        ONLY : nat
  USE scf,              ONLY : scf_type, rho, dna, drhoa, vna, dvna
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : omega
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE ener,             ONLY : ehh, edd, ehd, ehps
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(IN)::etot,ehart
  !
  REAL(DP)::charge
  INTEGER::ngrid
  INTEGER::is
  !get difference of NA and scf rho
  drhoa%of_r = rho%of_r - dna%of_r
  CALL vh_psolver(drhoa%of_r, edd, charge, dvna%of_r)
  !calculate host-defect interaction
  ngrid=dfftp%nr1*dfftp%nr2*dfftp%nr3
  do is=1,nspin
     ehd=SUM( vna%of_r(:,is)*drhoa%of_r(:,is) )
  enddo
  ehd=ehd*omega/ ( ngrid )
  CALL mp_sum(  ehd , intra_pool_comm )
  ehps=ehh+ehd+edd
  
  !add vna and dvna
  !dvna%of_r =  vna%of_r + dvna%of_r 
  !write(6,*) "ehh, ehd, edd, ehps =", ehh, ehd, edd, ehps
  !
  write(6,'(A,5F18.8)') "ehps, charge, ehh, ehd, edd =", ehps, charge, ehh, ehd, edd
  write(6,'(A,F18.8)') "Total psolver energy =", etot-ehart+ehps
end subroutine eh_psolver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine vh_psolver(rho, ehart, charge, v)
  USE kinds, ONLY : DP
  USE grid_dimensions, ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY: alat,at,omega
  USE fft_base,   ONLY : grid_scatter, grid_gather, dfftp
  USE mp_global, ONLY : me_pool,nproc_pool,inter_pool_comm,root_pool
  USE mp,        ONLY : mp_bcast
  use Poisson_Solver
  
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::rho(nrxx,nspin)
  REAL(DP), INTENT(INOUT)::v(nrxx,nspin)
  REAL(DP), INTENT(INOUT)::ehart,charge

  !Internal Variables
  REAL(DP)::ax,ay,az,dx,dy,dz
  REAL(DP)::eh,offset,chg
  CHARACTER(len=2) :: geocode
  INTEGER::itype_scf,ixc,ispin,ix,iy,iz
  INTEGER :: n3d,n3p,n3pi,i3xcsh,i3s
  REAL(DP), pointer :: karray(:) 
  REAL(DP), allocatable :: rhopot(:),xc_pot(:),aux(:)
  !Parallel vars
  INTEGER::iproc, nproc
  REAL(DP)::btoa,cell(3,3)
  !interpolation variables
  REAL(DP), allocatable ::rhou(:)
  integer::mr1,mr2,mr3
  REAL(DP)::hx,hy,hz

  allocate(rhopot(nr1*nr2*nr3),aux(nr1*nr2*nr3))
  rhopot(:)=1.0d-19
#ifdef __PARA
  do ispin=1,nspin
     aux=0.0d0
     CALL grid_gather(rho(:,ispin),aux)
     rhopot(:)=rhopot(:)+aux(:)
  enddo
  CALL mp_bcast(rhopot,root_pool)
#else
  do ispin=1,nspin
     rhopot(:)=rhopot(:)+rho(:,ispin)
  enddo
#endif
  deallocate(aux)

  !Espresso grid
  ax= SQRT(DOT_PRODUCT(at(:,1),at(:,1)))*alat
  ay= SQRT(DOT_PRODUCT(at(:,2),at(:,2)))*alat
  az= SQRT(DOT_PRODUCT(at(:,3),at(:,3)))*alat
  cell(:,:)=at(:,:)*alat

  dx= ax/real(nr1)
  dy= ay/real(nr2)
  dz= az/real(nr3)

  hx=min(dx,dy,dz)
  hy=hx
  hz=hx

  if(dx.ne.dy .or. dy.ne.dz .or. dz.ne.dx) then
     mr1=ceiling(ax/hx)
     mr2=ceiling(ay/hy)
     mr3=ceiling(az/hz)
     allocate(rhou(mr1*mr2*mr3))
     rhou=1.0d-19
     call fwdinterp(nr1,nr2,nr3,dx,dy,dz,rhopot,mr1,mr2,mr3,hx,rhou)
  else
     mr1=nr1
     mr2=nr2
     mr3=nr3
     rhou=1.0d-19
     allocate(rhou(mr1*mr2*mr3))
     rhou(:)=rhopot(:)
  endif
  deallocate(rhopot)

  write(*,*) "ax, ay, az, hx, hy, hz =",ax, ay, az, hx, hy, hz 

  iproc=me_pool
  nproc=nproc_pool
  geocode='F'
  itype_scf=16
  ixc=0

  call PS_dim4allocation(geocode,'G',iproc,nproc,mr1,mr2,mr3,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s)
  write(*,*) "n3d etc=",n3d,n3p,n3pi,i3xcsh,i3s
  call createKernel(iproc,nproc,geocode,mr1,mr2,mr3,hx,hy,hz,itype_scf,karray,.true.)

!!$  if(nr1*nr2*nr3 .ne. nr1x*nr2x*nr3x) STOP "vh_psolver: nr1*nr2*nr3 .ne. nr1x*nr2x*nr3x !!" ;


  allocate(xc_pot(n3pi))

  chg=0.
  do iz=1,mr3
     do iy=1,mr2
        do ix=1,mr1
           chg=chg+rhou(ix+(iy-1)*mr1+(iz-1)*mr1*mr2)
        enddo
     enddo
  enddo

  if(me_pool .eq. 0) call cube(cell,mr1,mr2,mr3,hx,hy,hz,1,mr1*mr2*mr3,rhou,111)

  offset=0
  call H_potential(geocode,'G',iproc,nproc,mr1,mr2,mr3,hx,hy,hz,&
       rhou,karray,xc_pot,eh,offset,.false.) !optional argument

  rhou(:)=rhou(:)*2.0d0 !Rydberg units

!  if(me_pool .eq. 0) call cube(cell,mr1,mr2,mr3,hx,hy,hz,1,mr1*mr2*mr3,rhou,222)

  if(dx .eq. dy .and. dx.eq.dz) then
     do ispin=1,nspin
#ifdef __PARA
        CALL grid_scatter(rhou,v(:,ispin))
#else
        v(:,ispin)=rhou(:)
#endif
     enddo
  else

  endif

  ehart=eh*2.0d0
  charge=chg*hx*hy*hz
  write(*,*) "EH,CHARGE=",ehart,charge
  deallocate(xc_pot,rhou)
end subroutine vh_psolver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fwdinterp(nr1,nr2,nr3,d1,d2,d3,rhopot,mr1,mr2,mr3,h1,rhou)
  use bspline

  integer::nr1,nr2,nr3,mr1,mr2,mr3
  double precision::d1,d2,d3,h1
  double precision::rhopot(nr1,nr2,nr3),rhou(mr1*mr2*mr3)

  integer::k1,k2,k3
  parameter (k1=4,k2=4,k3=4) !interpolation order
!  double precision::x1(nr1),x2(nr2),x3(nr3),kn1(nr1+k1),kn2(nr2+k2),kn3(nr3+k3),bsp(nr1*nr2*nr3)
  double precision, allocatable::x1(:),x2(:),x3(:),kn1(:),kn2(:),kn3(:),bsp(:),aux(:,:,:)
  double precision::u1(mr1),u2(mr2),u3(mr3),x,y,z

  integer::np1,np2,np3
  integer::i,j,k,ind,ii,jj,kk

  np1=nr1+1
  np2=nr2+1
  np3=nr3+1
  allocate(x1(np1),x2(np2),x3(np3),kn1(np1+k1),kn2(np2+k2),kn3(np3+k3),bsp(np1*np2*np3),aux(np1,np2,np3))

  do k=1,np3
     kk=mod(k-1,nr3)+1
     do j=1,np2
        jj=mod(j-1,nr2)+1
        do i=1,np1
           ii=mod(i-1,nr1)+1
           aux(i,j,k)=rhopot(ii,jj,kk)
        enddo
     enddo
  enddo

  do i=1,np1
     x1(i)=(i-1)*d1
  enddo
  do i=1,np2
     x2(i)=(i-1)*d2
  enddo
  do i=1,np3
     x3(i)=(i-1)*d3
  enddo

  call dbsnak(np1,x1,k1,kn1)
  call dbsnak(np2,x2,k2,kn2)
  call dbsnak(np3,x3,k3,kn3)
  call dbs3in(np1,x1,np2,x2,np3,x3,aux,np1,np2,k1,k2,k3,kn1,kn2,kn3,bsp)

  deallocate(aux)

  do i=1,mr1
     u1(i)=(i-1)*h1
  enddo
  do i=1,mr2
     u2(i)=(i-1)*h1
  enddo
  do i=1,mr3
     u3(i)=(i-1)*h1
  enddo

  do k=1,mr3
     z=min(u3(k),x3(np3))
     do j=1,mr2
        y=min(u2(j),x2(np2))
        do i=1,mr1
           x=min(u1(i),x1(np1))
           ind=(i+(j-1)*mr1+(k-1)*mr1*mr2)
           rhou(ind)=dbs3vl(x,y,z,k1,k2,k3,kn1,kn2,kn3,np1,np2,np3,bsp)
!           rhou(i,j,k)=dbs3vl(x,y,z,k1,k2,k3,kn1,kn2,kn3,nr1,nr2,nr3,bsp)
        enddo
     enddo
  enddo
  
end subroutine fwdinterp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine cube(acell,n01,n02,n03,hx,hy,hz,na,dim,rhopot,fh)
  double precision::acell(3,3),hx,hy,hz
  integer::n01,n02,n03,na,dim,fh
  double precision,intent(IN)::rhopot(dim)
  integer::ix,iy,iz

  write(fh,'(A)') "COM1"
  write(fh,'(A)') "COM2"
  write(fh,'(I5,3F12.6)') 1, 0.0, 0.0, 0.0
  write(fh,'(I5,3F12.6)') n01, hx, 0, 0
  write(fh,'(I5,3F12.6)') n02,  0,hy, 0
  write(fh,'(I5,3F12.6)') n03,  0, 0, hz
  write(fh,'(I5,4F12.6)') 1, 0.0, 0.0, 0.0, 0.0
  do iz=1,n03
     do iy=1,n02
        do ix=1,n01
           write(fh,'(F12.6)',ADVANCE="NO") rhopot(ix+(iy-1)*n01+(iz-1)*n01*n02)
           if(MOD(ix,6).eq.0) write(fh,'(A)') " "
        enddo
        write(fh,'(A)') " "
     enddo
  enddo
  
end subroutine cube
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$program PS_Test
!!$  use Poisson_Solver
!!$
!!$  implicit none
!!$  integer :: iproc,nproc
!!$  integer :: n01,n02,n03
!!$  integer :: itype_scf,ixc
!!$  integer :: n3d,n3p,n3pi,i3xcsh,i3s
!!$  integer :: ix,iy,iz
!!$  double precision, pointer :: karray(:)  
!!$  double precision, allocatable :: rhopot(:),xc_pot(:)
!!$  character(len=2) :: geocode
!!$  double precision::hx,hy,hz,acell    
!!$  double precision::eh,offset
!!$
!!$  iproc=0
!!$  nproc=1
!!$  geocode='F'
!!$  n01=32
!!$  n02=32
!!$  n03=32
!!$  acell=10.0d0
!!$  hx=acell/real(n01)
!!$  hy=acell/real(n02)
!!$  hz=acell/real(n03)
!!$  itype_scf=16
!!$  ixc=0
!!$
!!$  call createKernel(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,karray,.true.)
!!$  call PS_dim4allocation(geocode,'G',iproc,nproc,n01,n02,n03,ixc,&
!!$     n3d,n3p,n3pi,i3xcsh,i3s)
!!$  write(*,*) n3d,n3p,n3pi,i3xcsh,i3s
!!$  allocate(rhopot(n01*n02*n03))
!!$  allocate(xc_pot(1))
!!$  rhopot(:)=1.0d-19
!!$  call gaussian(n01,n02,n03,hx,hy,hz,n01*n02*n03,rhopot)
!!$  call cube(acell,n01,n02,n03,1,n01*n02*n03,rhopot,111)
!!$  offset=0
!!$  call H_potential(geocode,'G',iproc,nproc,n01,n02,n03,hx,hy,hz,&
!!$       rhopot,karray,xc_pot,eh,offset,.false.) !optional argument
!!$  call cube(acell,n01,n02,n03,1,n01*n02*n03,rhopot,222)
!!$end program PS_Test
!!$
!!$subroutine gaussian(n01,n02,n03,hx,hy,hz,dim,rhopot)
!!$  integer::n01,n02,n03,dim
!!$  double precision::hx,hy,hz,r,sum,gau,neu
!!$  double precision,intent(inout) :: rhopot(dim)
!!$  integer::ix,iy,iz,og(3)
!!$  double precision::se,sn,pi
!!$
!!$  pi=4.0d0*atan(1.0d0)
!!$  og(1)=floor(n01*1.0d0/2.)
!!$  og(2)=floor(n02*1.0d0/2.)
!!$  og(3)=floor(n03*1.0d0/2.)
!!$  sum=0.
!!$  se=0.75
!!$  sn=0.50
!!$  do iz=1,n03
!!$     do iy=1,n02
!!$        do ix=1,n01
!!$           r=sqrt( ((ix-og(1))*hx)**2 + ((iy-og(2))*hy)**2 + ((iz-og(3))*hz)**2 )
!!$           gau=(1.0d0/(se*sqrt(pi))**3)*exp(-(r/se)**2)
!!$           neu=(1.0d0/(sn*sqrt(pi))**3)*exp(-(r/sn)**2)
!!$           rhopot(ix+(iy-1)*n01+(iz-1)*n01*n02)=gau-neu
!!$           sum=sum+gau-neu
!!$        enddo
!!$     enddo
!!$  enddo
!!$  write(6,*) "sum =",sum*hx*hy*hz
!!$
!!$end subroutine gaussian
!!$
