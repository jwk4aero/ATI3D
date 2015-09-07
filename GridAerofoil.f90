!*****
!***** 3D FLAT-PLATE GRID GENERATION
!*****

 module gridgen

 use subroutineso
 implicit none

 integer(kind=int32),parameter :: lnaca=1000
 integer(kind=int32) :: lxi0,lxi1,lxi2,let0,lze0
 integer(kind=int32) :: lxit,lett,lxie0,lxis1,lxie1,lxis2,lete0,lets1,lxisz,im,jm

 real(kind=dp),dimension(0:lnaca,2) :: xnaca,ynaca

 real(kind=dp),dimension(:,:),allocatable :: xx,yy,zz
 real(kind=dp),dimension(:),allocatable :: zs
 real(kind=dp),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(kind=dp),dimension(:),allocatable :: pxi,qet

 real(kind=dp) :: rs,re,rp,ts,te,shs,she,shswle
 real(kind=dp) :: xa,xb,xc,xd,xe,xo,ya,yb,yc,yd,yo,sho,pp,qq
 real(kind=dp) :: am,err,tmp,tmpa,tmpb,gf

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb)

 integer(kind=int32),intent(in) :: ngridv,nthick,litr
 real(kind=dp),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb

    lxit=lxi0+lxi1+lxi2+2; lett=2*let0+1
    lxie0=lxi0; lxis1=lxie0+1; lxie1=lxis1+lxi1; lxis2=lxie1+1
    lete0=let0; lets1=lete0+1

    rs=nthick*domlen/(domlen+3); re=0
    shs=smgrid; she=shs

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lze0))
    allocate(xp(0:lxit,0:3),yp(0:lxit,0:3),xq(0:lett,0:3),yq(0:lett,0:3))
    allocate(pxi(0:lxit),qet(0:lett))

!----- WAVY LEADING-EDGE PROFILE

 if(myid==mo(mb)) then
    open(9,file=cgrid); close(9,status='delete')
    open(9,file=cgrid,access='stream',form='unformatted')

 do k=0,lze0
    zs(k)=span*(real(lze0-k,kind=dp)/lze0-half)

    xa=-domlen; xb=-half+wlea*sin(2*pi*(zs(k)-zs(0))/wlew); xc=half; xd=domlen-szth1; xe=domlen+szxt
    ya=-domlen; yb=0; yc=0; yd=domlen

    fctr=2*pi/wlew; shswle=shs*sqrt(1+0*(fctr*wlea*cos(fctr*(zs(k)-zs(0))))**2)

!----- AEROFOIL SURFACE GRID POINTS

 if(rs==0) then
    tmp=(xc-xb)/lnaca
 do n=1,2; do i=0,lnaca
    xnaca(i,n)=i*tmp+xb; ynaca(i,n)=0
 end do; end do
 else
    tmp=xc-xb
    open(8,file='aerofoil.dat',shared)
 do n=1,2; do i=0,lnaca
    read(8,*) xnaca(i,n),ynaca(i,n)
    xnaca(i,n)=tmp*xnaca(i,n)+xb; ynaca(i,n)=tmp*ynaca(i,n)
 end do; end do
    close(8)
 end if

 do n=1,2
    yp(lxis1,n)=yb
 if(rs==0) then
    ll=0; xo=xb; sho=shswle
 else
    ll=8 ! "ll" must be equal to or larger than 4.
    xp(lxis1,n)=xb
 do i=lxis1+1,lxis1+ll
    xp(i,n)=xp(i-1,n)+half*shswle; err=1
 do while(abs(err)>sml)
    yp(i,n)=ylagi(i,n)
    err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shswle-1; xp(i,n)=xp(i,n)-half*err*shswle
 end do
 end do
    xo=xp(lxis1+ll,n); sho=sum(xp(lxis1+ll-4:lxis1+ll,n)*(/3,-16,36,-48,25/))/12
 end if
    ip=lxis1+ll; im=lxi1-ll; call gridf(xp(:,n),pxi,xo,xc,sho,she,lxit,im,ip)
 do i=lxis1+ll+1,lxie1-1
    yp(i,n)=ylagi(i,n)
 end do
    yp(lxie1,n)=yc
 end do

!----- HORIZONTAL INTERFACE POINTS IN XI-DIRECTION

    n=1
    sho=tla/litr; ll=lxi0-2*(lxi0*sho-(-half-xa))/(sho-shswle+sml); xo=xa+ll*sho
    ip=0; im=ll; call gridf(xp(:,n),pxi,xa,xo,sho,sho,lxit,im,ip)
 if(ll/=lxi0) then
    ip=ll; im=lxi0-ll; call gridf(xp(:,n),pxi,xo,xb,sho,shswle,lxit,im,ip)
 end if
 if(k==0) then
    lxisz=lxi2*(minloc(abs(xa+szth1-xp(0:lxi0,n)),1)-1)/lxi0; lp=ll
 end if
    ip=lxis2; im=lxi2-lxisz; call gridf(xp(:,n),pxi,xc,xd,she,free,lxit,im,ip)
    ip=ip+im; im=lxisz; call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
 do m=1,2
 select case(m)
 case(1); is=0; ie=lxie0; ii=is; xo=xb-xa; yo=yb; case(2); is=lxis2; ie=lxit; ii=ie; xo=xd-xc; yo=yc
 end select
    am=yo/sin(half*pi)**2
 do i=is,ie
    gf=(-1)**(m+1)*(xp(i,n)-xp(ii,n))/xo; yp(i,n)=am*sin(half*pi*gf)**2
 end do
 do i=is,ie
    xp(i,n+1)=xp(i,n); yp(i,n+1)=yp(i,n)
 end do
 end do

!----- TOP & BOTTOM BOUNDARY POINTS IN XI-DIRECTION

    n=0
 if(rs==0) then
 do i=0,lxit
    xp(i,n)=xp(i,n+1)
 end do
 else
    tmpa=xb-rs; tmpb=xc+re
    ip=0; im=lxi0; call gridf(xp(:,n),pxi,xa,tmpa,sho,2*shs,lxit,im,ip)
    ip=lxis1; im=lxi1; call gridf(xp(:,n),pxi,tmpa,tmpb,2*shs,she,lxit,im,ip)
    ip=lxis2; im=lxi2-lxisz; call gridf(xp(:,n),pxi,tmpb,xd,she,free,lxit,im,ip)
    ip=ip+im; im=lxisz; call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
 end if
 do i=0,lxit
    yp(i,n)=ya; xp(i,n+3)=xp(i,n); yp(i,n+3)=yd
 end do

!----- VERTICAL INTERFACE POINTS IN ETA-DIRECTION

    yo=tlb*(2.5-cutlb); ll=yo/(half*(shs+sho))
 do m=1,2
    if(rs==0) then; fctr=1; else; fctr=sqrt(m*half); end if
    jp=lets1; jm=ll; call gridf(yq(:,m),qet,zero,yo,fctr*shs,sho,lett,jm,jp)
    jp=jp+jm; jm=let0-ll; call gridf(yq(:,m),qet,yo,yd,sho,free,lett,jm,jp)
 do j=0,lete0
    yq(j,m)=-yq(lett-j,m)
 end do
 end do
 if(rs==0) then
 do j=0,lett
    xq(j,1)=xb; xq(j,2)=xc
 end do
 else
 do n=1,2
 select case(n); case(1); js=0; je=lete0; nn=0; case(2); js=lets1; je=lett; nn=3; end select
 do m=1,2
 select case(m)
 case(1); i=lxis1
    tmpa=sqrt((xp(i,n)-xp(i-2,n))**2+(yp(i,n)-yp(i-2,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+1,n))**2+(yp(i,n)-yp(i+1,n))**2)
    tmp=pi-half*acos(half*(tmpa**2+tmpb**2-(xp(i+1,n)-xp(i-2,n))**2-(yp(i+1,n)-yp(i-2,n))**2)/(tmpa*tmpb))
 case(2); i=lxie1
    tmpa=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+2,n))**2+(yp(i,n)-yp(i+2,n))**2)
    tmp=half*acos(half*(tmpa**2+tmpb**2-(xp(i+2,n)-xp(i-1,n))**2-(yp(i+2,n)-yp(i-1,n))**2)/(tmpa*tmpb))
 end select
    res=tan(tmp); tmpa=abs(xp(i,nn)-xp(i,n)); tmpb=abs(yp(i,nn)-yp(i,n))
    qq=(half*res*tmpa)**2/abs(tmpb-res*tmpa); pp=2*sqrt(qq)/res
 do j=js,je
    xq(j,m)=(-1)**m*pp*(sqrt(abs(yq(j,m)-yp(i,n))+qq)-sqrt(qq))+xp(i,n)
 end do
 end do
 end do
 end if

!----- LEFT & RIGHT BOUNDARY POINTS IN ETA-DIRECTION

    yo=tlb*(2.5-cutlb); fctr=0.8; ll=yo/(half*(fctr*sho+sho))
 do m=0,3,3
    jp=lets1; jm=ll; call gridf(yq(:,m),qet,zero,yo,fctr*sho,sho,lett,jm,jp)
    jp=jp+jm; jm=let0-ll; call gridf(yq(:,m),qet,yo,yd,sho,free,lett,jm,jp)
 do j=0,lete0
    yq(j,m)=-yq(lett-j,m)
 end do
 end do
 do j=0,lett
    xq(j,0)=xa; xq(j,3)=xd
 end do
 if(k==0) then
    lq=ll
 end if

!----- GRID OUTPUT

 do n=0,2,2
 select case(n); case(0); js=0; je=lete0; case(2); js=lets1; je=lett; end select
 do m=0,2
 select case(m)
 case(0); is=0; ie=lxie0; nn=1; tmpa=0; tmpb=int(rs/abs(rs-sml))
 case(1); is=lxis1; ie=lxie1; nn=0; tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
 case(2); is=lxis2; ie=lxit; nn=1; tmpa=int(re/abs(re-sml)); tmpb=0
 end select
    gf=1
 do j=js,je; do i=is,ie
    pp=real(i-is,kind=dp)/(ie-is); qq=real(j-js,kind=dp)/(je-js)
    tmp=sin(half*pi*pp); ra0=((2-n)*real(je-j,kind=dp)/(je-js)+n*qq)/2; ra1=gf+(2-gf)*ra0**2
    pxi(i)=(1-nn)*tmp**ra1+nn*tmp**2
    ts=tmpa*tmpb*(1-pxi(i))+1-tmpb; te=1-ts
    xx(i,j)=(xp(i,n+1)-xp(i,n))*(ts*(xq(j,m)-xq(js,m)+qq*(1-tmpa))/(xq(je,m)-xq(js,m)+1-tmpa)&
           +te*(xq(j,m+1)-xq(js,m+1)+qq*(1-tmpb))/(xq(je,m+1)-xq(js,m+1)+1-tmpb))+xp(i,n)
    ts=1-pxi(i); te=1-ts
    yy(i,j)=(yp(i,n+1)-yp(i,n))*(ts*(yq(j,m)-yq(js,m))/(yq(je,m)-yq(js,m))&
           +te*(yq(j,m+1)-yq(js,m+1))/(yq(je,m+1)-yq(js,m+1)))+yp(i,n)
    zz(i,j)=zs(k)
 end do; end do
 end do
 end do
 select case(mb); case(0,1,2); js=0; je=lete0; case(3,4,5); js=lets1; je=lett; end select
 select case(mb)
 case(0,3); is=0; ie=lxie0; case(1,4); is=lxis1; ie=lxie1; case(2,5); is=lxis2; ie=lxit
 end select
    np=8*(je-js+1)*(ie-is+1)
    write(9,pos=np*k+1) ((xx(i,j),i=is,ie),j=js,je)
    write(9,pos=np*(k+lze0+1)+1) ((yy(i,j),i=is,ie),j=js,je)
    write(9,pos=np*(k+2*lze0+2)+1) ((zz(i,j),i=is,ie),j=js,je)
 end do

    close(9)

!----- GRID CHECKING

 if(ngridv==1) then
    open(8,file='misc/gridview'//cnzone//'.dat')
    write(8,*) 'variables=v1,v2'; write(8,"('zone i=',i4,' j=',i4)") ie-is+1,je-js+1
 do j=js,je; do i=is,ie
    write(8,'(2es15.7)') xx(i,j),yy(i,j)
 end do; end do
    close(8)
 else
    open(8,file='misc/gridview'//cnzone//'.dat'); close(8,status='delete')
 end if
 if(myid==0) then
    write(*,"('Grid generation is complete.')")
    write(*,"('Number of cells to resolve inflow turbulence in x:',i4)") lp
    write(*,"('Number of cells to resolve inflow turbulence in y:',i4)") lq
    write(*,"('Number of cells across sponge:',i4)") minloc(abs(xa+szth1-xp(:,1)),1)-1
 end if

 end if
    deallocate(xx,yy,zz,xp,yp,xq,yq,pxi,qet)

 end subroutine gridaerofoil

!===== LAGRANGIAN INTERPOLATION FOR Y-COORDINATES

 function ylagi(i,n) result(yi)

 integer(kind=int32),intent(in) :: i,n
 real(kind=dp) :: yi

    is=0; ie=lnaca; ii=minloc(abs(xnaca(:,n)-xp(i,n)),1)-1; ip=ii
    if(ii-is<=1) then; ip=is+2; end if
    if(ie-ii<=1) then; ip=ie-2; end if
    yi=0; alag(:)=xp(i,n)-xnaca(ip-2:ip+2,n)
 do jj=-2,2
    blag(:)=xnaca(ip+jj,n)-xnaca(ip-2:ip+2,n); ao=1; bo=1
 do ii=-2,2; if(ii/=jj) then
    ao=ao*alag(ii); bo=bo*blag(ii)
 end if; end do
    yi=yi+ao*ynaca(ip+jj,n)/bo
 end do

 end function ylagi

!=====

 end module gridgen

!*****
