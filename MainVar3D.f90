!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS

 integer,parameter :: int16=selected_int_kind(4),int32=selected_int_kind(9),int64=selected_int_kind(18)
 integer,parameter :: ieee32=selected_real_kind(6,37),ieee64=selected_real_kind(15,307)
 integer,parameter :: sp=ieee32,dp=ieee64

 integer(kind=int32),parameter :: ntdrv=0,ntflt=1,nrall=0,nrone=1,n45no=0,n45go=1
 integer(kind=int32),parameter :: lmd=11,lmf=8,lmp=max(lmd,lmf),mfbi=2,mbci=3
 integer(kind=int32),parameter :: liofs=16,liofl=24

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(kind=dp),parameter :: pi=3.141592653589793
 real(kind=dp),parameter :: zero=0,one=1,half=0.5,sqrt2=sqrt(2.0),sqrt2i=1/sqrt2,sml=1.0d-6,free=1.0d+6
 real(kind=dp),parameter :: gam=1.4,gamm1=gam-1,ham=1/gam,hamm1=1/gamm1
 real(kind=dp),parameter :: prndtl=0.71,gamm1prndtli=1/(gamm1*prndtl)

 real(kind=dp),parameter :: alpha=0.5862704032801503
 real(kind=dp),parameter :: beta=0.09549533555017055
 real(kind=dp),parameter :: aa=0.6431406736919156
 real(kind=dp),parameter :: ab=0.2586011023495066
 real(kind=dp),parameter :: ac=0.007140953479797375

 real(kind=dp),parameter :: beta20=0.03250008295108466
 real(kind=dp),parameter :: alpha21=0.3998040493524358
 real(kind=dp),parameter :: alpha23=0.7719261277615860
 real(kind=dp),parameter :: beta24=0.1626635931256900
 real(kind=dp),parameter :: a20=-0.1219006056449124
 real(kind=dp),parameter :: a21=-0.6301651351188667
 real(kind=dp),parameter :: a23=0.6521195063966084
 real(kind=dp),parameter :: a24=0.3938843551210350
 real(kind=dp),parameter :: a25=0.01904944407973912
 real(kind=dp),parameter :: a26=-0.001027260523947668

 real(kind=dp),parameter :: alpha10=0.08360703307833438
 real(kind=dp),parameter :: alpha12=2.058102869495757
 real(kind=dp),parameter :: beta13=0.9704052014790193
 real(kind=dp),parameter :: a10=-0.3177447290722621
 real(kind=dp),parameter :: a12=-0.02807631929593225
 real(kind=dp),parameter :: a13=1.593461635747659
 real(kind=dp),parameter :: a14=0.2533027046976367
 real(kind=dp),parameter :: a15=-0.03619652460174756
 real(kind=dp),parameter :: a16=0.004080281419108407

 real(kind=dp),parameter :: alpha01=5.912678614078549
 real(kind=dp),parameter :: beta02=3.775623951744012
 real(kind=dp),parameter :: a01=-3.456878182643609
 real(kind=dp),parameter :: a02=5.839043358834730
 real(kind=dp),parameter :: a03=1.015886726041007
 real(kind=dp),parameter :: a04=-0.2246526470654333
 real(kind=dp),parameter :: a05=0.08564940889936562
 real(kind=dp),parameter :: a06=-0.01836710059356763

 real(kind=dp),parameter,dimension(3) :: fex=(/45,-9,1/)/(mfbi*30.0)

!===== ALLOCATABLE MAIN ARRAYS

 integer(kind=int32),dimension(:,:),allocatable :: npc,lio
 integer(kind=int32),dimension(:),allocatable :: li,lcsz,lctz,mxc
 integer(kind=int32),dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer(kind=int32),dimension(:),allocatable :: lximb,letmb,lzemb,mo
 integer(kind=int64),dimension(:),allocatable :: lhmb

 real(kind=dp),dimension(:,:),allocatable :: qo,qa,qb,de
 real(kind=dp),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(kind=dp),dimension(:,:),allocatable :: xim,etm,zem
 real(kind=dp),dimension(:),allocatable :: p,yaco
 real(kind=dp),dimension(:),allocatable :: sbcc

 real(kind=dp),dimension(:,:),allocatable :: rr,ss

 real(kind=dp),dimension(:,:),allocatable :: xu,yu
 real(kind=dp),dimension(:,:),allocatable :: xl,yl
 real(kind=dp),dimension(:),allocatable :: sa,sb

 real(kind=dp),dimension(:,:),allocatable :: ran,sit,ait
 real(kind=dp),dimension(:),allocatable :: xit,yit,zit
 real(kind=dp),dimension(:),allocatable :: asz,bsz,csz
 real(kind=dp),dimension(:),allocatable :: times,vmpi

 real(kind=dp),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm

 real(kind=dp),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
 real(kind=dp),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
 real(kind=dp),dimension(:,:,:),allocatable,target :: send1,send2,send3
 real(kind=dp),dimension(:,:,:),allocatable,target :: recv1,recv2,recv3
 real(kind=dp),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

 real(kind=dp),dimension(:,:),allocatable :: varm
 real(kind=dp),dimension(:),allocatable :: varmin,varmax
 real(kind=sp),dimension(:),allocatable :: varr,vart,vara,varb,vmean

 character(13),dimension(:),allocatable :: ctecplt,cthead
 character(4),dimension(:),allocatable :: cfilet
 character(4),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

 integer(kind=int32),dimension(0:1,0:1,3) :: ndf
 integer(kind=int32),dimension(3,3) :: ijk
 integer(kind=int32),dimension(0:1,3) :: nbc,ncd,nsz
 integer(kind=int32),dimension(0:4) :: no
 integer(kind=int32),dimension(-2:2) :: ilag
 integer(kind=int32),dimension(3) :: nbsize,nbcs,nbce,ncds,ncde,ms,me
 integer(kind=int32) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,nrec
 integer(kind=int32) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,kp,l,lh,ll,lp,lq,ltomb
 integer(kind=int32) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nt,nz,ndati,nsigi,nout,nfile
 integer(kind=int32) :: nts,nscrn,nsgnl,ndata,nviscous,nkrk,nsmf,nfskp,nrestart
 integer(kind=int64) :: nlmx,llmb,llmo,lis,lie,ljs,lje

 real(kind=dp),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(kind=dp),dimension(-2:2,0:2,0:1) :: albed,albef
 real(kind=dp),dimension(mbci,mbci) :: cbca,cbcs
 real(kind=dp),dimension(5,5) :: xt
 real(kind=dp),dimension(0:5,0:2) :: abc
 real(kind=dp),dimension(0:1,0:1) :: pbcot
 real(kind=dp),dimension(0:lmp) :: sap
 real(kind=dp),dimension(mbci) :: rbci,sbci
 real(kind=dp),dimension(5) :: cha,dha
 real(kind=dp),dimension(-2:2) :: alag,blag,tlag
 real(kind=dp),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
 real(kind=dp),dimension(0:2) :: fam,fbm,fcm
 real(kind=dp) :: alphf,betf,fa,fb,fc
 real(kind=dp) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
 real(kind=dp) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltkbc,dto
 real(kind=dp) :: rhooo,poo,aoo,amachoo,srefoo,srefp1dre
 real(kind=dp) :: dt,dts,dte,dtk,dtko,dtsum,timo,tsam,wts,wte,wtime
 real(kind=dp) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

 character(1),dimension(0:4) :: cno
 character(3) :: cnzone,cndata
 character(5) :: cnnode
 character(7) :: czone
 character(13) :: coutput
 character(16) :: cinput
 character(16) :: cgrid
 character(18) :: cdata,cturb
 character(19) :: crestart

!===== INTEGER VARIABLES FOR MPI COMMANDS

 integer(kind=int32),dimension(:,:),allocatable :: ista
 integer(kind=int32),dimension(:),allocatable :: ireq
 integer(kind=int32) :: ir,mpro,npro,myid,itag,info,icom,ierr

!=====

 end module mainvar3d

!*****
