      subroutine bonder
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,lb(20,20)
      character*4 mnam,munit,ban,base*1
      integer*4 break,spline
      dimension rn(10,10),rx(10,10),id(55)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/lsf/bref(10,3,5),th1,th2,dis,rs2,ibref(5),iequ(9),
     1 ibd(20,2),ban(9),base(9)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      data id/1,                        0,
     1       11, 0,      0, 2, 3, 4, 7, 0,
     1       12,16,      0, 0, 5, 6, 8, 0,
     1       13,17,10*0, 0, 0, 0, 0, 9, 0,
     1       14, 0,10*0, 0, 0, 0, 0,10, 0,
     1       15/
      data lb/
     1 .false.,.true.,.true.,.true.,.true.,.true.,.false.,.false.,
     1 .false.,.false.,10*.false.,
     1 .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,
     1 .true.,.true.,10*.false.,
     1 .false.,.false.,.false.,.true.,.false.,.false.,.false.,.false.,
     1 .false.,.false.,10*.false.,
     1 .false.,.false.,.false.,.false.,.true.,.true.,.false.,.false.,
     1 .false.,.false.,10*.false.,
     1 .false.,.false.,.false.,.false.,.false.,.true.,.false.,.false.,
     1 .false.,.false.,10*.false.,
     1 .false.,.false.,.false.,.false.,.false.,.true.,.false.,.false.,
     1 .false.,.false.,10*.false.,280*.false./
      data rn/
c         h    c    n    o    p    s    f    cl   br
     1   0.00,0.98,0.94,0.90,1.38,1.28,0.00,0.00,0.00,0.00,
     1   0.00,1.25,1.20,1.10,0.00,1.66,1.28,1.65,1.80,2.00,
     1   0.00,0.00,0.00,1.15,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,1.35,1.40,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,1.81,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,2.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/
      data rx/
c         h    c    n    o    p    s    f    cl   br
     1   0.00,1.20,1.16,1.12,1.50,1.38,0.00,0.00,0.00,0.00,
     1   0.00,1.68,1.60,2.58,0.00,1.87,1.42,1.83,1.99,2.19,
     1   0.00,0.00,0.00,1.41,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,2.10,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,1.91,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,2.10,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/
c-----------------------------------------------------------setup data
      do i=1,20
      do j=1,i
      lb(j,i)=lb(i,j)
      enddo
      enddo
      do i=1,10
      do j=1,i
      rn(i,j)=rn(i,j)**2
      rn(j,i)=rn(i,j)
      rx(i,j)=rx(i,j)**2
      rx(j,i)=rx(i,j)
      enddo
      enddo
c--------------------------------------------------find chemical bonds
      do i=1,kam-1
      m=id(imch(i))
      if(m.eq.0) goto 10
      x0=corm(i,1)
      y0=corm(i,2)
      z0=corm(i,3)
      do j=i+1,kam
      n=id(imch(j))
      if(n.eq.0) goto 20
      if(lb(m,n)) then
      rmin=rn(m,n)
      rmax=rx(m,n)
      dx=corm(j,1)-x0
      dy=corm(j,2)-y0
      dz=corm(j,3)-z0
      r2=dx*dx+dy*dy+dz*dz
      if(r2.ge.rmin.and.r2.le.rmax) then
      ki=matd(i,7)
      kj=matd(j,7)
      if(ki.lt.6.and.kj.lt.6) then
      ki=ki+1
      matd(i,7)=ki
      matd(i,ki)=j
      kj=kj+1
      matd(j,7)=kj
      matd(j,kj)=i
      endif
      endif
      endif
 20   enddo
 10   enddo
c--------------------------------------------specially requested bonds
      if(ibond.gt.0) then
      do k=1,ibond
      i1=ibd(k,1)
      i2=ibd(k,2)
      l=matd(i1,7)
      l=l+1
      matd(i1,l)=i2
      matd(i1,7)=l
      l=matd(i2,7)
      l=l+1
      matd(i2,l)=i1
      matd(i2,7)=l
      enddo
      endif
      return
      end
