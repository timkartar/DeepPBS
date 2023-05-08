      subroutine calc
      include 'curves_data.inc'
      character*1 na
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly
      integer*4 break,spline
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/der/udx(n3,4),udy(n3,4),udz(n3,4),pdx(n3,4),pdy(n3,4),
     1 pdz(n3,4),upx(n3),upy(n3),upz(n3),usx(n3),usy(n3),usz(n3),
     1 umx(n3),umy(n3),umz(n3),qx(n3),qy(n3),qz(n3)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      call up
      scp(1)=0.
      scp(2)=0.
      scp(3)=0.
      scp(4)=0.
      do i=iste,iene
      dif(i)=0.
      enddo
      do k=1,ns
      no=nl
      if(comb) no=k
c--------------------------------------------------------helical terms
      if(.not.dinu) then
      do i=iste+1,iene
      if(i.eq.break.or.li(i-1,no).lt.0.or.li(i,no).lt.0) goto 100
      do j=1,3
      if(i.gt.ist) then
      rx1=rex(i-1,j,no)
      ry1=rey(i-1,j,no)
      rz1=rez(i-1,j,no)
      d1=rx1*ux(i-1)   +ry1*uy(i-1)   +rz1*uz(i-1)
      c1=rx1*bx(i-1,no)+ry1*by(i-1,no)+rz1*bz(i-1,no)
      else
      d1=efd(j,1,no)
      c1=efc(j,1,no)
      endif
      if(i.le.ien) then
      rx2=rex(i,j,no)
      ry2=rey(i,j,no)
      rz2=rez(i,j,no)
      d2=rx2*ux(i)   +ry2*uy(i)   +rz2*uz(i)
      c2=rx2*bx(i,no)+ry2*by(i,no)+rz2*bz(i,no)
      else
      d2=efd(j,2,no)
      c2=efc(j,2,no)
      endif
      qr(i,j,no)=d2-d1
      qp(i,j,no)=c2-c1
      si1=(d2-d1)**2
      dif(i-1)=dif(i-1)+si1*10
      scp(1)=scp(1)+si1
      si2=(c2-c1)**2
      dif(i-1)=dif(i-1)+si2
      scp(2)=scp(2)+si2
      enddo
100   enddo
c---------------------------------------------------------dinucleotide
      else
      do i=iste+2,iene
      if(i.eq.break.or.i-1.eq.break.or.
     1 li(i-2,no).lt.0.or.li(i,no).lt.0) goto 150
      do j=1,3
      if(i.gt.ist+1) then
      rx1=rex(i-2,j,no)
      ry1=rey(i-2,j,no)
      rz1=rez(i-2,j,no)
      d1=rx1*ux(i-2)   +ry1*uy(i-2)   +rz1*uz(i-2)
      c1=rx1*bx(i-2,no)+ry1*by(i-2,no)+rz1*bz(i-2,no)
      else
      d1=efd(j,1,no)
      c1=efc(j,1,no)
      endif
      if(i.le.ien) then
      rx2=rex(i,j,no)
      ry2=rey(i,j,no)
      rz2=rez(i,j,no)
      d2=rx2*ux(i)   +ry2*uy(i)   +rz2*uz(i)
      c2=rx2*bx(i,no)+ry2*by(i,no)+rz2*bz(i,no)
      else
      d2=efd(j,2,no)
      c2=efc(j,2,no)
      endif
      qr(i,j,no)=d2-d1
      qp(i,j,no)=c2-c1
      si1=(d2-d1)**2
      dif(i-2)=dif(i-2)+si1*5
      dif(i-1)=dif(i-1)+si1*5
      scp(1)=scp(1)+si1
      si2=(c2-c1)**2
      dif(i-2)=dif(i-2)+si2/2
      dif(i-1)=dif(i-1)+si2/2
      scp(2)=scp(2)+si2
      enddo
150   enddo
      endif
c-----------------------------------------------------------kink terms
      if(.not.line) then
      do i=ist+1,ien
      if(i.ne.break) then
      upx(i)=ux(i)-ux(i-1)
      upy(i)=uy(i)-uy(i-1)
      upz(i)=uz(i)-uz(i-1)
      qi1=upx(i)**2+upy(i)**2+upz(i)**2
      dif(i-1)=dif(i-1)+qi1*10
      scp(3)=scp(3)+qi1
      usx(i)=ux(i)+ux(i-1)
      usy(i)=uy(i)+uy(i-1)
      usz(i)=uz(i)+uz(i-1)
      dot=usx(i)**2+usy(i)**2+usz(i)**2
      umx(i)=usx(i)/dot
      umy(i)=usy(i)/dot
      umz(i)=usz(i)/dot
      dot=usx(i)*sx(i)+usy(i)*sy(i)+usz(i)*sz(i)
      qx(i)=sx(i)-umx(i)*dot
      qy(i)=sy(i)-umy(i)*dot
      qz(i)=sz(i)-umz(i)*dot
      qi2=qx(i)**2+qy(i)**2+qz(i)**2
      dif(i-1)=dif(i-1)+qi2
      scp(4)=scp(4)+qi2
      endif
      enddo
      endif
      enddo
      sum=10*(scp(1)+scp(3))+scp(2)+scp(4)
      return
      end
