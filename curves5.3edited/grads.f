      subroutine grads
      include 'curves_data.inc'
      character*1 na
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,even
      integer*4 break,spline
      dimension du1(4),du2(4),dv1(4),dv2(4),due1(4),due2(4),
     1 duo1(4),duo2(4)
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
      do m=1,4*n+4
      gra(m)=0.
      enddo
      do k=1,ns
      no=nl
      if(comb) no=k
      do j=1,3
c--------------------------------------------------linear helical axis
      if(line) then
      lof=0
      ido=ist
      iup=ien
      if(break.gt.0) iup=break-1
5     if(.not.dinu) then
      do i=ido+1,iup
      if(li(i-1,no).gt.0.and.li(i,no).gt.0) then
      rx2=rex(i,j,no)
      ry2=rey(i,j,no)
      rz2=rez(i,j,no)
      rx1=rex(i-1,j,no)
      ry1=rey(i-1,j,no)
      rz1=rez(i-1,j,no)
      do l=1,4
      dot1=udx(i-1,l)*rx1+udy(i-1,l)*ry1+udz(i-1,l)*rz1
      dot2=udx(i,l)*rx2+udy(i,l)*ry2+udz(i,l)*rz2
      ds1=20*qr(i,j,no)*(dot2-dot1)
      dot1=pdx(i-1,l)*rx1+pdy(i-1,l)*ry1+pdz(i-1,l)*rz1
      dot2=pdx(i,l)*rx2+pdy(i,l)*ry2+pdz(i,l)*rz2
      ds2=2*qp(i,j,no)*(dot2-dot1)
      gra(l+lof)=gra(l+lof)+ds1+ds2
      enddo
      endif
      enddo
         else
         do i=ido+2,iup
         if(li(i-2,no).gt.0.and.li(i,no).gt.0) then
         rx2=rex(i,j,no)
         ry2=rey(i,j,no)
         rz2=rez(i,j,no)
         rx1=rex(i-2,j,no)
         ry1=rey(i-2,j,no)
         rz1=rez(i-2,j,no)
         do l=1,4
         dot1=udx(i-2,l)*rx1+udy(i-2,l)*ry1+udz(i-2,l)*rz1
         dot2=udx(i,l)*rx2+udy(i,l)*ry2+udz(i,l)*rz2
         ds1=20*qr(i,j,no)*(dot2-dot1)
         dot1=pdx(i-2,l)*rx1+pdy(i-2,l)*ry1+pdz(i-2,l)*rz1
         dot2=pdx(i,l)*rx2+pdy(i,l)*ry2+pdz(i,l)*rz2
         ds2=2*qp(i,j,no)*(dot2-dot1)
         gra(l+lof)=gra(l+lof)+ds1+ds2
         enddo
         endif
      enddo
      endif
c-------------------------------------------------line+break loop back
         if(line.and.ido.eq.ist.and.break.gt.0) then
         ido=break
         iup=ien
         lof=4
         goto 5
         endif
c---------------------------------------------------------------------
      goto 50
      endif
c---------------------------------------------------------------------
      if(.not.dinu) then
      do l=1,4
      du1(l)=0.
      du2(l)=0.
      dv1(l)=0.
      dv2(l)=0.
      enddo
      else
      do l=1,4
      due1(l)=0.
      due2(l)=0.
      duo1(l)=0.
      duo2(l)=0.
      dv1(l)=0.
      dv2(l)=0.
      enddo
      endif
c--------------------------------------------------------helical terms
      if(.not.dinu) then
      do i=ist+1,ien
      il=i-ist+1
      if(i.eq.break.or.li(i-1,no).lt.0.or.li(i,no).lt.0) then
      do l=1,4
      m=(il-2)*4+l
      gra(m)=gra(m)+du1(l)+du2(l)
      du1(l)=0.
      du2(l)=0.
      enddo
      goto 100
      endif
      rx2=rex(i,j,no)
      ry2=rey(i,j,no)
      rz2=rez(i,j,no)
      rx1=rex(i-1,j,no)
      ry1=rey(i-1,j,no)
      rz1=rez(i-1,j,no)
      do l=1,4
      m=(il-2)*4+l
      dot=udx(i-1,l)*rx1+udy(i-1,l)*ry1+udz(i-1,l)*rz1
      ds1=-20*qr(i,j,no)*dot
      if(ends.and.i.eq.ist+1) ds1=ds1+20*qr(1,j,no)*dot
      dot=pdx(i-1,l)*rx1+pdy(i-1,l)*ry1+pdz(i-1,l)*rz1
      ds2=-2*qp(i,j,no)*dot
      if(ends.and.i.eq.ist+1) ds2=ds2+2*qp(1,j,no)*dot
      gra(m)=gra(m)+ds1+ds2+du1(l)+du2(l)
      dot=udx(i,l)*rx2+udy(i,l)*ry2+udz(i,l)*rz2
      du1(l)=20*qr(i,j,no)*dot
      if(ends.and.i.eq.ien) du1(l)=du1(l)-20*qr(ien+1,j,no)*dot
      dot=pdx(i,l)*rx2+pdy(i,l)*ry2+pdz(i,l)*rz2
      du2(l)=2*qp(i,j,no)*dot
      if(ends.and.i.eq.ien) du2(l)=du2(l)-2*qp(ien+1,j,no)*dot
      enddo
100   enddo
      else
      do i=ist+2,ien
      il=i-ist+1
      even=.false.
      if(mod(i,2).eq.0) even=.true.
      if(i.eq.break.or.i-1.eq.break.or.
     1 li(i-2,no).lt.0.or.li(i,no).lt.0) then
      if(even) then
      do l=1,4
      m=(il-3)*4+l
      gra(m)=gra(m)+due1(l)+due2(l)
      due1(l)=0.
      due2(l)=0.
      enddo
      else
      do l=1,4
      m=(il-3)*4+l
      gra(m)=gra(m)+duo1(l)+duo2(l)
      duo1(l)=0.
      duo2(l)=0.
      enddo
      endif
      goto 200
      endif
      rx2=rex(i,j,no)
      ry2=rey(i,j,no)
      rz2=rez(i,j,no)
      rx1=rex(i-2,j,no)
      ry1=rey(i-2,j,no)
      rz1=rez(i-2,j,no)
      do l=1,4
      m=(il-3)*4+l
      dot=udx(i-2,l)*rx1+udy(i-2,l)*ry1+udz(i-2,l)*rz1
      ds1=-20*qr(i,j,no)*dot
      if(ends.and.i.eq.ist+3) ds1=ds1+20*qr(ist+1,j,no)*dot
      dot=pdx(i-2,l)*rx1+pdy(i-2,l)*ry1+pdz(i-2,l)*rz1
      ds2=-2*qp(i,j,no)*dot
      if(ends.and.i.eq.ist+3) ds2=ds2+2*qp(ist+1,j,no)*dot
      if(even) then
      gra(m)=gra(m)+ds1+ds2+due1(l)+due2(l)
      else
      gra(m)=gra(m)+ds1+ds2+duo1(l)+duo2(l)
      endif
      if(even) then
      dot=udx(i,l)*rx2+udy(i,l)*ry2+udz(i,l)*rz2
      due1(l)=20*qr(i,j,no)*dot
      if(ends.and.i.eq.ien-1) due1(l)=due1(l)-20*qr(n+1,j,no)*dot
      dot=pdx(i,l)*rx2+pdy(i,l)*ry2+pdz(i,l)*rz2
      due2(l)=2*qp(i,j,no)*dot
      if(ends.and.i.eq.ien-1) due2(l)=due2(l)-2*qp(n+1,j,no)*dot
      else
      dot=udx(i,l)*rx2+udy(i,l)*ry2+udz(i,l)*rz2
      duo1(l)=20*qr(i,j,no)*dot
      if(ends.and.i.eq.ien-1) duo1(l)=duo1(l)-20*qr(n+1,j,no)*dot
      dot=pdx(i,l)*rx2+pdy(i,l)*ry2+pdz(i,l)*rz2
      duo2(l)=2*qp(i,j,no)*dot
      if(ends.and.i.eq.ien-1) duo2(l)=duo2(l)-2*qp(n+1,j,no)*dot
      endif
      enddo
200   enddo
      endif
c-----------------------------------------------------------kink terms
      if(j.eq.1) then
      do i=ist+1,ien
      il=i-ist+1
      if(i.eq.break) then
c      if(i.eq.break.or.li(i-1,no).lt.0.or.li(i,no).lt.0) then
      do l=1,4
      m=(il-2)*4+l
      gra(m)=gra(m)+dv1(l)+dv2(l)
      dv1(l)=0.
      dv2(l)=0.
      enddo
      goto 250
      endif
      do l=1,4
      m=(il-2)*4+l
      dq1=-20*(upx(i)*udx(i-1,l)+upy(i)*udy(i-1,l)+upz(i)*udz(i-1,l))
      dot1=udx(i-1,l)*sx(i)+udy(i-1,l)*sy(i)+udz(i-1,l)*sz(i)
      dot2=usx(i)*pdx(i-1,l)+usy(i)*pdy(i-1,l)+usz(i)*pdz(i-1,l)
      dot3=umx(i)*sx(i)+umy(i)*sy(i)+umz(i)*sz(i)
      dot4=usx(i)*udx(i-1,l)+usy(i)*udy(i-1,l)+usz(i)*udz(i-1,l)
      dqx=-pdx(i-1,l)-umx(i)*(dot1-dot2)-udx(i-1,l)*dot3+2*
     1 umx(i)*dot3*dot4
      dqy=-pdy(i-1,l)-umy(i)*(dot1-dot2)-udy(i-1,l)*dot3+2*
     1 umy(i)*dot3*dot4
      dqz=-pdz(i-1,l)-umz(i)*(dot1-dot2)-udz(i-1,l)*dot3+2*
     1 umz(i)*dot3*dot4
      dq2=2*(qx(i)*dqx+qy(i)*dqy+qz(i)*dqz)
      gra(m)=gra(m)+dq1+dq2+dv1(l)+dv2(l)
      dv1(l)=20*(upx(i)*udx(i,l)+upy(i)*udy(i,l)+upz(i)*udz(i,l))
      dot1=udx(i,l)*sx(i)+udy(i,l)*sy(i)+udz(i,l)*sz(i)
      dot2=usx(i)*pdx(i,l)+usy(i)*pdy(i,l)+usz(i)*pdz(i,l)
      dot4=usx(i)*udx(i,l)+usy(i)*udy(i,l)+usz(i)*udz(i,l)
      dqx=pdx(i,l)-umx(i)*(dot1+dot2)-udx(i,l)*dot3+2*umx(i)*dot3*dot4
      dqy=pdy(i,l)-umy(i)*(dot1+dot2)-udy(i,l)*dot3+2*umy(i)*dot3*dot4
      dqz=pdz(i,l)-umz(i)*(dot1+dot2)-udz(i,l)*dot3+2*umz(i)*dot3*dot4
      dv2(l)=2*(qx(i)*dqx+qy(i)*dqy+qz(i)*dqz)
      enddo
250   enddo
      endif
c---------------------------------------------------------------------
      il=ien-ist+1
      m=(il-1)*4
      even=.false.
      if(mod(ien,2).eq.0) even=.true.
      if(.not.dinu) then
      do l=1,4
      gra(m+l)=gra(m+l)+du1(l)+du2(l)+dv1(l)+dv2(l)
      enddo
      else
      if(even) then
      do l=1,4
      gra(m-4+l)=gra(m-4+l)+duo1(l)+duo2(l)
      gra(m+l)=gra(m+l)+due1(l)+due2(l)+dv1(l)+dv2(l)
      enddo
      else
      do l=1,4
      gra(m-4+l)=gra(m-4+l)+due1(l)+due2(l)
      gra(m+l)=gra(m+l)+duo1(l)+duo2(l)+dv1(l)+dv2(l)
      enddo
      endif
      endif
50    enddo
      enddo
      return
      end
