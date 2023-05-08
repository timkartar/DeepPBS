      subroutine up
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      integer*4 break,spline
      character*1 na
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      iup=ien
      ido=ist
      if(line) iup=1
5     do i=ido,iup
      nch=nl
      if(comb) nch=iact(i)
      xdi=hel(i,1,nch)
      ydi=hel(i,2,nch)
      cln=hel(i,4,nch)
      tip=hel(i,5,nch)
      ca=cos(cdr*(-tip))
      sa=sin(cdr*(-tip))
      rx=rex(i,2,nch)
      ry=rey(i,2,nch)
      rz=rez(i,2,nch)
      xx=rex(i,3,nch)
      yy=rey(i,3,nch)
      zz=rez(i,3,nch)
      tx=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1  (rx*rz*(1-ca)+ry*sa)*zz
      ty=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1  (ry*rz*(1-ca)-rx*sa)*zz
      tz=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1  (rz*rz+(1-rz*rz)*ca)*zz
      dx=ry*tz-rz*ty
      dy=rz*tx-rx*tz
      dz=rx*ty-ry*tx
      ca=cos(cdr*(-cln))
      sa=sin(cdr*(-cln))
      rx=dx
      ry=dy
      rz=dz
      xx=tx
      yy=ty
      zz=tz
      ux(i)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1     (rx*rz*(1-ca)+ry*sa)*zz
      uy(i)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1     (ry*rz*(1-ca)-rx*sa)*zz
      uz(i)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1     (rz*rz+(1-rz*rz)*ca)*zz
      wx=uy(i)*dz-uz(i)*dy
      wy=uz(i)*dx-ux(i)*dz
      wz=ux(i)*dy-uy(i)*dx
      hx(i)=rex(i,4,nch)-dx*xdi-wx*ydi
      hy(i)=rey(i,4,nch)-dy*xdi-wy*ydi
      hz(i)=rez(i,4,nch)-dz*xdi-wz*ydi
      ox(i)=hx(i)
      oy(i)=hy(i)
      oz(i)=hz(i)
      if(nch.ne.nl) then
      ux(i)=-ux(i)
      uy(i)=-uy(i)
      uz(i)=-uz(i)
      endif
c-----------------------------------------------------combined o point
      if(comb) then
      dot=0.
      m=1
      do k=1,ns
      if(k.ne.nch.and.(li(i,k).gt.0.or.(li(i,nch).eq.-1
     1 .and.li(i,k).eq.-1))) then
      m=m+1
      x=rex(i,4,k)-hx(i)
      y=rey(i,4,k)-hy(i)
      z=rez(i,4,k)-hz(i)
      dot=dot+ux(i)*x+uy(i)*y+uz(i)*z
      endif
      enddo
      dot=dot/m
      ox(i)=hx(i)+ux(i)*dot
      oy(i)=hy(i)+uy(i)*dot
      oz(i)=hz(i)+uz(i)*dot
      endif
      do k=1,ns
      no=nl
      if(comb) no=k
      if(li(i,no).gt.0.or.(li(i,nch).eq.-1
     1 .and.li(i,k).eq.-1)) then
      bx(i,no)=ox(i)-rex(i,4,no)
      by(i,no)=oy(i)-rey(i,4,no)
      bz(i,no)=oz(i)-rez(i,4,no)
      endif
      enddo
c---------------------------------------------------------------------
      if(i.gt.ido) then
      sx(i)=ox(i)-ox(i-1)
      sy(i)=oy(i)-oy(i-1)
      sz(i)=oz(i)-oz(i-1)
      endif
      enddo
c-------------------------------------------------line+break loop back
         if(line.and.ido.eq.ist.and.break.gt.0) then
         ido=break
         iup=break
         goto 5
         endif
c--------------------------------------------------linear helical axis
      if(line) then
      ido=ist
      iup=ien
      if(break.gt.0) iup=break-1
10    do i=ido+1,iup
      nch=nl
      if(comb) nch=iact(i)
      ux(i)=ux(ido)
      uy(i)=uy(ido)
      uz(i)=uz(ido)
      dx=rex(i,4,nch)-hx(ido)
      dy=rey(i,4,nch)-hy(ido)
      dz=rez(i,4,nch)-hz(ido)
      dot=dx*ux(ido)+dy*uy(ido)+dz*uz(ido)
      hx(i)=hx(ido)+ux(ido)*dot
      hy(i)=hy(ido)+uy(ido)*dot
      hz(i)=hz(ido)+uz(ido)*dot
      ox(i)=hx(i)
      oy(i)=hy(i)
      oz(i)=hz(i)
c-----------------------------------------------------combined o point
      if(comb) then
      dot=0.
      m=1
      do k=1,ns
      if(k.ne.nch.and.(li(i,k).gt.0.or.(li(i,nch).eq.-1
     1 .and.li(i,k).eq.-1))) then
      m=m+1
      x=rex(i,4,k)-hx(i)
      y=rey(i,4,k)-hy(i)
      z=rez(i,4,k)-hz(i)
      dot=dot+ux(ido)*x+uy(ido)*y+uz(ido)*z
      endif
      enddo
      dot=dot/m
      ox(i)=hx(i)+ux(ido)*dot
      oy(i)=hy(i)+uy(ido)*dot
      oz(i)=hz(i)+uz(ido)*dot
      endif
      do k=1,ns
      no=nl
      if(comb) no=k
      if(li(i,no).gt.0.or.(li(i,nch).eq.-1.and.li(i,k).eq.-1)) then
      bx(i,no)=ox(i)-rex(i,4,no)
      by(i,no)=oy(i)-rey(i,4,no)
      bz(i,no)=oz(i)-rez(i,4,no)
      endif
      enddo
      sx(i)=ox(i)-ox(i-1)
      sy(i)=oy(i)-oy(i-1)
      sz(i)=oz(i)-oz(i-1)
      enddo
c-------------------------------------------------line+break loop back
         if(line.and.ido.eq.ist.and.break.gt.0) then
         ido=break
         iup=ien
         goto 10
         endif
c---------------------------------------------------------------------
      endif
      return
      end
