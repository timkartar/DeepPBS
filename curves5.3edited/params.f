      subroutine params
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,kcopy
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,na*1
      integer*4 break,spline
      real*8 nx,ny,nz
      dimension ulx(0:n3,4),uly(0:n3,4),ulz(0:n3,4),
     1          px(0:n3,4),py(0:n3,4),pz(0:n3,4),
     1          wx(0:n3,4),wy(0:n3,4),wz(0:n3,4),
     1          vx(0:n3,4),vy(0:n3,4),vz(0:n3,4),
     1          vax(0:n3) ,vay(0:n3) ,vaz(0:n3)
      common/cha/file,lis,dna,axin,axout,daf,pdb,mcode
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/hol/uho(0:n3,3,4),hho(0:n3,3,4),vkin(0:n3,7,4),bend(n3,4),
     1 hold(0:n3,2,4),vold(0:n3,4),pal(0:n3,6,4),pab(0:n3,6,4),inv(4)
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      kcopy=.false.
      if(comb) then
      do i=iste,iene
      if(iact(i).ne.1) kcopy=.true.
      enddo
      endif
c-------------------------------------------set axis system directions
      do k=1,max0(ns,nst)
      id=idr(k)
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=ist-1
      iene=ien+1
      endif
      endif
      inv(k)=1
      rise=0.
      do i=iste+1,iene
      xu=uho(i,1,k)+uho(i-1,1,k)
      yu=uho(i,2,k)+uho(i-1,2,k)
      zu=uho(i,3,k)+uho(i-1,3,k)
      xp=hho(i,1,k)-hho(i-1,1,k)
      yp=hho(i,2,k)-hho(i-1,2,k)
      zp=hho(i,3,k)-hho(i-1,3,k)
      rise=rise+sign(1.d0,xu*xp+yu*yp+zu*zp)
      enddo
      if(id.eq. 1.and.rise.lt.0) inv(k)=-1
      if(id.eq.-1.and.rise.gt.0) inv(k)=-1
      enddo
c---------------------------------------------------------------------
      do k=1,max0(ns,nst)
      if(comb) then
      is=1
      else
      is=k
      ist=ng(k)
      ien=nr(k)
      iste=ist
      iene=ien
      if(ends) then
      iste=ist-1
      iene=ien+1
      endif
      endif
      hel(ist,3,k)=0.
      hel(ist,6,k)=0.
c====================================================local axis system
      do i=iste,iene
      if(comb.and.k.gt.1) then
      ulx(i,k)=-ulx(i,1)
      uly(i,k)=-uly(i,1)
      ulz(i,k)=-ulz(i,1)
      else
      ulx(i,k)=uho(i,1,is)*inv(k)
      uly(i,k)=uho(i,2,is)*inv(k)
      ulz(i,k)=uho(i,3,is)*inv(k)
      endif
      if(li(i,k).ge.-1) then
      x=rex(i,4,k)-hho(i,1,is)
      y=rey(i,4,k)-hho(i,2,is)
      z=rez(i,4,k)-hho(i,3,is)
      dot=x*ulx(i,k)+y*uly(i,k)+z*ulz(i,k)
      px(i,k)=hho(i,1,is)+ulx(i,k)*dot
      py(i,k)=hho(i,2,is)+uly(i,k)*dot
      pz(i,k)=hho(i,3,is)+ulz(i,k)*dot
      dax=rex(i,2,k)
      day=rey(i,2,k)
      daz=rez(i,2,k)
      dot=ulx(i,k)*dax+uly(i,k)*day+ulz(i,k)*daz
      x=dax-ulx(i,k)*dot
      y=day-uly(i,k)*dot
      z=daz-ulz(i,k)*dot
      r=sqrt(x*x+y*y+z*z)
      wx(i,k)=x/r
      wy(i,k)=y/r
      wz(i,k)=z/r
         if(i.gt.iste.and.li(i,k).eq.1.and.li(i-1,k).eq.1) then
         dit=wx(i-1,k)*wx(i,k)+wy(i-1,k)*wy(i,k)+wz(i-1,k)*wz(i,k)
         if(dit.lt.0.) then
         wx(i,k)=-wx(i,k)
         wy(i,k)=-wy(i,k)
         wz(i,k)=-wz(i,k)
         endif
         endif
      vx(i,k)=wy(i,k)*ulz(i,k)-wz(i,k)*uly(i,k)
      vy(i,k)=wz(i,k)*ulx(i,k)-wx(i,k)*ulz(i,k)
      vz(i,k)=wx(i,k)*uly(i,k)-wy(i,k)*ulx(i,k)
      kc=1
      if(li(i,k).lt.-1) kc=iact(i)
      if(comb.and.kc.gt.1.and.kc.eq.k) then
      dot=ulx(i,1)*ulx(i,kc)+uly(i,1)*uly(i,kc)+ulz(i,1)*ulz(i,kc)
      vx(i,1)=vx(i,k)*sign(1.d0,-dot)
      vy(i,1)=vy(i,k)*sign(1.d0,-dot)
      vz(i,1)=vz(i,k)*sign(1.d0,-dot)
      wx(i,1)=wx(i,k)*sign(1.d0,-dot)
      wy(i,1)=wy(i,k)*sign(1.d0,-dot)
      wz(i,1)=wz(i,k)*sign(1.d0,-dot)
      endif
      ex=rex(i,4,k)-px(i,k)
      ey=rey(i,4,k)-py(i,k)
      ez=rez(i,4,k)-pz(i,k)
      fx=rex(i,3,k)
      fy=rey(i,3,k)
      fz=rez(i,3,k)
c------------------------------------------------------xdi,ydi,cln,tip
      hel(i,1,k)=ex*vx(i,k)+ey*vy(i,k)+ez*vz(i,k)
      hel(i,2,k)=ex*wx(i,k)+ey*wy(i,k)+ez*wz(i,k)
      dot=wx(i,k)*dax+wy(i,k)*day+wz(i,k)*daz
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      hel(i,4,k)=acos(dot)*crd
      tx=wy(i,k)*daz-wz(i,k)*day
      ty=wz(i,k)*dax-wx(i,k)*daz
      tz=wx(i,k)*day-wy(i,k)*dax
      if(tx*vx(i,k)+ty*vy(i,k)+tz*vz(i,k).lt.0) hel(i,4,k)=-hel(i,4,k)
      qx=vy(i,k)*daz-vz(i,k)*day
      qy=vz(i,k)*dax-vx(i,k)*daz
      qz=vx(i,k)*day-vy(i,k)*dax
      rq=sqrt(qx*qx+qy*qy+qz*qz)
      dot=(qx*fx+qy*fy+qz*fz)/rq
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      hel(i,5,k)=acos(dot)*crd
      tx=qy*fz-qz*fy
      ty=qz*fx-qx*fz
      tz=qx*fy-qy*fx
      if(tx*dax+ty*day+tz*daz.lt.0) hel(i,5,k)=-hel(i,5,k)
      endif
      enddo
c=========================================================overall bend
      if(comb.and.k.gt.1) goto 150
      do ic=iste+1,iene-1
      if(.not.comb) then
      dx=hx(iene)-hx(iste)
      dy=hy(iene)-hy(iste)
      dz=hz(iene)-hz(iste)
      tx=hx(ic)-hx(iste)
      ty=hy(ic)-hy(iste)
      tz=hz(ic)-hz(iste)
         else
         dx=ox(iene)-ox(iste)
         dy=oy(iene)-oy(iste)
         dz=oz(iene)-oz(iste)
         tx=ox(ic)-ox(iste)
         ty=oy(ic)-oy(iste)
         tz=oz(ic)-oz(iste)
         endif
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      dx=dx/rd
      dy=dy/rd
      dz=dz/rd
      rx=ulx(ic,k)
      ry=uly(ic,k)
      rz=ulz(ic,k)
      drd=rx*dx+ry*dy+rz*dz
      drt=rx*tx+ry*ty+rz*tz
      cx=dx*drt/drd-tx
      cy=dy*drt/drd-ty
      cz=dz*drt/drd-tz
      rc=sqrt(cx*cx+cy*cy+cz*cz)
      xv=vx(ic,k)
      yv=vy(ic,k)
      zv=vz(ic,k)
      dot=(cx*xv+cy*yv+cz*zv)/rc
      if(abs(dot).gt.1.0) dot=sign(1.d0,dot)
      bend(ic,k)=acos(dot)*crd
      xn=yv*cz-zv*cy
      yn=zv*cx-xv*cz
      zn=xv*cy-yv*cx
      if(rx*xn+ry*yn+rz*zn.lt.0) bend(ic,k)=-bend(ic,k)
      enddo
      dot=ulx(iste,k)*ulx(iene,k)+uly(iste,k)*uly(iene,k)
     & +ulz(iste,k)*ulz(iene,k)
      if(abs(dot).gt.1.0) dot=sign(1.d0,dot)
      bend(iste,k)=acos(dot)*crd
      if(.not.comb) then
      dx=hx(iene)-hx(iene-1)
      dy=hy(iene)-hy(iste-1)
      dz=hz(iene)-hz(iste-1)
      tx=hx(iste+1)-hx(iste)
      ty=hy(iste+1)-hy(iste)
      tz=hz(iste+1)-hz(iste)
         else
         dx=ox(iene)-ox(iene-1)
         dy=oy(iene)-oy(iene-1)
         dz=oz(iene)-oz(iene-1)
         tx=ox(iste+1)-ox(iste)
         ty=oy(iste+1)-oy(iste)
         tz=oz(iste+1)-oz(iste)
         endif
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      rt=sqrt(tx*tx+ty*ty+tz*tz)
      dot=(dx*tx+dy*ty+dz*tz)/(rd*rt)
      if(abs(dot).gt.1.0) dot=sign(1.d0,dot)
      bend(iene,k)=acos(dot)*crd
150   enddo
c-------------------------------------------------------mean v vectors
      if(comb) then
      do i=iste,iene
      x=0.
      y=0.
      z=0.
      do k=1,ns
      if(li(i,k).ge.-1) then
      if(k.eq.1) then
      dot=-1.
      else
      dot=ulx(i,1)*ulx(i,k)+uly(i,1)*uly(i,k)+ulz(i,1)*ulz(i,k)
      endif
      x=x+vx(i,k)*sign(1.d0,-dot)
      y=y+vy(i,k)
      z=z+vz(i,k)
      endif
      enddo
      r=sqrt(x*x+y*y+z*z)
      vax(i)=x/r
      vay(i)=y/r
      vaz(i)=z/r
      enddo
      endif
c================================================================kinks
      do k=1,max0(ns,nst)
      if(.not.comb) then
      ist=ng(k)
      ien=nr(k)
      iste=ist
      iene=ien
      endif
      if(ends) then
      iste=ist-1
      iene=ien+1
      endif
      do i=iste+1,iene
c-----------------------------------------------mean plane axis system
      if(.not.comb) then
      nx=ulx(i-1,k)+ulx(i,k)
      ny=uly(i-1,k)+uly(i,k)
      nz=ulz(i-1,k)+ulz(i,k)
      rn=sqrt(nx*nx+ny*ny+nz*nz)
      nx=nx/rn
      ny=ny/rn
      nz=nz/rn
      x=px(i-1,k)-px(i,k)
      y=py(i-1,k)-py(i,k)
      z=pz(i-1,k)-pz(i,k)
      vkin(i,7,k)=sqrt(x*x+y*y+z*z)
      qx=(px(i-1,k)+px(i,k))/2.
      qy=(py(i-1,k)+py(i,k))/2.
      qz=(pz(i-1,k)+pz(i,k))/2.
      x=vx(i-1,k)+vx(i,k)
      y=vy(i-1,k)+vy(i,k)
      z=vz(i-1,k)+vz(i,k)
      else
      nx=ulx(i-1,1)+ulx(i,1)
      ny=uly(i-1,1)+uly(i,1)
      nz=ulz(i-1,1)+ulz(i,1)
      rn=sqrt(nx*nx+ny*ny+nz*nz)
      nx=nx/rn
      ny=ny/rn
      nz=nz/rn
      x=ox(i-1)-ox(i)
      y=oy(i-1)-oy(i)
      z=oz(i-1)-oz(i)
      vkin(i,7,1)=sqrt(x*x+y*y+z*z)
      qx=(ox(i-1)+ox(i))/2.
      qy=(oy(i-1)+oy(i))/2.
      qz=(oz(i-1)+oz(i))/2.
      x=vax(i-1)+vax(i)
      y=vay(i-1)+vay(i)
      z=vaz(i-1)+vaz(i)
      endif
      dot=nx*x+ny*y+nz*z
      dx=x-nx*dot
      dy=y-ny*dot
      dz=z-nz*dot
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      dx=dx/rd
      dy=dy/rd
      dz=dz/rd
      fx=ny*dz-nz*dy
      fy=nz*dx-nx*dz
      fz=nx*dy-ny*dx
      ka=k
      if(comb) ka=1
      km=k
      if(comb.and.li(i-1,k).lt.-1) km=iact(i-1)
      kp=k
      if(comb.and.li(i,k).lt.-1)   kp=iact(i)
c-----------------------------------------------u vector intersections
      x=qx-px(i-1,km)
      y=qy-py(i-1,km)
      z=qz-pz(i-1,km)
      dot1=nx*x+ny*y+nz*z
      dot2=nx*ulx(i-1,ka)+ny*uly(i-1,ka)+nz*ulz(i-1,ka)
      dl=dot1/dot2
      plx=px(i-1,km)+ulx(i-1,ka)*dl
      ply=py(i-1,km)+uly(i-1,ka)*dl
      plz=pz(i-1,km)+ulz(i-1,ka)*dl
      x=px(i,kp)-qx
      y=py(i,kp)-qy
      z=pz(i,kp)-qz
      dot1=nx*x+ny*y+nz*z
      dot2=nx*ulx(i,ka)+ny*uly(i,ka)+nz*ulz(i,ka)
      du=dot1/dot2
      pux=px(i,kp)-ulx(i,ka)*du
      puy=py(i,kp)-uly(i,ka)*du
      puz=pz(i,kp)-ulz(i,ka)*du
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) hel(i,3,k)=dl+du
c-----------------------------------------------------shift parameters
      x=pux-plx
      y=puy-ply
      z=puz-plz
      vkin(i,1,k)= dx*x+dy*y+dz*z
      vkin(i,2,k)=(fx*x+fy*y+fz*z)*idr(k)
      vkin(i,5,k)=sqrt(x*x+y*y+z*z)
      dot=ulx(i-1,ka)*ulx(i,ka)+uly(i-1,ka)*uly(i,ka)+
     1 ulz(i-1,ka)*ulz(i,ka)
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      vkin(i,6,k)=acos(dot)*crd
c------------------------------------------------------kink parameters
      tx=uly(i,ka)*dz-ulz(i,ka)*dy
      ty=ulz(i,ka)*dx-ulx(i,ka)*dz
      tz=ulx(i,ka)*dy-uly(i,ka)*dx
      rt=sqrt(tx*tx+ty*ty+tz*tz)
      dot=(fx*tx+fy*ty+fz*tz)/rt
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      cln=acos(dot)*crd
      x=fy*tz-fz*ty
      y=fz*tx-fx*tz
      z=fx*ty-fy*tx
      if(x*dx+y*dy+z*dz.lt.0.) cln=-cln
      rx=dy*tz-dz*ty
      ry=dz*tx-dx*tz
      rz=dx*ty-dy*tx
      rr=sqrt(rx*rx+ry*ry+rz*rz)
      dot=(ulx(i,ka)*rx+uly(i,ka)*ry+ulz(i,ka)*rz)/rr
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      tip=acos(dot)*crd
      x=ry*ulz(i,ka)-rz*uly(i,ka)
      y=rz*ulx(i,ka)-rx*ulz(i,ka)
      z=rx*uly(i,ka)-ry*ulx(i,ka)
      if(x*tx+y*ty+z*tz.lt.0.) tip=-tip
      vkin(i,3,k)=2*cln
      vkin(i,4,k)=2*tip*idr(k)
c--------------------------------------------------------------winding
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      hel(i,6,k)=0.
      do l=i-1,i
      sa=sin(cdr*(-cln))
      ca=cos(cdr*(cln))
      if(l.eq.i) sa=sin(cdr*(cln))
      fpx=(dx*dx+(1-dx*dx)*ca)*fx+(dx*dy*(1-ca)-dz*sa)*fy+
     1   (dx*dz*(1-ca)+dy*sa)*fz
      fpy=(dx*dy*(1-ca)+dz*sa)*fx+(dy*dy+(1-dy*dy)*ca)*fy+
     1   (dy*dz*(1-ca)-dx*sa)*fz
      fpz=(dx*dz*(1-ca)-dy*sa)*fx+(dy*dz*(1-ca)+dx*sa)*fy+
     1   (dz*dz+(1-dz*dz)*ca)*fz
      dot=fpx*wx(l,k)+fpy*wy(l,k)+fpz*wz(l,k)
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      x=fpy*wz(l,k)-fpz*wy(l,k)
      y=fpz*wx(l,k)-fpx*wz(l,k)
      z=fpx*wy(l,k)-fpy*wx(l,k)
      dot=x*ulx(l,ka)+y*uly(l,ka)+z*ulz(l,ka)
      if(l.eq.i-1.and.dot.gt.0.) wdg=-wdg
      if(l.eq.i  .and.dot.lt.0.) wdg=-wdg
      hel(i,6,k)=hel(i,6,k)+wdg
      enddo
      h=mod(hel(i,6,k),360.d0)
      if(abs(h).gt.180.) hel(i,6,k)=h-sign(360.d0,h)
      endif
c----------------------------------------------jumna style kink params
      if(k.eq.1.and..not.kcopy) then
      dx=vx(i-1,k)
      dy=vy(i-1,k)
      dz=vz(i-1,k)
      fx=wx(i-1,k)
      fy=wy(i-1,k)
      fz=wz(i-1,k)
      x=px(i,k)-px(i-1,k)
      y=py(i,k)-py(i-1,k)
      z=pz(i,k)-pz(i-1,k)
      dot1=ulx(i-1,k)*x+uly(i-1,k)*y+ulz(i-1,k)*z
      dot2=ulx(i-1,k)*ulx(i,k)+uly(i-1,k)*uly(i,k)+ulz(i-1,k)*ulz(i,k)
      dd=dot1/dot2
      hold(i,1,k)=dd
      psx=px(i,k)-ulx(i,k)*dd
      psy=py(i,k)-uly(i,k)*dd
      psz=pz(i,k)-ulz(i,k)*dd
      x=psx-px(i-1,k)
      y=psy-py(i-1,k)
      z=psz-pz(i-1,k)
      vold(i,1)= dx*x+dy*y+dz*z
      vold(i,2)=(fx*x+fy*y+fz*z)*idr(k)
      tx=uly(i,k)*dz-ulz(i,k)*dy
      ty=ulz(i,k)*dx-ulx(i,k)*dz
      tz=ulx(i,k)*dy-uly(i,k)*dx
      rt=sqrt(tx*tx+ty*ty+tz*tz)
      dot=(fx*tx+fy*ty+fz*tz)/rt
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      cln=acos(dot)*crd
      x=fy*tz-fz*ty
      y=fz*tx-fx*tz
      z=fx*ty-fy*tx
      if(x*dx+y*dy+z*dz.lt.0.) cln=-cln
      rx=dy*tz-dz*ty
      ry=dz*tx-dx*tz
      rz=dx*ty-dy*tx
      rr=sqrt(rx*rx+ry*ry+rz*rz)
      dot=(ulx(i,k)*rx+uly(i,k)*ry+ulz(i,k)*rz)/rr
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      tip=acos(dot)*crd
      x=ry*ulz(i,k)-rz*uly(i,k)
      y=rz*ulx(i,k)-rx*ulz(i,k)
      z=rx*uly(i,k)-ry*ulx(i,k)
      if(x*tx+y*ty+z*tz.lt.0.) tip=-tip
      vold(i,3)=cln
      vold(i,4)=tip*idr(k)
c-------------------------------------------jumna style winding params
      ca=cos(cdr*(cln))
      sa=sin(cdr*(cln))
      fpx=(dx*dx+(1-dx*dx)*ca)*fx+(dx*dy*(1-ca)-dz*sa)*fy+
     1   (dx*dz*(1-ca)+dy*sa)*fz
      fpy=(dx*dy*(1-ca)+dz*sa)*fx+(dy*dy+(1-dy*dy)*ca)*fy+
     1   (dy*dz*(1-ca)-dx*sa)*fz
      fpz=(dx*dz*(1-ca)-dy*sa)*fx+(dy*dz*(1-ca)+dx*sa)*fy+
     1   (dz*dz+(1-dz*dz)*ca)*fz
      dot=fpx*wx(i,k)+fpy*wy(i,k)+fpz*wz(i,k)
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      x=fpy*wz(i,k)-fpz*wy(i,k)
      y=fpz*wx(i,k)-fpx*wz(i,k)
      z=fpx*wy(i,k)-fpy*wx(i,k)
      if(x*ulx(i,k)+y*uly(i,k)+z*ulz(i,k).lt.0) wdg=-wdg
      hold(i,2,k)=wdg
      endif
c---------------------------------------------------------------------
      enddo
      enddo
c=================================================other strand zsh/wdg
      if(comb) then
      do k=2,ns
      dzl=0.
      dwl=0.
      do i=iste,iene
      if(li(i,k).ge.-1.and.li(i,1).ge.-1) then
      x=px(i,k)-px(i,1)
      y=py(i,k)-py(i,1)
      z=pz(i,k)-pz(i,1)
      dzu=sqrt(x*x+y*y+z*z)
      if(x*ulx(i,1)+y*uly(i,1)+z*ulz(i,1).lt.0.) dzu=-dzu
      dot=-(wx(i,1)*wx(i,k)+wy(i,1)*wy(i,k)+wz(i,1)*wz(i,k))
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      dwu=acos(dot)*crd
      x=wy(i,1)*wz(i,k)-wz(i,1)*wy(i,k)
      y=wz(i,1)*wx(i,k)-wx(i,1)*wz(i,k)
      z=wx(i,1)*wy(i,k)-wy(i,1)*wx(i,k)
      if(x*ulx(i,1)+y*uly(i,1)+z*ulz(i,1).gt.0.) dwu=-dwu
      hold(i,1,k)=hold(i,1,1)+dzu-dzl
      hold(i,2,k)=hold(i,2,1)+dwu-dwl
      if(i.eq.ng(k).and.i.ge.ng(1)) then
      hel(i,3,k)=hel(i,3,1)+dzu
      hel(i,6,k)=hel(i,6,1)+dwu
      else if(i.eq.ng(1).and.i.ge.ng(k)) then
      hel(i,3,1)=hel(i,3,k)-dzu
      hel(i,6,1)=hel(i,6,k)-dwu
      endif
      dzl=dzu
      dwl=dwu
      endif
      enddo
      enddo
      endif
      return
      end
