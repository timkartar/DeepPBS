      subroutine locpar
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,na*1
      integer*4 break,spline
      real*8 nx,ny,nz
      dimension ulx(0:n3),uly(0:n3),ulz(0:n3),
     1          px(0:n3),py(0:n3),pz(0:n3),
     1          wx(0:n3),wy(0:n3),wz(0:n3),
     1          vx(0:n3),vy(0:n3),vz(0:n3)
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
c==========================================local inter base parameters
      do k=1,max(ns,nst)
      id=idr(k)
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
c------------------------------------------------axis system direction
      inv(k)=1
      rise=0.
      do i=iste+1,iene
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      xu=rex(i,3,k)+rex(i-1,3,k)
      yu=rey(i,3,k)+rey(i-1,3,k)
      zu=rez(i,3,k)+rez(i-1,3,k)
      xp=rex(i,4,k)-rex(i-1,4,k)
      yp=rey(i,4,k)-rey(i-1,4,k)
      zp=rez(i,4,k)-rez(i-1,4,k)
      rise=rise+sign(1.d0,xu*xp+yu*yp+zu*zp)
      endif
      enddo
      if(id.eq. 1.and.rise.lt.0.) inv(k)=-1
      if(id.eq.-1.and.rise.gt.0.) inv(k)=-1
c-----------------------------------------------mean plane axis system
      if(comb.and.k.gt.1) then
      lu=inv(1)*idr(1)*idr(k)
      lv=-lu
      lw=-1
      else
      lu=inv(k)
      lv=lu
      lw=1
      endif
      do i=iste+1,iene
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      nx=rex(i-1,3,k)+rex(i,3,k)
      ny=rey(i-1,3,k)+rey(i,3,k)
      nz=rez(i-1,3,k)+rez(i,3,k)
      rn=sqrt(nx*nx+ny*ny+nz*nz)
      nx=lu*nx/rn
      ny=lu*ny/rn
      nz=lu*nz/rn
      qx=(rex(i-1,4,k)+rex(i,4,k))/2.
      qy=(rey(i-1,4,k)+rey(i,4,k))/2.
      qz=(rez(i-1,4,k)+rez(i,4,k))/2.
      x=lv*(rex(i-1,1,k)+rex(i,1,k))
      y=lv*(rey(i-1,1,k)+rey(i,1,k))
      z=lv*(rez(i-1,1,k)+rez(i,1,k))
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
c-----------------------------------------------u vector intersections
      x=qx-rex(i-1,4,k)
      y=qy-rey(i-1,4,k)
      z=qz-rez(i-1,4,k)
      dot1=nx*x+ny*y+nz*z
      dot2=lu*(nx*rex(i-1,3,k)+ny*rey(i-1,3,k)+nz*rez(i-1,3,k))
      dl=dot1/dot2
      plx=rex(i-1,4,k)+lu*rex(i-1,3,k)*dl
      ply=rey(i-1,4,k)+lu*rey(i-1,3,k)*dl
      plz=rez(i-1,4,k)+lu*rez(i-1,3,k)*dl
      x=rex(i,4,k)-qx
      y=rey(i,4,k)-qy
      z=rez(i,4,k)-qz
      dot1=nx*x+ny*y+nz*z
      dot2=lu*(nx*rex(i,3,k)+ny*rey(i,3,k)+nz*rez(i,3,k))
      du=dot1/dot2
      pux=rex(i,4,k)-lu*rex(i,3,k)*du
      puy=rey(i,4,k)-lu*rey(i,3,k)*du
      puz=rez(i,4,k)-lu*rez(i,3,k)*du
c-----------------------------------------------------shift parameters
      pal(i,3,k)=dl+du
      x=pux-plx
      y=puy-ply
      z=puz-plz
      pal(i,1,k)=  dx*x+dy*y+dz*z
      pal(i,2,k)= (fx*x+fy*y+fz*z)*idr(1)
c------------------------------------------------------kink parameters
      tx=lu*(rey(i,3,k)*dz-rez(i,3,k)*dy)
      ty=lu*(rez(i,3,k)*dx-rex(i,3,k)*dz)
      tz=lu*(rex(i,3,k)*dy-rey(i,3,k)*dx)
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
      dot=lu*(rex(i,3,k)*rx+rey(i,3,k)*ry+rez(i,3,k)*rz)/rr
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      tip=acos(dot)*crd
      x=lu*(ry*rez(i,3,k)-rz*rey(i,3,k))
      y=lu*(rz*rex(i,3,k)-rx*rez(i,3,k))
      z=lu*(rx*rey(i,3,k)-ry*rex(i,3,k))
      if(x*tx+y*ty+z*tz.lt.0.) tip=-tip
      pal(i,4,k)=2*cln
      pal(i,5,k)=2*tip*idr(1)
c--------------------------------------------------------------winding
      pal(i,6,k)=0.
      do l=i-1,i
      sa=sin(cdr*(-cln))
      ca=cos(cdr*(cln))
      if(l.eq.i) sa=sin(cdr*(cln))
      fpx=(dx*dx+(1-dx*dx)*ca)*fx+(dx*dy*(1-ca)-dz*sa)*fy+
     1 (dx*dz*(1-ca)+dy*sa)*fz
      fpy=(dx*dy*(1-ca)+dz*sa)*fx+(dy*dy+(1-dy*dy)*ca)*fy+
     1 (dy*dz*(1-ca)-dx*sa)*fz
      fpz=(dx*dz*(1-ca)-dy*sa)*fx+(dy*dz*(1-ca)+dx*sa)*fy+
     1 (dz*dz+(1-dz*dz)*ca)*fz
      dot=lw*(fpx*rex(l,2,k)+fpy*rey(l,2,k)+fpz*rez(l,2,k))
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      x=lw*(fpy*rez(l,2,k)-fpz*rey(l,2,k))
      y=lw*(fpz*rex(l,2,k)-fpx*rez(l,2,k))
      z=lw*(fpx*rey(l,2,k)-fpy*rex(l,2,k))
      dot=lu*(x*rex(l,3,k)+y*rey(l,3,k)+z*rez(l,3,k))
      if(l.eq.i-1.and.dot.gt.0.) wdg=-wdg
      if(l.eq.i .and.dot.lt.0.) wdg=-wdg
      pal(i,6,k)=pal(i,6,k)+wdg
      enddo
      h=mod(pal(i,6,k),360.d0)
      if(abs(h).gt.180.) pal(i,6,k)=h-sign(360.d0,h)
      endif
      enddo
      enddo
c====================================local inter base level parameters
      if(comb) then
      lu1=inv(1)
      lv1=lu1
      lw1=1
      do k=2,ns
      lu=inv(1)*idr(1)*idr(k)
      lv=-lu
      lw=-1
      do i=iste,iene
      if(li(i,k).ge.-1) then
      x=lu1*rex(i,3,1)+lu*rex(i,3,k)
      y=lu1*rey(i,3,1)+lu*rey(i,3,k)
      z=lu1*rez(i,3,1)+lu*rez(i,3,k)
      r=sqrt(x*x+y*y+z*z)
      ulx(i)=x/r
      uly(i)=y/r
      ulz(i)=z/r
      x=lv1*rex(i,1,1)+lv*rex(i,1,k)
      y=lv1*rey(i,1,1)+lv*rey(i,1,k)
      z=lv1*rez(i,1,1)+lv*rez(i,1,k)
      r=sqrt(x*x+y*y+z*z)
      vx(i)=x/r
      vy(i)=y/r
      vz(i)=z/r
      wx(i)=uly(i)*vz(i)-ulz(i)*vy(i)
      wy(i)=ulz(i)*vx(i)-ulx(i)*vz(i)
      wz(i)=ulx(i)*vy(i)-uly(i)*vx(i)
      px(i)=(rex(i,4,1)+rex(i,4,k))/2.
      py(i)=(rey(i,4,1)+rey(i,4,k))/2.
      pz(i)=(rez(i,4,1)+rez(i,4,k))/2.
      endif
      enddo
c-----------------------------------------------mean plane axis system
      do i=iste+1,iene
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      nx=ulx(i-1)+ulx(i)
      ny=uly(i-1)+uly(i)
      nz=ulz(i-1)+ulz(i)
      rn=sqrt(nx*nx+ny*ny+nz*nz)
      nx=nx/rn
      ny=ny/rn
      nz=nz/rn
      qx=(px(i-1)+px(i))/2.
      qy=(py(i-1)+py(i))/2.
      qz=(pz(i-1)+pz(i))/2.
      x=vx(i-1)+vx(i)
      y=vy(i-1)+vy(i)
      z=vz(i-1)+vz(i)
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
c-----------------------------------------------u vector intersections
      x=qx-px(i-1)
      y=qy-py(i-1)
      z=qz-pz(i-1)
      dot1=nx*x+ny*y+nz*z
      dot2=nx*ulx(i-1)+ny*uly(i-1)+nz*ulz(i-1)
      dl=dot1/dot2
      plx=px(i-1)+ulx(i-1)*dl
      ply=py(i-1)+uly(i-1)*dl
      plz=pz(i-1)+ulz(i-1)*dl
      x=px(i)-qx
      y=py(i)-qy
      z=pz(i)-qz
      dot1=nx*x+ny*y+nz*z
      dot2=nx*ulx(i)+ny*uly(i)+nz*ulz(i)
      du=dot1/dot2
      pux=px(i)-ulx(i)*du
      puy=py(i)-uly(i)*du
      puz=pz(i)-ulz(i)*du
c-----------------------------------------------------shift parameters
      pab(i,3,k)=dl+du
      x=pux-plx
      y=puy-ply
      z=puz-plz
      pab(i,1,k)=  dx*x+dy*y+dz*z
      pab(i,2,k)= (fx*x+fy*y+fz*z)*idr(1)
c------------------------------------------------------kink parameters
      tx=uly(i)*dz-ulz(i)*dy
      ty=ulz(i)*dx-ulx(i)*dz
      tz=ulx(i)*dy-uly(i)*dx
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
      dot=(ulx(i)*rx+uly(i)*ry+ulz(i)*rz)/rr
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      tip=acos(dot)*crd
      x=ry*ulz(i)-rz*uly(i)
      y=rz*ulx(i)-rx*ulz(i)
      z=rx*uly(i)-ry*ulx(i)
      if(x*tx+y*ty+z*tz.lt.0.) tip=-tip
      pab(i,4,k)=2*cln
      pab(i,5,k)=2*tip*idr(1)
c--------------------------------------------------------------winding
      pab(i,6,k)=0.
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
      dot=fpx*wx(l)+fpy*wy(l)+fpz*wz(l)
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      x=fpy*wz(l)-fpz*wy(l)
      y=fpz*wx(l)-fpx*wz(l)
      z=fpx*wy(l)-fpy*wx(l)
      dot=x*ulx(l)+y*uly(l)+z*ulz(l)
      if(l.eq.i-1.and.dot.gt.0.) wdg=-wdg
      if(l.eq.i .and.dot.lt.0.) wdg=-wdg
      pab(i,6,k)=pab(i,6,k)+wdg
      enddo
      h=mod(pab(i,6,k),360.d0)
      if(abs(h).gt.180.) pab(i,6,k)=h-sign(360.d0,h)
      endif
      enddo
      enddo
      endif
      return
      end
