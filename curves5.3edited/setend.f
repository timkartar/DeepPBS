      subroutine setend
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      integer*4 break,spline
      character na*1
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/hol/uho(0:n3,3,4),hho(0:n3,3,4),vkin(0:n3,7,4),bend(n3,4),
     1 hold(0:n3,2,4),vold(0:n3,4),pal(0:n3,6,4),pab(0:n3,6,4),inv(4)
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      do m=0,n+1,n+1
      id=1
      is=idr(nl)
      if(m.eq.n+1) then
      id=-1
      is=-idr(nl)
      endif
      xdi= hel(m,1,nl)
      ydi= hel(m,2,nl)
      rise=hel(m,3,nl)
      cln= hel(m,4,nl)
      tip= hel(m,5,nl)
      twis=hel(m,6,nl)
      rx=ux(m+id)
      ry=uy(m+id)
      rz=uz(m+id)
      uho(m,1,nl)=rx
      uho(m,2,nl)=ry
      uho(m,3,nl)=rz
c----------------------------------------------introduce end rise/twis
      ox(m)=ox(m+id)-rx*rise*is
      oy(m)=oy(m+id)-ry*rise*is
      oz(m)=oz(m+id)-rz*rise*is
      hho(m,1,nl)=ox(m)
      hho(m,2,nl)=oy(m)
      hho(m,3,nl)=oz(m)
      dax=rex(m+id,2,nl)
      day=rey(m+id,2,nl)
      daz=rez(m+id,2,nl)
      dot=rx*dax+ry*day+rz*daz
      xx=dax-rx*dot
      yy=day-ry*dot
      zz=daz-rz*dot
      r=sqrt(xx*xx+yy*yy+zz*zz)
      xx=xx/r
      yy=yy/r
      zz=zz/r
      ca=cos(cdr*(    twis))
      sa=sin(cdr*(-is*twis))
      wx=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1  (rx*rz*(1-ca)+ry*sa)*zz
      wy=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1  (ry*rz*(1-ca)-rx*sa)*zz
      wz=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1  (rz*rz+(1-rz*rz)*ca)*zz
      vx=wy*rz-wz*ry
      vy=wz*rx-wx*rz
      vz=wx*ry-wy*rx
      rex(m,4,nl)=ox(m)+vx*xdi+wx*ydi
      rey(m,4,nl)=oy(m)+vy*xdi+wy*ydi
      rez(m,4,nl)=oz(m)+vz*xdi+wz*ydi
      ca=cos(cdr*(cln))
      sa=sin(cdr*(cln))
      tx=(vx*vx+(1-vx*vx)*ca)*wx+(vx*vy*(1-ca)-vz*sa)*wy+
     1  (vx*vz*(1-ca)+vy*sa)*wz
      ty=(vx*vy*(1-ca)+vz*sa)*wx+(vy*vy+(1-vy*vy)*ca)*wy+
     1  (vy*vz*(1-ca)-vx*sa)*wz
      tz=(vx*vz*(1-ca)-vy*sa)*wx+(vy*vz*(1-ca)+vx*sa)*wy+
     1  (vz*vz+(1-vz*vz)*ca)*wz
      rex(m,2,nl)=tx
      rey(m,2,nl)=ty
      rez(m,2,nl)=tz
      qx=(vx*vx+(1-vx*vx)*ca)*rx+(vx*vy*(1-ca)-vz*sa)*ry+
     1  (vx*vz*(1-ca)+vy*sa)*rz
      qy=(vx*vy*(1-ca)+vz*sa)*rx+(vy*vy+(1-vy*vy)*ca)*ry+
     1  (vy*vz*(1-ca)-vx*sa)*rz
      qz=(vx*vz*(1-ca)-vy*sa)*rx+(vy*vz*(1-ca)+vx*sa)*ry+
     1  (vz*vz+(1-vz*vz)*ca)*rz
      ca=cos(cdr*(tip))
      sa=sin(cdr*(tip))
      rx=(tx*tx+(1-tx*tx)*ca)*qx+(tx*ty*(1-ca)-tz*sa)*qy+
     1  (tx*tz*(1-ca)+ty*sa)*qz
      ry=(tx*ty*(1-ca)+tz*sa)*qx+(ty*ty+(1-ty*ty)*ca)*qy+
     1  (ty*tz*(1-ca)-tx*sa)*qz
      rz=(tx*tz*(1-ca)-ty*sa)*qx+(ty*tz*(1-ca)+tx*sa)*qy+
     1  (tz*tz+(1-tz*tz)*ca)*qz
      rex(m,3,nl)=rx
      rey(m,3,nl)=ry
      rez(m,3,nl)=rz
      rex(m,1,nl)=ty*rz-tz*ry
      rey(m,1,nl)=tz*rx-tx*rz
      rez(m,1,nl)=tx*ry-ty*rx
c---------------------------------------------------------------------
      if(comb) then
      rex(m,4,2)=ox(m)+vx*xdi-wx*ydi
      rey(m,4,2)=oy(m)+vy*xdi-wy*ydi
      rez(m,4,2)=oz(m)+vz*xdi-wz*ydi
      tx=-tx
      ty=-ty
      tz=-tz
      qx=-qx
      qy=-qy
      qz=-qz
      rex(m,2,2)=tx
      rey(m,2,2)=ty
      rez(m,2,2)=tz
      ca=cos(cdr*(tip))
      sa=sin(cdr*(tip))
      rx=(tx*tx+(1-tx*tx)*ca)*qx+(tx*ty*(1-ca)-tz*sa)*qy+
     1  (tx*tz*(1-ca)+ty*sa)*qz
      ry=(tx*ty*(1-ca)+tz*sa)*qx+(ty*ty+(1-ty*ty)*ca)*qy+
     1  (ty*tz*(1-ca)-tx*sa)*qz
      rz=(tx*tz*(1-ca)-ty*sa)*qx+(ty*tz*(1-ca)+tx*sa)*qy+
     1  (tz*tz+(1-tz*tz)*ca)*qz
      rex(m,3,2)=rx
      rey(m,3,2)=ry
      rez(m,3,2)=rz
      rex(m,1,2)=ty*rz-tz*ry
      rey(m,1,2)=tz*rx-tx*rz
      rez(m,1,2)=tx*ry-ty*rx
      endif
      enddo
      return
      end
