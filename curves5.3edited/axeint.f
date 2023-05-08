      subroutine axeint
      include 'curves_data.inc'
      integer*4 spline,break
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly
      character*1 na
      dimension p(3,0:2),u(3,2),g(4,3),dx(0:1),dy(0:1),dz(0:1)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/gro/uxb(5000,0:n3),uyb(5000,0:n3),uzb(5000,0:n3),
     1 cor(n1,3),dya(n1,3),box(0:n1,4),boy(0:n1,4),
     1 boz(0:n1,4),ind(n3,0:n3),nma(n3),num,numa,nsu(4)
      common/hol/uho(0:n3,3,4),hho(0:n3,3,4),vkin(0:n3,7,4),bend(n3,4),
     1 hold(0:n3,2,4),vold(0:n3,4),pal(0:n3,6,4),pab(0:n3,6,4),inv(4)
c------------------------------------------u vectors at base pair levels
      do i=num,numa
      uxb(i,0)=uho(i,1,1)
      uyb(i,0)=uho(i,2,1)
      uzb(i,0)=uho(i,3,1)
      enddo
      do i=num+1,numa
      co=uxb(i-1,0)*uxb(i,0)+uyb(i-1,0)*uyb(i,0)+uzb(i-1,0)*uzb(i,0)
      if(co.lt.0.) then
      uxb(i,0)=-uxb(i,0)
      uyb(i,0)=-uyb(i,0)
      uzb(i,0)=-uzb(i,0)
      endif
      enddo
c------------------------------------------------------------smooth axis
      k=1
      do i=num,numa-1
      if(i+1.eq.break) then
      nma(i)=1
      ind(i,0)=k
      do j=1,3
      cor(k,j)=hho(i,j,1)
      enddo
      kbreak=k
      k=k+1
      goto 100
      endif
      do j=1,2
      ij=i+j-1
      p(1,j)=hho(ij,1,1)
      p(2,j)=hho(ij,2,1)
      p(3,j)=hho(ij,3,1)
      u(1,j)=uxb(ij,0)
      u(2,j)=uyb(ij,0)
      u(3,j)=uzb(ij,0)
      enddo
      gl=sqrt((p(1,2)-p(1,1))**2+(p(2,2)-p(2,1))**2+(p(3,2)-p(3,1))**2)
      do j=1,3
      g(1,j)=p(j,1)
      g(2,j)=u(j,1)
      g(3,j)= 3*(p(j,2)-p(j,1))/gl**2-(u(j,2)+2*u(j,1))/gl
      g(4,j)=-2*(p(j,2)-p(j,1))/gl**3+(u(j,2)+  u(j,1))/gl**2
      enddo
      rise=(hel(i+1,3,1)+hel(i+1,3,2))/2.
      nma(i)=int(abs(rise)*(nlevel+1)/3.4+0.5)
      if(nma(i).eq.0) nma(i)=1
      nmi=nma(i)
      do n=0,nmi-1
      s=gl*n/nmi
      ind(i,n)=k
      do j=1,3
      cor(k,j)=g(1,j)+g(2,j)*s+g(3,j)*s**2+g(4,j)*s**3
      enddo
      k=k+1
      enddo
100   enddo
      do j=1,3
      cor(k,j)=hho(numa,j,1)
      enddo
      ind(numa,0)=k
      nma(numa)=1
c----------------------------------------------------base pair dyad axes
      do i=num,numa
      k=ind(i,0)
      rx=rex(i,1,1)+rex(i,1,2)
      ry=rey(i,1,1)+rey(i,1,2)
      rz=rez(i,1,1)+rez(i,1,2)
      dot=uxb(i,0)*rx+uyb(i,0)*ry+uzb(i,0)*rz
      vx=rx-uxb(i,0)*dot
      vy=ry-uyb(i,0)*dot
      vz=rz-uzb(i,0)*dot
      v=sqrt(vx*vx+vy*vy+vz*vz)
      dya(k,1)=vx/v
      dya(k,2)=vy/v
      dya(k,3)=vz/v
      enddo
c-------------------------------------------------intermediate u vectors
      do i=num,numa-1
      nmi=nma(i)
      if(nmi.eq.0.or.nmi.eq.1) goto 101
      do n=1,nmi-1
      x=uxb(i,0)*(nmi-n)/nmi+uxb(i+1,0)*n/nmi
      y=uyb(i,0)*(nmi-n)/nmi+uyb(i+1,0)*n/nmi
      z=uzb(i,0)*(nmi-n)/nmi+uzb(i+1,0)*n/nmi
      um=sqrt(x*x+y*y+z*z)
      uxb(i,n)=x/um
      uyb(i,n)=y/um
      uzb(i,n)=z/um
c----------------------------------------------intermediate dyad vectors
      do l=0,1
      ik=i+l
      k=ind(ik,0)
      co=uxb(ik,0)*uxb(i,n)+uyb(ik,0)*uyb(i,n)+uzb(ik,0)*uzb(i,n)
      rx=uyb(i,n)*uzb(ik,0)-uyb(ik,0)*uzb(i,n)
      ry=uzb(i,n)*uxb(ik,0)-uzb(ik,0)*uxb(i,n)
      rz=uxb(i,n)*uyb(ik,0)-uxb(ik,0)*uyb(i,n)
      dp=rx*dya(k,1)+ry*dya(k,2)+rz*dya(k,3)
      dx(l)=dya(k,1)*co+dp*rx/(1+co)+ry*dya(k,3)-rz*dya(k,2)
      dy(l)=dya(k,2)*co+dp*ry/(1+co)+rz*dya(k,1)-rx*dya(k,3)
      dz(l)=dya(k,3)*co+dp*rz/(1+co)+rx*dya(k,2)-ry*dya(k,1)
      enddo
      co=dx(0)*dx(1)+dy(0)*dy(1)+dz(0)*dz(1)
      wx=dx(1)-dx(0)*co
      wy=dy(1)-dy(0)*co
      wz=dz(1)-dz(0)*co
      wm=sqrt(wx*wx+wy*wy+wz*wz)
      wx=wx/wm
      wy=wy/wm
      wz=wz/wm
      k=ind(i,n)
      an=acos(co)
      co=cos(n*an/nmi)
      si=sin(n*an/nmi)
      dya(k,1)=dx(0)*co+wx*si
      dya(k,2)=dy(0)*co+wy*si
      dya(k,3)=dz(0)*co+wz*si
      enddo
101   enddo
      return
      end
