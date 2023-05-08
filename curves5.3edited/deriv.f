      subroutine deriv
      include 'curves_data.inc'
      character*1 na
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly
      integer*4 break,spline
      dimension dot3(4)
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
c---------------------------------------------------------u & p derivs
      ido=ist
      iup=ien
      if(line) iup=1
5     do i=ido,iup
      nch=nl
      if(comb) nch=iact(i)
      id=1
      if(nch.ne.nl) id=-1
      dx=rex(i,2,nch)
      dy=rey(i,2,nch)
      dz=rez(i,2,nch)
      ex=rex(i,4,nch)-hx(i)
      ey=rey(i,4,nch)-hy(i)
      ez=rez(i,4,nch)-hz(i)
      vx=dy*uz(i)-dz*uy(i)
      vy=dz*ux(i)-dx*uz(i)
      vz=dx*uy(i)-dy*ux(i)
      rv=sqrt(vx*vx+vy*vy+vz*vz)
      udx(i,4)=-vx
      udy(i,4)=-vy
      udz(i,4)=-vz
      vx=vx/rv
      vy=vy/rv
      vz=vz/rv
      udx(i,3)=(uy(i)*vz-uz(i)*vy)*id
      udy(i,3)=(uz(i)*vx-ux(i)*vz)*id
      udz(i,3)=(ux(i)*vy-uy(i)*vx)*id
      udx(i,2)=0.
      udy(i,2)=0.
      udz(i,2)=0.
      udx(i,1)=0.
      udy(i,1)=0.
      udz(i,1)=0.
      pdx(i,1)=-vx*id
      pdy(i,1)=-vy*id
      pdz(i,1)=-vz*id
      pdx(i,2)=-udx(i,3)*id
      pdy(i,2)=-udy(i,3)*id
      pdz(i,2)=-udz(i,3)*id
      pdx(i,3)=(vy*ez-vz*ey)*id
      pdy(i,3)=(vz*ex-vx*ez)*id
      pdz(i,3)=(vx*ey-vy*ex)*id
      pdx(i,4)=dy*ez-dz*ey
      pdy(i,4)=dz*ex-dx*ez
      pdz(i,4)=dx*ey-dy*ex
c----------------------------------------------combined p point derivs
      if(comb) then
      do l=1,4
      dot3(l)=ux(i)*pdx(i,l)+uy(i)*pdy(i,l)+uz(i)*pdz(i,l)
      enddo
      m=0
      do k=1,ns
      if(li(i,k).gt.0.or.(li(i,nch).eq.-1.and.li(i,k).eq.-1)) m=m+1
      enddo
      do k=1,ns
      if(k.ne.nch.and.li(i,k).gt.0) then
      x=rex(i,4,k)-hx(i)
      y=rey(i,4,k)-hy(i)
      z=rez(i,4,k)-hz(i)
      dot1=(ux(i)*x+uy(i)*y+uz(i)*z)*id
      do l=1,4
      dot2=udx(i,l)*x+udy(i,l)*y+udz(i,l)*z
      pdx(i,l)=pdx(i,l)+(udx(i,l)*dot1+ux(i)*(dot2-dot3(l)))/m
      pdy(i,l)=pdy(i,l)+(udy(i,l)*dot1+uy(i)*(dot2-dot3(l)))/m
      pdz(i,l)=pdz(i,l)+(udz(i,l)*dot1+uz(i)*(dot2-dot3(l)))/m
      enddo
      endif
      enddo
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
      if(.not.comb) then
      dx1=rex(i,4,nl)-ox(ido)
      dy1=rey(i,4,nl)-oy(ido)
      dz1=rez(i,4,nl)-oz(ido)
      else
      x1=0.
      y1=0.
      z1=0.
      m=0
      do k=1,ns
      if(li(i,k).gt.0.or.(li(i,nch).eq.-1.and.li(i,k).eq.-1)) then
      m=m+1
      x1=x1+rex(i,4,k)
      y1=y1+rey(i,4,k)
      z1=z1+rez(i,4,k)
      endif
      enddo
      dx1=x1/m-ox(ido)
      dy1=y1/m-oy(ido)
      dz1=z1/m-oz(ido)
      endif
      dot1=dx1*ux(ido)+dy1*uy(ido)+dz1*uz(ido)
      do j=1,4
      udx(i,j)=udx(ido,j)
      udy(i,j)=udy(ido,j)
      udz(i,j)=udz(ido,j)
      dot2=dx1*udx(ido,j)+dy1*udy(ido,j)+dz1*udz(ido,j)
      dop3=pdx(ido,j)*ux(ido)+pdy(ido,j)*uy(ido)+pdz(ido,j)*uz(ido)
      pdx(i,j)=pdx(ido,j)+udx(ido,j)*dot1+ux(ido)*(dot2-dop3)
      pdy(i,j)=pdy(ido,j)+udy(ido,j)*dot1+uy(ido)*(dot2-dop3)
      pdz(i,j)=pdz(ido,j)+udz(ido,j)*dot1+uz(ido)*(dot2-dop3)
      enddo
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
