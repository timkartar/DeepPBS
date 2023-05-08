      subroutine lsfit(lequ,n,nind,dcor,rms,key)
      include 'curves_data.inc'
      character*4 mnam,munit,ban,base*1
      real*8 h(3,3),k(3,3)
      dimension cg(3),u(3,3),w(21),v(6,6),dcor(10,3),nind(10)
      common/lsf/bref(10,3,5),th1,th2,dis,rs2,ibref(5),iequ(9),
     1 ibd(20,2),ban(9),base(9)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      sq2=sqrt(2.d0)
      key=0
c----------------------------------select atoms of real base & find cg
      cg(1)=0.
      cg(2)=0.
      cg(3)=0.
      do i=1,n
      m=nind(i)
      do j=1,3
      cg(j)=cg(j)+corm(m,j)
      enddo
      enddo
      cg(1)=cg(1)/n
      cg(2)=cg(2)/n
      cg(3)=cg(3)/n
c-----------------------------------------------------find ls rotation
      do i=1,3
      do j=1,3
      u(i,j)=0.
      enddo
      enddo
      do l=1,n
      m=nind(l)
      do i=1,3
      do j=1,3
      u(i,j)=u(i,j)+bref(l,i,lequ)*(corm(m,j)-cg(j))/n
      enddo
      enddo
      enddo
      det=u(1,1)*(u(2,2)*u(3,3)-u(2,3)*u(3,2))
     1  -u(1,2)*(u(2,1)*u(3,3)-u(2,3)*u(3,1))
     1  +u(1,3)*(u(2,1)*u(3,2)-u(2,2)*u(3,1))
      if(abs(det).lt.1.d-9) then
      key=1
      return
      endif
      m=0
      do i=1,6
      do j=1,i
      m=m+1
      if(i.gt.3.and.j.le.3) then
      w(m)=u(j,i-3)
      else
      w(m)=0.
      endif
      enddo
      enddo
      call eigen(w,v,6)
      if(det.lt.0.0.and.abs(w(3)-w(6)).lt.1.e-6) then
      key=1
      return
      endif
      do i=1,3
      do j=1,3
      h(i,j)=sq2*v(i,j)
      k(i,j)=sq2*v(i+3,j)
      enddo
      enddo
      sn=h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))
     1 +h(2,3)*(h(3,1)*h(1,2)-h(1,1)*h(3,2))
     1 +h(3,3)*(h(1,1)*h(2,2)-h(2,1)*h(1,2))
      if(sn.lt.0.) then
      do i=1,3
      h(i,3)=-h(i,3)
      k(i,3)=-k(i,3)
      enddo
      endif
      do i=1,3
      do j=1,3
      u(i,j)=k(i,1)*h(j,1)+k(i,2)*h(j,2)+sign(1.d0,det)*k(i,3)*h(j,3)
      enddo
      enddo
c-----------------------------------------------------make ls rotation
      rms=0.
      do l=1,n
      m=nind(l)
      x0=bref(l,1,lequ)
      y0=bref(l,2,lequ)
      z0=bref(l,3,lequ)
      dcor(l,1)=u(1,1)*x0+u(1,2)*y0+u(1,3)*z0+cg(1)
      dcor(l,2)=u(2,1)*x0+u(2,2)*y0+u(2,3)*z0+cg(2)
      dcor(l,3)=u(3,1)*x0+u(3,2)*y0+u(3,3)*z0+cg(3)
      dx=corm(m,1)-dcor(l,1)
      dy=corm(m,2)-dcor(l,2)
      dz=corm(m,3)-dcor(l,3)
      rms=rms+(dx*dx+dy*dy+dz*dz)
      enddo
      rms=sqrt(rms/n)
      return
      end
