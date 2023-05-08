      subroutine plate(ns)
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,flag
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,base(10)*1,mnam,munit,snam,sunit,
     1 name,na*1
      integer*4 ioc(n3),icp(2*n3),break,spline
      dimension a(3),b(3),cl(8,3),iba(11),ibm(12),ibr(13),
     1 sugp(n3,3),sugv(n3,3),phop(n3,3),phov(n3,3)
      common/axe/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),matc(n1,7),kas,khs,kces
      common/bak/tor(0:n3,4,13),val(0:n3,4,2),suga(0:n3,4,2),
     1 nat(n3,4,15),flag(n3,4)
      common/cha/file,lis,dna,axin,axout,daf,pdb,mcode
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/hlx/pts(n3,3),ss(3),c(3),rep(3),rad,pit,pha,tbn,dan,isen
      common/hol/uho(0:n3,3,4),hho(0:n3,3,4),vkin(0:n3,7,4),bend(n3,4),
     1 hold(0:n3,2,4),vold(0:n3,4),pal(0:n3,6,4),pab(0:n3,6,4),inv(4)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      data base/'G','A','I','Y','L','C','T','U','P','M'/
      data iba/1,2,3,4,5,20000,3,10000,2,6,10000/
      data ibm/1,2,3,4,5,6,20000,2,10000,4,7,10000/
      data ibr/6,1,2,0,10000,2,3,4,5,10000,4,20000,6/
c--------------------------------------------------setup ior reference
      if(ior.gt.0) then
      id=1
      do kc=1,max(ns,nst)
      if(li(ior,kc).ge.-1) goto 100
      enddo
         do kc=1,max(ns,nst)
         if(li(ior,kc).eq.-1) goto 100
         enddo
100   if(kc.gt.1) id=-1
      dx=rex(ior,2,kc)*id
      dy=rey(ior,2,kc)*id
      dz=rez(ior,2,kc)*id
      usx=uho(ior,1,1)
      usy=uho(ior,2,1)
      usz=uho(ior,3,1)
      dot=usx*dx+usy*dy+usz*dz
      wx=dx-usx*dot
      wy=dy-usy*dot
      wz=dz-usz*dot
      r=sqrt(wx*wx+wy*wy+wz*wz)
      dax=(wy*usz-wz*usy)/r
      day=(wz*usx-wx*usz)/r
      daz=(wx*usy-wy*usx)/r
      psx=hho(ior,1,1)
      psy=hho(ior,2,1)
      psz=hho(ior,3,1)
      endif
c----------------------------------------------------------spline axes
      kces=0
      kas=0
      khs=0
      if(spline.gt.0) then
      call cubspl
      goto 500
      endif
c----------------------------------------------------------normal axes
      do kk=1,nst
      if(comb) then
      ist=1
      ien=nux
      else
      ist=ng(kk)
      ien=nr(kk)
      endif
      if(ends) then
      ist=ist-1
      ien=ien+1
      endif
      do i=ist,ien
      kc=kk
      id=1
      if(comb) then
      kc=iact(i)
      id=-1
      endif
      do j=1,11
      khs=khs+1
      if(iba(j).lt.10000) then
      mats(khs)=iba(j)+kas
      else
      mats(khs)=iba(j)
      endif
      enddo
      ks=kas
      kces=kces+1
      do l=1,6
      kas=kas+1
      snam(kas)='A'
      sunit(kas)='AXIS'
      nunis(kas)=kces
      enddo
      ux=uho(i,1,kk)
      uy=uho(i,2,kk)
      uz=uho(i,3,kk)
      dx=rex(i,2,kc)*id
      dy=rey(i,2,kc)*id
      dz=rez(i,2,kc)*id
      dot=ux*dx+uy*dy+uz*dz
      wx=dx-ux*dot
      wy=dy-uy*dot
      wz=dz-uz*dot
      r=sqrt(wx*wx+wy*wy+wz*wz)
      wx=wx/r
      wy=wy/r
      wz=wz/r
      vx=wy*uz-wz*uy
      vy=wz*ux-wx*uz
      vz=wx*uy-wy*ux
      x0=hho(i,1,kk)
      y0=hho(i,2,kk)
      z0=hho(i,3,kk)
      cors(ks+1,1)=x0-ux*1.2
      cors(ks+1,2)=y0-uy*1.2
      cors(ks+1,3)=z0-uz*1.2
      cors(ks+2,1)=x0
      cors(ks+2,2)=y0
      cors(ks+2,3)=z0
      cors(ks+3,1)=x0+ux*1.2
      cors(ks+3,2)=y0+uy*1.2
      cors(ks+3,3)=z0+uz*1.2
      cors(ks+4,1)=x0+ux*0.90-wx*0.15
      cors(ks+4,2)=y0+uy*0.90-wy*0.15
      cors(ks+4,3)=z0+uz*0.90-wz*0.15
      cors(ks+5,1)=x0+ux*0.90+wx*0.15
      cors(ks+5,2)=y0+uy*0.90+wy*0.15
      cors(ks+5,3)=z0+uz*0.90+wz*0.15
      cors(ks+6,1)=x0+vx*0.85
      cors(ks+6,2)=y0+vy*0.85
      cors(ks+6,3)=z0+vz*0.85
      enddo
      enddo
c---------------------------------------------------------ref and base
500   if(axonly) goto 600
      kh=0
      kl=0
      do kk=1,max0(ns,nst)
      ist=ng(kk)
      ien=nr(kk)
      if(ends) then
      ist=ist-1
      ien=ien+1
      endif
         io=ist-1
         do i=ist,ien
         if(li(i,kk).ge.-2) then
         io=io+1
         ioc(io)=i
         endif
         enddo
         ien=io
      ic=0
      khold=kh
      do i=ist,ien
      ii=ioc(i)
      if(li(ii,kk).lt.-1) goto 50
      kh=kh+1
      if(ii.ge.ng(kk).and.ii.le.nr(kk)) then
      kl=kl+1
      name=na(ii,kk)
      else
      name='L'
      if(kk.eq.2) name='M'
      endif
      numb=kh
      rk=0.5
      rl=4.5
      do j=1,10
      if(name.eq.base(j)) goto 200
      enddo
200   if(j.gt.5) then
      rk=2.25
      rl=2.75
      endif
      a(1)= rex(ii,1,kk)
      a(2)= rey(ii,1,kk)
      a(3)= rez(ii,1,kk)
      b(1)=-rex(ii,2,kk)
      b(2)=-rey(ii,2,kk)
      b(3)=-rez(ii,2,kk)
      cl(1,1)=rex(ii,4,kk)
      cl(1,2)=rey(ii,4,kk)
      cl(1,3)=rez(ii,4,kk)
      do j=1,3
      cl(2,j)=cl(1,j)-rk*b(j)
      cl(3,j)=cl(2,j)-1.25*a(j)
      cl(4,j)=cl(3,j)-rl*b(j)
      cl(5,j)=cl(4,j)+2.75*a(j)
      cl(6,j)=cl(5,j)+rl*b(j)
      cl(7,j)=cl(4,j)-0.327*b(j)-0.491*a(j)
      enddo
c---------------------------------------------------------------------
      do j=1,12
      khs=khs+1
      if(ibm(j).lt.10000) then
      mats(khs)=ibm(j)+kas
      else
      mats(khs)=ibm(j)
      endif
      enddo
      do l=1,7
      kas=kas+1
      snam(kas)='P'
      sunit(kas)=name
      nunis(kas)=numb
      do j=1,3
      cors(kas,j)=cl(l,j)
      enddo
      enddo
      ic=ic+1
      if(ii.ge.ng(kk).and.ii.le.nr(kk)) icp(ic)=kas
 50   enddo
c==============================================================ribbons
      k=kas
      ic=0
c------------------------------------------------------phosphorus data
      do i=ist,ien-1
      iof=ioc(i)
      if(idr(kk).eq.-1) iof=ioc(i+1)
      i3=nat(iof,kk,8)
      ip=nat(iof,kk,9)
      i5=nat(iof,kk,10)
      ia=0
      ib=0
      do m=1,matd(ip,7)
      l=abs(matd(ip,m))
      if(l.ne.i3.and.l.ne.i5) then
      if(ia.eq.0) then
      ia=l
      else
      ib=l
      endif
      endif
      enddo
      x=corm(ia,1)-corm(ib,1)
      y=corm(ia,2)-corm(ib,2)
      z=corm(ia,3)-corm(ib,3)
      r=sqrt(x*x+y*y+z*z)
      x=x*wid/r
      y=y*wid/r
      z=z*wid/r
      if(i.gt.ist) then
      dot=x*xr+y*yr+z*zr
      xr=x*sign(1.d0,dot)
      yr=y*sign(1.d0,dot)
      zr=z*sign(1.d0,dot)
      else
      xr=x
      yr=y
      zr=z
      endif
      phop(i,1)=corm(ip,1)
      phop(i,2)=corm(ip,2)
      phop(i,3)=corm(ip,3)
      phov(i,1)=xr
      phov(i,2)=yr
      phov(i,3)=zr
      enddo
c-----------------------------------------------------------sugar data
      do i=ist,ien
      ii=ioc(i)
      ia=nat(ii,kk,5)
      ib=nat(ii,kk,6)
      sugp(i,1)=(corm(ia,1)+corm(ib,1))/2
      sugp(i,2)=(corm(ia,2)+corm(ib,2))/2
      sugp(i,3)=(corm(ia,3)+corm(ib,3))/2
      if(i.eq.ist) then
      sugv(i,1)=phov(i,1)
      sugv(i,2)=phov(i,2)
      sugv(i,3)=phov(i,3)
      else if(i.eq.ien) then
      sugv(i,1)=phov(i-1,1)
      sugv(i,2)=phov(i-1,2)
      sugv(i,3)=phov(i-1,3)
      else
      x=phov(i-1,1)+phov(i,1)
      y=phov(i-1,2)+phov(i,2)
      z=phov(i-1,3)+phov(i,3)
      r=sqrt(x*x+y*y+z*z)
      sugv(i,1)=x*wid/r
      sugv(i,2)=y*wid/r
      sugv(i,3)=z*wid/r
      endif
      enddo
c---------------------------------------------------------build ribbon
      kh=khold
      do i=ist,ien
      ii=ioc(i)
      kh=kh+1
      if(i.lt.ien) then
      do j=1,6
      k=k+1
      snam(k)='R'
      sunit(k)='RIB'
      nunis(k)=kh
      if(j.eq.1) then
      cors(k,1)=sugp(i,1)+sugv(i,1)
      cors(k,2)=sugp(i,2)+sugv(i,2)
      cors(k,3)=sugp(i,3)+sugv(i,3)
      else if(j.eq.2) then
      cors(k,1)=sugp(i,1)
      cors(k,2)=sugp(i,2)
      cors(k,3)=sugp(i,3)
      else if(j.eq.3) then
      cors(k,1)=sugp(i,1)-sugv(i,1)
      cors(k,2)=sugp(i,2)-sugv(i,2)
      cors(k,3)=sugp(i,3)-sugv(i,3)
      else if(j.eq.4) then
      cors(k,1)=phop(i,1)-phov(i,1)
      cors(k,2)=phop(i,2)-phov(i,2)
      cors(k,3)=phop(i,3)-phov(i,3)
      else if(j.eq.5) then
      cors(k,1)=sugp(i+1,1)-sugv(i+1,1)
      cors(k,2)=sugp(i+1,2)-sugv(i+1,2)
      cors(k,3)=sugp(i+1,3)-sugv(i+1,3)
      else
      cors(k,1)=phop(i,1)+phov(i,1)
      cors(k,2)=phop(i,2)+phov(i,2)
      cors(k,3)=phop(i,3)+phov(i,3)
      endif
      enddo
      if(li(ii,kk).ge.-1) ic=ic+1
      do j=1,13
      if(j.eq.4.and.li(ii,kk).lt.-1) goto 55
      khs=khs+1
      if(j.eq.4) then
      mats(khs)=icp(ic)
      else if(ibr(j).lt.10000) then
      mats(khs)=ibr(j)+kas
      else
      mats(khs)=ibr(j)
      endif
 55   enddo
      khs=khs+1
      mats(khs)=kas+7
      khs=khs+1
      mats(khs)=10000
      else
      do j=1,3
      k=k+1
      snam(k)='R'
      sunit(k)='RIB'
      nunis(k)=kh
      if(j.eq.1) then
      cors(k,1)=sugp(i,1)+sugv(i,1)
      cors(k,2)=sugp(i,2)+sugv(i,2)
      cors(k,3)=sugp(i,3)+sugv(i,3)
      else if(j.eq.2) then
      cors(k,1)=sugp(i,1)
      cors(k,2)=sugp(i,2)
      cors(k,3)=sugp(i,3)
      else
      cors(k,1)=sugp(i,1)-sugv(i,1)
      cors(k,2)=sugp(i,2)-sugv(i,2)
      cors(k,3)=sugp(i,3)-sugv(i,3)
      endif
      enddo
      if(li(ii,kk).ge.-1) ic=ic+1
      do j=2,7
      if(j.eq.4.and.li(ii,kk).lt.-1) goto 60
      khs=khs+1
      if(j.eq.4) then
      mats(khs)=icp(ic)
      else if(ibr(j).lt.10000) then
      mats(khs)=ibr(j)+kas
      else
      mats(khs)=ibr(j)
      endif
 60   enddo
      khs=khs+1
      mats(khs)=10000
      endif
      kas=k
      enddo
      enddo
      khs=khs-1
c------------------------------------------------------------transform
600   call setd
      if(ior.gt.0) then
      ca=usz
      sa=sin(acos(ca))
      r=sqrt(usx*usx+usy*usy)
      if(r.gt.0.) then
      rx= usy/r
      ry=-usx/r
      rz= 0.
      else
      rx=1.
      ry=0.
      rz=0.
      endif
c----------------------------------------------------update local dyad
      dbx=(rx*rx+(1-rx*rx)*ca)*dax+(rx*ry*(1-ca)-rz*sa)*day+
     1    (rx*rz*(1-ca)+ry*sa)*daz
      dby=(rx*ry*(1-ca)+rz*sa)*dax+(ry*ry+(1-ry*ry)*ca)*day+
     1    (ry*rz*(1-ca)-rx*sa)*daz
c--------------------rotate u vector to z axis and center on ref point
      do i=1,kas
      xx=cors(i,1)-psx
      yy=cors(i,2)-psy
      zz=cors(i,3)-psz
      cors(i,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1          (rx*rz*(1-ca)+ry*sa)*zz
      cors(i,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1          (ry*rz*(1-ca)-rx*sa)*zz
      cors(i,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1          (rz*rz+(1-rz*rz)*ca)*zz
      enddo
      do i=1,kam
      xx=corm(i,1)-psx
      yy=corm(i,2)-psy
      zz=corm(i,3)-psz
      corm(i,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1          (rx*rz*(1-ca)+ry*sa)*zz
      corm(i,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1          (ry*rz*(1-ca)-rx*sa)*zz
      corm(i,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1          (rz*rz+(1-rz*rz)*ca)*zz
      enddo
c-----------------------------------------------rotate dyad to -x axis
      th=pi-acos(dbx)
      if(dby.lt.0.) th=-th
      ca=cos(th)
      sa=sin(th)
      do i=1,kas
      xx=cors(i,1)
      yy=cors(i,2)
      cors(i,1)=xx*ca-yy*sa
      cors(i,2)=xx*sa+yy*ca
      enddo
      do i=1,kam
      xx=corm(i,1)
      yy=corm(i,2)
      corm(i,1)=xx*ca-yy*sa
      corm(i,2)=xx*sa+yy*ca
      enddo
      endif
c------------------------------------------------find subunit linkages
      name=sunit(1)
      numb=nunis(1)
      kces=1
      do i=2,kas
      if(name.ne.sunit(i).or.numb.ne.nunis(i)) then
      kces=kces+1
      name=sunit(i)
      numb=nunis(i)
      endif
      enddo
c--------------------------------------------------write plate macfile
      if(dna.ne.' ') then
      kfi=index(dna,' ')-1
      open(unit=1,file=dna(:kfi)//'.dna',status='new')
      write(1,5) '# Curve result from '//file
      write(1,6) kas,kces
      write(1,7) (snam(i),cors(i,1),cors(i,2),cors(i,3),2,
     1 20,0.0,sunit(i),nunis(i),0,i,i=1,kas)
      write(1,8) (i,(matc(i,j),j=1,6),i=1,kas)
5     format(a)
6     format(2i5)
7     format(a4,3f10.5,2i3,f8.4,1x,a4,2i4,i5)
8     format(i5,') ',6i5,1x,i5,') ',6i5)
      close(unit=1)
      endif
      if(pdb.ne.' ') call pdbout(pdb)
c---------------------------------------macfile reoriented coordinates
      if(ior.ne.0) then
      kfi=index(file,'.')-1
      open(unit=1,file=file(:kfi)//'_.mac',status='new')
      write(1,5) '# '//mcode
      write(1,6) kam,kcen
      write(1,7) (mnam(i),corm(i,1),corm(i,2),corm(i,3),imch(i),
     1 imty(i),dmon(i),munit(i),nunit(i),icm(i),i,i=1,kam)
      write(1,8) (i,(matd(i,j),j=1,6),i=1,kam)
      close(unit=1)
      endif
      return
      end
