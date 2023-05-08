      subroutine locate
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,purine,flag
      character*4 mnam,munit,na*1,name,ban,base*1
      integer*4 break,spline
      dimension inc(n3),dcor(10,3),nind(10)
      common/bak/tor(0:n3,4,13),val(0:n3,4,2),suga(0:n3,4,2),
     1 nat(n3,4,15),flag(n3,4)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/lsf/bref(10,3,5),th1,th2,dis,rs2,ibref(5),iequ(9),
     1 ibd(20,2),ban(9),base(9)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),nn,nl,ns
c--------------------------------------------------------------find c1
      lc=0
      do i=1,kam
      if(mnam(i).eq.'C1'''.or.mnam(i).eq.'C1*') then
      lc=lc+1
      inc(lc)=i
      endif
      enddo
c----------------------------------------------find base reference pts
      do n=1,nst
      do m=1,nux
      na(m,n)='-'
      enddo
      do m=ng(n),nr(n)
      iu=ni(m,n)
         if(iu.eq.0) then
         name='-'
         goto 20
         endif
      is=ncen(iu-1)+1
      ie=ncen(iu)
      name=munit(is)
      do l=1,9
      if(name(:1).eq.base(l)) goto 10
      enddo
         if(li(m,n).le.0) then
         name='*'
         if(li(m,n).eq.-1) li(m,n)=-2
         if(li(m,n).eq.-2) write(6,*) '.... no base for residue ',m,
     1                     ' in strand ',n
         do i=is,ie
         if(mnam(i).eq.'C1'''.or.mnam(i).eq.'C1*') nat(m,n,3)=i
         enddo
         goto 20
         endif
      write(6,*) '  ---- Base not located, check numbers ----'
      stop
10    lsav=l
      lequ=iequ(l)
      na(m,n)=name
      mn=9
      purine=.true.
      if(lsav.gt.4) then
      mn=6
      purine=.false.
      endif
      kb=1
      i1=0
      i2=0
      i3=0
      i4=0
      do i=is,ie
      name=mnam(i)
      if(name.eq.'C1'''.or.name.eq.'C1*') i1=i
      if(purine) then
      if(name.eq.'N9') i2=i
      if(name.eq.'C4') i3=i
      if(name.eq.'C8') i4=i
      else
      if(lsav.lt.9) then
      if(name.eq.'N1') i2=i
      if(name.eq.'C2') i3=i
      if(name.eq.'C6') i4=i
      else
      if(name.eq.'C5') i2=i
      if(name.eq.'C4') i3=i
      if(name.eq.'N1') i4=i
      endif
      endif
c------------------------------------------base atom search for ls fit
      if(fit) then
      do mm=1,mn
      if(name.eq.ban(mm)) then
      kb=kb+1
      nind(mm+1)=i
      endif
      enddo
      endif
      enddo
      if(fit.and.kb.ne.mn+1) then
         if(li(m,n).le.-1.and.i1.ne.0) then
         if(li(m,n).eq.-1) li(m,n)=-2
         write(6,*) '.... no base for residue ',m,' in strand ',n
         nat(m,n,3)=i1
         goto 20
         else
         write(6,*) '  ---- Not all base atoms for LS-fit found ----'
         stop
         endif
      endif
c----------------------------------------------------------sugar retry
      if(i1.eq.0.and.i2.ne.0) then
      x0=corm(i2,1)
      y0=corm(i2,2)
      z0=corm(i2,3)
      do l=1,lc
      j=inc(l)
      x=x0-corm(j,1)
      y=y0-corm(j,2)
      z=z0-corm(j,3)
      r2=x*x+y*y+z*z
      if(r2.lt.rs2) then
      i1=j
      goto 30
      endif
      enddo
      endif
c----------------------------------------------------------------check
30    name=' '
      if(i1.eq.0) name='C1'''
      if(purine) then
      if(i2.eq.0) name='N9'
      if(i3.eq.0) name='C4'
      if(i4.eq.0) name='C8'
      else
      if(lsav.lt.9) then
      if(i2.eq.0) name='N1'
      if(i3.eq.0) name='C2'
      if(i4.eq.0) name='C6'
      else
      if(i2.eq.0) name='C5'
      if(i3.eq.0) name='C4'
      if(i4.eq.0) name='N1'
      endif
      endif
      if(name.ne.' ') then
         if(li(m,n).eq.-1.and.i1.ne.0) then
         write(6,*) '.... no base for residue ',m,' in strand ',n
         li(m,n)=-2
         nat(m,n,3)=i1
         goto 20
         endif
      write(6,40) n,m,munit(is),iu,name
40    format(/2x,'Strand: ',i2,' Unit: ',i3,' (',a4,i3,
     1 ') atom ',a4,' absent')
      stop
      endif
      nat(m,n,1) =i4
      nat(m,n,2) =i2
      nat(m,n,3) =i1
      nat(m,n,15)=i3
c----------------------------------------------ls fit of standard base
      if(fit.and.li(m,n).ge.-1) then
      nind(1)=i1
      call lsfit(lequ,kb,nind,dcor,rms,key)
      if(n.eq.1.and.m.eq.ng(1)) write(6,50)
50    format(/2x,'Least squares fitting of standard bases ...',
     1     //3x,'Str',3x,'Pos',2x,'Base',10x,'Rms (ang)'/)
      if(key.eq.0) then
      write(6,60) n,m,munit(is),nunit(is),rms
60    format(2x,i3,' : ',i3,')  ',a4,1x,i3,4x,f7.3)
         else
         write(6,62) n,m,munit(is),nunit(is)
62       format(2x,i3,' : ',i3,')  ',a4,1x,i3,4x,
     1   'error ... fit impossible')
         goto 64
         endif
      i1=1
      if(purine) then
      i2=10
      i3=5
      else
      if(lsav.lt.9) then
      i2=2
      i3=3
      else
      i2=6
      i3=5
      endif
      endif
      x0=dcor(i2,1)
      y0=dcor(i2,2)
      z0=dcor(i2,3)
      ax=dcor(i1,1)-x0
      ay=dcor(i1,2)-y0
      az=dcor(i1,3)-z0
      cx=dcor(i3,1)-x0
      cy=dcor(i3,2)-y0
      cz=dcor(i3,3)-z0
      else
64    x0=corm(i2,1)
      y0=corm(i2,2)
      z0=corm(i2,3)
      ax=corm(i1,1)-x0
      ay=corm(i1,2)-y0
      az=corm(i1,3)-z0
      cx=corm(i3,1)-x0
      cy=corm(i3,2)-y0
      cz=corm(i3,3)-z0
      endif
c----------------------------------------------------------find normal
      rx=ay*cz-az*cy
      ry=az*cx-ax*cz
      rz=ax*cy-ay*cx
      r=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/r
      ry=ry/r
      rz=rz/r
      rex(m,3,n)=rx
      rey(m,3,n)=ry
      rez(m,3,n)=rz
      ra=sqrt(ax*ax+ay*ay+az*az)
c--------------------------------------------construct reference point
      ca=cos(cdr*(th1))
      sa=sin(cdr*(th1))
      fac=dis/ra
      xx=ax*fac
      yy=ay*fac
      zz=az*fac
      rex(m,4,n)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1          (rx*rz*(1-ca)+ry*sa)*zz+x0
      rey(m,4,n)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1          (ry*rz*(1-ca)-rx*sa)*zz+y0
      rez(m,4,n)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1          (rz*rz+(1-rz*rz)*ca)*zz+z0
c--------------------------------------------------construct perp dyad
      cb=cos(cdr*(th2))
      sb=sin(cdr*(th2))
      xx=ax/ra
      yy=ay/ra
      zz=az/ra
      dx=(rx*rx+(1-rx*rx)*cb)*xx+(rx*ry*(1-cb)-rz*sb)*yy+
     1  (rx*rz*(1-cb)+ry*sb)*zz
      dy=(rx*ry*(1-cb)+rz*sb)*xx+(ry*ry+(1-ry*ry)*cb)*yy+
     1  (ry*rz*(1-cb)-rx*sb)*zz
      dz=(rx*rz*(1-cb)-ry*sb)*xx+(ry*rz*(1-cb)+rx*sb)*yy+
     1  (rz*rz+(1-rz*rz)*cb)*zz
      rex(m,2,n)=dx
      rey(m,2,n)=dy
      rez(m,2,n)=dz
c-------------------------------------------------------construct dyad
      rex(m,1,n)=dy*rz-dz*ry
      rey(m,1,n)=dz*rx-dx*rz
      rez(m,1,n)=dx*ry-dy*rx
20    enddo
c--------------------------------------------------------end extension
      if(ends) then
      id=1
      if(n.gt.1) id=-1
      do m=1,2
      l=ng(n)-1
      if(m.eq.2) l=nr(n)+1
      xdi=hel(l,1,n)
      ydi=hel(l,2,n)
      cln=hel(l,4,n)
      tip=hel(l,5,n)
      ct=cos(cdr*(cln))
      st=sin(cdr*(cln))
      cp=cos(cdr*(   tip))
      sp=sin(cdr*(id*tip))
      rx3= sp*id
      ry3=-st*cp*id
      rz3= ct*cp*id
      rx2= 0.
      ry2=ct*id
      rz2=st*id
      rx1=ry2*rz3-rz2*ry3
      ry1=rz2*rx3-rx2*rz3
      rz1=rx2*ry3-ry2*rx3
      if(comb) then
      efd(1,m,n)=rz1
      efd(2,m,n)=rz2
      efd(3,m,n)=rz3
      else
      efd(1,m,n)=rz1*id
      efd(2,m,n)=rz2*id
      efd(3,m,n)=rz3*id
      endif
      efc(1,m,n)=-rx1*xdi-ry1*ydi
      efc(2,m,n)=-rx2*xdi-ry2*ydi
      efc(3,m,n)=-rx3*xdi-ry3*ydi
      enddo
      endif
c---------------------------------------------------------------------
      enddo
      return
      end
