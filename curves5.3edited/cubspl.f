      subroutine cubspl
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      character*4 snam,sunit,na*1
      integer*4 spline,break
      dimension f(4,3),p(3,2),u(3,2)
      common/axe/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),matc(n1,7),kas,khs,kces
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
      k=0
      l=0
c--------------------------------------------------------------strands
      do kk=1,nst
      kces=kces+1
      ks=k+1
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
c------------------------------------------------------------residues
      do i=ist,ien-1
c--------------------------------------------arrow head at break point
         if(i+1.eq.break) then
         l=l+1
         mats(l)=10000
         if(idr(kk).eq.1) then
         kbr=k
         i1=break-1
         i2=i1-1
            else
            kbr=k+1
            i1=break
            i2=i1+1
            endif
         ax=hho(i1,1,kk)-hho(i2,1,kk)
         ay=hho(i1,2,kk)-hho(i2,2,kk)
         az=hho(i1,3,kk)-hho(i2,3,kk)
         r=sqrt(ax*ax+ay*ay+az*az)
         ax=ax/r
         ay=ay/r
         az=az/r
         kc=kk
         if(comb) kc=iact(i1)
         dot=ax*rex(i1,2,kc)+ay*rey(i1,2,kc)+az*rez(i1,2,kc)
         bx=rex(i1,2,kc)-ax*dot
         by=rey(i1,2,kc)-ay*dot
         bz=rez(i1,2,kc)-az*dot
         r=sqrt(bx*bx+by*by+bz*bz)
         bx=bx/r
         by=by/r
         bz=bz/r
         do j=1,2
         k=k+1
         is=-1
         if(j.eq.2) is=1
         snam(k)='A'
         sunit(k)='AXIS'
         nunis(k)=kk
         cors(k,1)=hho(i1,1,kk)-ax*0.55+is*bx*0.18
         cors(k,2)=hho(i1,2,kk)-ay*0.55+is*by*0.18
         cors(k,3)=hho(i1,3,kk)-az*0.55+is*bz*0.18
         enddo
         mats(l+1)=10000
         mats(l+2)=kbr
         mats(l+3)=k-1
         mats(l+4)=k
         mats(l+5)=20000
         mats(l+6)=kbr
         mats(l+7)=10000
         l=l+7
      goto 100
      endif
c--------------------------------------------------------------------
      do j=1,2
      ij=i+j-1
      p(1,j)=hho(ij,1,kk)
      p(2,j)=hho(ij,2,kk)
      p(3,j)=hho(ij,3,kk)
      u(1,j)=uho(ij,1,kk)
      u(2,j)=uho(ij,2,kk)
      u(3,j)=uho(ij,3,kk)
      enddo
c---------------------------------------------------------------cubspl
      gl=sqrt((p(1,2)-p(1,1))**2+(p(2,2)-p(2,1))**2+(p(3,2)-p(3,1))**2)
      do j=1,3
      f(1,j)=p(j,1)
      f(2,j)=u(j,1)
      f(3,j)=( 3/gl**2)*(p(j,2)-p(j,1))-(1/gl   )*(u(j,2)+2*u(j,1))
      f(4,j)=(-2/gl**3)*(p(j,2)-p(j,1))+(1/gl**2)*(u(j,2)+  u(j,1))
      enddo
      do m=0,spline+1
      t=gl*m/(spline+1)
      k=k+1
      l=l+1
      mats(l)=k
      snam(k)='A'
      sunit(k)='AXIS'
      nunis(k)=kk
      do j=1,3
      cors(k,j)=f(1,j)+f(2,j)*t+f(3,j)*t**2+f(4,j)*t**3
      enddo
      enddo
100   enddo
      khs=k
c--------------------------------------------------------3' arrow head
      if(idr(kk).eq.1) then
      i1=ien
      i2=ien-1
      kbr=khs
      else
      i1=ist
      i2=ist+1
      kbr=ks
      endif
      ax=hho(i1,1,kk)-hho(i2,1,kk)
      ay=hho(i1,2,kk)-hho(i2,2,kk)
      az=hho(i1,3,kk)-hho(i2,3,kk)
      r=sqrt(ax*ax+ay*ay+az*az)
      ax=ax/r
      ay=ay/r
      az=az/r
      kc=kk
      if(comb) kc=iact(i1)
      dot=ax*rex(i1,2,kc)+ay*rey(i1,2,kc)+az*rez(i1,2,kc)
      bx=rex(i1,2,kc)-ax*dot
      by=rey(i1,2,kc)-ay*dot
      bz=rez(i1,2,kc)-az*dot
      r=sqrt(bx*bx+by*by+bz*bz)
      bx=bx/r
      by=by/r
      bz=bz/r
      do j=1,2
      k=k+1
      is=-1
      if(j.eq.2) is=1
      snam(k)='A'
      sunit(k)='AXIS'
      nunis(k)=kk
      cors(k,1)=hho(i1,1,kk)-ax*0.55+is*bx*0.18
      cors(k,2)=hho(i1,2,kk)-ay*0.55+is*by*0.18
      cors(k,3)=hho(i1,3,kk)-az*0.55+is*bz*0.18
      enddo
      mats(l+1)=10000
      mats(l+2)=kbr
      mats(l+3)=k-1
      mats(l+4)=k
      mats(l+5)=20000
      mats(l+6)=kbr
      mats(l+7)=10000
      khs=l+7
      kas=k
      enddo
      return
      end
