      subroutine analy
      include 'curves_data.inc'
      character*1 na
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly
      integer*4 break,spline
      dimension w(2*n3*(4*n3+13)),var(4*n3),grm(4*n3),scale(4*n3)
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
      external move
      call calc
      call deriv
      call grads
      write(6,8) sum,scp(1)*10,scp(2),scp(3)*10,scp(4)
8     format(/2x,'FIRST SUM= ',f8.3,' CPTS: ',4f8.3)
      if(.not.mini) return
c---------------------------------------------------------minimisation
      write(6,5) acc,maxn,nvar
5     format(/2x,'MINIMISATION: ACC = ',e10.3,' MAXN= ',i4,' NVAR= ',i4)
      summ=sum
      m=0
      iup=ien
      if(line) iup=1
      do i=ist,iup
      nch=nl
      if(comb) nch=iact(i)
      do l=1,4
      m=m+1
      if(l.le.2) then
      var(m)=hel(i,l,nch)
      grm(m)=gra(m)
      scale(m)=0.5
      else
      var(m)=hel(i,l+1,nch)
      grm(m)=gra(m)*cdr
      scale(m)=1.5
      endif
      enddo
      enddo
         if(line.and.break.gt.0) then
         do l=1,4
         m=l+4
         if(l.le.2) then
         var(m)=hel(break,l,nch)
         grm(m)=gra(m)
         scale(m)=0.5
         else
         var(m)=hel(break,l+1,nch)
         grm(m)=gra(m)*cdr
         scale(m)=1.5
         endif
         enddo
         endif
      nd=1+(nvar*(nvar+1))/2
      nw=nd+nvar
      nxa=nw+nvar
      nga=nxa+nvar
      nxb=nga+nvar
      ngb=nxb+nvar
      call minfor(move,nvar,var,summ,grm,scale,acc,w,w(nd),w(nw),
     1 w(nxa),w(nga),w(nxb),w(ngb),maxn)
      write(6,10) sum,scp(1)*10,scp(2),scp(3)*10,scp(4)
10    format(/2x,'FINAL SUM= ',f8.3,' CPTS: ',4f8.3/)
      write(6,20) (grm(i),i=1,nvar)
20    format(2x,'GRA=',8e9.2)
      call title('A','Global axis parameters',22)
      write(6,*) ' '
      do i=ist,ien
      x=hx(i)
      y=hy(i)
      z=hz(i)
      if(comb) then
      x=ox(i)
      y=oy(i)
      z=oz(i)
      endif
      if(i.lt.ien) then
      write(6,101) i,ux(i),uy(i),uz(i),x,y,z,dif(i)
101   format(2x,i3,') U: ',3f8.3,'  P: ',3f8.3,:,'  D: ',f8.3)
      else
      write(6,101) i,ux(i),uy(i),uz(i),x,y,z
      endif
      enddo
      return
      end
