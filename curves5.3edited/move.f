      subroutine move(var,summ,grm)
      include 'curves_data.inc'
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,na*1
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      integer*4 break,spline
      dimension var(4*n3),grm(4*n3)
      common/cha/file,lis,dna,axin,axout,daf,pdb,mcode
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
      mof=4*(break-1)
      m=0
      iup=ien
      if(line) iup=1
      do i=ist,iup
      nch=nl
      if(comb) nch=iact(i)
      hel(i,1,nch)=var(m+1)
      hel(i,2,nch)=var(m+2)
      hel(i,4,nch)=var(m+3)
      hel(i,5,nch)=var(m+4)
      m=m+4
      enddo
         if(line.and.break.gt.0) then
         hel(break,1,nch)=var(5)
         hel(break,2,nch)=var(6)
         hel(break,4,nch)=var(7)
         hel(break,5,nch)=var(8)
         endif
      call calc
      call deriv
      call grads
      if(icyc.eq.0) then
      ref=sum
      write(6,*) ' '
      endif
      icyc=icyc+1
      delta=sum-ref
      ref =sum
      summ=sum
      m=0
      do i=ist,ien
      do l=1,4
      m=m+1
      if(l.le.2) then
      grm(m)=gra(m)
      else
      grm(m)=gra(m)*cdr
      endif
      enddo
      enddo
c---------------------------------------------------------------output
      if(supp) then
      write(6,10) icyc,sum,delta
      else
      write(6,*) ' '
      write(6,10) icyc,sum,delta
      write(6,*) ' '
      write(6,30) (grm(i),i=1,nvar)
      write(6,*) ' '
      endif
10    format(2x,'STEP ',i4,' SUM= ',f8.3,' DEL= ',e10.3)
30    format(2x,'GRA=',8e9.2)
      return
      end
