      subroutine gradt
      include 'curves_data.inc'
      character*1 na
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      integer*4 break,spline
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
      write(6,5)
5     format(/10x,'---- GRADTEST ----'/)
      call calc
      call deriv
      call grads
      ref=sum
      write(6,6) sum
6     format(/2x,'REFERENCE SUM= ',e10.3,
     1 //19x,'Num ',5x,'Anal',5x,'Diff')
c---------------------------------------------------test var gradients
      m=0
      ido=ist
      iup=ien
      if(line) iup=1
15    do i=ido,iup
      nch=nl
      if(comb) nch=iact(i)
      do l=1,5
      if(l.eq.3) goto 200
      hel(i,l,nch)=hel(i,l,nch)+deltv
      call calc
      der=(sum-ref)/deltv
      m=m+1
      grd=gra(m)
      if(l.gt.3) grd=grd*cdr
      write(6,10) i,l,der,grd,der-grd
10    format(2x,'VAR ',i3,i3,')  ',3f9.4)
      hel(i,l,nch)=hel(i,l,nch)-deltv
200   enddo
      enddo
c-------------------------------------------------line+break loop back
         if(line.and.ido.eq.ist.and.break.gt.0) then
         ido=break
         iup=break
         goto 15
         endif
c---------------------------------------------------------------------
      return
      end
