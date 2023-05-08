      subroutine backbo
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 flag,space,grv,old,axonly
      character*4 mnam,munit,na*1,nam,atom(15)
      integer*4 break,spline
      dimension set(8),pset(2),ita(16,4),store(16),nroot(15)
      common/bak/tor(0:n3,4,13),val(0:n3,4,2),suga(0:n3,4,2),
     1 nat(n3,4,15),flag(n3,4)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      data set/106.2540,100.7339,102.3817,78.1441,39.4467,-34.6199,
     1 -133.1470,-156.8818/,pset/39.70,154.76/
      data atom/'C6C8','C2N9','C1''','C2''','C3''','C4''','O4''',
     1 'O3''','P','O5''','C5''','C4''','C5''','O5''','C2C4'/
      data nroot/0,0,0,3,4,5,3,5,8,9,10,11,6,13,0/
      data ita/7,3,4,1,7,3,4,5,6,6, 5,14,13, 8, 9,15,
     1         3,4,5,2,3,4,5,6,7,5, 8,13, 6, 9,10, 2,
     1         4,5,6,3,4,5,6,7,3,8, 9, 6, 5,10,11, 3,
     1         0,0,0,7,5,6,7,3,4,9,10, 5, 8,11,12, 7/
      space=.true.
      do k=1,nst
      do i=ng(k),nr(k)
      if(li(i,k).ge.-2) then
      flag(i,k)=.false.
c-------------------------------------------------find necessary atoms
      n=nat(i,k,3)
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'O1'''.or.nam.eq.'O1*'.or.
     1  nam.eq.'O4'''.or.nam.eq.'O4*') nat(i,k,7)=m
      if(nam.eq.'C2'''.or.nam.eq.'C2*') nat(i,k,4)=m
      enddo
      n=nat(i,k,4)
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'C3'''.or.nam.eq.'C3*') nat(i,k,5)=m
      enddo
      n=nat(i,k,5)
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'C4'''.or.nam.eq.'C4*') nat(i,k,6)=m
      if(nam.eq.'O3'''.or.nam.eq.'O3*') nat(i,k,8)=m
      enddo
      n=nat(i,k,6)
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'C5'''.or.nam.eq.'C5*') nat(i,k,13)=m
      enddo
      n=nat(i,k,13)
      if(n.gt.0) then
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'O5'''.or.nam.eq.'O5*') nat(i,k,14)=m
      enddo
      endif
      if((idr(k).eq. 1.and.i.lt.nr(k)).or.
     1  (idr(k).eq.-1.and.i.gt.ng(k))) then
      n=nat(i,k,8)
      if(n.gt.0) then
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'P') nat(i,k,9)=m
      enddo
      endif
      n=nat(i,k,9)
      if(n.gt.0) then
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'O5'''.or.nam.eq.'O5*') nat(i,k,10)=m
      enddo
      endif
      n=nat(i,k,10)
      if(n.gt.0) then
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'C5'''.or.nam.eq.'C5*') nat(i,k,11)=m
      enddo
      endif
      n=nat(i,k,11)
      if(n.gt.0) then
      do l=1,matd(n,7)
      m=abs(matd(n,l))
      nam=mnam(m)
      if(nam.eq.'C4'''.or.nam.eq.'C4*') nat(i,k,12)=m
      enddo
      endif
      endif
c---------------------------------------------------------------check
      do l=1,15
      if(li(i,k).eq.-2.and.(l.lt.3.or.l.eq.15)) goto 100
      if(idr(k).eq. 1.and.i.eq.ng(k).and.l.eq.14) goto 100
      if(idr(k).eq.-1.and.i.eq.nr(k).and.l.eq.14) goto 100
      if(idr(k).eq. 1.and.i.eq.nr(k).and.l.ge.9.and.l.le.12) goto 100
      if(idr(k).eq.-1.and.i.eq.ng(k).and.l.ge.9.and.l.le.12) goto 100
      if(nat(i,k,l).eq.0) then
      if(space) then
      write(6,*) ' '
      space=.false.
      endif
      m=ni(i,k)
      j=ncen(m)
      lr=nroot(l)
      if(lr.gt.0) then
      jr=nat(i,k,lr)
      if(jr.gt.0) write(6,10) k,i,munit(j),nunit(j),atom(lr),jr,atom(l)
10    format(2x,'Strand: ',i2,' Unit: ',i3,' (',a4,i3,
     1 ') linkage from atom ',a4,'(',i4,') to ',a4,' absent')
      else
      write(6,20) k,i,munit(j),nunit(j),atom(l)
20    format(2x,'Strand: ',i2,' Unit: ',i3,' (',a4,i3,
     1 ')  atom ',a4,' absent')
      endif
      endif
100   enddo
c---------------------------------------------find angles and torsions
      do l=1,16
      store(l)=999.
      enddo
      isg=0
      do l=1,16
      i1=nat(i,k,ita(l,1))
      i2=nat(i,k,ita(l,2))
      i3=nat(i,k,ita(l,3))
      i4=0
      if(ita(l,4).gt.0) i4=nat(i,k,ita(l,4))
         if(i1.gt.0.and.i2.gt.0.and.i3.gt.0) then
         if(l.eq.11) val(i,k,1)=torp(i1,i2,i3,0)
         if(l.eq.14) val(i,k,2)=torp(i1,i2,i3,0)
            if(l.le.3.or.i4.gt.0) then
            store(l)=torp(i1,i2,i3,i4)
            if(l.ge.5.and.l.le.9) isg=isg+1
            endif
         endif
      enddo
      if(isg.lt.5) flag(i,k)=.true.
      do l=1,6
      tor(i,k,l)=store(l)
      enddo
      do l=10,16
      tor(i,k,l-3)=store(l)
      enddo
c----------------------------------------------------find sugar pucker
      a=0.
      b=0.
      do l=1,5
      j=l+5
      if(l.eq.5) j=5
      a=a+store(j)*cos(cdr*(144.*(l-1)))
      b=b+store(j)*sin(cdr*(144.*(l-1)))
      enddo
      a= a*2./5.
      b=-b*2./5.
      amp=sqrt(a*a+b*b)
      if(amp.gt.0.) then
      cp=a/amp
      sp=b/amp
      if(abs(cp).gt.1.) cp=sign(1.d0,cp)
      pha=acos(cp)*crd
      if(sp.lt.0.) pha=360.-pha
      else
      pha=0.
      endif
      suga(i,k,1)=amp
      suga(i,k,2)=pha
      endif
      enddo
      enddo
c-----------------------------------------------------------------ends
      if(ends) then
      do k=1,nst
      m=ng(k)-1
      n=nr(k)+1
      do l=1,8
      tor(m,k,l)=set(l)
      tor(n,k,l)=set(l)
      enddo
      if(idr(k).eq.1) then
      do l=7,8
      tor(n,k,l)=0.
      tor(n-1,k,l)=set(l)
      enddo
      else
      do l=7,8
      tor(m,k,l)=0.
      tor(m+1,k,l)=set(l)
      enddo
      endif
      do l=1,2
      suga(m,k,l)=pset(l)
      suga(n,k,l)=pset(l)
      enddo
      enddo
      endif
      return
      end
