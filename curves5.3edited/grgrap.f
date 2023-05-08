      subroutine grgrap
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,seen(n1)
      character*4 na*1,atoms(9),mnam,munit,lis*20,chain(52)*1,mpdb(n1)*1
      integer*4 spline,break,matm(1000)
      dimension ohm(2),ogm(2),hpi(2),gqi(2),dm(2),inn(2),jnn(2),ini(2),
     1 gqn(2),hpn(2),jni(2),radius(9),ilis(n1),grox(8,48,2),
     1 groy(8,48,2),groz(8,48,2),ifr(n1),nfr(0:999),matt(7)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/gro/uxb(200,0:29),uyb(200,0:29),uzb(200,0:29),cor(900,3),
     1 dya(900,3),box(0:900,4),boy(0:900,4),boz(0:900,4),ind(40,0:30),
     1 nma(n3),num,numa,nsu(4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,nbac,ends,
     1 supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      data chain/
     1 'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
     1 'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data atoms/'C1''','C2''','C3''','C4''','O1''','O3''','P',
     1 'O5''','C5'''/,radius/1.6,1.6,1.6,1.6,1.4,1.4,2.9,1.4,1.6/
      return
c-------------------------------------------------matrix for graphic
      grox(1,k,l)=corx
      groy(1,k,l)=cory
      groz(1,k,l)=corz
      if(l.eq.1)of=oa0*(nmi-n)/nmi+oa1*n/nmi
      if(l.eq.2)of=ob0*(nmi-n)/nmi+ob1*n/nmi
      if(pq+dm(l).eq.198.)then
      oh=ohm(l)
      og=ogm(l)
      if(hp.ne.99.0.and.gq.eq.99.)then
      grox(2,k,l)=corx+oh*dyax
      groy(2,k,l)=cory+oh*dyay
      groz(2,k,l)=corz+oh*dyaz
      do j=3,4
      grox(j,k,l)=box(in,1)
      groy(j,k,l)=boy(in,1)
      groz(j,k,l)=boz(in,1)
      enddo
      elseif(hp.eq.99.0.and.gq.ne.99.)then
      grox(2,k,l)=corx+og*dyax
      groy(2,k,l)=cory+og*dyay
      groz(2,k,l)=corz+og*dyaz
      do j=3,4
      grox(j,k,l)=box(jn,2)
      groy(j,k,l)=boy(jn,2)
      groz(j,k,l)=boz(jn,2)
      enddo
      else
      if(l.eq.1)indx=-1
      if(l.eq.2)indx=+1
      do j=2,4
      grox(j,k,l)=corx+dyax*indx*2
      groy(j,k,l)=cory+dyay*indx*2
      groz(j,k,l)=corz+dyaz*indx*2
      enddo
      endif
      elseif(pq+dm(l).lt.198.)then
      if(pq.lt.99.)then
      im=inn(l)
      jm=jnn(l)
      r=hpn(l)/(hpn(l)+gqn(l))
      else
      im=ini(l)
      jm=jni(l)
      r=hpi(l)/(hpi(l)+gqi(l))
      endif
      grox(2,k,l)=r*box(jm,2)+(1-r)*box(im,1)
      groy(2,k,l)=r*boy(jm,2)+(1-r)*boy(im,1)
      groz(2,k,l)=r*boz(jm,2)+(1-r)*boz(im,1)
      grox(3,k,l)=box(im,1)
      groy(3,k,l)=boy(im,1)
      groz(3,k,l)=boz(im,1)
      grox(4,k,l)=box(jm,2)
      groy(4,k,l)=boy(jm,2)
      groz(4,k,l)=boz(jm,2)
      endif
      do j=5,7
      grox(j,k,l)=corx+of*dyax+(j-6)*tex
      groy(j,k,l)=cory+of*dyay+(j-6)*tey
      groz(j,k,l)=corz+of*dyaz+(j-6)*tez
      enddo
      grox(8,k,l)=grox(2,k,l)
      groy(8,k,l)=groy(2,k,l)
      groz(8,k,l)=groz(2,k,l)
c--------------------------------------graphic: strands and axis
      k=0
      do i=0,insu
      k=k+1
      matm(k)=k
      munit(k)='BACK'
      nunit(k)=1
      corm(k,1)=box(i,1)
      corm(k,2)=boy(i,1)
      corm(k,3)=boz(i,1)
      enddo
      l=k+1
      matm(l)=10000
      do i=0,jnsu
      l=l+1
      k=k+1
      matm(l)=k
      munit(k)='BACK'
      nunit(k)=2
      corm(k,1)=box(i,2)
      corm(k,2)=boy(i,2)
      corm(k,3)=boz(i,2)
      enddo
      l=l+1
      matm(l)=10000
      do i=0,nsu(3)
      l=l+1
      k=k+1
      matm(l)=k
      munit(k)='BACK'
      nunit(k)=3
      corm(k,1)=box(i,3)
      corm(k,2)=boy(i,3)
      corm(k,3)=boz(i,3)
      enddo
      l=l+1
      matm(l)=10000
      do i=1,ind(nux,0)
      l=l+1
      k=k+1
      matm(l)=k
 
      munit(k)='AXE'
      nunit(k)=4
      do j=1,3
      corm(k,j)=cor(i,j)
      enddo
      enddo
      kam=k
 
c----------------------------------------graphic: grooves
      kmax=ind(numa,0)
      k=0
      l=0
      kh1=k
 
      do i=1,kmax
      do m=1,4
      k=k+1
      munit(k)='MIN'
      nunit(k)=1
      corm(k,1)=grox(m,i,1)
      corm(k,2)=groy(m,i,1)
      corm(k,3)=groz(m,i,1)
      enddo
      enddo
      do i=1,kmax
      do m=1,4
      k=k+1
      munit(k)='MAJ'
      nunit(k)=1
      corm(k,1)=grox(m,i,2)
      corm(k,2)=groy(m,i,2)
      corm(k,3)=groz(m,i,2)
      enddo
      enddo
      k=kh1
      do i=1,2*kmax
      matm(l+1)=10000
      matm(l+2)=k+1
      matm(l+3)=k+2
      matm(l+4)=k+3
      matm(l+5)=10000
      matm(l+6)=k+2
      matm(l+7)=k+4
      l=l+7
      k=k+4
      enddo
 
      kh2=k
 
c      do i=1,kmax
      do m=5,8
      k=k+1
      munit(k)='DMIN'
      nunit(k)=1
      corm(k,1)=grox(m,i,1)
      corm(k,2)=groy(m,i,1)
      corm(k,3)=groz(m,i,1)
      enddo
c      enddo
c      do i=1,kmax
      do m=5,8
      k=k+1
      munit(k)='DMAJ'
      nunit(k)=1
      corm(k,1)=grox(m,i,2)
      corm(k,2)=groy(m,i,2)
      corm(k,3)=groz(m,i,2)
      enddo
c      enddo
      kam=k
      k=kh2
      do i=1,2
      matm(l+1)=10000
      matm(l+2)=k+1
      matm(l+3)=k+2
      matm(l+4)=k+3
      matm(l+5)=10000
      matm(l+6)=k+2
      matm(l+7)=k+4
      l=l+7
      k=k+4
      enddo
 
      khm=l
c      kam=k
 
         do i=1,kam
         do j=1,7
         matd(i,j)=0
         enddo
         enddo
      do i=1,khm-1
      k=matm(i)
      kn=abs(k)
      if(k.lt.10000) then
      ik=matd(kn,7)
      kp=matm(i+1)
      kpn=abs(kp)
      if(kp.lt.10000) then
      ikp=matd(kpn,7)
      ik=ik+1
      matd(kn,ik)=kp
      matd(kn,7)=ik
      ikp=ikp+1
      matd(kpn,ikp)=sign(k,kp)
      matd(kpn,7)=ikp
      endif
      if(kp.eq.20000) then
      kpp=matm(i+2)
      kppn=abs(kpp)
      ikpp=matd(kppn,7)
      ik=ik+1
      matd(kn,ik)=kpp
      matd(kn,7)=ik
      ikpp=ikpp+1
      matd(kppn,ikpp)=sign(k,kpp)
      matd(kppn,7)=ikpp
      endif
      endif
      enddo
c-------------------------------------------------------find fragments
      ifr(1)=1
      seen(1)=.true.
      do i=2,kam
      ifr(i)=0
      seen(i)=.false.
      enddo
      n=1
      kfr=1
      ilis(1)=kfr
      k=0
20    k=k+1
      i=ilis(k)
      do j=1,matd(i,7)
      in=abs(matd(i,j))
      if(.not.seen(in)) then
      ifr(in)=kfr
      seen(in)=.true.
      n=n+1
      ilis(n)=in
      endif
      enddo
      if(k.lt.n) goto 20
      do i=2,kam
      if(ifr(i).eq.0) then
      kfr=kfr+1
      ifr(i)=kfr
      seen(i)=.true.
      n=n+1
      ilis(n)=i
      goto 20
      endif
      enddo
c-------------------------------------------------------------------add chain id
      k=0
      nfr(0)=0
      ifrag=ifr(1)
      do i=1,kam
      ir=ifr(i)
      if(kfr.eq.1) then
      mpdb(i)=' '
      else
      mpdb(i)=chain(ir)
      endif
      if(ir.ne.ifrag) then
      k=k+1
      nfr(k)=i-1
      ifrag=ir
      endif
      enddo
      k=k+1
      nfr(k)=kam
c--------------------------------------------------------------output brk format
      kfi=index(lis,' ')-1
      open(unit=2,file='groove.pdb',status='unknown')
      write(2,5) 'HEADER Cur4 groove graphic'
5     format(a)
 
      kl=0
      do k=1,kfr
      do i=nfr(k-1)+1,nfr(k)
      kl=kl+1
      write(2,10) 'ATOM',kl,' ',munit(i),mpdb(i),
     1 nunit(i),corm(i,1),corm(i,2),corm(i,3)
10    format(a4,2x,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3)
      enddo
      kl=kl+1
      iend=nfr(k)
      write(2,10) 'TER ',kl,' ',munit(iend),mpdb(iend),
     1 nunit(iend)
      enddo
 
      do i=1,kam
      l=matd(i,7)
      ii=i+ifr(i)-1
      do j=1,l
      jj=abs(matd(i,j))
      matt(j)=jj+ifr(jj)-1
      enddo
      write(2,41) 'CONECT',ii,(matt(j),j=1,l)
41    format(a6,5i5)
      enddo
      close(2)
      return
      end
