      subroutine pdbout(pdb)
      include 'curves_data.inc'
      character*4 snam,sunit,pdb*32,chain(52)*1,mpdb(n1)*1
      logical*2 seen(n1)
      dimension ifr(n1),nfr(0:52),matt(7),ilis(n1)
      common/axe/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),matc(n1,7),kas,khs,kces
      data chain/
     1 'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
     1 'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
c-------------------------------------------------------find fragments
      ifr(1)=1
      seen(1)=.true.
      do i=2,kas
      ifr(i)=0
      seen(i)=.false.
      enddo
      n=1
      kfr=1
      ilis(1)=kfr
      k=0
20    k=k+1
      i=ilis(k)
      do j=1,matc(i,7)
      in=abs(matc(i,j))
      if(.not.seen(in)) then
      ifr(in)=kfr
      seen(in)=.true.
      n=n+1
      ilis(n)=in
      endif
      enddo
      if(k.lt.n) goto 20
      do i=2,kas
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
      do i=1,kas
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
      nfr(k)=kas
c--------------------------------------------------------------output brk format
      kfi=index(pdb,' ')-1
      open(unit=2,file=pdb(:kfi)//'.pdb',status='new')
      write(2,5) 'HEADER    '//pdb(:kfi)//' from Curves'
5     format(a)
      kl=0
      do k=1,kfr
      do i=nfr(k-1)+1,nfr(k)
      kl=kl+1
      write(2,10) 'ATOM',kl,snam(i),sunit(i),mpdb(i),
     1 nunis(i),cors(i,1),cors(i,2),cors(i,3)
10    format(a4,2x,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3)
      enddo
      kl=kl+1
      iend=nfr(k)
      write(2,10) 'TER ',kl,' ',sunit(iend),mpdb(iend),
     1 nunis(iend)
      enddo
      do i=1,kas
      l=matc(i,7)
      ii=i+ifr(i)-1
      do j=1,l
      jj=abs(matc(i,j))
      matt(j)=jj+ifr(jj)-1
      enddo
      write(2,11) 'CONECT',ii,(matt(j),j=1,l)
11    format(a6,5i5)
      enddo
      write(2,12) 1,0,0,0,0,0,0,0,kas,0,kas,0
12    format('MASTER',4x,12i5,/'END')
      close(2)
      return
      end
