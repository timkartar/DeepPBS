      program aacur
c*********************************************************************
c                                                                    *
c      --------------------- Curves 5.3 --------------------         *
c                                                                    *
c     helical and kink parameters for irregular nucleic acids        *
c                                                                    *
c     updated to take account of the cambridge definitions of        *
c     helicoidal parameters and to include the calculation of        *
c     local inter-base and inter-base pair parameters.               *
c                                                                    *
c     updated to include treatment of dinucleotide repeats           *
c                                                                    *
c     updated to include backbone ribbons                            *
c                                                                    *
c     updated to print dif (dinucleotide irregularity function)      *
c                                                                    *
c     updated to treat 1-4 strands of unequal lengths                *
c     and an optional break point                                    *
c                                                                    *
c     updated to use include file for variable dimensions            *
c     and to use fortran namelist interpreter                        *
c                                                                    *
c     updated to include PDB format graphic output                   *
c     updated to include groove geometry                             *
c     updated to include superhelical axis input                     *
c     updated to new .axe file format                                *
c     updated to treat mispairs, abasic sites and bulges             *
c     updated to allow for line=.t. with a break                     *
c     updated to include new bending defn. and helical axis fit      *
c     updated to avoid crash after problem with lsfit                *
c                                                                    *
c     Copyright:  R.Lavery & H.Sklenar              Jan 1998         *
c*********************************************************************
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 lkink,flag,grv,old,axonly
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,na*1,mnam,munit,ban,dir(-1:1)*5,
     1 base*1,snam,sunit,lini*80
      integer*4 spline,break
      dimension xdi(n3),ydi(n3),cln(n3),tip(n3)
      common/axe/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),matc(n1,7),kas,khs,kces
      common/bak/tor(0:n3,4,13),val(0:n3,4,2),suga(0:n3,4,2),
     1 nat(n3,4,15),flag(n3,4)
      common/cha/file,lis,dna,axin,axout,daf,pdb,mcode
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/der/udx(n3,4),udy(n3,4),udz(n3,4),pdx(n3,4),pdy(n3,4),
     1 pdz(n3,4),upx(n3),upy(n3),upz(n3),usx(n3),usy(n3),usz(n3),
     1 umx(n3),umy(n3),umz(n3),qx(n3),qy(n3),qz(n3)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/eds/dif(0:n3),efd(3,2,4),efc(3,2,4),ist,ien,iste,iene
      common/gro/uxb(200,0:29),uyb(200,0:29),uzb(200,0:29),cor(900,3),
     1 dya(900,3),box(0:900,4),boy(0:900,4),boz(0:900,4),ind(40,0:30),
     1 nma(n3),num,numa,nsu(4)
      common/hlx/pts(n3,3),ss(3),c(3),rep(3),rad,pit,pha,tbn,dan,isen
      common/hol/uho(0:n3,3,4),hho(0:n3,3,4),vkin(0:n3,7,4),bend(n3,4),
     1 hold(0:n3,2,4),vold(0:n3,4),pal(0:n3,6,4),pab(0:n3,6,4),inv(4)
      common/lsf/bref(10,3,5),th1,th2,dis,rs2,ibref(5),iequ(9),
     1 ibd(20,2),ban(9),base(9)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      common/min/gra(4*n3),sum,ref,scp(4),nvar,icyc
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      data dir/'3''-5''',' ','5''-3'''/
c---------------------------------------------------------initial data
      file=' '
      lis=' '
      dna=' '
      axin=' '
      axout=' '
      daf=' '
      pdb=' '
      acc=1.d-6
      wid=0.75
      maxn=500
      ior=0
      spline=3
      break=-1
      ibond=0
      nlevel=3
      nbac=7
      fit=.false.
      ends=.false.
      supp=.true.
      comb=.false.
      test=.false.
      mini=.true.
      rest=.false.
      line=.false.
      zaxe=.false.
      dinu=.false.
      grv=.false.
      old=.true.
      axonly=.false.
c-----------------------------------------------------------input data
      call nml
      if(ibond.gt.0) then
      do i=1,ibond
      read(5,*) ibd(i,1),ibd(i,2)
      enddo
      endif
      if(ends.and.line) then
      write(6,*) '  ---- Ends not allowed with LINE ----'
      stop
      endif
      if(line.and..not.mini) mini=.true.
      if(zaxe.and.mini) mini=.false.
      if(ends.and.zaxe) then
      write(6,*) '  ---- Ends not allowed with ZAXE ----'
      stop
      endif
      do i=1,4
      nu(i)=0
      enddo
      read(5,*) nst,nu(1),nu(2),nu(3),nu(4)
      if(ends.and.nst.gt.2) then
      write(6,*) '  ---- Ends not allowed with more than 2 strands ----'
      stop
      endif
      nux=0
      l=0
      do i=1,4
      if(nu(i).eq.0) goto 6
      if(abs(nu(i)).gt.nux) nux=abs(nu(i))
      l=l+1
      enddo
6     if(l.ne.nst) then
      write(6,*) '  ---- Error in number of strands ----'
      stop
      endif
      if(nst.eq.1.and.comb) comb=.false.
      nt=0
      do k=1,nst
      idr(k)=sign(1,nu(k))
      read(5,*) (ni(j,k),j=1,nux)
      ng(k)=0
      nr(k)=0
      l=0
      do j=1,nux
      m=ni(j,k)
      li(j,k)=1
         if(m.lt.0) then
         m=-m
         ni(j,k)=m
         li(j,k)=-1
         endif
      if(m.gt.0) then
      l=l+1
      nr(k)=j
      if(ng(k).eq.0) ng(k)=j
      endif
      enddo
         do j=1,nux
         if(ni(j,k).eq.0) then
         li(j,k)=-4
         if(j.lt.ng(k).or.j.gt.nr(k)) li(j,k)=-3
            if(li(j,k).eq.-4.and..not.comb) then
            write(6,*) '  ---- No internal gaps unless comb=.t. ----'
            stop
            endif
         endif
         enddo
      nu(k)=l
      nt=nt+nu(k)
      if(l.ne.nux.and.ends) then
      write(6,*)
     1 '  ---- ENDS only for strands of equal length ----'
      stop
      endif
      enddo
      if(nu(1).lt.nux.and.axout.ne.' ') then
      write(6,*)
     1 '  ---- AXOUT only if 1st strand has no zeros ----'
      stop
      endif
      if(nst.gt.1.and..not.comb.and.axout.ne.' ') then
      write(6,*)
     1 '  ---- AXOUT only if COMB for multiple strands ----'
      stop
      endif
         do i=1,nux
         do k=1,nst
         if(li(j,k).ge.-1) goto 120
         enddo
         write(6,*) '  ---- Level with no bases not allowed ----'
         stop
 120     enddo
c-----------------------------------------------------------input xytp
      if(axin.ne.' ') rest=.true.
      inpv=1
      if(.not.rest) then
      if(nst.gt.1.and..not.comb) inpv=nst
      else
      inpv=nt
      if(comb) inpv=nux
      endif
      if(axin.eq.' ') then
      do i=1,inpv
      read(5,*,err=150) xdi(i),ydi(i),cln(i),tip(i)
      enddo
      goto 160
150   write(6,*) '  ---- Insufficient XYTP I/P lines ----'
      stop
      else
      kfi=index(axin,' ')-1
      open(unit=1,file=axin(:kfi)//'.axe',status='old')
      read(1,5) lini
5     format(a)
      do i=1,inpv
      read(1,35) xdi(i),ydi(i),cln(i),tip(i)
      read(1,*)
      enddo
35    format(2x,2f9.4,9x,2f9.4)
      if(index(lini,'#').eq.0) then
      do i=1,inpv
      xdi(i)=-xdi(i)
      tip(i)=-tip(i)
      enddo
      endif
      close(unit=1)
      endif
c------------------------------------------------------------end input
160   if(ends) then
      read(5,*) (hel(0,j,1),j=1,6)
      read(5,*) (hel(nux+1,j,1),j=1,6)
      munit(0)='VIRT'
      nunit(0)=0
      if(nst.gt.1) then
      do k=2,nst
      do j=1,6
      hel(0,j,k)=hel(0,j,1)
      hel(nux+1,j,k)=hel(nux+1,j,1)
      enddo
      enddo
      endif
      do k=1,nst
      li(0,k)=1
      li(nux+1,k)=1
      enddo
      endif
c------------------------------------------------------input structure
      call input
      call locate
      call backbo
      if(break.gt.0) write(6,18) break-1,break
18    format(/2x,'Break point between levels ',i2,' and ',i2)
      write(6,20) nst,nt,kam,kcen
20    format(
     1 /2x,'Strand= ',i4,' Nucleo= ',i4,' Atoms = ',i4,' Units = ',i4)
c----------------------------------------------------------choose case
      ns=1
      if(comb) then
      ns=nst
      nst=1
      endif
c------------------------------------------------------------z as axis
      if(zaxe) then
      if(comb) then
      write(6,32) nux
      do i=1,nux
      do j=1,ns
      if(li(i,j).gt.0) then
      iact(i)=j
      goto 110
      endif
      enddo
        do j=1,ns
        if(li(i,j).eq.-1) then
        iact(i)=j
        goto 110
        endif
        enddo
      write(6,*) '  ---- Level with no bases not allowed ----'
      stop
110   enddo
      endif
      do nl=1,nst
      if(.not.comb) then
      ist=ng(nl)
      ien=nr(nl)
      iste=ist
      iene=ien
      else
      ist=1
      ien=nux
      iste=ist
      iene=ien
      endif
      do i=ist,ien
      uho(i,1,nl)=0.
      uho(i,2,nl)=0.
      uho(i,3,nl)=idr(nl)
      hho(i,1,nl)=0.
      hho(i,2,nl)=0.
      nn=nl
      if(li(i,nl).lt.0) nn=iact(i)
      hho(i,3,nl)=rez(i,4,nn)
      ox(i)=0.
      oy(i)=0.
      oz(i)=rez(i,4,nn)
      enddo
      enddo
      do nl=1,max(ns,nst)
      write(6,30) nl,nu(nl),dir(idr(nl)),(na(i,nl),i=1,nux)
      enddo
      goto 200
      endif
c----------------------------------------------write helical variables
      write(6,*) ' '
      do k=1,inpv
      write(6,25) k,xdi(k),ydi(k),cln(k),tip(k)
      enddo
25    format(2x,'Input ',i3,') Xdisp= ',f7.2,' Ydisp= ',f7.2,
     1 ' Inclin= ',f7.2,' Tip= ',f7.2)
c-----------------------------------------------------------write ends
      if(ends) then
      write(6,*) ' '
      write(6,26) 0,(hel(0,j,1),j=1,6)
      write(6,26) nux+1,(hel(nux+1,j,1),j=1,6)
26    format(2x,'ENDS  ',i3,') Xd=',f7.2,' Yd=',f7.2,' Ri=',f7.2,
     1 ' In=',f7.2,' Tp=',f7.2,' Tw=',f7.2)
      endif
c------------------------------------------------------setup variables
      m=0
      do k=1,max(ns,nst)
      if(comb) then
      ist=1
      ien=nux
      else
      ist=ng(k)
      ien=nr(k)
      endif
      do i=ist,ien
      m=m+1
      if(comb) then
      l=1
      if(rest) l=i
      else
      l=k
      if(rest) l=m
      endif
      hel(i,1,k)=xdi(l)
      hel(i,2,k)=ydi(l)
      hel(i,4,k)=cln(l)
      hel(i,5,k)=tip(l)
      enddo
      enddo
c------------------------------------------------------for each strand
      do nl=1,nst
      n=nu(nl)
      if(comb) n=nux
      nvar=4*n
      if(line) nvar=4
      if(line.and.break.gt.0) nvar=8
      if(.not.comb) then
      write(6,*)
      write(6,30) nl,nu(nl),dir(idr(nl)),(na(i,nl),i=1,nux)
30    format(2x,'Strand ',i2,' has ',i3,' bases (',a5,'): ',150a1,
     1 :,/35x,150a1)
      ist=ng(nl)
      ien=nr(nl)
      iste=ist
      iene=ien
      if(ends) then
      iste=ist-1
      iene=ien+1
      endif
      do i=iste,iene
      iact(i)=nl
      enddo
      else
      write(6,32) n
32    format(/2x,'Combined strands have ',i4,' levels ...'/)
      do k=1,ns
      write(6,30) k,nu(k),dir(idr(k)),(na(i,k),i=1,nux)
      enddo
      ist=1
      ien=nux
      iste=ist
      iene=ien
      if(ends) then
      iste=ist-1
      iene=ien+1
      endif
      do i=iste,iene
      do j=1,ns
      if(li(i,j).gt.0) then
      iact(i)=j
      goto 100
      endif
      enddo
        do j=1,ns
        if(li(i,j).eq.-1) then
        iact(i)=j
        goto 100
        endif
        enddo
      write(6,*) '  ---- Level with no bases not allowed ----'
      stop
100   enddo
      endif
      icyc=0
      call analy
      if(ends) call setend
      if(test) call gradt
c----------------------------------------------------------save params
      do i=ist,ien
      uho(i,1,nl)=ux(i)
      uho(i,2,nl)=uy(i)
      uho(i,3,nl)=uz(i)
      hho(i,1,nl)=ox(i)
      hho(i,2,nl)=oy(i)
      hho(i,3,nl)=oz(i)
      enddo
      enddo
200   call params
      call locpar
      call outaxe
      if(dna.ne.' '.or.pdb.ne.' ') call plate(ns)
      if(grv) call groove
      call grgrap
c---------------------------------------------------------axe file o/p
      if(axout.ne.' ') then
      kfi=index(axout,' ')-1
      open(unit=1,file=axout(:kfi)//'.axe',status='new')
      nss=max(ns,nst)
      if(old) then
      write(1,50) nss,nu(1)*idr(1),nu(2)*idr(2),nu(3)*idr(3),
     1 nu(4)*idr(4)
50    format(5i4,6x,'#')
      else
      write(1,51) nss,nu(1)*idr(1),nu(2)*idr(2),nu(3)*idr(3),
     1 nu(4)*idr(4),0.,0.
51    format(5i4,2f8.3)
      endif
      nbto=0
      do k=1,nss
      ist=ng(k)
      ien=nr(k)
      if(ends) then
      ist=ng(k)-1
      ien=nr(k)+1
      endif
      nbto=nbto+ien-ist+1
      do i=ist,ien
      if(old) then
      write(1,52) hel(i,1,k),hel(i,2,k),hold(i,1,k),hel(i,4,k),
     1 hel(i,5,k),hold(i,2,k)
      else
      write(1,52) hel(i,1,k),hel(i,2,k),hold(i,1,k),hel(i,4,k),
     1 hel(i,5,k),hold(i,2,k),(val(i,k,j),j=1,2)
      endif
      if((idr(k).eq.1.and.i.eq.ien).or.
     1  (idr(k).eq.-1.and.i.eq.ist)) then
      write(1,52) (tor(i,k,j),j=1,6),0.0,0.0
      else
      write(1,52) (tor(i,k,j),j=1,8)
      endif
52    format(2x,8f9.4)
      enddo
      enddo
      write(1,54) (0.0,j=1,nbto)
      write(1,54) (0.0,j=1,nbto)
54    format(2x,7f9.4)
         ist=ng(1)
         ien=nr(1)
         if(ends) then
         ist=ng(1)-1
         ien=nr(1)+1
         endif
      do i=ist+1,ien
         lkink=.true.
         if(zaxe.or.line) lkink=.false.
         if(line.and.break.eq.i) lkink=.true.
      write(1,56) vold(i,1),vold(i,2),vold(i,3),vold(i,4),lkink
56    format(2x,4f9.4,l2)
      enddo
      close(unit=1)
      endif
      end
