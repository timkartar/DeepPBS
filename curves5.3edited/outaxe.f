      subroutine outaxe
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,flag,dtf,skip
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,na*1,mnam,munit,strand(4)*3,sym*1,
     1 base(9)*1,sugt(10)*8,puck*8,clin*56
      real*8 ktl,kpr,nx,ny,nz,hela(6)
      integer*4 bcod(0:n3,4),dcod(0:n3,4),tcod(0:n3,4),nbase(9),
     1 dimer(4,4),trimer(4,4,4),break,spline,jtran(7)
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
      common/vec/ux(n3),uy(n3),uz(n3),hx(n3),hy(n3),hz(n3),sx(n3),
     1 sy(n3),sz(n3),qr(n3,3,4),qp(n3,3,4),bx(n3,4),by(n3,4),
     1 bz(n3,4),ox(0:n3),oy(0:n3),oz(0:n3),n,nl,ns
      data sugt/'C3''-endo','C4''-exo ','O1''-endo','C1''-exo ',
     1'C2''-endo','C3''-exo ','C4''-endo','O1''-exo ',
     1'C1''-endo','C2''-exo '/
      data jtran/13,9,10,7,8,11,12/
      data strand/'1st','2nd','3rd','4th'/,
     1    base/'G','A','C','T','I','U','P','Y','R'/,
     1   nbase/ 1 , 2 , 3 , 4 , 1 , 4 , 4 , 2 , 3 /
     1 ((dimer(i,j),j=1,4),i=1,4)   /  1,  2,  5,  6,
     1                                 3,  4, -6,  7,
     1                                 8,  9, -1, -3,
     1                                -9, 10, -2, -4/,
     1 ((trimer(i,1,j),j=1,4),i=1,4)/  1,  2,  5,  6,
     1                                 3,  4,  7,  8,
     1                                 9, 10, 13, 14,
     1                                11, 12, 15, 16/,
     1 ((trimer(i,2,j),j=1,4),i=1,4)/ 17, 18, 21, 22,
     1                                19, 20, 23, 24,
     1                                25, 26, 29, 30,
     1                                27, 28, 31, 32/,
     1 ((trimer(i,3,j),j=1,4),i=1,4)/-13,-15, -5, -7,
     1                               -14,-16, -6, -8,
     1                                -9,-11, -1, -3,
     1                               -10,-12, -2, -4/,
     1 ((trimer(i,4,j),j=1,4),i=1,4)/-29,-31,-21,-23,
     1                               -30,-32,-22,-24,
     1                               -25,-27,-17,-19,
     1                               -26,-28,-18,-20/
c--------------------------------------setup optional data file output
      dtf=.false.
      if(daf.ne.' ') then
      dtf=.true.
      kfi=index(daf,' ')-1
      open(unit=2,file=daf(:kfi)//'.dat',status='new')
      write(2,8) daf,file,max0(ns,nst),(nu(j),j=1,4),break,comb,zaxe,
     1 line,ends,dinu
8     format(2x,a8,1x,a20,6i4,5l2/)
      endif
c------------------------------------------setup base & junction codes
      do k=1,max0(ns,nst)
      do i=ng(k),nr(k)
      bcod(i,k)=0
      do l=1,9
      if(na(i,k).eq.base(l)) bcod(i,k)=nbase(l)
      enddo
      enddo
      tcod(1,k)=0
      tcod(nr(k),k)=0
      do i=ng(k)+1,nr(k)
      ib=bcod(i-1,k)
      ii=bcod(i,k)
      if(ib.gt.0.and.ii.gt.0) then
      dcod(i,k)=dimer(ib,ii)
      if(idr(k).lt.0) dcod(i,k)=dimer(ii,ib)
         else
         dcod(i,k)=0
         endif
      if(i.lt.nr(k)) then
      ia=bcod(i+1,k)
      if(ib.gt.0.and.ii.gt.0.and.ia.gt.0) then
      tcod(i,k)=trimer(ib,ii,ia)
      if(idr(k).lt.0) tcod(i,k)=trimer(ia,ii,ib)
         else
         tcod(i,k)=0
         endif
      endif
      enddo
      if(ends) then
      bcod(ng(k)-1,k)=0
      bcod(nr(k)+1,k)=0
      tcod(ng(k)-1,k)=0
      tcod(nr(k)+1,k)=0
      endif
      enddo
c-------------------------------------------------------helical params
      call title('B','Global Base-Axis Parameters',27)
      if(dtf) write(2,*) ' Section B:'
      do k=1,max0(ns,nst)
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
      write(6,5) strand(k)
5     format(/4x,a3,' strand',6x,
     1 'Xdisp',4x,'Ydisp',3x,'Inclin',5x,'Tip',4x,'Bc',2x,'Tc',
     1 /20x,'(dx) ',4x,'(dy) ',3x,'(eta) ',3x,'(theta)'/)
      do i=iste,iene
      if(li(i,k).ge.-1) then
      iu=ni(i,k)
      ia=ncen(iu)
      sym=' '
      if(li(i,k).eq.-1) sym='*'
      write(6,10) i,sym,munit(ia),nunit(ia),hel(i,1,k),hel(i,2,k),
     1 hel(i,4,k),hel(i,5,k),bcod(i,k),tcod(i,k)
10    format(2x,i3,')',a1,a4,i3,2x,4f9.2,1x,i3,1x,i3)
      if(dtf) write(2,90) k,i,hel(i,1,k),hel(i,2,k),hel(i,4,k),
     1 hel(i,5,k),bcod(i,k),tcod(i,k)
90    format(2i3,4f10.3,2i4)
      else
      write(6,10) i,' ','----'
      endif
      enddo
      if(dtf) write(2,*) ' '
      enddo
c----------------------------------------------------multi-strand data
      if(comb) then
      call title('C','Global Base pair-Axis Parameters',32)
      if(dtf) write(2,*) ' Section C:'
      do k=2,ns
      nav=0
      write(6,60) k
60    format(/4x,'Strand 1 with strand ',i1,' ...')
      write(6,6)
6     format(/4x,'Duplex',10x,
     1 'Xdisp',4x,'Ydisp',3x,'Inclin',5x,'Tip',4x,'Bc',2x,'Tc',
     1 /20x,'(dx) ',4x,'(dy) ',3x,'(eta) ',3x,'(theta)'/)
      do j=1,6
      hela(j)=0.0
      enddo
      do i=iste,iene
      if(li(i,1).ge.-1.and.li(i,k).ge.-1) then
      nav=nav+1
      ia=ncen(ni(i,1))
      ic=ncen(ni(i,k))
      sym=' '
      if(li(i,1).eq.-1.or.li(i,k).eq.-1) sym='*'
      xdi=(hel(i,1,1)+hel(i,1,k))/2.
      ydi=(hel(i,2,1)-hel(i,2,k))/2.
      cln=aver(hel(i,4,1),hel(i,4,k))
      tip=aver(hel(i,5,1),-hel(i,5,k))
      if(abs(tip).gt.180.0) tip=tip-sign(360.d0,tip)
      if(abs(cln).gt.180.0) cln=cln-sign(360.d0,cln)
      if(idr(1).lt.idr(k)) then
      ydi=-ydi
      tip=-tip
      endif
      hela(1)=hela(1)+xdi
      hela(2)=hela(2)+ydi
      hela(4)=hela(4)+cln
      hela(5)=hela(5)+tip
      write(6,12) i,sym,munit(ia)(:1),nunit(ia),munit(ic)(:1),
     1 nunit(ic),xdi,ydi,cln,tip,bcod(i,1),tcod(i,1)
12    format(2x,i3,')',a1,a1,i3,'-',a1,i3,4f9.2,1x,i3,1x,i3)
      if(dtf) write(2,90) 0,i,xdi,ydi,cln,tip,bcod(i,1),tcod(i,1)
      else
      write(6,12) i,' ','-'
      endif
      enddo
      write(6,112) hela(1)/nav,hela(2)/nav,hela(4)/nav,hela(5)/nav
112   format(/4x,'Average: ',3x,6f9.2)
      if(dtf) write(2,*) ' '
      enddo
c---------------------------------------------------------------------
      call title('D','Global Base-Base Parameters',27)
      if(dtf) write(2,*) ' Section D:'
      do k=2,ns
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
      nav=0
      write(6,60) k
      write(6,7)
7     format(/4x,'Duplex',10x,'Shear',4x,'Stretch',2x,
     1 'Stagger',2x,'Buckle',3x,'Propel',2x,'Opening',2x,'Bc',2x,'Tc',
     1 /20x,'(Sx) ',4x,' (Sy)  ',2x,' (Sz)  ',2x,'(kappa)',
     1 2x,'(omega)',1x,'(sigma)'/)
      stg=0.
      opn=0.
      do j=1,6
      hela(j)=0.0
      enddo
      skip=.false.
      do i=iste,iene
      if(li(i,1).ge.-1.and.li(i,k).ge.-1) then
      nav=nav+1
      ia=ncen(ni(i,1))
      ic=ncen(ni(i,k))
      sym=' '
      if(li(i,1).eq.-1.or.li(i,k).eq.-1) sym='*'
      str=hel(i,2,1)+hel(i,2,k)
      pro=diff(hel(i,5,1),-hel(i,5,k))
      if(abs(pro).gt.180.0) pro=pro-sign(360.d0,pro)
      stg=stg+hel(i,3,1)-hel(i,3,k)
      opn=opn+hel(i,6,1)-hel(i,6,k)
         if(i.gt.iste) then
         if(li(i-1,1).lt.-1.or.li(i-1,k).lt.-1) then
         stg=0.
         opn=0.
         endif
         endif
      if(idr(1).ge.idr(k)) then
      shr=hel(i,1,1)-hel(i,1,k)
      buc=diff(hel(i,4,1),hel(i,4,k))
      else
      shr=hel(i,1,k)-hel(i,1,1)
      buc=diff(hel(i,4,k),hel(i,4,1))
      endif
      hela(1)=hela(1)+shr
      hela(2)=hela(2)+str
      hela(3)=hela(3)+stg
      hela(4)=hela(4)+buc
      hela(5)=hela(5)+pro
      hela(6)=hela(6)+opn
      write(6,14) i,sym,munit(ia)(:1),nunit(ia),munit(ic)(:1),
     1 nunit(ic),shr,str,stg,buc,pro,opn,bcod(i,1),tcod(i,1)
14    format(2x,i3,')',a1,a1,i3,'-',a1,i3,6f9.2,1x,i3,1x,i3)
      if(dtf) write(2,92) 0,i,shr,str,stg,buc,pro,opn,bcod(i,1),
     1 tcod(i,1)
92    format(2i3,6f10.3,2i4)
         else
         write(6,14) i,' ','-'
         endif
      enddo
      write(6,112) (hela(j)/nav,j=1,6)
      if(dtf) write(2,*) ' '
      enddo
      endif
c-------------------------------------------global junction parameters
      call title('E','Global Inter-Base Parameters',28)
      if(dtf) write(2,*) ' Section E:'
      do k=1,max0(ns,nst)
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
      write(6,11) strand(k)
11    format(/4x,a3,' strand',6x,'Shift',4x,'Slide',5x,'Rise',5x,
     1 'Tilt',5x,'Roll',4x,'Twist',3x,'Dc',
     1/20x,'(Dx) ',4x,'(Dy) ',5x,'(Dz)',5x,'(tau)',4x,
     1 '(rho)',2x,'(Omega)'/)
      is=k
      if(comb) is=1
      do i=iste+1,iene
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      ia=ncen(ni(i,k))
      ib=ncen(ni(i-1,k))
      sym=' '
      if(li(i-1,k).eq.-1.or.li(i,k).eq.-1) sym='*'
      shif=hel(i,1,k)+vkin(i,1,is)-hel(i-1,1,k)
      slid=hel(i,2,k)+vkin(i,2,is)*idr(k)-hel(i-1,2,k)
      rise=hel(i,3,k)
      tilt=hel(i,4,k)+vkin(i,3,is)-hel(i-1,4,k)
      roll=hel(i,5,k)+vkin(i,4,is)*idr(k)-hel(i-1,5,k)
      twis=hel(i,6,k)
      if(abs(roll).gt.180.) roll=roll-sign(360.d0,roll)
      write(6,16) i,sym,munit(ib)(:1),nunit(ib),munit(ia)(:1),
     1 nunit(ia),shif,slid,rise,tilt,roll,twis,dcod(i,k)
16    format(2x,i3,')',a1,a1,i3,'/',a1,i3,6f9.2,1x,i3)
      if(dtf) write(2,92) k,i,shif,slid,rise,tilt,roll,twis,dcod(i,k)
      else
      write(6,16) i,' ','-'
      endif
      enddo
      if(dtf) write(2,*) ' '
      enddo
c-----------------------------------------------ds junction parameters
      if(comb) then
      call title('F','Global Inter-Base pair Parameters',33)
      if(dtf) write(2,*) ' Section F:'
      do k=2,ns
      nav=0
      write(6,60) k
      write(6,15)
15    format(/4x,'Duplex',10x,'Shift',4x,'Slide',5x,'Rise',5x,
     1 'Tilt',5x,'Roll',4x,'Twist',3x,'Dc',
     1 /20x,'(Dx) ',4x,'(Dy) ',5x,'(Dz)',5x,'(tau)',
     1 4x,'(rho)',2x,'(Omega)'/)
      do j=1,6
      hela(j)=0.
      enddo
      do i=iste+1,iene
      if(li(i-1,1).ge.-1.and.li(i,1).ge.-1.and.
     1   li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      nav=nav+1
      ia=ncen(ni(i,1))
      ib=ncen(ni(i-1,1))
      sym=' '
      if(li(i,1).eq.-1.or.li(i-1,1).eq.-1.or.
     1   li(i,k).eq.-1.or.li(i-1,k).eq.-1) sym='*'
      xs=(hel(i,1,1)+hel(i,1,k))/2.
      xm=(hel(i-1,1,1)+hel(i-1,1,k))/2.
      ys=(hel(i,2,1)-hel(i,2,k))/2.
      ym=(hel(i-1,2,1)-hel(i-1,2,k))/2.
      ts=aver(hel(i,4,1),hel(i,4,k))
      tm=aver(hel(i-1,4,1),hel(i-1,4,k))
      ps=aver(hel(i,5,1),-hel(i,5,k))
      pm=aver(hel(i-1,5,1),-hel(i-1,5,k))
      if(idr(1).lt.idr(k)) then
      ys=-ys
      ym=-ym
      ps=-ps
      pm=-pm
      endif
      shif=xs+vkin(i,1,1)-xm
      slid=ys+vkin(i,2,1)-ym
      rise=(hel(i,3,1)+hel(i,3,k))/2.
      tilt=ts+vkin(i,3,1)-tm
      roll=ps+vkin(i,4,1)-pm
      twis=(hel(i,6,1)+hel(i,6,k))/2.
      if(abs(roll).gt.180.) roll=roll-sign(360.d0,roll)
      hela(1)=hela(1)+shif
      hela(2)=hela(2)+slid
      hela(3)=hela(3)+rise
      hela(4)=hela(4)+tilt
      hela(5)=hela(5)+roll
      hela(6)=hela(6)+twis
      write(6,16) i,sym,munit(ib)(:1),nunit(ib),munit(ia)(:1),
     1 nunit(ia),shif,slid,rise,tilt,roll,twis,dcod(i,1)
      if(dtf) write(2,92) 0,i,shif,slid,rise,tilt,roll,twis,dcod(i,1)
      else
      write(6,16) i,' ','-'
      endif
      enddo
      write(6,112) (hela(j)/nav,j=1,6)
      if(dtf) write(2,*) ' '
      enddo
      endif
c--------------------------------------------local junction parameters
      call title('G','Local Inter-Base Parameters',27)
      if(dtf) write(2,*) ' Section G:'
      do k=1,max0(ns,nst)
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
      write(6,11) strand(k)
      do i=iste+1,iene
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      ia=ncen(ni(i,k))
      ib=ncen(ni(i-1,k))
      write(6,16) i,' ',munit(ib)(:1),nunit(ib),munit(ia)(:1),
     1 nunit(ia),(pal(i,j,k),j=1,6),dcod(i,k)
      if(dtf) write(2,92) k,i,(pal(i,j,k),j=1,6),dcod(i,k)
      else
      write(6,16) i,' ','-'
      endif
      enddo
      if(dtf) write(2,*) ' '
      enddo
c-----------------------------------------------ds junction parameters
      if(comb) then
      call title('H','Local Inter-Base pair Parameters',32)
      if(dtf) write(2,*) ' Section H:'
      do k=2,ns
      nav=0
      write(6,60) k
      write(6,15)
      do j=1,6
      hela(j)=0.
      enddo
      do i=iste+1,iene
      if(li(i-1,1).ge.-1.and.li(i,1).ge.-1.and.
     1   li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      nav=nav+1
      ia=ncen(ni(i,1))
      ib=ncen(ni(i-1,1))
      hela(1)=hela(1)+pab(i,1,k)
      hela(2)=hela(2)+pab(i,2,k)
      hela(3)=hela(3)+pab(i,3,k)
      hela(4)=hela(4)+pab(i,4,k)
      hela(5)=hela(5)+pab(i,5,k)
      hela(6)=hela(6)+pab(i,6,k)
      write(6,16) i,' ',munit(ib)(:1),nunit(ib),munit(ia)(:1),
     1 nunit(ia),(pab(i,j,k),j=1,6),dcod(i,1)
      if(dtf) write(2,92) 0,i,(pab(i,j,k),j=1,6),dcod(i,1)
      else
      write(6,16) i,' ','-'
      endif
      enddo
      write(6,112) (hela(j)/nav,j=1,6)
      if(dtf) write(2,*) ' '
      enddo
      endif
c------------------------------------------------------axis parameters
      call title('I','Global Axis Curvature',21)
      if(dtf) write(2,*) ' Section I:'
      do k=1,nst
      if(.not.comb) then
      iste=ng(k)
      iene=nr(k)
      if(ends) then
      iste=iste-1
      iene=iene+1
      endif
      endif
      if(comb) then
      write(6,9)
9     format(/4x,'Duplex',10x,'Ax',6x,'Ay',5x,'Ainc',4x,'Atip',
     1       4x,'Adis',4x,'Angle',3x,'Path',3x,'Dc'/)
      else
      write(6,91) strand(k)
91    format(/4x,a3,' strand',6x,'Ax',6x,'Ay',5x,'Ainc',4x,'Atip',
     1       4x,'Adis',4x,'Angle',3x,'Path',3x,'Dc'/)
      endif
      do i=iste+1,iene
      dsx=vkin(i,1,k)
      dsy=vkin(i,2,k)
      ktl=vkin(i,3,k)
      kpr=vkin(i,4,k)
      dis=vkin(i,5,k)
      the=vkin(i,6,k)
      pat=vkin(i,7,k)
      if(li(i-1,k).ge.-1.and.li(i,k).ge.-1) then
      ia=ncen(ni(i,k))
      ib=ncen(ni(i-1,k))
      write(6,18) i,munit(ib)(:1),nunit(ib),munit(ia)(:1),nunit(ia),
     1 dsx,dsy,ktl,kpr,dis,the,pat,dcod(i,k)
18    format(2x,i2,') ',a1,i3,'/',a1,i3,7f8.2,1x,i3)
      else
      write(6,19) i,'-','-',dsx,dsy,ktl,kpr,dis,the,pat,dcod(i,k)
19    format(2x,i2,') ',a1,3x,'/',a1,3x,7f8.2,1x,i3)
      endif
      if(dtf) write(2,92) k,i,dsx,dsy,ktl,kpr,the,pat,dcod(i,k)
      enddo
      if(dtf) write(2,*) ' '
      enddo
      if(dtf) close(unit=2)
c----------------------------------------------------------helical fit
      do k=1,nst
         if(.not.comb) then
         iste=ng(k)
         iene=nr(k)
         if(ends) then
         iste=iste-1
         iene=iene+1
         endif
         endif
      buu=bend(iste,k)
      bpp=bend(iene,k)
      bend(iste,k)=0.
      bend(iene,k)=0.
      write(6,101) buu,bpp
 101  format(/2x,'Overall axis bend ... UU= ',f7.2,' PP= ',f7.2/)
c---------------------------------------------------bending parameters
      if(comb) then
      write(6,17)
17    format(/4x,'Duplex',7x,'Offset',3x,
     & 'L.Dir  ... wrt end-to-end vector'/)
      else
      write(6,20) strand(k)
20    format(/4x,a3,' strand',3x,'Offset',3x,
     & 'L.Dir  ... wrt end-to-end vector'/)
      endif
      x0=hho(iste,1,k)
      y0=hho(iste,2,k)
      z0=hho(iste,3,k)
      nx=hho(iene,1,k)-x0
      ny=hho(iene,2,k)-y0
      nz=hho(iene,3,k)-z0
      rn=sqrt(nx*nx+ny*ny+nz*nz)
      nx=nx/rn
      ny=ny/rn
      nz=nz/rn
      pl=0.
      do i=iste,iene
      if(i.gt.iste) pl=pl+vkin(i,7,k)
      x=hho(i,1,k)-x0
      y=hho(i,2,k)-y0
      z=hho(i,3,k)-z0
      dot=x*nx+y*ny+z*nz
      dx=x-nx*dot
      dy=y-ny*dot
      dz=z-nz*dot
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      if(li(i,k).ge.-1) then
      ia=ncen(ni(i,k))
      write(6,22) i,munit(ia)(:1),nunit(ia),rd,bend(i,k)
22    format(2x,i2,') ',a1,i3,5x,2f8.2)
      else
      write(6,23) i,'-',rd
23    format(2x,i2,') ',a1,8x,f8.2)
      endif
      enddo
      write(6,24) pl,rn,100.*(1.-rn/pl)
24    format(/4x,'Path length= ',f8.2,'  End-to-end= ',f8.2,
     1 '  Shortening= ',f8.2,' %')
      enddo
c------------------------------------------------------backbone params
      call title('J','Backbone Parameters',19)
      do k=1,max0(ns,nst)
      if(.not.comb) then
      ist=ng(k)
      ien=nr(k)
      endif
c-----------------------------------------------------sugar parameters
      write(6,35) strand(k)
35    format(/4x,a3,' strand',2x,' C1''-C2''',' C2''-C3''','  Phase ',
     1 '  Ampli ',2x,'Pucker  ','  C1'' ','  C2'' ','  C3'' '/)
      do i=ist,ien
      if(li(i,k).ge.-2) then
      iu=ni(i,k)
      ia=ncen(iu)
      if(flag(i,k)) then
      write(6,45) i,munit(ia),nunit(ia)
      else
      am=suga(i,k,1)
      ph=suga(i,k,2)
      puck=sugt(1+int(ph/36.))
      write(6,45) i,munit(ia),nunit(ia),tor(i,k,5),tor(i,k,6),ph,am,
     1 puck,tor(i,k,1),tor(i,k,2),tor(i,k,3)
45    format(2x,i3,')',a4,i3,2x,4f8.2,2x,a8,1x,3f6.1)
      endif
      else
      write(6,45) i,'----'
      endif
      enddo
c--------------------------------------------------backbone parameters
      write(6,47)
47    format( /4x,'Torsions',4x,'   Chi  ','  Gamma ','  Delta ',
     1              '  Epsil ','  Zeta  ','  Alpha ','  Beta  ',
     1       /16x,'  C1''-N ',' C5''-C4''',' C4''-C3''',' C3''-O3''',
     1            '  O3''-P ','  P-O5'' ',' O5''-C5'''/)
      do i=ist,ien
      if(li(i,k).ge.-2) then
      iu=ni(i,k)
      ia=ncen(iu)
         do jj=1,7
         j=jtran(jj)
         js=(jj-1)*8+1
         if(tor(i,k,j).lt.990) then
         write(clin(js:),'(f8.2)') tor(i,k,j)
         else
         write(clin(js:),'(''  ......'')')
         endif
         enddo
      write(6,31) i,munit(ia),nunit(ia),clin
31    format(2x,i3,')',a4,i3,2x,a56)
      else
      write(6,31) i,'----'
      endif
      enddo
      write(6,*) ' '
      enddo
      return
      end
