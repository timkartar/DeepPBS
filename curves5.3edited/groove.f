      subroutine groove
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly,inters,sgn1,sgn2,st1,st2
      character*4 mnam,munit,na*1,atoms(9),bato,c1*1,c5(2)*6,
     1 c6(2)*2,c7(2)*6,c8(2)*4,clear*21,ln(500)*72,c9*6
      integer*4 spline,break
      dimension hpm(2),gqm(2),ohm(2),ogm(2),hpn(2),gqn(2),ohn(2),ogn(2),
     1 pqn(2),hpi(2),gqi(2),ohi(2),ogi(2),pqi(2),dm(2),dep(2),inm(0:2),
     1 jnm(0:2),inn(2),jnn(2),ini(2),jni(2),ian(2),radius(9)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/gro/uxb(5000,0:n3),uyb(5000,0:n3),uzb(5000,0:n3),
     1 cor(n1,3),dya(n1,3),box(0:n1,4),boy(0:n1,4),
     1 boz(0:n1,4),ind(n3,0:n3),nma(n3),num,numa,nsu(4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      data clear/'   --       --     --'/
      data atoms/'C1''','C2''','C3''','C4''','O1''','O3''','P',
     1 'O5''','C5'''/,radius/1.6,1.6,1.6,1.6,1.4,1.4,2.9,1.4,1.6/
      ox(in,l)=box(in,l)-corx
      oy(in,l)=boy(in,l)-cory
      oz(in,l)=boz(in,l)-corz
      om(in,l)=sqrt(ox(in,l)**2+oy(in,l)**2+oz(in,l)**2)
      scal(in,l)=dyax*ox(in,l)+dyay*oy(in,l)+dyaz*oz(in,l)
      scau(in,l)=uux *ox(in,l)+uuy *oy(in,l)+uuz *oz(in,l)
      scat(in,l)=tex *ox(in,l)+tey *oy(in,l)+tez *oz(in,l)
      if(nbac.gt.9) stop '  ---- nbac too big, cannot exceed 9 ----'
      bato=atoms(nbac)
      vdw=radius(nbac)
      nbac=nbac+2
      if(nbac.eq.11) nbac=13
      call title('K','Groove parameters',17)
c--------------------------------------------------------find back-atoms
      write(6,11) bato,abs(nu(1)),nlevel
11    format(/2x,'Atom defining backbone: ',a4,2x,i3,' levels, ',
     1 i2,' sub-levels')
      le=0
      numa=nux
      do i=1,nux
      l=ni(i,1)*ni(i,2)
      if(le.eq.0.and.l.ne.0) num=i
      if(le.ne.0.and.l.eq.0) numa=i-1
      le=l
      enddo
      call axeint
      call bacint(bato)
c---------------------------------------------------first data for depth
      write(6,14)
14    format(/3x,'Levels',10x,'Minor groove',16x,'Major groove',
     1 /5x,'i',2x,'n',6x,'Width',4x,'Depth',2x,'Angle',7x,'Width',
     1  4x,'Depth',2x,'Angle',4x,'Diam'/)
      ea=-3.
      eb= 4.
      xdi=(hel(num,1,1)+hel(num,1,2))/2.
      tip=aver(hel(num,5,1),-hel(num,5,2))
      if(abs(tip).gt.180.) tip=tip-sign(360.d0,tip)
      if(cos(tip*cdr).lt.0.) xdi=-xdi
      oa1=xdi+ea
      ob1=xdi+eb
      kli=0
      do i=num,numa
      oa0=oa1
      ob0=ob1
      xdi=(hel(i+1,1,1)+hel(i+1,1,2))/2.
      tip=aver(hel(i+1,5,1),-hel(i+1,5,2))
      if(abs(tip).gt.180.) tip=tip-sign(360.d0,tip)
      if(cos(tip*cdr).lt.0.) xdi=-xdi
      oa1=xdi+ea
      ob1=xdi+eb
c-----------------------------------------intersections with the strands
      insu=nsu(1)
      jnsu=nsu(2)
      nmi=nma(i)
      do n=0,nmi-1
      k=ind(i,n)
      dyax=dya(k,1)
      dyay=dya(k,2)
      dyaz=dya(k,3)
      corx=cor(k,1)
      cory=cor(k,2)
      corz=cor(k,3)
      uux=uxb(i,n)
      uuy=uyb(i,n)
      uuz=uzb(i,n)
      tex=uuy*dyaz-uuz*dyay
      tey=uuz*dyax-uux*dyaz
      tez=uux*dyay-uuy*dyax
      inin=20*(i-5)-30
      inma=20*(i+5)-30
      inma=min0(insu,inma)
      inin=max0(0,inin)
      jnin=20*(i-5)-30
      jnma=20*(i+5)-30
      jnin=max0(0,jnin)
      jnma=min0(jnsu,jnma)
      hpm(1)=99.
      hpm(2)=99.
      gqm(1)=99.
      gqm(2)=99.
      pqn(1)=99.
      pqn(2)=99.
      dm(1)=99.
      dm(2)=99.
      ca=99.
      bd=99.
      dep(1)=99.
      dep(2)=99.
c--------------------------------------------------reference on strand 1
      pqee=1.
      pqe=2.
      dee=1.
      de=2.
      hpee=1.
      hpe=2.
      lee=0
      le=0
      jee=0
      je=0
      pqex=99.
      inex=-1
      scaue=scau(inin,1)
      rad1=0.
      x1=0.
      do in=inin,inma,2
      scu=scau(in,1)
      if(scu*scaue.lt.0.)x1=(om(in-2,1)*scu-om(in,1)*scaue)/(scu-scaue)
      if(x1.gt.rad1) rad1=x1
      inters=.false.
      sgn1=.false.
      sgn2=.false.
      j=0
      oh=scal(in,1)
      if(oh.lt.0.)l=1
      if(oh.ge.0.)l=2
      if(l.eq.le.and.le.eq.lee)sgn1=.true.
      hp=sqrt(om(in,1)**2-oh**2)
      if(sgn1.and.hp.gt.hpe.and.hpe.lt.hpee.and.hpe.lt.hpm(l))then
      hpm(l)=hpe
      ohm(l)=ohe
      inm(l)=in-2
      endif
c---------------------------------------------intersection with strand 2
      vx=dyay*oz(in,1)-dyaz*oy(in,1)
      vy=dyaz*ox(in,1)-dyax*oz(in,1)
      vz=dyax*oy(in,1)-dyay*ox(in,1)
      s=vx*ox(jnin,2)+vy*oy(jnin,2)+vz*oz(jnin,2)
      do jm=jnin,jnma-10,10
      r=s
      s=vx*ox(jm+10,2)+vy*oy(jm+10,2)+vz*oz(jm+10,2)
      if(r*s.lt.0.)then
      s=r
      jn=jm
      do while(r*s.gt.0.0.and.jn.lt.jm+10)
      jn=jn+1
      r=s
      s=vx*ox(jn,2)+vy*oy(jn,2)+vz*oz(jn,2)
      enddo
      a=scau(in,1)*scau(jn,2)
      b=scat(in,1)*scat(jn,2)
      if(a.le.0.0.and.b.le.0.)then
      inters=.true.
      jnt=jn
      og=scal(jn,2)
      if(og.lt.0.)j=1
      if(og.ge.0.)j=2
      if(j.eq.je.and.je.eq.jee)sgn2=.true.
      gq=sqrt(om(jn,2)**2-og**2)
      pq=sqrt((box(in,1)-box(jn,2))**2+(boy(in,1)-boy(jn,2))**2
     1         +(boz(in,1)-boz(jn,2))**2)
      if(in.eq.inex.and.pq.gt.pqex)then
      jnt=jnex
      og=ogex
      j=jex
      gq=gqex
      pq=pqex
      endif
      inex=in
      jnex=jn
      ogex=og
      jex=j
      gqex=gq
      pqex=pq
      endif
      endif
      enddo
c--------------------------------search for minimum I1I2 (or inflection)
      if(inters)then
      d=pq-pqee
      if(l.eq.j.and.sgn1.and.sgn2) then
      if(pq.gt.pqe.and.pqe.lt.pqee.and.pqe.lt.pqn(j)) then
      pqn(j)=pqe
      inn(j)=in-2
      jnn(j)=jne
      ohn(j)=ohe
      ogn(j)=oge
      hpn(j)=hpe
      gqn(j)=gqe
      endif
      if(d*de.gt.0.0.and.d*dee.gt.0.0.and.abs(d).gt.abs(de).and.abs(de)
     1 .lt.abs(dee).and.abs(de).lt.dm(j)) then
c     1 .lt.abs(dee).and.abs(de).lt.dm(j).and.abs(de+dee).lt.0.9) then
      pqi(j)=pqee
      ini(j)=in-4
      jni(j)=jnee
      ohi(j)=ohee
      ogi(j)=ogee
      hpi(j)=hpee
      gqi(j)=gqee
      dm(j)=abs(de)
      endif
      endif
      pqee=pqe
      pqe=pq
      jnee=jne
      jne=jnt
      ogee=oge
      oge=og
      gqee=gqe
      gqe=gq
      dee=de
      de=d
      endif
      hpee=hpe
      hpe=hp
      ohee=ohe
      ohe=oh
      lee=le
      le=l
      jee=je
      je=j
      scaue=scu
      enddo
c--------------------------------------------------reference on strand 2
      pqee=1.
      pqe=2.
      dee=1.
      de=2.
      gqee=1.
      gqe=2.
      lee=0
      le=0
      jee=0
      je=0
      pqex=99.
      jnex=-1
      scaue=scau(jnin,2)
      rad2=0.
      x2=0.
      do jn=jnin,jnma,2
      scu=scau(jn,2)
      if(scu*scaue.lt.0)x2=(om(jn-2,2)*scu-om(jn,2)*scaue)/(scu-scaue)
      if(x2.gt.rad2) rad2=x2
      inters=.false.
      sgn1=.false.
      sgn2=.false.
      l=0
      og=scal(jn,2)
      if(og.lt.0.)j=1
      if(og.ge.0.)j=2
      if(j.eq.je.and.je.eq.jee)sgn2=.true.
      gq=sqrt(om(jn,2)**2-og**2)
      if(sgn2.and.gq.gt.gqe.and.gqe.lt.gqee.and.gqe.lt.gqm(j))then
      gqm(j)=gqe
      ogm(j)=oge
      jnm(j)=jn-2
      endif
c---------------------------------------------intersection with strand 1
      vx=dyay*oz(jn,2)-dyaz*oy(jn,2)
      vy=dyaz*ox(jn,2)-dyax*oz(jn,2)
      vz=dyax*oy(jn,2)-dyay*ox(jn,2)
      s=vx*ox(inin,1)+vy*oy(inin,1)+vz*oz(inin,1)
      do im=inin,inma-10,10
      r=s
      s=vx*ox(im+10,1)+vy*oy(im+10,1)+vz*oz(im+10,1)
      if(r*s.lt.0.) then
      s=r
      in=im
      do while(r*s.gt.0.0.and.in.lt.im+10)
      in=in+1
      r=s
      s=vx*ox(in,1)+vy*oy(in,1)+vz*oz(in,1)
      enddo
      a=scau(in,1)*scau(jn,2)
      b=scat(in,1)*scat(jn,2)
      if(a.le.0.0.and.b.le.0.)then
      inters=.true.
      intr=in
      oh=scal(in,1)
      if(oh.lt.0.)l=1
      if(oh.ge.0.)l=2
      if(l.eq.le.and.le.eq.lee)sgn1=.true.
      hp=sqrt(om(in,1)**2-oh**2)
      pq=sqrt((box(in,1)-box(jn,2))**2+(boy(in,1)-boy(jn,2))**2
     1         +(boz(in,1)-boz(jn,2))**2)
      if(jn.eq.jnex.and.pq.gt.pqex)then
      intr=inex
      oh=ohex
      l=lex
      hp=hpex
      pq=pqex
      endif
      inex=in
      jnex=jn
      ohex=oh
      lex=l
      hpex=hp
      pqex=pq
      endif
      endif
      enddo
c--------------------------------search for minimum I1I2 (or inflection)
      if(inters)then
      d=pq-pqee
      if(l.eq.j.and.sgn1.and.sgn2)then
      if(pq.gt.pqe.and.pqe.lt.pqee.and.pqe.lt.pqn(j))then
      pqn(j)=pqe
      inn(j)=ine
      jnn(j)=jn-2
      ohn(j)=ohe
      ogn(j)=oge
      hpn(j)=hpe
      gqn(j)=gqe
      endif
      if(d*de.gt.0.0.and.d*dee.gt.0.0.and.abs(d).gt.abs(de).and.abs(de)
     1 .lt.abs(dee).and.abs(de).lt.dm(j)) then
c     1 .lt.abs(dee).and.abs(de).lt.dm(j).and.abs(de+dee).lt.0.9) then
      pqi(j)=pqee
      ini(j)=inee
      jni(j)=jn-4
      ohi(j)=ohee
      ogi(j)=ogee
      hpi(j)=hpee
      gqi(j)=gqee
      dm(j)=abs(de)
      endif
      endif
      pqee=pqe
      pqe=pq
      inee=ine
      ine=intr
      ohee=ohe
      ohe=oh
      hpee=hpe
      hpe=hp
      dee=de
      de=d
      endif
      ogee=oge
      oge=og
      gqee=gqe
      gqe=gq
      lee=le
      le=l
      jee=je
      je=j
      scaue=scu
      enddo
c---------------------------------------------------selections for width
      do l=1,2
      pq=pqn(l)
      if(pq.lt.99.)then
      im=inn(l)
      jm=jnn(l)
      if(im.eq.0.or.im.eq.insu)then
      pq=99.
      pqn(l)=99.
      endif
      if(jm.eq.0.or.jm.eq.jnsu)then
      pq=99.
      pqn(l)=99.
      endif
      endif
      if(pq.eq.99.)then
      in=inm(l)
      jn=jnm(l)
      hp=hpm(l)
      gq=gqm(l)
      if(hp.ne.99.0.and.(in.eq.0.or.in.eq.insu))then
      hp=99.
      hpm(l)=99.
      endif
      if(gq.ne.99.0.and.(jn.eq.0.or.jn.eq.jnsu))then
      gq=99.
      gqm(l)=99.
      endif
      if(hp.ne.99.0.and.gq.ne.99.)then
      if(hp.le.gq)then
      gq=99.
      gqm(l)=99.
      elseif(hp.gt.gq)then
      hp=99.
      hpm(l)=99.
      endif
      endif
      endif
c-----------------------------------------------------------------depths
      if(l.eq.1)oa=oa0*(nmi-n)/nmi+oa1*n/nmi
      if(l.eq.2)oa=ob0*(nmi-n)/nmi+ob1*n/nmi
      if(pq+dm(l).lt.198.)then
      if(pq.lt.99.)then
      im=inn(l)
      jm=jnn(l)
      r=hpn(l)/(hpn(l)+gqn(l))
      else
      im=ini(l)
      jm=jni(l)
      r=hpi(l)/(hpi(l)+gqi(l))
      endif
      ohx=r*box(jm,2)+(1-r)*box(im,1)-corx
      ohy=r*boy(jm,2)+(1-r)*boy(im,1)-cory
      ohz=r*boz(jm,2)+(1-r)*boz(im,1)-corz
      oh=ohx*dyax+ohy*dyay+ohz*dyaz
      if(l.eq.1) ca=oa-oh
      if(l.eq.1) dep(1)=oa-oh
      if(l.eq.2) bd=oh-oa
      if(l.eq.2) dep(2)=oh-oa
c-----------------------------------------------------------------angles
      vx=oy(im,1)*dyaz-oz(im,1)*dyay
      vy=oz(im,1)*dyax-ox(im,1)*dyaz
      vz=ox(im,1)*dyay-oy(im,1)*dyax
      vm=sqrt(vx*vx+vy*vy+vz*vz)
      vx=vx/vm
      vy=vy/vm
      vz=vz/vm
      co=tex*vx+tey*vy+tez*vz
      si=uux*vx+uuy*vy+uuz*vz
      if(co.lt.0.) co=-co
      an=dacos(co)*crd
      if(co.gt.1.) an=0.
      if(co*si.gt.0.) an=-an
      ian(l)=nint(an)
      endif
c---------------------------------------------------------prepare output
      c6(l)=' '
      write(c8(l),71) ian(l)
      if(pqn(l).lt.99.) then
      write(c5(l),73) pqn(l)-2*vdw
      write(c7(l),73) dep(l)
      else if(pqn(l).eq.99..and.dm(l).lt.99.) then
      write(c5(l),73) pqi(l)-2*vdw
      c6(l)='*'
      write(c7(l),73) dep(l)
      else if(pqn(l)+dm(l).eq.198.) then
      c5(l)='   -- '
      c7(l)='   -- '
      c8(l)='  --'
      endif
      enddo
      ia=ncen(ni(i,1))
      c1=' '
      if(n.eq.0) c1=munit(ia)
      c9='   -- '
      if(rad1*rad2.ne.0.)write(c9,73) rad1+rad2
      kli=kli+1
      write(ln(kli),74) c1,i,n,c5(1),c6(1),c7(1),c8(1),c1,
     1 c5(2),c6(2),c7(2),c8(2),c9
71    format(i4)
73    format(f6.2)
74    format(2x,a1,2i3,4x,a6,a1,2x,a6,2x,a4,4x,a1,2x,a6,a1,2x,a6,2x,a4,
     1 4x,a6)
      enddo
      enddo
c--------------------------------------------------------clip and output
      st1=.true.
      st2=.true.
      do i=1,kli
      if(ln(i)(17:18).ne.'--') then
      if(ln(i)(20:20).ne.'*') st1=.false.
      if(st1.and.ln(i)(20:20).eq.'*') ln(i)(14:34)=clear
      endif
      if(ln(i)(45:46).ne.'--') then
      if(ln(i)(48:48).ne.'*') st2=.false.
      if(st2.and.ln(i)(48:48).eq.'*') ln(i)(42:62)=clear
      endif
      enddo
      st1=.true.
      st2=.true.
      do i=kli,1,-1
      if(ln(i)(17:18).ne.'--') then
      if(ln(i)(20:20).ne.'*') st1=.false.
      if(st1.and.ln(i)(20:20).eq.'*') ln(i)(14:34)=clear
      endif
      if(ln(i)(45:46).ne.'--') then
      if(ln(i)(48:48).ne.'*') st2=.false.
      if(st2.and.ln(i)(48:48).eq.'*') ln(i)(42:62)=clear
      endif
      enddo
      do i=1,kli
      write(6,'(a)') ln(i)
      enddo
      return
      end
