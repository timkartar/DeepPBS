      subroutine input
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,
     1 grv,old,axonly
      character*4 file*32,lis*32,dna*32,axin*32,axout*32,daf*32,
     1 pdb*32,mcode*32,name,mnam,munit,anam,auni*3,acha*1,
     1 lini*80,ext*3,nam(55)*2
      integer*4 break,spline
      dimension ind(55)
      common/cha/file,lis,dna,axin,axout,daf,pdb,mcode
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,
     1 old,axonly
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      data nam/'H ','HE',
     1 'LI','BE',       'B ','C ','N ','O ','F ','NE',
     1 'NA','MG',       'AL','SI','P ','S ','CL','AR',
     1 'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     1                  'GA','GE','AS','SE','BR','KR',
     1 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',
     1                  'IN','SN','SB','TE','I ','XE',
     1 'CS'/
      data ind/1,2,2,2,1,1,1,1,1,2,2,2,2,2,1,1,2,2,
     1         1,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,
     1         2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
      kfi=index(file,' ')-1
         do i=kfi,1,-1
         if(file(i:i).eq.'.') then
         ext=file(i+1:)
         goto 11
         endif
         enddo
11    open(unit=1,file=file(:kfi),status='old')
      if(ext.eq.'mac'.or.ext.eq.'MAC') then
c------------------------------------------------------read mac format
      read(1,5) mcode
      if(mcode(:1).ne.'#') then
      write(6,*) '----  Wrong .MAC format ----'
      stop
      else
      mcode=mcode(2:)
9     read(1,5) lini
      if(lini(:1).eq.'#') goto 9
      read(lini,6) kam,kcen
      endif
      read(1,7) (mnam(i),corm(i,1),corm(i,2),corm(i,3),imch(i),
     1 imty(i),dmon(i),munit(i),nunit(i),icm(i),i=1,kam)
      read(1,8) ((matd(i,j),j=1,6),i=1,kam)
      do i=1,kam
      m=0
      do j=1,6
      if(matd(i,j).ne.0) m=m+1
      enddo
      matd(i,7)=m
      enddo
5     format(a)
6     format(2i5)
7     format(a4,3f10.5,2i3,f8.4,1x,a4,2i4)
8     format(7x,6i5,8x,6i5)
c------------------------------------------------------read pdb format
      else if(ext.eq.'pdb'.or.ext.eq.'PDB'.or.
     1        ext.eq.'brk'.or.ext.eq.'BRK') then
      i=0
      read(1,105) lini
105   format(a80)
         if(lini(:4).eq.'ATOM') then
         mcode='No Title'
         goto 101
         endif
      mcode=lini(8:50)
108   if(mcode(:1).eq.' ') then
      mcode=mcode(2:)
      goto 108
      endif
100   read(1,105,end=110) lini
101   if(lini(:4).ne.'ATOM') goto 100
      read(lini,106) name,anam,auni,acha,nuni,x,y,z
106   format(a4,8x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
      i=i+1
      mnam(i)=  anam
      munit(i)= auni//acha
      nunit(i)= nuni
      corm(i,1)=x
      corm(i,2)=y
      corm(i,3)=z
      goto 100
110   kam=i
c---------------------------------------------read general data format
      else if(ext.eq.'dat'.or.ext.eq.'DAT') then
      read(1,5) mcode
      i=0
20    i=i+1
      read(1,15,end=50) mnam(i),corm(i,1),corm(i,2),corm(i,3),
     1 munit(i),nunit(i)
15    format(a4,3f10.5,1x,a4,i4)
      goto 20
50    kam=i-1
c-----------------------------------------------------------------trap
      else
      write(6,*) '  ---- Unknown geometry file type ----'
      stop
      endif
      close(unit=1)
c----------------------------------left justify atom and subunit names
      do i=1,kam
10    if(mnam(i)(:1).eq.' '.or.(ichar(mnam(i)(:1)).ge.48
     1 .and.ichar(mnam(i)(:1)).le.57)) then
      mnam(i)=mnam(i)(2:)
      goto 10
      endif
12    if(munit(i)(:1).eq.' ') then
      munit(i)=munit(i)(2:)
      goto 12
      endif
      enddo
c--------------------------------------------------------find subunits
      k=1
      ncen(0)=0
      name=munit(1)
      numb=nunit(1)
      do i=2,kam
      if(munit(i).ne.name.or.nunit(i).ne.numb) then
      ncen(k)=i-1
      k=k+1
      name=munit(i)
      numb=nunit(i)
      endif
      enddo
      ncen(k)=kam
      kcen=k
c-----------------find atomic numbers and bonding if pdb or dat format
      if(ext.eq.'pdb'.or.ext.eq.'PDB'.or.
     1   ext.eq.'brk'.or.ext.eq.'BRK'.or.
     1   ext.eq.'dat'.or.ext.eq.'DAT') then
      do i=1,kam
      name=mnam(i)
      dmon(i)=0.
      imty(i)=0
      icm(i)=0
      do j=1,55
      k=ind(j)
      if(name(:k).eq.nam(j)(:k)) then
      imch(i)=j
      goto 200
      endif
      enddo
      if(name(:1).eq.'M') then
      imch(i)=6
      goto 200
      else if(name(:1).eq.'W') then
      imch(i)=8
      goto 200
      endif
      write(6,201) i,name
201   format(2x,'---- Non-standard atom name ',i3,': ',a4,' ----')
      stop
200   enddo
      call bonder
      endif
      return
      end
