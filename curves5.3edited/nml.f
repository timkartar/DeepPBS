      subroutine nml
      include 'curves_data.inc'
      parameter (n_real=2,n_int=7,n_log=13,n_cha=7,n_tot=29)
      character*32 vc,input(n_tot)*10,line*1000,day*9,date*24,lj*1
      logical*2 vo,iflag(n_tot),first,last,start
      integer*4 vi,nmls(n_tot)
      common/cha/vc(n_cha)
      common/drl/vr(n_real),vi(n_int),vo(n_log)
      data input/
     1 'acc','wid','maxn','ior','ibond','spline',
     1 'break','nlevel','nbac','ends','supp','comb','dinu',
     1 'mini','rest','line','zaxe','fit','test','grv','old',
     1 'axonly','file','lis','dna','axin','axout','daf','pdb'/
      nint=n_int+n_real
      nlog=n_log+nint
         do i=1,n_tot
         iflag(i)=.false.
         nmls(i)=index(input(i),' ')-1
         enddo
      first=.true.
      last=.false.
 10   read(5,5) line
 5    format(a)
      im=index(line,'&')
      if(im.gt.0) then
      if(.not.first) last=.true.
      if(first.and.index(line(im+1:),'&').ne.0) last=.true.
      endif
      do k=1,100
      if(line(k:k).eq.'=') then
      kl=k
      start=.true.
      do j=k-1,1,-1
      lj=line(j:j)
      if(start.and.lj.ne.' ') then
      start=.false.
      jh=j
      else if(.not.start.and.(lj.eq.' '.or.lj.eq.',')) then
      jl=j+1
      goto 15
      endif
      enddo
      goto 50
 15      do i=1,n_tot
         if(line(jl:jh).eq.input(i)) then
         iflag(i)=.true.
 17         kl=kl+1
            if(line(kl:kl).eq.' ') goto 17
            do j=kl,100
            lj=line(j:j)
            if(lj.eq.' '.or.lj.eq.','.or.lj.eq.'&') then
            kh=j-1
            goto 19
            endif
            enddo
 19      if(i.le.n_real) then
         read(line(kl:kh),*,err=50) vr(i)
         goto 25
         else if(i.le.nint) then
         read(line(kl:kh),*,err=50) vi(i-n_real)
         goto 25
         else if(i.le.nlog) then
         read(line(kl:kh),*,err=50) vo(i-nint)
         goto 25
         else
         if(line(kl:kl).eq.'''') kl=kl+1
         if(line(kh:kh).eq.'''') kh=kh-1
         if(kh.ge.kl) read(line(kl:kh),5,err=50) vc(i-nlog)
         goto 25
         endif
         endif
         enddo
         goto 50
      endif
 25   enddo
      first=.false.
      if(.not.last) goto 10
c-------------------------------------------------------------------------output
      if(vc(2).ne.' ') then
      kfi=index(vc(2),' ')-1
      open(unit=6,file=vc(2)(:kfi)//'.lis',status='new')
      endif
      call fdate(date)
      day=date(9:11)//date(5:8)//date(23:)
      write(6,200) day
200   format(
     1/5x,'***********************************',15x,'**************',
     1/5x,'******  CURVES 5.3 R.L. 1998  *****',15x,'*  ',a9,   ' *',
     1/5x,'***********************************',15x,'**************',
     1 //)
      do i=1,n_tot
      nm=nmls(i)
      if(iflag(i)) then
      do j=1,nm
      ic=ichar(input(i)(j:j))
      input(i)(j:j)=char(ic-32)
      enddo
      endif
      enddo
      write(6,8) (input(nlog+j),vc(j),j=1,n_cha)
8     format(2x,a5,': ',a32,:,2x,a5,': ',a32)
      write(6,*)
      write(6,12) (input(j),vr(j),j=1,n_real)
12    format(4(:,2x,a5,': ',f8.3))
      write(6,*)
      write(6,14) (input(n_real+j),vi(j),j=1,n_int)
14    format(5(:,2x,a5,': ',i5))
      write(6,*)
      write(6,16) (input(nint+j),vo(j),j=1,n_log)
16    format(5(:,2x,a5,': ',l5))
c-------------------------------------------------------------------------------
      return
50    write(6,55) line(jl:jh)
55    format(/2x,'---- Error in namelist input for ',a,' ----'/)
      stop
      end
