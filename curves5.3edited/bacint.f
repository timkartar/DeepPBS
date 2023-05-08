      subroutine bacint(bato)
      include 'curves_data.inc'
      logical*2 ends,supp,comb,dinu,mini,rest,line,zaxe,fit,
     1 flag,test,grv,old
      integer*4 spline,break
      character*4 mnam,munit,na*1,bato
      dimension px(0:n3),py(0:n3),pz(0:n3),tx(n3),ty(n3),tz(n3)
      common/bak/tor(0:n3,4,13),val(0:n3,4,2),suga(0:n3,4,2),
     1 nat(n3,4,15),flag(n3,4)
      common/dat/rex(0:n3,4,4),rey(0:n3,4,4),rez(0:n3,4,4),
     1 hel(0:n3,6,4),idr(4),ni(0:n3,4),ng(4),nr(4),nu(4),nux,
     1 nt,nst,iact(0:n3),li(0:n3,4),na(n3,4)
      common/drl/acc,wid,maxn,ior,ibond,spline,break,nlevel,
     1 nbac,ends,supp,comb,dinu,mini,rest,line,zaxe,fit,test,grv,old
      common/gro/uxb(5000,0:n3),uyb(5000,0:n3),uzb(5000,0:n3),
     1 cor(n1,3),dya(n1,3),box(0:n1,4),boy(0:n1,4),
     1 boz(0:n1,4),ind(n3,0:n3),nma(n3),num,numa,nsu(4)
      common/mac/corm(n1,3),dmon(n1),mnam(n1),munit(0:n1),
     1 imch(n1),imty(n1),icm(n1),matd(n1,7),nunit(0:n1),
     1 ncen(0:n2),kam,khm,lkm,kcen
      do k=1,2
      ix=numa
      le=0
      do i=num,numa
      l=nat(i,k,nbac)
      if(le.eq.0.and.l.ne.0) im=i
      if(le.ne.0.and.l.eq.0) ix=i-1
      le=l
      enddo
      nsu(k)=20*(ix-im)
      do i=im,ix
      l=nat(i,k,nbac)
      px(i)=corm(l,1)
      py(i)=corm(l,2)
      pz(i)=corm(l,3)
      enddo
c-------------------------------------------tangents on inner back-atoms
      do i=im+1,ix-1
      a=(px(i+1)-px(i))**2+(py(i+1)-py(i))**2+(pz(i+1)
     1 -pz(i))**2
      b=(px(i-1)-px(i))**2+(py(i-1)-py(i))**2+(pz(i-1)
     1 -pz(i))**2
      tgx=a*(px(i)-px(i-1))+b*(px(i+1)-px(i))
      tgy=a*(py(i)-py(i-1))+b*(py(i+1)-py(i))
      tgz=a*(pz(i)-pz(i-1))+b*(pz(i+1)-pz(i))
      tm=sqrt(tgx**2+tgy**2+tgz**2)
      tx(i)=tgx/tm
      ty(i)=tgy/tm
      tz(i)=tgz/tm
      enddo
c----------------------------------------tangents on terminal back-atoms
      do i=im,ix,ix-im
      if(i.eq.im) l=1
      if(i.eq.ix) l=-1
      if(i.eq.ix) nua=numa-1
      if(dinu) then
      if(i.eq.im) nua=num+2
      a=hel(nua,6,k)+hel(nua-l,6,k)
      uax=uxb(nua,0)+uxb(nua-1,0)
      uay=uyb(nua,0)+uyb(nua-1,0)
      uaz=uzb(nua,0)+uzb(nua-1,0)
      um=sqrt(uax**2+uay**2+uaz**2)
      uax=uax/um
      uay=uay/um
      uaz=uaz/um
      ax=tx(i+2*l)
      ay=ty(i+2*l)
      az=tz(i+2*l)
      else
      if(i.eq.im)nua=num+1
      a=(hel(nua,6,k)+hel(nua+1,6,k))/2.
      uax=uxb(nua,0)
      uay=uyb(nua,0)
      uaz=uzb(nua,0)
      ax=tx(i+l)
      ay=ty(i+l)
      az=tz(i+l)
      endif
      ap=ax*uax+ay*uay+az*uaz
      ca=cos(cdr*(a))
      sa=sin(cdr*(a))
      tx(i)=ap*uax+(ax-ap*uax)*ca+(ay*uaz-az*uay)*sa*l
      ty(i)=ap*uay+(ay-ap*uay)*ca+(az*uax-ax*uaz)*sa*l
      tz(i)=ap*uaz+(az-ap*uaz)*ca+(ax*uay-ay*uax)*sa*l
      enddo
c--------------------------------------------------backbone space curves
      do i=im,ix-1
      pxi=px(i)
      pyi=py(i)
      pzi=pz(i)
      pxs=px(i+1)
      pys=py(i+1)
      pzs=pz(i+1)
      txi=tx(i)
      tyi=ty(i)
      tzi=tz(i)
      txs=tx(i+1)
      tys=ty(i+1)
      tzs=tz(i+1)
      d=sqrt((pxs-pxi)**2+(pys-pyi)**2+(pzs-pzi)**2)
      fx= 3*(pxs-pxi)-(txs+2*txi)*d
      fy= 3*(pys-pyi)-(tys+2*tyi)*d
      fz= 3*(pzs-pzi)-(tzs+2*tzi)*d
      gx=-2*(pxs-pxi)+(txs+  txi)*d
      gy=-2*(pys-pyi)+(tys+  tyi)*d
      gz=-2*(pzs-pzi)+(tzs+  tzi)*d
      do m=0,19
      in=20*(i-im)+m
      r=m/20.
      box(in,k)=pxi+txi*d*r+fx*r**2+gx*r**3
      boy(in,k)=pyi+tyi*d*r+fy*r**2+gy*r**3
      boz(in,k)=pzi+tzi*d*r+fz*r**2+gz*r**3
c---------------------------------------------------------last back-atom
      enddo
      enddo
      in=20*(ix-im)
      box(in,k)=px(ix)
      boy(in,k)=py(ix)
      boz(in,k)=pz(ix)
      enddo
      return
      end
