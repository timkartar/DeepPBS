      subroutine setd
      include 'curves_data.inc'
      character*4 snam,sunit
      common/axe/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),matc(n1,7),kas,khs,kces
      do i=1,kas
      do j=1,7
      matc(i,j)=0
      enddo
      enddo
      do i=1,khs-1
      k=mats(i)
      kn=abs(k)
      if(k.lt.10000) then
      ik=matc(kn,7)
      kp=mats(i+1)
      kpn=abs(kp)
      if(kp.lt.10000) then
      ikp=matc(kpn,7)
      ik=ik+1
      matc(kn,ik)=kp
      matc(kn,7)=ik
      ikp=ikp+1
      matc(kpn,ikp)=sign(k,kp)
      matc(kpn,7)=ikp
      endif
      if(kp.eq.20000) then
      kpp=mats(i+2)
      kppn=abs(kpp)
      ikpp=matc(kppn,7)
      ik=ik+1
      matc(kn,ik)=kpp
      matc(kn,7)=ik
      ikpp=ikpp+1
      matc(kppn,ikpp)=sign(k,kpp)
      matc(kppn,7)=ikpp
      endif
      endif
      enddo
      return
      end
