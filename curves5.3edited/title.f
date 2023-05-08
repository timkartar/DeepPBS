      subroutine title(l,text,n)
      integer*4 n
      character*60 text*(*),line,l*1
      data line/
     1 '------------------------------------------------------------'/
      write(6,10) line(:n+6),l,text(:n),line(:n+6)
10    format(/2x,a,/2x,'|',a1,'| ',a,' |',/2x,a)
      return
      end
