      subroutine minfor(move,n,x,f,g,scale,acc,h,d,w,xa,ga,xb,gb,maxn)
      include 'curves_data.inc'
      dimension x(*),g(*),scale(*),h(*),d(*),w(*),xa(*),
     1 ga(*),xb(*),gb(*)
      external move
      nfun=0
      call move(x,f,g)
      nfun=1
      itr=0
      np=n+1
c     set the hessian to a diagonal matrix depending on scale(.)
      c=0.
      do 30 i=1,n
30    c=max(c,abs(g(i)*scale(i)))
      if (c.le.0.) c=1.
      k=(n*np)/2
      do 40 i=1,k
40    h(i)=0.
      k=1
      do 50 i=1,n
      h(k)=0.01*c/scale(i)**2
50    k=k+np-i
c     set some variables for the first iteration
      dff=0.
110   fa=f
      isfv=1
      do 120 i=1,n
      xa(i)=x(i)
120   ga(i)=g(i)
c     begin the iteration by giving the required printing
130   itr=itr+1
c     calculate the search direction of the iteration
      do 150 i=1,n
150   d(i)=-ga(i)
      call mc11e (h,n,d,w,n)
c     calculate a lower bound on the step-length
c     and the initial directional derivative
      c=0.
      dga=0.
      do 160 i=1,n
      c=max(c,abs(d(i)/scale(i)))
160   dga=dga+ga(i)*d(i)
c     test if the search direction is downhill
      if (dga.ge.0.) go to 240
c     set the initial step-length of the line search
      stmin=0.
      stepbd=0.
      steplb=acc/c
      fmin=fa
      gmin=dga
      step=1.
      if (dff.le.0.) step=min(step,1./c)
      if (dff.gt.0.) step=min(step,(dff+dff)/(-dga))
170   c=stmin+step
c     test whether func has been called maxn times
      if (nfun.ge.maxn) go to 250
      nfun=nfun+1
c     calculate another function value and gradient
      do 180 i=1,n
180   xb(i)=xa(i)+c*d(i)
      call move(xb,fb,gb)
c     store this function value if it is the smallest so far
      isfv=min(2,isfv)
      if (fb.gt.f) go to 220
      if (fb.lt.f) go to 200
      gl1=0.
      gl2=0.
      do 190 i=1,n
      gl1=gl1+(scale(i)*g(i))**2
190   gl2=gl2+(scale(i)*gb(i))**2
      if (gl2.ge.gl1) go to 220
200   isfv=3
      f=fb
      do 210 i=1,n
      x(i)=xb(i)
210   g(i)=gb(i)
c     calculate the directional derivative at the new point
220   dgb=0.
      do 230 i=1,n
230   dgb=dgb+gb(i)*d(i)
c     branch if we have found a new lower bound on the step-length
      if (fb-fa.le.0.1*c*dga) go to 280
c     finish the iteration if the current step is steplb
      if (step.gt.steplb) go to 270
240   if (isfv.ge.2) go to 110
c     at this stage the whole calculation is complete
250   if(nfun.lt.maxn) call move(x,f,g)
      return
c     calculate a new step-length by cubic interpolation
270   stepbd=step
      c=gmin+dgb-3.*(fb-fmin)/step
      c=gmin/(c+gmin-sqrt(c*c-gmin*dgb))
      step=step*max(0.1d0,c)
      go to 170
c     set the new bounds on the step-length
280   stepbd=stepbd-step
      stmin=c
      fmin=fb
      gmin=dgb
c     calculate a new step-length by extrapolation
      step=9.*stmin
      if (stepbd.gt.0.) step=0.5*stepbd
      c=dga+3.*dgb-4.*(fb-fa)/stmin
      if (c.gt.0.) step=min(step,stmin*max(1.d0,-dgb/c))
      if (dgb.lt.0.7*dga) go to 170
c     test for convergence of the iterations
      isfv=4-isfv
      if (stmin+step.le.steplb) go to 240
c     revise the second derivative matrix
      ir=-n
      do 290 i=1,n
      xa(i)=xb(i)
      xb(i)=ga(i)
      d(i)=gb(i)-ga(i)
290   ga(i)=gb(i)
      call mc11a (h,n,xb,1./dga,w,ir,1,0.d0)
      ir=-ir
      call mc11a (h,n,d,1./(stmin*(dgb-dga)),d,ir,0,0.d0)
c     branch if the rank of the new matrix is deficient
      if (ir.lt.n) go to 250
c     begin another iteration
      dff=fa-fb
      fa=fb
      go to 130
      end
