      function aver(a,b)
      include 'curves_data.inc'
      ca=cos(cdr*(a))
      sa=sin(cdr*(a))
      cb=cos(cdr*(b))
      sb=sin(cdr*(b))
      diff=acos(ca*cb+sa*sb)*crd
      if(cb*sa-sb*ca.lt.0.) diff=-diff
      aver=b+diff/2.
      return
      end
