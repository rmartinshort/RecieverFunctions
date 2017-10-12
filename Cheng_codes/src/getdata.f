      subroutine getdata(x,k1,k2,npts)
      dimension x(1)
      do 25 i=k1,k2
      k0=i-k1+1
      x(k0)=x(i)
   25 continue
      do 10 i=k2-k1+2,npts
      x(i)=0.
   10 continue
      return
      end
