      subroutine setdata(x,n1,n2,n)
      dimension x(n),y(6000)
      do 10 i=1,n1-1
      y(i)=0.0
   10 continue
      do 25 i=n1,n2
      k0=i-n1+1
      y(i)=x(k0)
   25 continue
      if(n2.eq.n) go to 35
      do 30 i=n2+1,n
      y(i)=0.0
   30 continue
   35 do 20 i=1,n
      x(i)=y(i)
   20 continue
      return
      end

       

