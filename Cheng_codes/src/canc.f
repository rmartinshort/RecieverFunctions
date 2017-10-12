      subroutine canc(x,npts,n,ndelay,iter0,coff)
      dimension x(1),rx(5000),sx(5000),ex(5000),w(100)
      do 1 i=1,5000
      rx(i)=0.0
      sx(i)=0.0
      ex(i)=0.0
 1    continue
      do 2 i=1,n
      w(i)=0.0
 2    continue
      k11=1
      do 3 i=ndelay,npts
      rx(i-ndelay+1)=x(i)
 3    continue
      dmean1=0.
      dmean2=0.
      do 12 i=1,npts-ndelay+1
      dmean1=dmean1+x(i)
      dmean2=dmean2+rx(i)
 12   continue
      dmean1=dmean1/float(npts)
      dmean2=dmean2/float(npts)
      do 13 i=1,npts-ndelay+1
      x(i)=x(i)-dmean1
      rx(i)=rx(i)-dmean2
 13   continue
      iter=0
 5    continue
      iter=iter+1
      do 11 j=n,npts-ndelay+1
      sx(j)=0.0
      do 10 i=1,n
      sx(j)=sx(j)+w(i)*rx(j-i+1)
 10   continue
      ex(j)=x(j)-sx(j)
 11   continue      
      do 41 i=1,n
      do 4 j=n,npts-ndelay+1
      w(i)=w(i)+2*coff*ex(j)*rx(j-i+1)
 4    continue
 41   continue
      if(iter.le.iter0) go to 5
      do 42 i=1,npts-ndelay+1
      x(i)=sx(i)
 42   continue
      do 43 i=npts-ndelay+2,npts
      x(i)=0.0
 43   continue
      end
