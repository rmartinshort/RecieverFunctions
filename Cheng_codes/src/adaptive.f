      dimension x(5000),rx(5000),sx(5000),ex(5000),w(100)
      character*32 infil,outfil,cutfil
      integer ounit
      common /innout/inunit,ounit
      inunit=5
      ounit=6
      write(*,*) 'input adaptive filter length: '
      read(*,*) n
      write(*,*) 'input time delay : '
      read(*,*) ndelay
      write(*,*) 'input iterative times : '
      read(*,*) iter0
      write(*,*) 'input convergence coffcient : '
      read(*,*) coff
      do 1 i=1,5000
      x(i)=0.0
      rx(i)=0.0
      sx(i)=0.0
      ex(i)=0.0
 1    continue
      do 2 i=1,n
      w(i)=0.0
 2    continue
      k11=1
      call asktxt('specify input file :  ',infil)
      call asktxt('specify cut file :  ',cutfil)
      call asktxt('specify output file: ',outfil)
      call rsac1(infil,x,npts,beg,dt,5000,nerr)
      call getfhv('a',a0,nerr)
      beg=beg-a0
      do 3 i=ndelay,npts
      rx(i-ndelay+1)=x(i)
 3    continue
      k12=(a0+120.-beg)/dt+1
      call getdata(x,k11,k12,npts)
      call getdata(rx,k11,k12,npts)
      npts=k12-k11+1
      dmean1=0.
      dmean2=0.
      do 12 i=1,npts
      dmean1=dmean1+x(i)
      dmean2=dmean2+rx(i)
 12   continue
      dmean1=dmean1/float(npts)
      dmean2=dmean2/float(npts)
      do 13 i=1,npts
      x(i)=x(i)-dmean1
      rx(i)=rx(i)-dmean2
 13   continue
      call wsac1(cutfil,x,npts,beg,dt,nerr)
      call wsac1('refsig',rx,npts,beg,dt,nerr)
      iter=0
 5    continue
      iter=iter+1
      do 11 j=n,npts
      sx(j)=0.0
      do 10 i=1,n
      sx(j)=sx(j)+w(i)*rx(j-i+1)
 10   continue
      ex(j)=x(j)-sx(j)
 11   continue      
      do 41 i=1,n
      do 4 j=n,npts
      w(i)=w(i)+2*coff*ex(j)*rx(j-i+1)
 4    continue
 41   continue
      if(iter.le.iter0) go to 5
      call wsac1(outfil,sx,npts,beg,dt,nerr)
      call wsac1('wcoff',w,n,beg,dt,nerr)
      end
