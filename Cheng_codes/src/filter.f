       subroutine filter(x,mm,nn,npts,n,beg,dt)
       dimension x(mm,nn),a(3,3),d(3),v(3,3),b(3),z(3),x1(8000,3)
       dimension dc(8000,3),coff(8000),dc1(3)
       logical ev
       ev=.true.
       do 1 i=1,8000
       coff(i)=0.
       do 2 j=1,3
       dc(i,j)=0.
       x1(i,j)=x(i,j)
 2     continue
 1     continue
       do 50 k=1,npts-2*n
       rmean=0.
       tmean=0.
       zmean=0.
       do 20 j=k,k+2*n
       zmean=zmean+x1(j,1)
       rmean=rmean+x1(j,2)
       tmean=tmean+x1(j,3)
  20   continue
       rmean=rmean/float(2*n+1)
       tmean=tmean/float(2*n+1)
       zmean=zmean/float(2*n+1)
       do 10 i=1,3
       do 10 j=1,3
       a(i,j)=0.0
  10   continue
       do 30 j=k,k+2*n
       a(1,1)=a(1,1)+(x1(j,2)-rmean)**2
       a(2,2)=a(2,2)+(x1(j,3)-tmean)**2
       a(3,3)=a(3,3)+(x1(j,1)-zmean)**2
       a(1,2)=a(1,2)+(x1(j,3)-tmean)*(x1(j,2)-rmean)
       a(1,3)=a(1,3)+(x1(j,1)-zmean)*(x1(j,2)-rmean)
       a(2,3)=a(2,3)+(x1(j,1)-zmean)*(x1(j,3)-tmean)
 30    continue
       a(2,1)=a(1,2)
       a(3,1)=a(1,3)
       a(3,2)=a(2,3)
       do 40 i=1,3
       do 40 j=1,3
 40    a(i,j)=a(i,j)/float(2*n+1)
       call jacobi(3,ev,a,d,v,b,z,krt)
       amax1=max(d(1),d(2),d(3))
       amin=min(d(1),d(2),d(3))
       do 41 i=1,3
       if(d(i).eq.amax1) k1=i
       if(d(i).eq.amin) k3=i
       if(d(i).gt.amin.and.d(i).lt.amax1) k2=i
 41    continue
       amax2=d(k2)
       if(abs(amax1).le.0.0000001) then
       coff(k+n)=0.
       else
       coff(k+n)=1-amax2/amax1
       endif
       do 42 i=1,3
       dc(k+n,i)=v(i,k1)
 42    continue
 50    continue
       call wsac1('fcoff',coff,npts,beg,dt,nerr)
       call wsac1('rcoff',dc(1,1),npts,beg,dt,nerr)
       call wsac1('tcoff',dc(1,2),npts,beg,dt,nerr)
       call wsac1('zcoff',dc(1,3),npts,beg,dt,nerr)
       nn0=n/2.
       do 60 k=1,npts-3*n
       coff1=0.
       dc1(1)=0.0
       dc1(2)=0.0
       dc1(3)=0.0
       do 70 k0=k,k+n
       coff1=coff1+coff(k0+n)
       do 80 j0=1,3
       dc1(j0)=dc1(j0)+dc(k0+n,j0)
 80    continue
 70    continue
       coff1=coff1/float(n+1)
       do 90 j0=1,3
       dc1(j0)=dc1(j0)/float(n+1)
 90    continue
       x(k+n+nn0,1)=x1(k+n+nn0,1)*coff1*dc1(3)
       x(k+n+nn0,2)=x1(k+n+nn0,2)*coff1*dc1(1)
       x(k+n+nn0,3)=x1(k+n+nn0,3)*coff1*dc1(2)
 60    continue
       do 61 k=1,n+nn0
       x(k,1)=0.
       x(k,2)=0.
       x(k,3)=0.
 61    continue
       do 71 k=npts-3*n+1,mm
       x(k,1)=0.
       x(k,2)=0.
       x(k,3)=0.
 71    continue
       return
       end
