      program wienerdecon
c
c  ********************************************************************
c  *                                                                  *
c  *  fortran 77 program to perform deconvolution on surface waves    *
c  *  recorded on two stations to obtain interstation Green function  *
c  *  - using the method of Wiener filter based on Levinson algorithm *
c  *                                                                  *
c  ********************************************************************
c
      dimension rzz(4000),rzr(4000)
      dimension a(4000),rr(4000),zr(4000),rzz1(4000)
      complex data(4500,2)
      character*32 evnam1,evnam2,outfil
      logical yes,yesno
      integer ounit
c **********************************************************************
c
c
      common /win/ data
      common /innout/ inunit,ounit
      inunit=5
      ounit=6
      pi=3.1415926
    2 call asktxt('specify station 1 name: ',evnam1)
      call asktxt('specify station 2 name: ',evnam2)
      call asktxt('specify output file: ',outfil)
      agauss=ask('Gaussian scale, a= ')
      t0=ask('time delay, t0= ')
         call zero(data(1,1),1,9000)
         call zero(data(1,2),1,9000)
         call rdata(evnam1,data(1,1),npts,dt)
         call rdata(evnam2,data(1,2),npts,dt)
      yes=yesno('window data (y or n)?')
      if(.not.yes) go to 4
        blen=ask('length of window (in Secs): ')
        call window(blen,npts,dt,2)
 4    nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 18 i=1,2
      call dfftr(data(1,i),nft,'forward',dt)
      do 19 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      gauss=-1.*w*w/(4.*agauss*agauss)
      data(j,i)=data(j,i)*exp(gauss)
   19 continue
      call dfftr(data(1,i),nft,'inverse',cdelf)
   18 continue
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)
c
c    *******************************************************************
c
      call crossfun1(data(1,1),data(1,1),rzz,npts,npts)
      call crossfun1(data(1,1),data(1,2),rzr,npts,npts)
      do 101 i=1,npts
      rzz1(i)=rzz(i)
 101  continue
      a(1)=1./rzz(1)
      rr(1)=rzr(1)/rzz(1)
      zr(1)=rzz1(1)/rzz(1)
      do 10 i=2,npts
      a(i)=0.
      rr(i)=0.
      zr(i)=0.
      coff=0.
      err=0.
      ezr=0.
      do 20 j=2,i
      coff=coff+a(i-j+1)*rzz(j)
      err=err+rr(i-j+1)*rzz(j)
      ezr=ezr+zr(i-j+1)*rzz(j)
  20  continue
      bot=1-coff*coff
      ih=(i+1)/2
      do 30 k=1,ih
      top=(a(i-k+1)-coff*a(k))/bot
      a(k)=(a(k)-coff*a(i-k+1))/bot
      a(i-k+1)=top
  30  continue
      do 50 k=1,i
      rr(k)=rr(k)+(rzr(i)-err)*a(i-k+1)
      zr(k)=zr(k)+(rzz1(i)-ezr)*a(i-k+1)
  50  continue
  10  continue
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 70 i=1,2
      call zero(data(1,i),1,9000)
70    continue
      do 60 i=1,nfpts
      cc1=zr(2*(i-1)+1)
      cc2=zr(2*i)
      cc3=rr(2*(i-1)+1)
      cc4=rr(2*i)
      data(i,1)=cmplx(cc1,cc2)
      data(i,2)=cmplx(cc3,cc4)
  60  continue
      do 8 i=1,2
      call dfftr(data(1,i),nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      gauss=-1.*w*w/(4.*agauss*agauss)
      data(j,i)=data(j,i)*exp(gauss)*exp(cmplx(0.,w*t0))
   9  continue
      call dfftr(data(1,i),nft,'inverse',cdelf)
 8    continue
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 80 j=1,npts
      data(j,2)=data(j,2)/dmax
 80   continue
        call wdata(outfil,data(1,2),npts,dt)
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
      subroutine window(b,npts,dt,ndat)
      dimension data(9000,2)
      common  /win/data
      data pi/3.1415926/
      bb=pi/b
      nend=ifix(b/dt+.5)+1
      do 1 i=1,nend
      t=float(i-1)*dt
      windo=0.5*(1.+cos(bb*t+pi))
      do 2 j=1,ndat
      data(i,j)=data(i,j)*windo
      data(npts+1-i,j)=data(npts+1-i,j)*windo
    2 continue
    1 continue
      return
      end
