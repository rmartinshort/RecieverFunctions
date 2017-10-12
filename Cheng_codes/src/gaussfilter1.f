c      subroutine gaussfilter(data,npts,dt,agauss)
c *********************************************************
c *    gaussfilter apply Gaussian filter to data          *
c *********************************************************
      dimension data(2000)
      complex d(4500)
      intrinsic aimag
 
      write(*,*) 'input gauss:'
      read(*,*) agauss
      dt=0.1
      beg = -10.
      npts=1024
      pi=3.1415926
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      write(*,*) dt,nft,fny,delf,cdelf
c      call dfftr(data,nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      gauss=-w*w/(4.*agauss*agauss)
      d(j)=exp(gauss)*exp(cmplx(0.,w*beg))
    9 continue
      call dfftr(d,nft,'inverse',cdelf)
      call minmax(d(1),npts,dmin,dmax,dmean)
      do 111 j=1,nfpts
      d(j)=d(j)/dmax
 111  continue
      call wsac1('filter.r',d(1),npts,beg,dt,nerr)
c      do 111 j=1,nfpts
c      data(2*(j-1)+1)=real(d(j))
c      data(2*j)=aimag(d(j))
c 111  continue
      end
