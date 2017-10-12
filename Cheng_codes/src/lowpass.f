      subroutine lowpass(data,npts,dt,f0)
      dimension data(9000,3)
      complex d(9000,3)
      intrinsic aimag
      pi=3.1415926
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      nlowf=f0*nfpts/fny+1
c      write(*,*) 'nlowf=', nlowf,'nfpts=',nfpts,'fny=',fny,
c     * 'cdelf=',cdelf,'delf=',delf,'nft=',nft
      do 8 i=1,3
c      write(*,*) 'i=',i
      call dfftr(data(1,i),nft,'forward',dt)
c      write(*,*) 'i=',i
      do 9 j=1,nlowf
c      freq=float(j-1)*delf
c      w=2.*pi*freq
      cc1=data(2*(j-1)+1,i)
      cc2=data(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
   9  continue
      do 10 j=nlowf+1,nfpts
      d(j,i)=cmplx(0.,0.)
  10  continue
      call dfftr(d(1,i),nft,'inverse',cdelf)
      do 111 j=1,nfpts
      data(2*(j-1)+1,i)=real(d(j,i))
      data(2*j,i)=aimag(d(j,i))
 111  continue
   8  continue
      return
      end
