      program theoryqn
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of langston (1979)
c
c    *** written by t.j. owens - july 1982
c    *** modified by G Randall for gaussian normalization
c    *** modified by C Ammon for water-level normaliztion 880711
c
c  *************************************************************************
c
      dimension d2(4000)
      character eqfile*32,outfil*32,comp(2)*2,outc(2)*4
      complex data(2000,2)
      double precision gnorm
      logical yes,yesno
      integer blank,ounit
c **********************************************************************
c
c   parameter definitions may be found in sacio comments
c
      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.r'/
     *    ,outc/'.eqz','.eqr'/
      inunit=5
      ounit=6
      pi=3.141592654
  2   outfil='                                '
      eqfile='                                '
      call asktxt('Specify event name: ',eqfile)
      iblank=blank(eqfile)
      do 1 i=1,2
         call zero(data(1,i),1,4000)
         eqfile(1:iblank + 3)=eqfile(1:iblank)//comp(i)
         call rsac1(eqfile,data(1,i),npts,beg,dt,5000,nerr)
    1 continue
c
c    *******************************************************************
c
      yes=yesno('Window data (y or n)? ')
      if(.not.yes) go to 4
         blen=ask('Length of window (in secs): ')
         call window(blen,npts,dt,3)
    4 nft=npowr2(npts)
      nfpts=nft/2 + 1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)
c
c  change delf to normalize deconvolution so dfftr works
c  dt from forward step is cancelled in decon, so no delf
c  on inverse, just 1/nft
c
      cdelf = 1. / float( nft )
      do 5 i=1,2
         call dfftr(data(1,i),nft,'forward',dt)
         if(i.ne.1) go to 5
         d2max=0.
         do 7 j=1,nfpts
            d2(j)=real(data(j,i)*conjg(data(j,i)))
            if(d2(j).gt.d2max) d2max=d2(j)
    7    continue
    5 continue
      c=ask('Trough filler, c =  ')
      agauss=ask('Gaussian scale, a = ')
      phi1=c*d2max
      gnorm = 0.0d0
         do 9 j=1,nfpts
            freq=float(j-1)*delf
            w=2.*pi*freq
            phi=phi1
            if(d2(j).gt.phi) phi=d2(j)
            gauss=-w*w/(4.*agauss*agauss)
	    gnorm = gnorm + exp( gauss )
            data(j,2)=data(j,2)*conjg(data(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
            data(j,2)=data(j,2)*exp(cmplx(0.,w*beg))
    9    continue
         call dfftr(data(1,2),nft,'inverse',cdelf)
c
c     deconvolve the vertical from itself using the
c         specified water-level parameter and gaussian
c
      do 19 j=1,nfpts
	 freq=float(j-1)*delf
	 w=2.*pi*freq
	 phi=phi1
	 gauss=-w*w/(4.*agauss*agauss)
	 if(d2(j).gt.phi) phi=d2(j)
	 data(j,1)=data(j,1)*conjg(data(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
	 data(j,1)=data(j,1)*exp(cmplx(0.,w*beg))
19    continue

c
c*************************************************************
c
c     inverse transform the equalized vertical component
c
      call dfftr(data(1,1),nft,'inverse',cdelf)

c
c     compute the maximum value of the vertical component
c        to include it in the normalization factor gnorm
c
      call minmax(data(1,1),npts,dmin,dmax,dmean)
c
c     normalize the dmax for the transforms and gaussian
c
c     dmax = dmax *  float(nfpts) / gnorm
c
c************************************************************* 
c
c     gnorm = dmax * gnorm / nfpts
c
c     note that (dmax * gnorm / nfpts) = the unormalized dmax
c
      gnorm = dmax
      do 111 j = 1,npts
	 data(j,2) = data(j,2) / gnorm
 111  continue
           outfil(1:iblank + 7)=eqfile(1:iblank)//'_qn'//outc(2)
           call wsac1(outfil,data(1,2),npts,beg,dt,nerr)
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
      subroutine window(b0,npts,dt,ndat)
      dimension data(4000,2)
      common /win/ data
      data pi/3.1415926/
      bb=pi/b0
      nend=ifix((b/dt)+.5) +1
      do 1 i=1,nend
      t=float(i-1)*dt
      windo=.5*(1. + cos(bb*t+pi))
         do 2 j=1,ndat
         data(i,j)=data(i,j)*windo
         data(npts+1-i,j)=data(npts+1-i,j)*windo
    2    continue
    1 continue
      return
      end
