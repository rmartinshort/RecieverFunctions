c      program zdeco.f
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of maximum 
c     entrop deconvolution 
c
c    *** modifying by zty - June 2001 ***
c
c  *************************************************************************
c      bptime=('pwave beginning time before A (in Sec positive): ')
c      eptime=('pwave ending time after A (in Sec): ')
c      tpw= ('pwave source time window: ')
c      iow>0 or=0 ('window data (y or n)?')
c      blen= ('length of window (in Sec),blen = ')
c      iom>0 or =0 ('remove mean value (y or n)?' )
c      iof>0 or =0 ('filter raw data (y or n)?' ) 
c      agauss= ('Gaussian scale, a= ')
c      c= ('trough filter,c = ')
C	BEG=BEGIN POINT OF TIME
C	A0=ARRIVE TIME OF DIRECT P WAVE
	
      dimension dat(9000,3),caz(3),amax(3)
      complex d(4500,3)
      dimension ef(4000),eb(4000),sz(4000),sr(4000),st(4000),hz(4000)
      dimension hr(4000),ht(4000),a1(4000),a2(4000)
c      character eqfile*80,comp(3,2)*2,outc(3)*4
      character*32 outfil,eqfile,IAJ

c      logical yes,yesno,rel,wrot, rmean, water,filter
c      integer blank,ounit
c      intrinsic aimag
c      real mean1,mean2,mean3
c
c
c **********************************************************************
c
c
      common /win/ dat
      common /innout/ inunit,ounit
c      data
c     *     comp/'.z','.n','.e','.z','.r','.t'/
c     *    ,outc/'.eqz','.eqr','.eqt'/
      inunit=5
      ounit=6
      pi=3.1415926

      open (5,file='ddeco',status='old')
2	continue
	read(*,500) IAJ
	eqfile='F'//IAJ
	outfil='RF'//IAJ
500	format(a12)	

	read(5,*) bptime,eptime,blen
	read(5,*) agauss,c
	read(5,*) (caz(k),k=1,3)
	read(5,*) beg,a0
	read(5,*) iow,iof,ior,iom
	read(5,*) iowp,iowf,iowr
	close(5)

c         call zero(dat(1,i),1,9000)
	do 61 k=1,3
	do 61 i=1,9000
61	dat(i,k)=0.	

c         call rsac1(eqfile,data(1,i),npts,beg,dt,9000,nerr)
c         call getfhv('cmpaz',caz(i),nerr)
	open(7,file=eqfile,status='old')
	read(7,*) nzs
	do 62 k=1,nzs
	read(7,*) npts,dt,baz,amax(K)
	read(7,*) (dat(i,k),i=1,npts)
62	continue	
	close(7)	
	do 1 j=1,npts
	yz=dat(j,3)
	yn=dat(j,1)
	ye=-dat(j,2)
	dat(j,1)=yz
	dat(j,2)=yn
	dat(j,3)=ye
  1   continue

	if (iowp.ne.0) then
	open(8,file='PRE',status='unknown')
	write(8,*) nzs
       do 411 i=1,nzs
	ymax=0.
	do 465 j=1,npts
	xx=abs(dat(j,i))
	if (xx.gt.ymax) ymax=xx
465	continue
	write(8,*) npts,dt,ymax
	write(8,*) (dat(j,i),j=1,npts)
  411   continue
	close(8) 		
	end if
c
c     remove the mean value from the raw data
c
       if(iom.eq.0) goto  1000

      do 22 i=1,3
c      do 70 j=1,npts
c70	dc(j)=dat(j,i)      
      call minmax(dat(1,i),npts,dmin,dmax,dmean)
      do 23 j=1,npts
      dat(j,i)=dat(j,i)-dmean
23    continue
22    continue

1000  continue

c         call getfhv('a',a0,nerr)
c         call getfhv('baz',baz,nerr)
      if(iow.eq.0) go to 4
         call window(blen,npts,dt,3)
 4    continue

c 
c     filter the high frequency from the raw data with Gaussian filter
c
      if(iof.eq.0) go to 5
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 8 i=1,3
c      do 81 j=1,npts
c81	dc(j)=dat(j,i)      
      call dfftr(dat(1,i),nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=dat(2*(j-1)+1,i)
      cc2=dat(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(gauss)
    9 continue
c      do 82 j=1,npts
c82	dcc(j)=d(j,i)      
      call dfftr(d(1,i),nft,'inverse',delf)
      do 111 j=1,nfpts
      dat(2*(j-1)+1,i)=real(d(j,i))
      dat(2*j,i)=aimag(d(j,i))
 111  continue
    8 continue

c      write(*,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)

	if (iowf.ne.0) then
	open(8,file='FRE',status='unknown')
	write(8,*) nzs
       do 211 i=1,nzs
	ymax=0.
	do 265 j=1,npts
	xx=abs(dat(j,i))
	if (xx.gt.ymax) ymax=xx
265	continue
	write(8,*) npts,dt,ymax
	write(8,*) (dat(j,i),j=1,npts)
  211   continue
	close(8) 		
	end if
c
c      cut off the expected P-wave data from the raw data
c
5      k11=1
c	k11=(a0-bptime)/dt+1
c	k12=(eptime+a0)/dt+1
       if(bptime.gt.0.) k11=(a0-bptime-beg)/dt+1
       k12=(eptime+a0-beg)/dt+1
       beg=beg-a0
c	beg=-a0
       if(bptime.gt.0.) beg=-bptime
       do 3 i=1,3
       call getdata(dat(1,i),k11,k12,npts)
  3    continue
       npts=k12-k11+1

c
c      rotate the raw data in ZNE direction to ZRT direction froming
c      radial and tangent component
c
       if(ior.ne.0) call rotate(dat,9000,3,baz,caz,npts)

	if (iowr.ne.0) then
 	open(8,file='RRE',status='unknown')
	write(8,*) nzs
        do 311 i=1,nzs
	ymax=0.
	do 365 j=1,npts
	xx=abs(dat(j,i))
	if (xx.gt.ymax) ymax=xx
365	continue
	write(8,*) npts,dt,ymax
	write(8,*) (dat(j,i),j=1,npts)
311     continue
	close(8) 		
	end if
c       if(.not.wrot) go to 133
c       do 103 i=1,3
c       outfil(1:iblank1+2)=evnam(1:iblank1)//comp(i,2)
c       call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
c 103   continue
c 133   continue
c
c    *******************************************************************
c
c     now to calculate predective error filter factor with Levision 
c     iterative algorithm

c
c     Initalize
c
      a1(1)=1.
      auto0=0.
      crossr0=0.
      crosst0=0.
      do 10 i=1,npts
      hz(i)=0.
      hr(i)=0.
      ht(i)=0.
      a2(i)=0.
c
c     autocorrelation of vertical component and crosscorrelation 
c     between vertical, radial, and tangent components
c
      auto0=auto0+dat(i,1)*dat(i,1)
      crossr0=crossr0+dat(i,1)*dat(i,2)
      crosst0=crosst0+dat(i,1)*dat(i,3)
c   
c     preliminary forward and backward predective error
c
      ef(i)=dat(i,1)
      eb(i)=dat(i,1)
   10 continue
      stable=c*auto0
	write(*,*) 'cr',auto0,crossr0,crosst0
c
c     preliminary receiver funtion value (1th order)
c
      hz(1)=1.
      hr(1)=crossr0/auto0
      ht(1)=crosst0/auto0
	write(*,*)hz(1),hr(1),ht(1)
c
c     preliminary filter residual
c
      do 20 i=1,npts
      sz(i)=dat(i,1)-hz(1)*dat(i,1)
      sr(i)=dat(i,2)-hr(1)*dat(i,1)
      st(i)=dat(i,3)-ht(1)*dat(i,1) 
   20 continue
c
c     Begin recurive .........
c
      do 30 i=2,npts
c
c     calculate k+1 order reflection coefficient rk from k order
c
      top=0.
      bot=0.
      do 45 j=i,npts
      bot=bot+ef(j)*ef(j)+eb(j-i+1)*eb(j-i+1)
      top=top+eb(j-i+1)*ef(j)
   45 continue
      bot=bot+stable
      rk=2.*top/bot
c
c     calculate k+1 order predective error filter factor a2(j) form 
c      k orer
c
      a1(i)=0.
      do 60 j=1,i
      a2(j)=a1(j)-rk*a1(i-j+1)
   60 continue
c
c     calculate k+1 order forward and backward predective error ef(j) 
c      and  eb(j)
c
      do 80 j=i,npts
      tmp=ef(j)
      ef(j)=ef(j)-rk*eb(j-i+1)
      eb(j-i+1)=eb(j-i+1)-rk*tmp
   80 continue

      do 85 j=1,i
      a1(j)=a2(j)
   85 continue

c
c     computer k+1 order reflection coeffient(rz, rr, rt) about receiver 
c       function  from k order
c
      bot=0.
      ztop=0.
      rtop=0.
      ttop=0.
      do 90 j=i,npts
      bot=bot+eb(j-i+1)*eb(j-i+1)
      ztop=ztop+sz(j)*eb(j-i+1)
      rtop=rtop+sr(j)*eb(j-i+1)
      ttop=ttop+st(j)*eb(j-i+1)
   90 continue
      bot=bot+0.5*stable
      rz=ztop/bot
      rr=rtop/bot
      rt=ttop/bot

c
c     computer k+1 order receiver function hz(j), hr(j), ht(j) from
c       k order
c
      do 100 j=1,i
      hz(j)=hz(j)+rz*a1(i-j+1)
      hr(j)=hr(j)+rr*a1(i-j+1)
      ht(j)=ht(j)+rt*a1(i-j+1)
  100 continue
c
c      computer k+1 order filter residual sz(j), sr(j), st(j) from k order
c
      do 110 j=i,npts
      sz(j)=sz(j)-rz*eb(j-i+1)
      sr(j)=sr(j)-rr*eb(j-i+1)
      st(j)=st(j)-rt*eb(j-i+1)
  110 continue
   30 continue

      do 120 j=1,npts
      dat(j,1)=hz(j)
      dat(j,2)=hr(j)
      dat(j,3)=ht(j)
  120 continue
 	open(8,file='PPRE',status='unknown')

	write(8,*) nzs
       do 511 i=1,nzs
	ymax=0.
	do 565 j=1,npts
	xx=abs(dat(j,i))
	if (xx.gt.ymax) ymax=xx
565	continue
	write(8,*) npts,dt,ymax
	write(8,*) (dat(j,i),j=1,npts)
511     continue
	close(8) 		
c
c     filter the receiver function with Gaussian filter to remove high 
c     frequency 
c
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 18 i=1,3
c      do 83 j=1,npts
c83	dc(j)=dat(j,i)      
      call dfftr(dat(1,i),nft,'forward',dt)
      do 19 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=dat(2*(j-1)+1,i)
      cc2=dat(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(cmplx(0.,w*beg))*exp(gauss)
   19 continue
c      do 84 j=1,npts
c84	dcc(j)=d(j,i)      
      call dfftr(d(1,i),nft,'inverse',cdelf)
      do 101 j=1,nfpts
      dat(2*(j-1)+1,i)=real(d(j,i))
      dat(2*j,i)=aimag(d(j,i))
 101  continue
  18  continue
  	write(*,*)'beg:',beg
c
c     normalize receiver function with the maximum of vertical component
c
c      do 71 j=1,npts
c71	dc(j)=dat(j,1)      
      call minmax(dat(1,1),npts,dmin,dmax,dmean)
      do 108 j=1,3
      do 108 i=1,npts
      dat(i,j)=dat(i,j)/dmax
108   continue
c
c      output the receiver function
c
	open(8,file=outfil,status='unknown')
	nz=2
	write(8,*) nz
       do 11 i=2,3
	ymax=0.
	do 65 j=1,npts
	xx=abs(dat(j,i))
	if (xx.gt.ymax) ymax=xx
65	continue
	write(8,*) npts,dt,ymax
	write(8,*) (dat(j,i),j=1,npts)
c        outfil(1:iblank1+4)=evnam(1:iblank1)//outc(i)
c        call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
   11  continue
	close(8) 		
c      yes=yesno('Try another (y or n)? ')
c      if(yes) go to 2
      stop
      end
c
c***********************************************************
c
      subroutine window(b,npts,dt,ndat)
      dimension dat(9000,3)
      common  /win/dat
      data pi/3.1415926/
      bb=pi/b
      nend=ifix(b/dt+.5)+1
      do 1 i=1,nend
      t=float(i-1)*dt
      windo=0.5*(1.+cos(bb*t+pi))
      do 2 j=1,ndat
      dat(i,j)=dat(i,j)*windo
      dat(npts+1-i,j)=dat(npts+1-i,j)*windo
    2 continue
    1 continue
      return
      end

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

      subroutine rotate(x,m,n,baz,az,npts)
c
c   rotates horz. components of x into radial & tangential components
c      given the back azimuth and their orientations
c
c   baz   = back azimuth from station to source in degrees
c   az(i) = + direction of each horz. comp.
c
c      conventions --
c
c          radial is positive away form source
c          tangential is positive clockwise from + radial direction
c
      dimension x(m,n),az(1)
      integer ounit
      common /innout/ inunit,ounit
      rad(deg)=deg/57.295779
      azck=az(2) + 90.
      diff=abs(azck - az(3))
      if(diff.lt..01) go to 1
      if(diff.lt.179.) write(*,100) az(2),az(3),azck,diff
         do 2 i=1,npts
    2     x(i,3)=-x(i,3)
    1 a=sin(rad(baz) - rad(az(2)))
      b=cos(rad(baz) - rad(az(2)))
      do 4 i=1,npts
         radial = -x(i,2)*b - x(i,3)*a
         trans  =  x(i,2)*a - x(i,3)*b
         x(i,2) = radial
         x(i,3) = trans
    4 continue
      return
  100 format(' p r o b l e m   i n   r o t a t e ',/,
     *       1x,'+ horz direction: 1= ',f8.4,' 2= ',f8.4,/,
     *       1x,'for proper rotation 2= ',f8.4,' difference = ',f8.4)
      end
      subroutine minmax(x,npts,min,max,mean)
      dimension x(1)
      real min,max,mean
      min=9.0e+19
      max=-9.0e+19
      mean=0.
      do 1 i=1,npts
           if(x(i).gt.max) max=x(i)
           if(x(i).lt.min) min=x(i)
           mean=mean + x(i)
    1 continue
      mean=mean/float(npts)
      return
      end

      function npowr2(n)
c
c finds the next power of 2 .ge.n
c
      ipowr=alog10(2.*float(n)-1.)/.301029996
      if(n.eq.1) ipowr=1
      npowr2=2**ipowr
      return
      end

      subroutine dfftr (x,nft,dirctn,delta)
c                                              a.shakal, 1/78, 15 jul 80
c           this subroutine does a fast fourier transform on a real
c        time series.  it requires 1/2 the storage and e1/2 the time
c        required by a complex fft.
c
c     forward transform, "call dfftr(x,nft,'forward',dt)":
c           input = x(1),x(2),..,x(nft) = real time series of nft points
c          output = x(1),x(2),..,x(nft+2) = nft/2+1 complex spectral poi
c        these spectral points are identical to the first nft/2+1 return
c        by subroutine fft (i.e., pos freq terms).  thus, the coefficien
c        at fj, the j-th frequency point (where fj = (j-1)*delf, j=1,nft
c        and delf = 1/(nft*dt)), is in x(i-1),x(i), where i=2j.  x(1) is
c        dc term, x(2) = 0 (because real time series), x(nft+1) is real
c        of nyquist coef, and x(nft+2) is imaginary part (0 because real
c        series).
c
c     inverse transform, "call dfftr(x,nft,'inverse',delf)":
c        input and output are interchanged.
c
c           if this subroutine is called with 'forward', and then with '
c        and delf of 1/(nft*dt), the original time series is recovered.
c        identical results (but for scaling) can be obtained by calling
c        fft(x,nft,isign), but in fft a real time series must be stored
c        complex array with zero imaginary parts, which requires 2*nft p
c        of array x.  also, the coefs returned by the fft will differ by
c        n-scaling, since fft's leave out the dt,delf of the approximate
c        integrations.  this subroutine calls fft.
c           this subroutine is a modification of the subroutine 'fftr',
c        written by c.frasier.  the principal modifications are:
c             1) the delt,delf of the integrations are included to make
c                a discrete approximation to the fourier transform.
c             2) the storage of the spectrum (on output if forward, or i
c                if inverse) has x(2) = zero, with the nyquist component
c                x(nft+1), with x(nft+2) = 0.
c
      logical forwrd, invrse
      character dirctn*7
      complex  csign, c1, c2, c3, speci, specj
      real x(*)
      integer ounit
      common /innout/ inunit,ounit
      pi = 3.1415927
c
      call locast(dirctn,invrse,forwrd)
c
      nftby2 = nft/2
      if (.not.(forwrd)) go to 20001
c            forward transform..
      call fft (x,nftby2,-1)
      x1 = x(1)
      x(1) = x1 + x(2)
      x(2) = x1 - x(2)
      sign = -1.
      go to 20002
20001 if (.not.(invrse)) go to 10001
c            adjust nyquist element storage for inverse transform
      x(2) = x(nft+1)
      x(nft+1) = 0.
      sign = +1.
      go to 20002
10001 stop 'dirctn bad to dfftr'
c
c           manipulate elements as appropropriate for a 1/2 length
c        complex fft, after the forward fft, or before the inverse.
20002 piovrn = pi*sign/float(nftby2)
      csign = cmplx(0.,sign)
      do 10 i = 3,nftby2,2
      j = nft-i+2
      c1 = cmplx(x(i)+x(j), x(i+1)-x(j+1))
      c2 = cmplx(x(i)-x(j), x(i+1)+x(j+1))
      w = piovrn*float(i/2)
      c3 = cmplx(cos(w),sin(w))*c2
      speci = c1 + csign*c3
      x(i) = real(speci)/2.
      x(i+1) = aimag(speci)/2.
      specj = conjg(c1) + csign*conjg(c3)
      x(j) = real(specj)/2.
      x(j+1) = aimag(specj)/2.
   10 continue
      x(nftby2+2) = -x(nftby2+2)
      if (.not.(forwrd)) go to 20004
c            include dt of integration, for forward transform...
      dt = delta
      do 9000  i = 1,nft
 9000 x(i) = x(i)*dt
c            adjust storage of the nyquist component...
      x(nft+1) = x(2)
      x(nft+2) = 0.
      x(2) = 0.
      go to 20005
20004 if (.not.(invrse)) go to 10002
      x1 = x(1)
      x(1) = (x1+x(2))/2.
      x(2) = (x1-x(2))/2.
c            do the inverse transform...
      call fft (x,nftby2,+1)
c            in the inverse transform, include the df of the integration
c            and a factor of 2 because only doing half the integration
c            (i.e., just over the positive freqs).
      twodf = 2.*delta
      do 9002  i = 1,nft
 9002 x(i) = x(i)*twodf
10002 continue
20005 return
      end

      subroutine fft(dat,nn,isign)
c                                              a.shakal, 1/78, 10 jul 80
c        cooley-tukey 'fast fourier trnasform' in ansi fortran 77.
c
c           transform(j) = sum {data(i)*w**u(i-1)*(j-1)e}, where i and
c        j run from 1 to nn, and w = exp(sign*twopi*sqrtu-1e/nn).
c        data is a one-dimensional complex array (i.e., the real and
c        imaginary parts of the data are located immediately adjacent
c        in storage, such as fortran places them) whose length nn is
c        a power of two.  isign is +1 or -1, giving the sign of the
c        transform.  transform values are returned in array data,
c        replacing the input data.  the time is proportional to
c        n*log2(n), rather than the non-fft n**2.  modified from the
c        fortran ii coding from n.brenner's mit-ll tech rept.
c
      real dat(1)
      pi = 3.1415926
c
      n = 2*nn
      j = 1
      do 5 i = 1,n,2
      if (.not.(i .lt. j)) go to 20001
      tempr = dat(j)
      tempi = dat(j+1)
      dat(j) = dat(i)
      dat(j+1) = dat(i+1)
      dat(i) = tempr
      dat(i+1) = tempi
20001 m = n/2
    3 if (.not.(j .gt. m)) go to 20004
      j = j-m
      m = m/2
      if (m .ge. 2) go to 3
20004 j = j+m
   5  continue
c
c
      mmax = 2
    6 if (.not.(mmax .ge. n)) go to 20007
      return
20007 if (.not.(mmax .lt. n)) go to 10001
      istep = 2*mmax
      pibymx = pi*float(isign)/float(mmax)
c
      do 8 m = 1,mmax,2
      theta = pibymx*float(m-1)
      wr = cos(theta)
      wi = sin(theta)
      do 8 i = m,n,istep
      j = i + mmax
      tempr = wr*dat(j) - wi*dat(j+1)
      tempi = wr*dat(j+1) + wi*dat(j)
      dat(j) = dat(i) - tempr
      dat(j+1) = dat(i+1) - tempi
      dat(i) = dat(i) + tempr
      dat(i+1) = dat(i+1) + tempi
   8  continue
      mmax = istep
      go to 6
10001 continue
20008 return
      end
      subroutine locast(dirctn,invrse,forwrd)
      character dirctn*7
      logical forwrd,invrse
      integer ounit
      common /innout/ inunit,ounit
      if(dirctn.eq.'forward') go to 1
      if(dirctn.eq.'inverse') go to 2
      write(ounit,100)dirctn
  100 format(1x,a7,2x,'is meaningless to dfftr, use forward or inverse
     *only')
      invrse=.false.
      forwrd=.false.
      return
    1 invrse=.false.
      forwrd=.true.
      return
    2 invrse=.true.
      forwrd=.false.
      return
      end
