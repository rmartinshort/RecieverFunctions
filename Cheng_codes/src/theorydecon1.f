      program burgdecon
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of maximum 
c     entrop deconvolution 
c
c    *** written by Wu qingju - October 1993 ***
c
c  *************************************************************************
c
      dimension data(2000,2)
      complex d(1000,2)
      dimension ef(2000),eb(2000),sz(2000),sr(2000),hz(2000)
      dimension hr(2000),a1(2000),a2(2000)
      character eqfile*32,comp(2)*2,outc(2)*4
      character*32 evnam,outfil
      logical yes,yesno,rel
      integer blank,ounit
      intrinsic aimag
c
c
c **********************************************************************
c
c
      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.r'/
     *    ,outc/'.eqz','.eqr'/
      inunit=5
      ounit=6
      pi=3.1415926
    2 call asktxt('specify event name: ',evnam)
      yes=yesno('window data (y or n)?')
      if(yes) blen=ask('length of window (in Sec),blen = ')
      agauss=ask('Gaussian scale, a= ')
      c=ask('trough filter,c = ')
      iblank=blank(evnam)
      outfil(1:32)='                                '
      if(rel) isyntp=1
      do 1 i=1,2
         eqfile(1:iblank+2)=evnam(1:iblank)//comp(i)
         call zero(data(1,i),1,2000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,4000,nerr)
  1   continue
      if(.not.yes) go to 4
         call window(blen,npts,dt,2)
 4    continue
c
c    *******************************************************************
c
      k11=1
      k12=(5-beg)/dt+1
      call getdata(data(1,1),k11,k12,npts)
      do 1800 i=1,2
      eqfile(1:iblank+3)=eqfile(1:iblank)//'_'//comp(i)
      write(*,*) eqfile(1:iblank+3)
      call wsac1(eqfile,data(1,i),npts,beg,dt,nerr)
1800  continue
      a1(1)=1.
      auto0=0.
      crossr0=0.
      do 10 i=1,npts
      hz(i)=0.
      hr(i)=0.
      a2(i)=0.
      auto0=auto0+data(i,1)*data(i,1)
      crossr0=crossr0+data(i,1)*data(i,2)
      ef(i)=data(i,1)
      eb(i)=data(i,1)
   10 continue
      stable=c*auto0
      hz(1)=auto0/auto0
      hr(1)=crossr0/auto0
      do 20 i=1,npts
      sz(i)=data(i,1)-hz(1)*data(i,1)
      sr(i)=data(i,2)-hr(1)*data(i,1)
   20 continue
      do 30 i=2,npts
      auto=0.0
      top=0.
      bot=0.
      do 45 j=i,npts
      bot=bot+ef(j)*ef(j)+eb(j-i+1)*eb(j-i+1)
      top=top+eb(j-i+1)*ef(j)
   45 continue
      bot=bot+stable
      rk=2.*top/bot
      a1(i)=0.
      do 60 j=1,i
      a2(j)=a1(j)-rk*a1(i-j+1)
   60 continue
      do 80 j=i,npts
      tmp=ef(j)
      ef(j)=ef(j)-rk*eb(j-i+1)
      eb(j-i+1)=eb(j-i+1)-rk*tmp
   80 continue
      do 85 j=1,i
      a1(j)=a2(j)
   85 continue
      bot=0.
      ztop=0.
      rtop=0.
      do 90 j=i,npts
      bot=bot+eb(j-i+1)*eb(j-i+1)
      ztop=ztop+sz(j)*eb(j-i+1)
      rtop=rtop+sr(j)*eb(j-i+1)
   90 continue
      bot=stable+bot
      rz=ztop/bot
      rr=rtop/bot
      do 100 j=1,i
      hz(j)=hz(j)+rz*a1(i-j+1)
      hr(j)=hr(j)+rr*a1(i-j+1)
  100 continue
      do 110 j=i,npts
      sz(j)=sz(j)-rz*eb(j-i+1)
      sr(j)=sr(j)-rr*eb(j-i+1)
  110 continue
   30 continue
      do 120 j=1,npts
      data(j,1)=hz(j)
      data(j,2)=hr(j)
  120 continue
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 18 i=1,2
      call dfftr(data(1,i),nft,'forward',dt)
      do 19 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=data(2*(j-1)+1,i)
      cc2=data(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(cmplx(0.,w*beg))*exp(gauss)
   19 continue
      call dfftr(d(1,i),nft,'inverse',cdelf)
      do 101 j=1,nfpts
      data(2*(j-1)+1,i)=real(d(j,i))
      data(2*j,i)=aimag(d(j,i))
 101  continue
  18  continue
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 108 j=1,2
      do 108 i=1,npts
      data(i,j)=data(i,j)/dmax
108   continue
       do 11 i=2,2
        outfil(1:iblank+4)=evnam(1:iblank)//outc(i)
        write(*,*) outfil
        call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
   11  continue
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
      subroutine window(b,npts,dt,ndat)
      dimension data(2000,2)
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
