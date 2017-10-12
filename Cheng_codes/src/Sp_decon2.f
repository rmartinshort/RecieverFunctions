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
      dimension data(20000,3),temp(20000),data1(2000,3)
     *    ,data2(2000),caz(3)
      complex d(10000,3)
      dimension ef(4000),eb(4000),sz(4000),sr(4000),st(4000),hz(4000)
      dimension hr(4000),ht(4000),a1(4000),a2(4000)
      character eqfile*80,comp(3,2)*2,outc(3)*4
      character*32 evnam,stanam,outfil
      logical yes,yesno,rel,rmean,filter,wdat
      integer blank,ounit
      intrinsic aimag
      real min, max
c
c
c **********************************************************************
c
c
      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.n','.e','.z','.r','.t'/
     *    ,outc/'.eqz','.eqr','.eqt'/
      inunit=5
      ounit=6
      max=-9.0e+19
      min=9.0e+19
      pi=3.1415926
    2 call asktxt('specify event name: ',evnam)
      call asktxt('specify station name: ',stanam)
      rel=yesno('Geographic(zne) data (y or n)? ')
      ks=ask('which chanel as decon: 2 (r) or 3 (t) ') 
      stime=ask('S-P time:' )
      btime=ask('pwave beginning time before S (in Sec positive): ')
      etime=ask('pwave ending time after S (in Sec): ')
      swin=ask('S wavelet window (in Sec): ')
      yes=yesno('window data (y or n)?')
      if(yes) blen=ask('length of window (in Sec),blen = ')
      rmean=yesno('remove mean value (y or n)?' )
      filter=yesno('filter raw data (y or n)?' ) 
      agauss=ask('Gaussian scale, a= ')
      c=ask('trough filter,c = ')
      iblank1=blank(evnam)
      iblank2=blank(stanam)
      iblank=iblank1+iblank2+19
      eqfile(1:16)='/mnt1/ding/tele/'
      outfil(1:32)='                                '
      isyntp=2
      if(rel) isyntp=1
      do 1 i=1,3
         eqfile(1:iblank)=eqfile(1:16)//evnam(1:iblank1)
     *   //'/'//stanam(1:iblank2)//comp(i,isyntp)
         call zero(data(1,i),1,20000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,20000,nerr)
         call getfhv('cmpaz',caz(i),nerr)
  1   continue

c
c     remove the mean value from the raw data
c
      if(.not.rmean) goto  1000
      do 22 i=1,3
      call minmax(data(1,i),npts,dmin,dmax,dmean)
      do 23 j=1,npts
      data(j,i)=data(j,i)-dmean
 23   continue
 22   continue

1000  continue

         call getfhv('a',a0,nerr)
         call getfhv('baz',baz,nerr)
      if(.not.yes) go to 4
         call window(blen,npts,dt,3)
 4    continue
c 
c     filter the high frequency from the raw data with Gaussian filter
c
      if(.not.filter) go to 5

      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)

      do 8 i=1,3
      call dfftr(data(1,i),nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=data(2*(j-1)+1,i)
      cc2=data(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(gauss)
    9 continue
      call dfftr(d(1,i),nft,'inverse',delf)
      do 111 j=1,nfpts
      data(2*(j-1)+1,i)=real(d(j,i))
      data(2*j,i)=aimag(d(j,i))
 111  continue
    8 continue
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)
c
c
 5     continue
            
c
c      rotate the raw data in ZNE direction to ZRT direction froming
c      radial and tangent component
c

       if(rel) call rotate(data,20000,3,baz,caz,npts)
c
c      cut off the expected S-wave data from the raw data
c
       k11=(stime-btime-beg)/dt+1
       k12=(stime+etime-beg)/dt+1
       beg=-blen

       do 2101 i=1,3
       call getdata(data(1,i),k11,k12,npts)
 2101  continue
       npts=k12-k11+1
       call window(blen,npts,dt,3)
       k11=btime/dt+1
       k12=npts
       npts1=k12-k11+1
       npts2=(btime+swin)/dt+1
       do 2103 i=1,3
       do 2102 j=1,npts2
       data1(npts2-j+1,i)=data(j,i)
 2102  continue
       call getdata(data(1,i),k11,k12,npts)       
       outfil(1:iblank1+4)=evnam(1:iblank1)//'.1'//comp(i,2)
       call wsac1(outfil,data(1,i),npts1,beg,dt,nerr)
       outfil(1:iblank1+4)=evnam(1:iblank1)//'.2'//comp(i,2)
       call wsac1(outfil,data1(1,i),npts2,beg,dt,nerr)
 2103  continue
       
        k=1
 3001   if(k.gt.2) goto 3000
        if(k.eq.1) then
         npts=npts1
         do 1101 i=1,npts
         temp(i)=data(i,1)
         data(i,1)=data(i,ks)
         data(i,ks)=temp(i)
 1101    continue
        else
         npts=npts2
        do 2105 i=1,npts
         data(i,1)=data1(i,ks)
         data(i,ks)=data1(i,1)
 2105   continue
        endif
c
c     now to calculate predective error filter factor with Levision 
c     iterative algorithm
c
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
      auto0=auto0+data(i,1)*data(i,1)
      crossr0=crossr0+data(i,1)*data(i,2)
      crosst0=crosst0+data(i,1)*data(i,3)
c
c     preliminary forward and backward predective error
c
      ef(i)=data(i,1)
      eb(i)=data(i,1)
   10 continue
      stable=c*auto0
c
c     preliminary receiver funtion value (1th order)
c
      hz(1)=auto0/auto0
      hr(1)=crossr0/auto0
      ht(1)=crosst0/auto0
c
c     preliminary filter residual
c
      do 20 i=1,npts
      sz(i)=data(i,1)-hz(1)*data(i,1)
      sr(i)=data(i,2)-hr(1)*data(i,1)
      st(i)=data(i,3)-ht(1)*data(i,1) 
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
      data(j,1)=hz(j)
      data(j,2)=hr(j)
      data(j,3)=ht(j)
  120 continue
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
      call dfftr(data(1,i),nft,'forward',dt)
      if(k.eq.1) delay=-blen
      if(k.eq.2) delay=-swin
      do 19 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=data(2*(j-1)+1,i)
      cc2=data(2*j,i)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(cmplx(0.,w*delay))*exp(gauss)
   19 continue
      call dfftr(d(1,i),nft,'inverse',cdelf)
      do 101 j=1,nfpts
      data(2*(j-1)+1,i)=real(d(j,i))
      data(2*j,i)=aimag(d(j,i))
 101  continue
  18  continue
c
c     normalize receiver function with the maximum of vertical component
c
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 108 i=1,npts
      data(i,ks)=data(i,ks)/dmax
108   continue
      if(k.eq.1) then
      do 4000 i=1,npts
       data2(i)=data(i,ks)
       if(abs(data2(i)).gt.max) then
        kmax1=i
        max=abs(data2(i))
       else
       endif
4000  continue
      else
      do 3100 i=1,npts
      temp(npts-i+1)=data(i,ks)
3100  continue
      max=-9.0e+19
      do 3102 i=1,npts
      data(i,ks)=temp(i)
       if(abs(temp(i)).gt.max) then
         max=abs(temp(i))
         kmax2=i
       else
       endif
3102  continue
      endif
c
c      output the receiver function
c
        if(k.eq.1) outfil(1:iblank1+6)=evnam(1:iblank1)//'.1'//outc(1)
        if(k.eq.2) outfil(1:iblank1+6)=evnam(1:iblank1)//'.2'//outc(1)
        call wsac1(outfil,data(1,ks),npts,beg,dt,nerr)
      k=k+1
      goto 3001
3000  continue
      npts=kmax2+npts1-kmax1
      do 5000 i=1,npts2
      temp(i)=data(i,ks)
5000  continue
      do 5001 i=-6,6
      temp(kmax2+i)=(data(kmax2+i,ks)+data2(kmax1+i))*0.5
5001  continue
      do 5002 i=kmax2+7,npts
      temp(i)=data2(kmax1+i-kmax2)
5002  continue
      outfil(1:32)='                                '
      outfil(1:iblank1+4)=evnam(1:iblank1)//outc(1)
      call wsac1(outfil,temp,npts,beg,dt,nerr)
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
c
c***********************************************************
c
      subroutine window(b,npts,dt,ndat)
      dimension data(20000,3)
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

      subroutine window_1(data,b,npts,dt)
      dimension data(npts)
      data pi/3.1415926/
      bb=pi/b
      nend=ifix(b/dt+.5)+1
      do 1 i=1,nend
      t=float(i-1)*dt
      windo=0.5*(1.+cos(bb*t+pi))
      data(i)=data(i)*windo
      data(npts+1-i)=data(npts+1-i)*windo
    1 continue
      return
      end
