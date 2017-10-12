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
      dimension data(9000,3),caz(3),dat(9000,3,30)
      complex d(4500,3)
      dimension ef(4000,30),eb(4000,30)
      dimension sz(4000,30),sr(4000,30),st(4000,30)
      dimension hz(4000)
      dimension hr(4000),ht(4000),a1(4000),a2(4000)
      character eqfile*80,comp(3,2)*2,outc(3)*4
      character*32 evnam,stanam,outfil
      logical yes,yesno,rel,wrot, rmean, water,filter
      integer blank,ounit
      intrinsic aimag
      real mean1,mean2,mean3
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
      pi=3.1415926
c    2 call asktxt('specify event name: ',evnam)
c    2 call asktxt('specify event name: ',eqfile)
c      call asktxt('specify station name: ',stanam)
      outfil(1:32)='                                '
      call asktxt('specify outfile name: ',outfil)
      rel=yesno('Geographic(zne) data (y or n)? ')
c      wrot=yesno('write the raw rotated data (y or n)? ')
      bptime=ask('pwave beginning time before A (in Sec positive): ')
      eptime=ask('pwave ending time after A (in Sec): ')
      yes=yesno('window data (y or n)?')
      if(yes) blen=ask('length of window (in Sec),blen = ')
      rmean=yesno('remove mean value (y or n)?' )
      filter=yesno('filter raw data (y or n)?' ) 
      agauss=ask('Gaussian scale, a= ')
      c=ask('trough filter,c = ')
c      iblank1=blank(evnam)
c      iblank2=blank(stanam)
c      iblank=iblank1+iblank2+19
c      eqfile(1:16)='/mnt/tibet/tele/'
c      outfil(1:32)='                                '
      isyntp=2
      if(rel) isyntp=1
      write(*,*) 'File num: '
      read(*,*) num_file
      do 5000 kkk = 1, num_file
      call asktxt('specify event name: ',eqfile)
      iblank = blank(eqfile)
      do 1 i=1,3
c         eqfile(1:iblank)=eqfile(1:16)//evnam(1:iblank1)
c     *   //'/'//stanam(1:iblank2)//comp(i,isyntp)
         eqfile(1:iblank+2)=eqfile(1:iblank)//comp(i,isyntp)
         write(*,*) eqfile
         call zero(dat(1,i,kkk),1,9000)
         call rsac1(eqfile,dat(1,i,kkk),npts,beg,dt,9000,nerr)
         if(rel) call getfhv('cmpaz',caz(i),nerr)
  1   continue

c
c     remove the mean value from the raw data
c
       if(.not.rmean) goto  1000

      do 22 i=1,3
      call minmax(dat(1,i,kkk),npts,dmin,dmax,dmean)
      do 23 j=1,npts
      dat(j,i,kkk)=dat(j,i,kkk)-dmean
23    continue
22    continue

1000  continue

         call getfhv('a',a0,nerr)
         if(rel) call getfhv('baz',baz,nerr)
c       if(yes) call window(blen,npts,dt,3)
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
      call dfftr(dat(1,i,kkk),nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      cc1=dat(2*(j-1)+1,i,kkk)
      cc2=dat(2*j,i,kkk)
      d(j,i)=cmplx(cc1,cc2)
      gauss=-w*w/(4.*agauss*agauss)
      d(j,i)=d(j,i)*exp(gauss)
    9 continue
      call dfftr(d(1,i),nft,'inverse',delf)
      do 111 j=1,nfpts
      dat(2*(j-1)+1,i,kkk)=real(d(j,i))
      dat(2*j,i,kkk)=aimag(d(j,i))
 111  continue
    8 continue
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)

c
c      cut off the expected P-wave data from the raw data
c
5      k11=1
       if(bptime.gt.0.) k11=(a0-bptime-beg)/dt+1
       k12=(eptime+a0-beg)/dt+1
       beg=beg-a0
       if(bptime.gt.0.) beg=-bptime
       do 3 i=1,3
       call getdata(dat(1,i,kkk),k11,k12,npts)
  3    continue
       npts=k12-k11+1

c
c      rotate the raw data in ZNE direction to ZRT direction froming
c      radial and tangent component
c
       do 3002 jjj = 1, npts
       do 3003 mmm = 1,3
          data(jjj,mmm) = dat(jjj,mmm,kkk)
3003   continue
3002   continue
       if(yes) call window(blen,npts,dt,3)
       if(rel) call rotate(data,9000,3,baz,caz,npts)
      do 3004 jjj = 1, npts
      do 3005 mmm = 1, 3
         dat(jjj,mmm,kkk) = data(jjj,mmm)
3005  continue
3004  continue
5000  continue
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
      do 2000 kkk = 1, num_file 
      auto0=auto0+dat(i,1,kkk)*dat(i,1,kkk)
      crossr0=crossr0+dat(i,1,kkk)*dat(i,2,kkk)
      crosst0=crosst0+dat(i,1,kkk)*dat(i,3,kkk)
c   
c     preliminary forward and backward predective error
c
      ef(i,kkk)=dat(i,1,kkk)
      eb(i,kkk)=dat(i,1,kkk)
2000  continue
   10 continue
      stable=c*auto0
c
c     preliminary receiver funtion value (1th order)
c
      hz(1)=1.
      hr(1)=crossr0/auto0
      ht(1)=crosst0/auto0
c
c     preliminary filter residual
c
      do 2001 kkk = 1, num_file
      do 20 i=1,npts
      sz(i,kkk)=dat(i,1,kkk)-hz(1)*dat(i,1,kkk)
      sr(i,kkk)=dat(i,2,kkk)-hr(1)*dat(i,1,kkk)
      st(i,kkk)=dat(i,3,kkk)-ht(1)*dat(i,1,kkk) 
   20 continue
2001  continue
c
c     Begin recurive .........
c
      do 30 i=2,npts
c
c     calculate k+1 order reflection coefficient rk from k order
c
      top=0.
      bot=0.
      do 2002 kkk = 1, num_file
      do 45 j=i,npts
      bot=bot+ef(j,kkk)*ef(j,kkk)+eb(j-i+1,kkk)*eb(j-i+1,kkk)
      top=top+eb(j-i+1,kkk)*ef(j,kkk)
   45 continue
2002  continue
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
      do 2003 kkk = 1, num_file
      do 80 j=i,npts
      tmp=ef(j,kkk)
      ef(j,kkk)=ef(j,kkk)-rk*eb(j-i+1,kkk)
      eb(j-i+1,kkk)=eb(j-i+1,kkk)-rk*tmp
   80 continue
2003  continue

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
      do 2004 kkk = 1, num_file
      do 90 j=i,npts
      bot=bot+eb(j-i+1,kkk)*eb(j-i+1,kkk)
      ztop=ztop+sz(j,kkk)*eb(j-i+1,kkk)
      rtop=rtop+sr(j,kkk)*eb(j-i+1,kkk)
      ttop=ttop+st(j,kkk)*eb(j-i+1,kkk)
   90 continue
2004  continue
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
      do 2005 kkk = 1, num_file
      do 110 j=i,npts
      sz(j,kkk)=sz(j,kkk)-rz*eb(j-i+1,kkk)
      sr(j,kkk)=sr(j,kkk)-rr*eb(j-i+1,kkk)
      st(j,kkk)=st(j,kkk)-rt*eb(j-i+1,kkk)
  110 continue
2005  continue
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
c
c     normalize receiver function with the maximum of vertical component
c
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 108 j=1,3
      do 108 i=1,npts
      data(i,j)=data(i,j)/dmax
108   continue
c
c      output the receiver function
c
      iblank = blank(outfil)
      write(*,*) iblank
       do 11 i=2,3
        outfil(1:iblank+4)=outfil(1:iblank)//outc(i)
        write(*,*) outfil
        call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
   11  continue
c      yes=yesno('Try another (y or n)? ')
c      if(yes) go to 2
      stop
      end
c
c***********************************************************
c
      subroutine window(b,npts,dt,ndat)
      dimension data(9000,3)
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