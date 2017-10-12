      program theorywiener
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of Wiener
c     filter based on Levinson algorithm
c
c    *** written by Wu Qinju - October 1993
c
c  *************************************************************************
c
      dimension rzz(4000),rzr(4000)
      dimension a(4000),rr(4000),zr(4000),rzz1(4000)
      complex data(2000,3), temp(2000)
      character comp(3)*2,outc(3)*4
      character*32 eqfile,outfil
      logical yes,yesno,rel
      integer blank,ounit
c **********************************************************************
c
c
      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.r','.t'/
     *    ,outc/'.eqz','.eqr','.eqt'/
      inunit=5
      ounit=6
      pi=3.1415926
  2   eqfile(1:32)='                                '
      outfil(1:32)='                                '
      call asktxt('specify event name: ',eqfile)
      ks=ask('which chanel to be decon: 2 or 3  ')
      sptime=ask('Sp-S time, sptime= ')
      agauss=ask('Gaussian scale, a= ')
      c=ask('trough filter, c= ')
      rel=yesno('filter primary data (y or n)? ')
      iblank=blank(eqfile)
      do 1 i=1,3
         eqfile(1:iblank+2)=eqfile(1:iblank)//comp(i)
         call zero(data(1,i),1,4000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,4000,nerr)
  1   continue
      if(rel) then
      yes=yesno('window data (y or n)?')
      if(.not.yes) go to 4
        blen=ask('length of window (in Secs): ')
        call window(blen,npts,dt,3)
 4    nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 18 i=2,3
      call dfftr(data(1,i),nft,'forward',dt)
      do 19 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
c      gauss=-w*w/(4.*agauss*agauss)
c      data(j,i)=data(j,i)*exp(gauss)*exp(cmplx(0.,w*sptime))
      data(j,i)=data(j,i)*exp(cmplx(0.,w*sptime))
   19 continue
      call dfftr(data(1,i),nft,'inverse',delf)
   18 continue
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)
      else
      endif
c
c    *******************************************************************
c
      do 1101 i=1,npts/2
      temp(i)=data(i,1)
      data(i,1)=data(i,ks)
      data(i,ks)=temp(i)
 1101 continue

      call crossfun(data(1,1),data(1,1),rzz,npts,npts)
      call crossfun(data(1,1),data(1,2),rzr,npts,npts)
      do 101 i=1,npts
      rzz1(i)=rzz(i)
 101  continue
      stable=c*rzz(1)
      a(1)=1./(rzz(1)+stable)
      rr(1)=rzr(1)/(rzz(1)+stable)
      zr(1)=1.
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
      call zero(data(1,i),1,4000)
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
      gauss=-w*w/(4.*agauss*agauss)
      data(j,i)=data(j,i)*exp(gauss)*exp(cmplx(0.,w*(beg+sptime)))
   9  continue
      call dfftr(data(1,i),nft,'inverse',cdelf)
 8    continue
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 80 i=1,2
      do 80 j=1,npts
      data(j,i)=data(j,i)/dmax
 80   continue
      do 111 i=2,3
        outfil(1:iblank+11)=eqfile(1:iblank)//'_wiener'//outc(i)
        call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
 111    continue
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
      subroutine window(b,npts,dt,ndat)
      dimension data(4000,3)
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
