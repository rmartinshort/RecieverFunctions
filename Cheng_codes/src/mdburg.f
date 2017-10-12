      program pwaveqn
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
      dimension data(8000,3),caz(3)
      double precision ef1(4000),eb1(4000),sr(4000),st(4000),hr(4000)
      double precision ht(4000),ef2(4000),eb2(4000),a1(4000),a2(4000)
      character eqfile*80,comp(3,2)*2,outc(3)*4
      character*32 evnam,stanam,outfil
      logical yes,yesno
      integer blank,ounit
      double precision den1,den2,den3,dent2,denr2,auto0,crosst0,crossr0
      double precision rk,rt1,rr1,den12,stable
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
      common /win/ data
      common /innout/ inunit,ounit
      include 'hdr.inc'
      data
     *     comp/'.z','.n','.e','.z','.r','.t'/
     *    ,outc/'.eqz','.eqr','.eqt'/
      inunit=5
      ounit=6
      pi=3.1415926
    2 call asktxt('specify event name: ',evnam)
      decon=1.
c      agauss=ask('Gaussian scale, a= ')
      iblank1=blank(evnam)
      isyntp=2
      do 1 i=1,3
         eqfile(1:iblank1+2)=evnam(1:iblank1)//comp(i,isyntp)
         call zero(data(1,i),1,8000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,8000,nerr)
  1   continue
c      yes=yesno('window data (y or n)?')
c      if(.not.yes) go to 4
c         blen=ask('length of window (in Secs): ')
c         call window(blen,npts,dt,3)
c    4    nft=npowr2(npts)
c      nfpts=nft/2+1
c      fny=1./(2.*dt)
c      delf=fny/float(nft/2)
c      cdelf=1./float(nft)
c      do 8 i=1,nbs
c      call dfftr(data(1,i),nft,'forward',dt)
c      do 9 j=1,nfpts
c      freq=float(j-1)*delf
c      w=2.*pi*freq
c      gauss=-w*w/(4.*agauss*agauss)
c      data(j,i)=data(j,i)*exp(gauss)
c    9 continue
c      call dfftr(data(1,i),nft,'inverse',cdelf)
c    8 continue
c      write(ounit,102) npts,nft,fny,delf,dt
c  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
c     *         'delf=',f8.4,1x,'dt=',f6.3)
c
c    *******************************************************************
c
      do 22 i=1,3
      call minmax(data(1,i),npts,dmin,dmax,dmean)
      auto0=0.0
      do 21 j=1,npts
      data(j,i)=data(j,i)-dmean
      auto0=auto0+data(j,i)*data(j,i)
   21 continue
   22 continue
      a1(1)=1.
      a2(1)=1.
      auto0=0d0
      crossr0=0d0
      crosst0=0d0
      write(*,*) 'filter length: '
      read(*,*) mm
      do 5 i=2,mm
      a1(i)=0.
      a2(i)=0.
    5 continue
      do 10 i=1,npts
      hr(i)=0.
      ht(i)=0.
      auto0=auto0+data(i,1)*data(i,1)
      crossr0=crossr0+data(i,1)*data(i,2)
      crosst0=crosst0+data(i,1)*data(i,3)
      ef1(i)=data(i,1)
      eb1(i)=data(i,1)
      ef2(i)=ef1(i)
      eb2(i)=eb1(i)
   10 continue
      cov=sqrt(auto0/npts)
      do 12 i=1,npts
c      data(i,1)=data(i,1)/cov
      ef1(i)=data(i,1)
      eb1(i)=data(i,1)
      ef2(i)=ef1(i)
      eb2(i)=eb1(i)
   12 continue
      hr(1)=crossr0/auto0
      ht(1)=crosst0/auto0
      do 20 i=1,npts
      sr(i)=data(i,2)-hr(1)*data(i,1)
      st(i)=data(i,3)-ht(1)*data(i,1) 
   20 continue
      do 30 i=2,mm
      den1=0d0
      den2=0d0
      den3=0d0
      do 45 j=i,npts
      den1=den1+ef1(j)*ef1(j)
      den2=den2+eb1(j-i+1)*eb1(j-i+1)
      den3=den3+eb1(j-i+1)*ef1(j)
   45 continue
      den12=den1+den2
      rk=2.*den3/den12
      do 60 j=1,i
      a2(j)=a1(j)-rk*a1(i-j+1)
   60 continue
      do 80 j=i,npts
      ef2(j)=ef1(j)-rk*eb1(j-i+1)
      eb2(j-i+1)=eb1(j-i+1)-rk*ef1(j)
   80 continue
      do 70 j=i,npts
      ef1(j)=ef2(j)
      eb1(j-i+1)=eb2(j-i+1)
   70 continue
      do 85 j=1,i
      a1(j)=a2(j)
   85 continue
      den1=0d0
      denr2=0d0
      dent2=0d0
      do 90 j=i,npts
      den1=den1+eb2(j-i+1)*eb2(j-i+1)
      denr2=denr2+sr(j)*eb2(j-i+1)
      dent2=dent2+st(j)*eb2(j-i+1)
   90 continue
      rr1=denr2/den1
      rt1=dent2/den1
      do 100 j=1,i
      hr(j)=hr(j)+rr1*a2(i-j+1)
      ht(j)=ht(j)+rt1*a2(i-j+1)
  100 continue
      do 110 j=i,npts
      sr(j)=sr(j)-rr1*eb2(j-i+1)
      st(j)=st(j)-rt1*eb2(j-i+1)
  110 continue
   30 continue
      do 120 j=1,mm
      data(j,1)=a1(j)
      data(j,2)=hr(j)
      data(j,3)=ht(j)
  120 continue
       iblank=blank(outfil)
       do 11 i=1,3
        outfil(1:iblank1+4)=evnam(1:iblank1)//outc(i)
        write(*,*) 'outfil=',outfil
        call minmax(data(1,i),npts,rheader(2),rheader(3),rheader(57))
            call wsac1(outfil,data(1,i),mm,beg,dt,nerr)
   11  continue
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
      subroutine window(b,npts,dt,ndat)
      dimension data(8000,3)
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
