      program wienerdecon
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
      dimension caz(3),rzz(4000),rzr(4000),rzt(4000)
      dimension a(4000),rr(4000),tr(4000),zr(4000),rzz1(4000)
      dimension tmp1(4000), tmp2(4000), tmp3(4000)
      complex data(4500,3)
      character eqfile*80,comp(3,2)*2,outc(3)*4
      character*32 evnam,stanam,outfil
      logical yes,yesno,rel,rel1
      integer blank,ounit
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
      call asktxt('specify outfil name: ',outfil)
      rel=yesno('Geographic(zne) data (y or n)? ')
c      rel1=yesno('write the raw rotated data (y or n)? ')
      b1=ask('pwave beginning time before A (in Sec positive): ')
      e1=ask('pwave ending time after a (in Sec): ')
      agauss=ask('Gaussian scale, a= ')
      c=ask('trough filter, c=')
      yes=yesno('window data (y or n)?')
      if(yes)  blen=ask('length of window (in Secs): ')
c      iblank1=blank(evnam)
c      iblank2=blank(stanam)
c      iblank=iblank1+iblank2+25
c      eqfile(1:22)='/var/data/field/pwave/'
      write(*,*) 'file num: '
      read(*,*) numfile
      do 6666 ifile = 1, numfile
      call asktxt('specify station name: ',eqfile)
      iblank = blank(eqfile)
      isyntp=2
      if(rel) isyntp=1
      do 1 i=1,3
c         eqfile(1:iblank)=eqfile(1:22)//evnam(1:iblank1)
c     *   //'/'//stanam(1:iblank2)//comp(i,isyntp)
         eqfile(1:iblank+2) = eqfile(1:iblank)//comp(i,isyntp)
         write(*,*) eqfile
         call zero(data(1,i),1,9000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,9000,nerr)
         if(rel) call getfhv('cmpaz',caz(i),nerr)
  1   continue
         call getfhv('a',a0,nerr)
         if(rel) call getfhv('baz',baz,nerr)

c 
      do 1111 i=1,3
      call minmax(data(1,i),npts,dmin,dmax,dmean)
      write(*,*) dmin,dmax,dmean
      do  1111 j=1,npts/2+1
      data(j,i)=data(j,i)-cmplx(dmean,dmean)
 1111 continue

c     if(.not.yes) go to 4
c        call window(blen,npts,dt,3)
 4    continue
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
      gauss=-w*w/(4.*agauss*agauss)
      data(j,i)=data(j,i)*exp(gauss)
   19 continue
      call dfftr(data(1,i),nft,'inverse',cdelf)
   18 continue
      write(ounit,102) npts,nft,fny,delf,dt
  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *         'delf=',f8.4,1x,'dt=',f6.3)
       k11=1
       if(b1.gt.0.) k11=(a0-b1-beg)/dt+1
       k12=(e1+a0-beg)/dt+1
       beg=beg-a0
       if(b1.gt.0.) beg=-b1
       do 3 i=1,3
       call getdata(data(1,i),k11,k12,npts)
  3    continue
       npts=k12-k11+1
       if(yes) call window(blen,npts,dt,3)
       if(rel) call rotate(data,9000,3,baz,caz,npts)
c       if(.not.rel1) go to 133
c       do 103 i=1,3
c       outfil(1:iblank1+2)=evnam(1:iblank1)//comp(i,2)
c       call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
c 103   continue
c 133   continue
c
c    *******************************************************************
c
      call crossfun(data(1,1),data(1,1),tmp1,npts,npts)
      call crossfun(data(1,1),data(1,2),tmp2,npts,npts)
      call crossfun(data(1,1),data(1,3),tmp3,npts,npts)
      do 5555 i = 1, npts
         rzz(i) = rzz(i) + tmp1(i)
         rzr(i) = rzr(i) + tmp2(i)
         rzt(i) = rzt(i) + tmp3(i)
5555  continue
6666  continue
      do 101 i=1,npts
      rzz1(i)=rzz(i)
 101  continue
      stable=c*rzz(1)
      a(1)=1./(rzz(1)+stable)
      rr(1)=rzr(1)/(rzz(1)+stable)
      tr(1)=rzt(1)/(rzz(1)+stable)
      zr(1)=rzz1(1)/(rzz(1)+stable)
      do 10 i=2,npts
      a(i)=0.
      rr(i)=0.
      tr(i)=0.
      zr(i)=0.
      coff=0.
      err=0.
      etr=0.
      ezr=0.
      do 20 j=2,i
      coff=coff+a(i-j+1)*rzz(j)
      err=err+rr(i-j+1)*rzz(j)
      etr=etr+tr(i-j+1)*rzz(j)
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
      tr(k)=tr(k)+(rzt(i)-etr)*a(i-k+1)
      zr(k)=zr(k)+(rzz1(i)-ezr)*a(i-k+1)
  50  continue
  10  continue
      nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      do 70 i=1,3
      call zero(data(1,i),1,9000)
70    continue
      do 60 i=1,nfpts
      cc1=zr(2*(i-1)+1)
      cc2=zr(2*i)
      cc3=rr(2*(i-1)+1)
      cc4=rr(2*i)
      cc5=tr(2*(i-1)+1)
      cc6=tr(2*i)
      data(i,1)=cmplx(cc1,cc2)
      data(i,2)=cmplx(cc3,cc4)
      data(i,3)=cmplx(cc5,cc6)
  60  continue
      do 8 i=1,3
      call dfftr(data(1,i),nft,'forward',dt)
      do 9 j=1,nfpts
      freq=float(j-1)*delf
      w=2.*pi*freq
      gauss=-w*w/(4.*agauss*agauss)
      data(j,i)=data(j,i)*exp(gauss)*exp(cmplx(0.,w*beg))
   9  continue
      call dfftr(data(1,i),nft,'inverse',cdelf)
 8    continue
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      do 80 i=2,3
      do 80 j=1,npts
      data(j,i)=data(j,i)/dmax
 80   continue
      iblank = blank(outfil)
      do 111 i=1,3
        outfil(1:iblank+4)=outfil(1:iblank)//outc(i)
        call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
 111    continue
c      yes=yesno('Try another (y or n)? ')
c      if(yes) go to 2
      stop
      end
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
