      dimension caz(3)
      character comp(3,2)*2
      integer blank,ounit
      dimension data(20000,3),data1(20000,3),data2(20000,3)
      character*32 evnam,stanam
      character*80 outfil,infile
      logical yes,yesno
      common /innout/inunit,ounit
      data comp/'.z','.n','.e','.z','.r','.t'/
      inunit=5
      ounit=6
      call asktxt('specify event name: ', evnam)
      call asktxt('specify station name: ',stanam)
      yes=yesno('window data (y or n)?')
      if(yes) blen=ask('length of window (in Sec): ')
      iblank1=blank(evnam)
      iblank2=blank(stanam)
      do 10 i=1,3
      infile(1:iblank1+iblank2+19)='/mnt1/ding/tele/'//
     * evnam(1:iblank1)//'/'//stanam(1:iblank2)//comp(i,1)
      write(*,*) 'infile=', infile
      call zero(data(1,i),1,20000)
      call zero(data1(1,i),1,20000)
      call rsac1(infile,data(1,i),npts,beg,dt,20000,nerr)
      call getfhv('cmpaz',caz(i),nerr)
  10  continue
      call getfhv('a',a0,nerr)
      call getfhv('baz',baz,nerr)
      call rotate(data,20000,3,baz,caz,npts)
      if(.not.yes) goto 4
      xxxx=1.
c      call window(blen,npts,dt,3)
  4   nft=npowr2(npts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
      cdelf=1./float(nft)
      write(ounit,102) npts,nft,fny,delf,dt
 102  format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
     *          'delf=',f8.4,1x,'dt=',f6.3)
      do 100 i=1,3
      do 105 j=1,npts
      data2(j,i)=data(j,i)*data(j,1)
 105  continue
      outfil(1:iblank2+2)=stanam(1:iblank2)//comp(i,2)
      write(*,*) 'outfil=',outfil
      call wsac1(outfil,data2(1,i),npts,beg,dt,nerr)
 100  continue
      end
