      character*32 stanam,evnam,dirdat,stanam1,filgrop
      character*80 file1,evgrop,file(500),file2(500)
      character*7 c1,c2
      dimension x(10000),y(200)
      integer ounit,blank
      common /innout/inunit,ounit
      ounit=6
      inunit=5
      ifil=7
      ifil0=8
      ngrop=0
      nstack=0
      call asktxt('specify station name: ',stanam)
      call asktxt('specify outfil name: ',filgrop)
      iblank1=blank(stanam)
      dirdat(1:25)='/var/data/home/wu/datdir/'
      evgrop(1:31)=dirdat(1:25)//'datfil'
      open(unit = ifil, file = evgrop)
      rewind(ifil)
      nsort = 1
  1   continue
      read(ifil,'(a)',end = 2) evnam
      iblank0=blank(evnam)
      file1(1:iblank0+25)=dirdat(1:25)//evnam(1:iblank0)
      open(unit=ifil0,file=file1)
      rewind(ifil0)
  3   continue
      read(ifil0,'(a)',end=4) stanam1
      if(stanam1.eq.stanam) then
      ngrop=ngrop+1
      file(ngrop)(1:iblank0+iblank1+25)='/var/data/field/pwave/'
     * //evnam(1:iblank0)//'/'//stanam(1:iblank1)//'.e'
      goto 5
      else
      endif
      goto 3
  4   continue
  5   continue
      nsort = nsort + 1
      goto 1
  2   continue
      nsort = nsort- 1
      if(nsort .eq. 0) then
         write(ounit,*) 'No filenames in list(?)'
         stop
      endif
      write(*,*) 'nsort=',nsort,'ngrop=',ngrop
      do 10 i=1,ngrop
      call rsac1(file(i),x,nlen,beg,del,10000,nerr)
      call getfhv('baz',baz,nerr)
      call getfhv('gcarc',gcarc,nerr)
      y(i)=baz
      write(c1,'(f7.3)') baz
      write(c2,'(f7.3)') gcarc
c      if(abs(baz1-baz).le.20.and.abs(gcarc1-garc).le.10.) then
c      nstack=nstack+1
      file2(i)(1:iblank0+18)=file(i)(23:22+iblank0)//'  '//c1//'  '//c2
c      else
c      endif
  10  continue
      do 30 j=1,ngrop-1
      do 20 i=j+1,ngrop
      if(y(i).ge.y(j))go to 20
      amin=y(i)
      y(i)=y(j)
      y(j)=amin
      evnam(1:iblank0+18)=file2(i)
      file2(i)=file2(j)
      file2(j)=evnam
  20  continue
  30  continue     
      i=1
      open(ifil,file=filgrop)
      rewind(ifil)
  12  continue
      write(ifil,'(a)') file2(i)
      i=i+1
      if(i.le.ngrop) goto 12
        end
        
