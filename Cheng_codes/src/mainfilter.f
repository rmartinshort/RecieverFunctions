      program mainfilter
c
c  ********************************************************************
c
c  *************************************************************************
c
      dimension caz(3)
      character eqfile*80,outfil*32,comp(3,2)*2,outc(2)*5,knm(2)*8
      character evnam*32,stanam*32
      complex data(4000,3)
      logical yesno,rel,rel1
      integer blank,ounit
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.n','.e','.z','.r','.t'/
     *    ,outc/'.eqr ','.eqt '/
     *     ,knm/'radial  ','tangentl'/
      outfil(1:32)='                                '
      inunit=5
      ounit=6
      call asktxt('Specify event name: ',evnam)
      call asktxt('specify station name: ',stanam)
      rel=yesno('Geographic(zne) data (y or n)? ')
      rel1=yesno('full or part  data to be filtered (y or n)? ')
      tw=ask('input polization filter window length (in Sec): ')
      eqfile(1:22)='/var/data/field/pwave/' 
      iblank1=blank(evnam)
      iblank2=blank(stanam)
      iblank=iblank1+iblank2+22
      isyntp=2
      if(rel) isyntp=1
c      call asktxt('specify outfil: ',outfil)
c      iblank0=blank(outfil)
      do 1 i=1,3
         call zero(data(1,i),1,8000)
         eqfile(1:iblank + 3)=eqfile(1:22)//evnam(1:iblank1)
     *   //'/'//stanam(1:iblank2)//comp(i,isyntp)
         call rsac1(eqfile,data(1,i),npts,beg,dt,8000,nerr)
         call getfhv('cmpaz',caz(i),nerr)
    1 continue
      call getfhv('a',a0,nerr)
      call getfhv('baz',baz,nerr)
      beg=beg-a0
      ktw=0.5*tw/dt
      if (.not.rel1) then
      b1=ask('p wave ending time(in Sec): ')
      k12=(b1+a0-beg)/dt+1
      k11=1
      do 3 i=1,3
      call getdata(data(1,i),k11,k12,npts)
  3   continue
      npts=k12-k11+1
      else
      endif
      if(rel) then
      call rotate(data,8000,3,baz,caz,npts)
      do 101 i=1,3
      outfil(1:iblank1+2)=evnam(1:iblank1)//comp(i,2)
      call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
101   continue
      else
      endif
c
c    *******************************************************************
c
      call filter(data,8000,3,npts,ktw,beg,dt)
      do 11 i=1,3
      outfil(1:iblank1+5)=evnam(1:iblank1)//'.pf'//comp(i,2)    
      call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
   11 continue
      stop
      end
